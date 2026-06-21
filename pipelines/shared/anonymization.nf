#!/usr/bin/env nextflow

nextflow.enable.dsl=2

shared_scripts_dir = "${projectDir}/pipelines/shared/scripts"

/** 
* This anonymizes reads and the gold standard assembly.
* Takes:
*   reads_ch Channel containing all read files for grouped for every sample.
*   seed_file_read_simulation_ch
* Emits: 
*   A channel containing the merged fastq file and the merged bam file over all genomes for every sample.
**/
workflow anonymization {

    take: reads_ch
    take: seed_file_read_simulation_ch
    take: seed_file_gsa_ch
    take: seed_file_pooled_gsa_ch
    take: samplewise_gsa_ch // tuple val(sim_type), val(sample_id), path(gsa)
    take: bam_file_list_per_sample_ch
    take: pooled_gsa_ch        // tuple val(sim_type), path(gsa_fasta)  – one per type
    take: merged_bam_ch        // tuple val(sim_type), path(bam)        – one per type
    take: genome_location_file_ch
    take: metadata_ch

    main:

        genome_location_file_val_ch = genome_location_file_ch.first()
        metadata_val_ch = metadata_ch.first()

        def types_list = (params.type instanceof List) ? params.type : [params.type]
        def has_se_type = types_list.any { it == "nanosim3" || it == "pbsim3" }
        def has_pe_type = types_list.any { it == "art" || it == "art_modern" || it == "wgsim" || it == "iss" }

        // get the seed for every sample
        seed_ch = seed_file_read_simulation_ch.splitCsv(sep:'\t', skip:2)

        // reads_ch carries tuple(sim_type, sample_id, …reads…)
        // join on sample_id only – seed_ch is tuple(sample_id, seed)
        reads_seed_ch = reads_ch
            .map { row -> tuple(row[1], row) }   // key = sample_id
            .combine(seed_ch, by: 0)             // attach seed
            .map { sample_id, row, seed -> row + [seed] }  // append seed to original tuple

        // Route each type's reads to the correct shuffle process
        se_reads_seed_ch = reads_seed_ch.filter { row ->
            def sim_type = row[0]
            sim_type == "nanosim3" || sim_type == "pbsim3"
        }.map { row ->
            // row: [sim_type, sample_id, reads_files..., seed]
            tuple(
                row[0],
                row[1],
                (row[2] instanceof List ? row[2].flatten() : row[2]),
                row[-1]
            )  // (sim_type, sample_id, read_files, seed)
        }

        pe_reads_seed_ch = reads_seed_ch.filter { row ->
            def sim_type = row[0]
            sim_type == "art" || sim_type == "art_modern" || sim_type == "wgsim" || sim_type == "iss"
        }.map { row ->
            // row: [sim_type, sample_id, r1_files, r2_files, seed]
            tuple(
                row[0],
                row[1],
                (row[2] instanceof List ? row[2].flatten() : row[2]),
                (row[3] instanceof List ? row[3].flatten() : row[3]),
                row[-1]
            )  // (sim_type, sample_id, r1, r2, seed)
        }

        read_mapping_input_ch = Channel.empty()
        if (has_se_type) {
            out_shuffle_se = shuffle(se_reads_seed_ch)
            // gs_read_mapping(out_shuffle_se[1].combine(genome_location_file_val_ch).combine(metadata_val_ch))
            read_mapping_input_ch = read_mapping_input_ch.mix(out_shuffle_se[1])
        }
        if (has_pe_type) {
            out_shuffle_pe = shuffle_paired_end(pe_reads_seed_ch)
            // gs_read_mapping(out_shuffle_pe[1].combine(genome_location_file_val_ch).combine(metadata_val_ch))
            read_mapping_input_ch = read_mapping_input_ch.mix(out_shuffle_pe[1])
        }
        gs_read_mapping(read_mapping_input_ch.combine(genome_location_file_val_ch).combine(metadata_val_ch))

        // anonymize assembly per type per sample
        // samplewise_gsa_ch: tuple(sim_type, sample_id, gsa)
        seed_gsa_ch = seed_file_gsa_ch.splitCsv(sep:'\t', skip:2)
        // join on sample_id (index 1 of the gsa tuple)
        gsa_seed_ch = samplewise_gsa_ch.map { sim_type, sample_id, gsa -> tuple(sample_id, sim_type, gsa) }.combine(seed_gsa_ch, by: 0)
            .map { sample_id, sim_type, gsa, seed -> tuple(sim_type, sample_id, gsa, seed) }
        shuffle_gsa_typed(gsa_seed_ch)
        bam_for_read_positions_ch = bam_file_list_per_sample_ch.map { sim_type, sample_id, bams -> tuple(sample_id, bams) }.groupTuple()
            .map { sample_id, bams_list -> tuple(sample_id, bams_list.flatten()) }
        read_start_positions_from_dir_of_bam(bam_for_read_positions_ch)

        // bam_file_list_per_sample_ch: tuple(sim_type, sample_id, [bam_paths])
        // join on (sim_type, sample_id) with shuffle_gsa_typed.out[1]: tuple(sim_type, sample_id, tmp_mapping)
        gs_contig_ch = shuffle_gsa_typed.out[1]
            .map { sim_type, sample_id, mapping -> tuple(sample_id, sim_type, mapping) }
            .combine(read_start_positions_from_dir_of_bam.out, by: 0)
            .map { sample_id, sim_type, mapping, read_pos -> tuple(sim_type, sample_id, mapping, read_pos) }
            .combine(genome_location_file_val_ch)
            .combine(metadata_val_ch)
        gs_contig_mapping(gs_contig_ch)

        // anonymize pooled gold standard assembly – one per type
        seed_pooled_gsa_ch = seed_file_pooled_gsa_ch.splitCsv(sep:'\t', skip:2)
        // Extract just the seed scalar and wrap it in a value channel so it can
        // be broadcast to all sim_types without the channel closing after one use.
        pooled_seed_val_ch = seed_pooled_gsa_ch.map { row -> row[0] }.first()
        pooled_gsa_with_seed_ch = pooled_gsa_ch.combine(pooled_seed_val_ch)
        shuffle_pooled_gsa(pooled_gsa_with_seed_ch)
        read_start_positions_from_merged_bam(merged_bam_ch)
        pooled_gs_contig_mapping(shuffle_pooled_gsa.out[1].join(read_start_positions_from_merged_bam.out), genome_location_file_val_ch, metadata_val_ch)
}

/*
* This process shuffles and anonymizes single-end reads.
* Takes:
*    A list with the paths to all read files gruoped by sample id and the generated seed.
* Output:
*    The anonymous read file for the given sample.
*    The temp reads mapping file for the given sample, containing the read id and the anonymous read id.
 */
process shuffle {

    conda "conda-forge::biopython=1.83 conda-forge::pigz"

    input:
    tuple val(sim_type), val(sample_id), path(read_files), val(seed)

    output:
    tuple val(sim_type), val(sample_id), path(anonymous_reads_file)
    tuple val(sim_type), val(sample_id), path(tmp_reads_mapping_file)

    script:
    anonymous_reads_file = "anonymous_reads_${sim_type}.fq"
    tmp_reads_mapping_file = "tmp_reads_mapping_${sim_type}.tsv"
    threads = Math.max(1, ((task.cpus ?: 1) as int))
    """
    set -euo pipefail

    touch ${anonymous_reads_file}
    touch ${tmp_reads_mapping_file}
    get_seeded_random() { seed="\$1"; openssl enc -aes-256-ctr -pass pass:"\$seed" -nosalt < /dev/zero 2>/dev/null; };
    cat ${read_files} \
      | paste -d "\\t" - - - - \
      | shuf --random-source=<(get_seeded_random ${seed}) \
      | tr "\\t" "\\n" \
      | python3 ${shared_scripts_dir}/anonymizer.py  -prefix S${sample_id}R -format fastq -map ${tmp_reads_mapping_file} -out ${anonymous_reads_file} -s
    mkdir --parents ${params.outdir}/sample_${sample_id}/reads/${sim_type}
    pigz -p ${threads} -k ${anonymous_reads_file}
    cp ${anonymous_reads_file}.gz ${params.outdir}/sample_${sample_id}/reads/${sim_type}/
    """
}

/*
* This process shuffles and anonymizes paired-end reads.
* Takes:
*    A list with the paths to all read files gruoped by sample id and the generated seed.
* Output:
*    The anonymous read file for the given sample.
*    The temp reads mapping file for the given sample, containing the read id and the anonymous read id.
 */
process shuffle_paired_end {

    conda "conda-forge::biopython=1.83 conda-forge::python=3.11.5 conda-forge::pigz"

    input:
    tuple val(sim_type), val(sample_id), path(first_read_files), path(second_read_files), val(seed)

    output:
    tuple val(sim_type), val(sample_id), path(anonymous_reads_file)
    tuple val(sim_type), val(sample_id), path(tmp_reads_mapping_file)

    script:
    anonymous_reads_file = "anonymous_reads_${sim_type}.fq"
    tmp_reads_mapping_file = "tmp_reads_mapping_${sim_type}.tsv"
    threads = Math.max(1, ((task.cpus ?: 1) as int))
    """
    set -euo pipefail

    touch ${anonymous_reads_file}
    touch ${tmp_reads_mapping_file}
    cat ${first_read_files} > first_reads.fq
    cat ${second_read_files} > second_reads.fq
    paste -d " " - - - - <first_reads.fq > first_reads_clustered.fq
    paste -d " " - - - - <second_reads.fq > second_reads_clustered.fq
    paste -d ' ' first_reads_clustered.fq second_reads_clustered.fq  > sample${sample_id}_${sim_type}_interweaved.fq
    get_seeded_random() { seed="\$1"; openssl enc -aes-256-ctr -pass pass:"\$seed" -nosalt < /dev/zero 2>/dev/null; };
    shuf --random-source=<(get_seeded_random ${seed}) sample${sample_id}_${sim_type}_interweaved.fq | tr " " "\\n" | tr -d '\\000' | python3 ${shared_scripts_dir}/anonymizer.py -prefix S${sample_id}R -format fastq -map ${tmp_reads_mapping_file} -out ${anonymous_reads_file}
    mkdir --parents ${params.outdir}/sample_${sample_id}/reads/${sim_type}
    pigz -p ${threads} -k ${anonymous_reads_file}
    cp ${anonymous_reads_file}.gz ${params.outdir}/sample_${sample_id}/reads/${sim_type}/
    """
}

/*
* This process created a gold standard read mapping file for one sample.
* Takes:
*   The temp reads mapping file for the given sample, containing the read id and the anonymous read id.
*   A file containing all reference genome locations.
*   The metadata file.
* Output:
*    The reads mapping file for the given sample,.
 */
process gs_read_mapping {

    conda "conda-forge::biopython=1.83 conda-forge::python=3.11.5 conda-forge::pigz"

    input:
    tuple val(sim_type), val(sample_id), path(tmp_reads_mapping_file), path(genome_locations_file), path(metadata_file)

    output:
    tuple val(sim_type), val(sample_id), path(reads_mapping_file)

    script:
    reads_mapping_file = "reads_mapping_${sim_type}.tsv"
    simulator = ""
    real_fastq = ""

    if (sim_type == "nanosim3") {
        if (params.pipeline.equals("metatranscriptomic")) {
            real_fastq = "-nanosim_real_fastq"
        } else if (params.nanosim3 && params.nanosim3.containsKey('simulate_fastq_directly') && params.nanosim3.simulate_fastq_directly) {
            real_fastq = "-nanosim_real_fastq"
        }
    } else if (sim_type == "wgsim") {
        simulator = "-simulator wgsim"
    } else if (sim_type == "art_modern") {
        simulator = "-simulator art_modern"
    } else if (sim_type == "iss") {
        simulator = "-simulator iss"
    }

    if(params.pipeline.equals("metatranscriptomic")){
        metatranscriptomic = "-metatranscriptomic"
    } else {
        metatranscriptomic = ""
    }
    threads = Math.max(1, ((task.cpus ?: 1) as int))
    """
    touch ${reads_mapping_file}
    python ${shared_scripts_dir}/goldstandardfileformat.py -input ${tmp_reads_mapping_file} -genomes ${genome_locations_file} -metadata ${metadata_file} -out ${reads_mapping_file} -projectDir ${projectDir} ${real_fastq} ${simulator} ${metatranscriptomic}
    mkdir --parents ${params.outdir}/sample_${sample_id}/reads/${sim_type}
    pigz -p ${threads} -k ${reads_mapping_file}
    cp ${reads_mapping_file}.gz ${params.outdir}/sample_${sample_id}/reads/${sim_type}/
    """
}

/*
* This process shuffles and anonymizes the gsa per sample.
* Takes:
*    A list with the paths to all read files gruoped by sample id and the generated seed.
* Output:
*    The anonymous read file for the given sample.
*    The temp reads mapping file for the given sample, containing the read id and the anonymous read id.
 */
process shuffle_gsa {

    conda "conda-forge::biopython=1.83 conda-forge::pigz"

    input:
    tuple val(sample_id), path(gsa_file), val(seed)

    output:
    tuple val(sample_id), path(anonymous_gsa_file)
    tuple val(sample_id), path(tmp_reads_mapping_file)

    script:
    anonymous_gsa_file = 'anonymous_gsa.fasta'
    tmp_reads_mapping_file = 'tmp_reads_mapping.tsv'
    threads = Math.max(1, ((task.cpus ?: 1) as int))
    """
    touch ${anonymous_gsa_file}
    touch ${tmp_reads_mapping_file}
    get_seeded_random() { seed="\$1"; openssl enc -aes-256-ctr -pass pass:"\$seed" -nosalt < /dev/zero 2>/dev/null; };
    cat ${gsa_file} \\
      | paste -d "\\t" - - \\
      | shuf --random-source=<(get_seeded_random ${seed}) \\
      | tr "\\t" "\\n" \\
      | tr -d '\\000' \\
      | python3 ${shared_scripts_dir}/anonymizer.py  -prefix S${sample_id}C -format fasta -map ${tmp_reads_mapping_file} -out ${anonymous_gsa_file} -s
    mkdir --parents ${params.outdir}/sample_${sample_id}/contigs
    pigz -p ${threads} -k ${anonymous_gsa_file}
    cp ${anonymous_gsa_file}.gz ${params.outdir}/sample_${sample_id}/contigs/
    """
}

/*
* Per-type variant of shuffle_gsa. Carries sim_type so outputs go into
* type-specific subdirectories and file names don't collide across types.
*/
process shuffle_gsa_typed {

    conda "conda-forge::biopython=1.83 conda-forge::pigz"

    input:
    tuple val(sim_type), val(sample_id), path(gsa_file), val(seed)

    output:
    tuple val(sim_type), val(sample_id), path(anonymous_gsa_file)
    tuple val(sim_type), val(sample_id), path(tmp_reads_mapping_file)

    script:
    anonymous_gsa_file = "anonymous_gsa_${sim_type}.fasta"
    tmp_reads_mapping_file = "tmp_reads_mapping_${sim_type}.tsv"
    threads = Math.max(1, ((task.cpus ?: 1) as int))
    """
    touch ${anonymous_gsa_file}
    touch ${tmp_reads_mapping_file}
    get_seeded_random() { seed="\$1"; openssl enc -aes-256-ctr -pass pass:"\$seed" -nosalt < /dev/zero 2>/dev/null; };
    cat ${gsa_file} \\
      | paste -d "\\t" - - \\
      | shuf --random-source=<(get_seeded_random ${seed}) \\
      | tr "\\t" "\\n" \\
      | tr -d '\\000' \\
      | python3 ${shared_scripts_dir}/anonymizer.py  -prefix S${sample_id}C -format fasta -map ${tmp_reads_mapping_file} -out ${anonymous_gsa_file} -s
    mkdir --parents ${params.outdir}/sample_${sample_id}/contigs/${sim_type}
    pigz -p ${threads} -k ${anonymous_gsa_file}
    cp ${anonymous_gsa_file}.gz ${params.outdir}/sample_${sample_id}/contigs/${sim_type}/
    """
}

/*
* This process shuffles and anonymizes the pooled gsa.
* Takes:
*    The path to the pooled gsa file and the generated seed.
* Output:
*    The anonymous gsa file.
*    The temp reads mapping file for the given sample, containing the read id and the anonymous read id.
 */
process shuffle_pooled_gsa {

    conda "conda-forge::biopython=1.83 conda-forge::pigz"

    input:
    tuple val(sim_type), path(gsa_file), val(seed)

    output:
    tuple val(sim_type), path(anonymous_gsa_pooled)
    tuple val(sim_type), path(tmp_reads_mapping_file)

    script:
    anonymous_gsa_pooled = "anonymous_gsa_pooled_${sim_type}.fasta"
    tmp_reads_mapping_file = "tmp_reads_mapping_pooled_${sim_type}.tsv"
    threads = Math.max(1, ((task.cpus ?: 1) as int))
    """
    touch ${anonymous_gsa_pooled}
    touch ${tmp_reads_mapping_file}
    get_seeded_random() { seed="\$1"; openssl enc -aes-256-ctr -pass pass:"\$seed" -nosalt < /dev/zero 2>/dev/null; };
    cat ${gsa_file} \\
      | paste -d "\\t" - - \\
      | shuf --random-source=<(get_seeded_random ${seed}) \\
      | tr "\\t" "\\n" \\
      | tr -d '\\000' \\
      | python3 ${shared_scripts_dir}/anonymizer.py -prefix PC -format fasta -map ${tmp_reads_mapping_file} -out ${anonymous_gsa_pooled} -s
    mkdir --parents ${params.outdir}/pooled_gsa
    pigz -p ${threads} -k ${anonymous_gsa_pooled}
    cp ${anonymous_gsa_pooled}.gz ${params.outdir}/pooled_gsa/
    """
}

/*
* This process parses 'read' start positions from bam files in a directory.
* Takes:
*   The list of bam files per sample id.
* Output:
*    A file containing the read start posotions for the given sample.
 */
process read_start_positions_from_dir_of_bam {

    conda 'bioconda::samtools=1.13'
    
    input:
    tuple val(sample_id), path(list_bam_files)

    output:
    tuple val(sample_id), path(filename)

    when: params.gsa  // Only execute this process when gsa is set to true

    script:
    filename = "read_start_positions"
    """
    set -o pipefail
    for bamfile in ${list_bam_files}; do
        samtools view "\$bamfile" | awk '{print \$3 "\\t" \$4}' >> ${filename}
    done
    """
}

/*
* This process parses 'read' start positions from bam files in a directory.
* Takes:
*   The list of bam files per sample id.
* Output:
*    A file containing the read start posotions for the given sample.
 */
process read_start_positions_from_merged_bam {

    conda 'bioconda::samtools=1.13'

    input:
    tuple val(sim_type), path(merged_bam_files)

    output:
    tuple val(sim_type), path(filename)

    script:
    filename = "read_start_positions_pooled_${sim_type}"
    """
    set -o pipefail
    samtools view ${merged_bam_files} | awk '{print \$3 "\\t" \$4}' >> ${filename}
    """
}

/*
* This process created a gold standard contig mapping file for one sample per type.
* Carries sim_type so outputs go into type-specific subdirectories.
* Takes:
*   The temp reads mapping file for the given sample, containing the read id and the anonymous read id.
*   A file containing all reference genome locations.
*   The metadata file.
* Output:
*    The reads mapping file for the given sample,.
 */
process gs_contig_mapping {

    conda "conda-forge::biopython=1.83 conda-forge::python=3.11.5 conda-forge::pigz"

    input:
    tuple val(sim_type), val(sample_id), path(tmp_contig_mapping_file), path(read_start_positions), path(genome_locations_file), path(metadata_file)

    output:
    tuple val(sim_type), val(sample_id), path(gsa_mapping_file)

    script:
    gsa_mapping_file = "gsa_mapping_${sim_type}.tsv"
    simulator = ""
    real_fastq = ""
    if (sim_type == "nanosim3") {
        if (params.nanosim3 && params.nanosim3.containsKey('simulate_fastq_directly') && params.nanosim3.simulate_fastq_directly) {
            real_fastq = "-nanosim_real_fastq"
        }
    } else if (sim_type == "wgsim") {
        simulator = "-simulator wgsim"
    } else if (sim_type == "art_modern") {
        simulator = "-simulator art_modern"
    } else if (sim_type == "iss") {
        simulator = "-simulator iss"
    }

    if(params.pipeline.equals("metatranscriptomic")){
        metatranscriptomic = "-metatranscriptomic"
    } else {
        metatranscriptomic = ""
    }
    threads = Math.max(1, ((task.cpus ?: 1) as int))
    """
    touch ${gsa_mapping_file}
    python ${shared_scripts_dir}/goldstandardfileformat.py -contig -input ${tmp_contig_mapping_file} -genomes ${genome_locations_file} -metadata ${metadata_file} -out ${gsa_mapping_file} -projectDir ${projectDir} ${real_fastq} ${simulator} ${metatranscriptomic} -read_positions ${read_start_positions}
    mkdir --parents ${params.outdir}/sample_${sample_id}/contigs/${sim_type}
    pigz -p ${threads} -k ${gsa_mapping_file}
    cp ${gsa_mapping_file}.gz ${params.outdir}/sample_${sample_id}/contigs/${sim_type}/
    """
}

/*
* This process created a gold standard read mapping file for one sample.
* Takes:
*   The temp reads mapping file for the given sample, containing the read id and the anonymous read id.
*   A file containing all reference genome locations.
*   The metadata file.
* Output:
*    The reads mapping file for the given sample,.
 */
process pooled_gs_contig_mapping {

    conda "conda-forge::biopython=1.83 conda-forge::python=3.11.5 conda-forge::pigz"

    input:
    tuple val(sim_type), path(tmp_contig_mapping_file), path(read_start_positions)
    path(genome_locations_file)
    path(metadata_file)

    output:
    tuple val(sim_type), path(gsa_mapping_file)

    script:
    gsa_mapping_file = "gsa_pooled_${sim_type}_mapping.tsv"
    simulator = ""
    real_fastq = ""
    if (sim_type == "nanosim3") {
        if (params.nanosim3 && params.nanosim3.containsKey('simulate_fastq_directly') && params.nanosim3.simulate_fastq_directly) {
            real_fastq = "-nanosim_real_fastq"
        }
    } else if (sim_type == "wgsim") {
        simulator = "-simulator wgsim"
    } else if (sim_type == "art_modern") {
        simulator = "-simulator art_modern"
    } else if (sim_type == "iss") {
        simulator = "-simulator iss"
    }

    if(params.pipeline.equals("metatranscriptomic")){
        metatranscriptomic = "-metatranscriptomic"
    } else {
        metatranscriptomic = ""
    }
    threads = Math.max(1, ((task.cpus ?: 1) as int))

    """
    touch ${gsa_mapping_file}
    python ${shared_scripts_dir}/goldstandardfileformat.py -contig -input ${tmp_contig_mapping_file} -genomes ${genome_locations_file} -metadata ${metadata_file} -out ${gsa_mapping_file} -projectDir ${projectDir} ${real_fastq} ${simulator} ${metatranscriptomic} -read_positions ${read_start_positions}
    mkdir --parents ${params.outdir}/pooled_gsa
    pigz -p ${threads} -k ${gsa_mapping_file}
    cp ${gsa_mapping_file}.gz ${params.outdir}/pooled_gsa/
    """
}

/**
* This sub-workflow anonymizes the merged gold standard assemblies (one per
* custom combination defined via `merged_gsa_combinations`). For each
* combination it shuffles the contigs of the merged GSA, anonymizes them and
* produces a `gsa_mapping.tsv.gz` file mapping the anonymized contig ids back
* to their original "labels".
*
* Takes:
*   merged_gsa_ch              tuple val(combination_id), val(sample_ids), path(gsa_fasta)
*   merged_bam_per_combination tuple val(combination_id), val(sample_ids), path(merged_bam)
*   seed_file_merged_gsa_ch    path  seed_merged_gsa_anonymisation.txt
*   genome_location_file_ch    path  genome_locations.tsv
*   metadata_ch                path  metadata.tsv
**/
workflow anonymize_merged_gsa {

    take: merged_gsa_ch
    take: merged_bam_per_combination
    take: seed_file_merged_gsa_ch
    take: genome_location_file_ch
    take: metadata_ch

    main:

        genome_location_file_val_ch = genome_location_file_ch.first()
        metadata_val_ch = metadata_ch.first()

        // parse seeds: header is 2 lines, then combination_id\tseed
        seed_merged_ch = seed_file_merged_gsa_ch.splitCsv(sep:'\t', skip:2).map { combination_id, seed -> tuple(combination_id.toString(), seed) }

        // index merged gsa channel by stringified combination_id for joining
        merged_gsa_keyed_ch = merged_gsa_ch.map { sim_label, combination_id, sample_ids, gsa -> tuple(combination_id.toString(), sim_label, sample_ids, gsa) }

        merged_bam_keyed_ch = merged_bam_per_combination.map { sim_label, combination_id, sample_ids, bam -> tuple(combination_id.toString(), sim_label, bam) }

        // shuffle/anonymize the contigs of the merged GSA per combination
        shuffle_input_ch = merged_gsa_keyed_ch
            .map { combination_id, sim_label, sample_ids, gsa -> tuple(combination_id, sim_label, sample_ids, gsa) }
            .combine(seed_merged_ch, by: 0)
            .map { combination_id, sim_label, sample_ids, gsa, seed -> tuple(sim_label, combination_id, sample_ids, gsa, seed) }
        shuffle_merged_gsa(shuffle_input_ch)

        // start positions from the corresponding merged BAM
        read_start_positions_from_merged_combination_bam(merged_bam_keyed_ch)

        // build the gsa_mapping.tsv per combination
        mapping_in_ch = shuffle_merged_gsa.out[1]
            .map { sim_label, combination_id, sample_ids, mapping -> tuple(tuple(sim_label, combination_id), sample_ids, mapping) }
            .join(read_start_positions_from_merged_combination_bam.out.map { sim_label, combination_id, read_pos -> tuple(tuple(sim_label, combination_id), read_pos) })
            .map { key, sample_ids, mapping, read_pos -> tuple(key[0], key[1], sample_ids, mapping, read_pos) }
            .combine(genome_location_file_val_ch)
            .combine(metadata_val_ch)
        merged_gs_contig_mapping(mapping_in_ch)
}

/*
* This process shuffles and anonymizes the merged gsa of one custom sample
* combination.
* Takes:
*    A tuple with the combination id, the list of sample ids, the merged gsa
*    file and the generated seed.
* Output:
*    The anonymous merged gsa file (per combination).
*    The temp contig mapping file (per combination), containing the original
*    contig id and the anonymous contig id.
*/
process shuffle_merged_gsa {

    conda "conda-forge::biopython=1.83 conda-forge::pigz"

    input:
    tuple val(sim_label), val(combination_id), val(sample_ids), path(gsa_file), val(seed)

    output:
    tuple val(sim_label), val(combination_id), val(sample_ids), path(anonymous_gsa_file)
    tuple val(sim_label), val(combination_id), val(sample_ids), path(tmp_reads_mapping_file)

    script:
    samples_str = sample_ids.join('_')
    anonymous_gsa_file = "anonymous_gsa_merged_samples_${samples_str}_${sim_label}.fasta"
    tmp_reads_mapping_file = "tmp_contig_mapping_merged_samples_${samples_str}_${sim_label}.tsv"
    threads = Math.max(1, ((task.cpus ?: 1) as int))
    """
    touch ${anonymous_gsa_file}
    touch ${tmp_reads_mapping_file}
    get_seeded_random() { seed="\$1"; openssl enc -aes-256-ctr -pass pass:"\$seed" -nosalt < /dev/zero 2>/dev/null; };
    cat ${gsa_file} \\
      | paste -d "\\t" - - \\
      | shuf --random-source=<(get_seeded_random ${seed}) \\
      | tr "\\t" "\\n" \\
      | tr -d '\\000' \\
      | python3 ${shared_scripts_dir}/anonymizer.py -prefix M${combination_id}C -format fasta -map ${tmp_reads_mapping_file} -out ${anonymous_gsa_file} -s
    mkdir --parents ${params.outdir}/merged_gsa
    pigz -p ${threads} -k ${anonymous_gsa_file}
    cp ${anonymous_gsa_file}.gz ${params.outdir}/merged_gsa/
    """
}

/*
* This process parses 'read' start positions from a merged BAM file produced
* for one custom sample combination.
* Takes:
*   A tuple with the combination id and the merged BAM file.
* Output:
*   A tuple with the combination id and the file containing the read start
*   positions.
*/
process read_start_positions_from_merged_combination_bam {

    conda 'bioconda::samtools=1.13'

    input:
    tuple val(combination_id), val(sim_label), path(merged_bam_file)

    output:
    tuple val(sim_label), val(combination_id), path(filename)

    script:
    filename = "read_start_positions_combination_${combination_id}_${sim_label}"
    """
    set -o pipefail
    samtools view ${merged_bam_file} | awk '{print \$3 "\\t" \$4}' >> ${filename}
    """
}

/*
* This process creates a gold standard contig mapping file for a merged gsa of
* one custom sample combination.
* Takes:
*   The temp contig mapping file, the read start positions of the merged BAM,
*   the genome locations file and the metadata file (all per combination).
* Output:
*   The gsa_mapping.tsv file for the given combination.
*/
process merged_gs_contig_mapping {

    conda "conda-forge::biopython=1.83 conda-forge::python=3.11.5 conda-forge::pigz"

    input:
    tuple val(sim_label), val(combination_id), val(sample_ids), path(tmp_contig_mapping_file), path(read_start_positions), path(genome_locations_file), path(metadata_file)

    output:
    tuple val(sim_label), val(combination_id), val(sample_ids), path(gsa_mapping_file)

    script:
    samples_str = sample_ids.join('_')
    gsa_mapping_file = "gsa_mapping.tsv"

    // Resolve target types
    // If this label is hybrid (starts with "hybrid_"), we parse the hybridized types
    def is_hybrid = sim_label.startsWith("hybrid_")
    def current_types = []
    if (is_hybrid) {
        current_types = sim_label.substring(7).split('_') as List
    } else {
        current_types = [sim_label]
    }

    def types_list = (params.type instanceof List) ? params.type : [params.type]
    simulator = ""
    real_fastq = ""
    if (current_types.contains("nanosim3")) {
        if (params.nanosim3 && params.nanosim3.containsKey('simulate_fastq_directly') && params.nanosim3.simulate_fastq_directly) {
            real_fastq = "-nanosim_real_fastq"
        }
    }
    if (current_types.size() == 1) {
        def sim_type = current_types[0]
        if (sim_type == "wgsim") {
            simulator = "-simulator wgsim"
        } else if (sim_type == "art_modern") {
            simulator = "-simulator art_modern"
        } else if (sim_type == "iss") {
            simulator = "-simulator iss"
        }
    }

    if(params.pipeline.equals("metatranscriptomic")){
        metatranscriptomic = "-metatranscriptomic"
    } else {
        metatranscriptomic = ""
    }
    threads = Math.max(1, ((task.cpus ?: 1) as int))
    """
    touch ${gsa_mapping_file}
    python ${shared_scripts_dir}/goldstandardfileformat.py -contig -input ${tmp_contig_mapping_file} -genomes ${genome_locations_file} -metadata ${metadata_file} -out ${gsa_mapping_file} -projectDir ${projectDir} ${real_fastq} ${simulator} ${metatranscriptomic} -read_positions ${read_start_positions}
    mkdir --parents ${params.outdir}/merged_gsa
    pigz -p ${threads} -k ${gsa_mapping_file}
    # Use the descriptive per-combination name in the output directory to avoid
    # collisions between combinations.
    cp ${gsa_mapping_file}.gz ${params.outdir}/merged_gsa/gsa_mapping_merged_samples_${samples_str}_${sim_label}.tsv.gz
    """
}

/**
* This sub-workflow anonymizes the hybrid gold standard assemblies (one per
* sample). For each sample it shuffles the contigs of the hybrid GSA,
* anonymizes them and produces a `gsa_mapping.tsv.gz` file mapping the
* anonymized contig ids back to their original labels.
*
* Takes:
*   hybrid_gsa_ch              tuple val(sample_id), path(hybrid_gsa_fasta)
*   hybrid_bam_ch              tuple val(sample_id), path(hybrid_bam)
*   seed_file_hybrid_gsa_ch    path  seed_hybrid_gsa_anonymisation.txt
*   genome_location_file_ch    path  genome_locations.tsv
*   metadata_ch                path  metadata.tsv
**/
workflow anonymize_hybrid_gsa {

    take: hybrid_gsa_ch
    take: hybrid_bam_ch
    take: seed_file_hybrid_gsa_ch
    take: genome_location_file_ch
    take: metadata_ch
    take: hybrid_types

    main:

        genome_location_file_val_ch = genome_location_file_ch.first()
        metadata_val_ch = metadata_ch.first()

        // build a stable label from the type names, e.g. "iss_nanosim3"
        types_label = (hybrid_types instanceof List ? hybrid_types : [hybrid_types]).join('_')

        // parse seeds: header is 2 lines, then sample_id\tseed
        seed_hybrid_ch = seed_file_hybrid_gsa_ch.splitCsv(sep:'\t', skip:2).map { sample_id, seed -> tuple(sample_id.toString(), seed) }

        // key by stringified sample_id for joining
        hybrid_gsa_keyed_ch = hybrid_gsa_ch.map { sample_id, gsa -> tuple(sample_id.toString(), gsa) }

        hybrid_bam_keyed_ch = hybrid_bam_ch.map { sample_id, bam -> tuple(sample_id.toString(), bam) }

        // shuffle/anonymize the contigs of the hybrid GSA per sample
        shuffle_input_ch = hybrid_gsa_keyed_ch.join(seed_hybrid_ch)
        shuffle_hybrid_gsa(shuffle_input_ch, types_label)

        // start positions from the hybrid BAM
        read_start_positions_from_hybrid_bam(hybrid_bam_keyed_ch)

        // build the gsa_mapping.tsv per sample
        mapping_in_ch = shuffle_hybrid_gsa.out[1]
            .join(read_start_positions_from_hybrid_bam.out)
            .combine(genome_location_file_val_ch)
            .combine(metadata_val_ch)
        hybrid_gs_contig_mapping(mapping_in_ch, types_label)
}

/*
* This process shuffles and anonymizes the hybrid gsa of one sample.
* Takes:
*    A tuple with the sample id, the hybrid gsa file and the generated seed.
* Output:
*    The anonymous hybrid gsa file (per sample).
*    The temp contig mapping file (per sample), containing the original
*    contig id and the anonymous contig id.
*/
process shuffle_hybrid_gsa {

    conda "conda-forge::biopython=1.83 conda-forge::pigz"

    input:
    tuple val(sample_id), path(gsa_file), val(seed)
    val types_label

    output:
    tuple val(sample_id), path(anonymous_gsa_file)
    tuple val(sample_id), path(tmp_contig_mapping_file)

    script:
    anonymous_gsa_file = "anonymous_gsa_hybrid_${types_label}_sample_${sample_id}.fasta"
    tmp_contig_mapping_file = "tmp_contig_mapping_hybrid_${types_label}_sample_${sample_id}.tsv"
    threads = Math.max(1, ((task.cpus ?: 1) as int))
    """
    touch ${anonymous_gsa_file}
    touch ${tmp_contig_mapping_file}
    get_seeded_random() { seed="\$1"; openssl enc -aes-256-ctr -pass pass:"\$seed" -nosalt < /dev/zero 2>/dev/null; };
    cat ${gsa_file} \\
      | paste -d "\\t" - - \\
      | shuf --random-source=<(get_seeded_random ${seed}) \\
      | tr "\\t" "\\n" \\
      | tr -d '\\000' \\
      | python3 ${shared_scripts_dir}/anonymizer.py -prefix H${sample_id}C -format fasta -map ${tmp_contig_mapping_file} -out ${anonymous_gsa_file} -s
    mkdir --parents ${params.outdir}/sample_${sample_id}/hybrid_gsa
    pigz -p ${threads} -k ${anonymous_gsa_file}
    cp ${anonymous_gsa_file}.gz ${params.outdir}/sample_${sample_id}/hybrid_gsa/
    """
}

/*
* This process parses 'read' start positions from the hybrid BAM file of one
* sample (combining all sequencing types).
* Takes:
*   A tuple with the sample id and the hybrid BAM file.
* Output:
*   A tuple with the sample id and the file containing the read start positions.
*/
process read_start_positions_from_hybrid_bam {

    conda 'bioconda::samtools=1.13'

    input:
    tuple val(sample_id), path(hybrid_bam_file)

    output:
    tuple val(sample_id), path(filename)

    script:
    filename = "read_start_positions_hybrid_${sample_id}"
    """
    set -o pipefail
    samtools view ${hybrid_bam_file} | awk '{print \$3 "\\t" \$4}' >> ${filename}
    """
}

/*
* This process creates a gold standard contig mapping file for the hybrid gsa
* of one sample.
* Takes:
*   The temp contig mapping file, the read start positions of the hybrid BAM,
*   the genome locations file and the metadata file (all per sample).
* Output:
*   The gsa_mapping.tsv file for the given sample.
*/
process hybrid_gs_contig_mapping {

    conda "conda-forge::biopython=1.83 conda-forge::python=3.11.5 conda-forge::pigz"

    input:
    tuple val(sample_id), path(tmp_contig_mapping_file), path(read_start_positions), path(genome_locations_file), path(metadata_file)
    val types_label

    output:
    tuple val(sample_id), path(gsa_mapping_file)

    script:
    gsa_mapping_file = "gsa_hybrid_${types_label}_mapping_sample_${sample_id}.tsv"
    def types_list = (params.type instanceof List) ? params.type : [params.type]
    def hybrid_param = params.containsKey('hybrid') && params.hybrid ? (params.hybrid instanceof Boolean ? types_list : (params.hybrid instanceof List ? params.hybrid : [params.hybrid])) : []
    def hybrid_types = hybrid_param.intersect(types_list)
    simulator = ""
    real_fastq = ""
    if (hybrid_types.contains("nanosim3")) {
        if (params.nanosim3 && params.nanosim3.containsKey('simulate_fastq_directly') && params.nanosim3.simulate_fastq_directly) {
            real_fastq = "-nanosim_real_fastq"
        }
    }
    // Since hybrid processes are only run when hybrid_types.size() > 1,
    // we never have a single simulator type, so simulator remains "".

    if(params.pipeline.equals("metatranscriptomic")){
        metatranscriptomic = "-metatranscriptomic"
    } else {
        metatranscriptomic = ""
    }
    threads = Math.max(1, ((task.cpus ?: 1) as int))
    """
    touch ${gsa_mapping_file}
    python ${shared_scripts_dir}/goldstandardfileformat.py -contig -input ${tmp_contig_mapping_file} -genomes ${genome_locations_file} -metadata ${metadata_file} -out ${gsa_mapping_file} -projectDir ${projectDir} ${real_fastq} ${simulator} ${metatranscriptomic} -read_positions ${read_start_positions}
    mkdir --parents ${params.outdir}/sample_${sample_id}/hybrid_gsa
    pigz -p ${threads} -k ${gsa_mapping_file}
    cp ${gsa_mapping_file}.gz ${params.outdir}/sample_${sample_id}/hybrid_gsa/
    """
}

/**
* This sub-workflow anonymizes the hybrid pooled gold standard assembly.
* It shuffles the contigs of the hybrid pooled GSA, anonymizes them and
* produces a `gsa_pooled_hybrid_<types>_mapping.tsv.gz` file mapping the
* anonymized contig ids back to their original labels.
**/
workflow anonymize_hybrid_pooled_gsa {

    take: hybrid_pooled_gsa_ch
    take: hybrid_pooled_bam_ch
    take: seed_file_pooled_gsa_ch
    take: genome_location_file_ch
    take: metadata_ch
    take: hybrid_types

    main:

        genome_location_file_val_ch = genome_location_file_ch.first()
        metadata_val_ch = metadata_ch.first()

        // build a stable label from the type names, e.g. "iss_nanosim3"
        types_label = (hybrid_types instanceof List ? hybrid_types : [hybrid_types]).join('_')

        // parse seeds: same seed as pooled gsa
        seed_pooled_gsa_ch = seed_file_pooled_gsa_ch.splitCsv(sep:'\t', skip:2)
        pooled_seed_val_ch = seed_pooled_gsa_ch.map { row -> row[0] }.first()

        // shuffle/anonymize the contigs of the hybrid pooled GSA
        shuffle_hybrid_pooled_gsa(hybrid_pooled_gsa_ch.combine(pooled_seed_val_ch), types_label)

        // start positions from the hybrid pooled BAM
        read_start_positions_from_hybrid_pooled_bam(hybrid_pooled_bam_ch, types_label)

        // build the gsa_mapping.tsv for hybrid pooled
        mapping_in_ch = shuffle_hybrid_pooled_gsa.out[1]
            .join(read_start_positions_from_hybrid_pooled_bam.out)
            .combine(genome_location_file_val_ch)
            .combine(metadata_val_ch)
        hybrid_pooled_gs_contig_mapping(mapping_in_ch, types_label)
}

/*
* This process shuffles and anonymizes the hybrid pooled gsa.
*/
process shuffle_hybrid_pooled_gsa {

    conda "conda-forge::biopython=1.83 conda-forge::pigz"

    input:
    tuple path(gsa_file), val(seed)
    val types_label

    output:
    path anonymous_gsa_file
    tuple val(types_label), path(tmp_contig_mapping_file)

    script:
    anonymous_gsa_file = "anonymous_gsa_pooled_hybrid_${types_label}.fasta"
    tmp_contig_mapping_file = "tmp_contig_mapping_pooled_hybrid_${types_label}.tsv"
    threads = Math.max(1, ((task.cpus ?: 1) as int))
    """
    touch ${anonymous_gsa_file}
    touch ${tmp_contig_mapping_file}
    get_seeded_random() { seed="\$1"; openssl enc -aes-256-ctr -pass pass:"\$seed" -nosalt < /dev/zero 2>/dev/null; };
    cat ${gsa_file} \\
      | paste -d "\\t" - - \\
      | shuf --random-source=<(get_seeded_random ${seed}) \\
      | tr "\\t" "\\n" \\
      | tr -d '\\000' \\
      | python3 ${shared_scripts_dir}/anonymizer.py -prefix HPC -format fasta -map ${tmp_contig_mapping_file} -out ${anonymous_gsa_file} -s
    mkdir --parents ${params.outdir}/hybrid_pooled_gsa
    pigz -p ${threads} -k ${anonymous_gsa_file}
    cp ${anonymous_gsa_file}.gz ${params.outdir}/hybrid_pooled_gsa/
    """
}

/*
* This process parses 'read' start positions from the hybrid pooled BAM file.
*/
process read_start_positions_from_hybrid_pooled_bam {

    conda 'bioconda::samtools=1.13'

    input:
    path hybrid_pooled_bam_file
    val types_label

    output:
    tuple val(types_label), path(filename)

    script:
    filename = "read_start_positions_hybrid_pooled_${types_label}"
    """
    set -o pipefail
    samtools view ${hybrid_pooled_bam_file} | awk '{print \$3 "\\t" \$4}' >> ${filename}
    """
}

/*
* This process creates a gold standard contig mapping file for the hybrid pooled gsa.
*/
process hybrid_pooled_gs_contig_mapping {

    conda "conda-forge::biopython=1.83 conda-forge::python=3.11.5 conda-forge::pigz"

    input:
    tuple val(types_label), path(tmp_contig_mapping_file), path(read_start_positions), path(genome_locations_file), path(metadata_file)
    val types_label_val

    output:
    path gsa_mapping_file

    script:
    gsa_mapping_file = "gsa_pooled_hybrid_${types_label}_mapping.tsv"
    def types_list = (params.type instanceof List) ? params.type : [params.type]
    def hybrid_param = params.containsKey('hybrid') && params.hybrid ? (params.hybrid instanceof Boolean ? types_list : (params.hybrid instanceof List ? params.hybrid : [params.hybrid])) : []
    def hybrid_types = hybrid_param.intersect(types_list)
    simulator = ""
    real_fastq = ""

    // Nanosim FastQ check is preserved if nanosim3 is a part of the hybridized types
    if (hybrid_types.contains("nanosim3")) {
        if (params.nanosim3 && params.nanosim3.containsKey('simulate_fastq_directly') && params.nanosim3.simulate_fastq_directly) {
            real_fastq = "-nanosim_real_fastq"
        }
    }
    // Since hybrid processes are only ran when hybrid_types.size() > 1,
    // we never have a single simulator type, so simulator remains "".

    if(params.pipeline.equals("metatranscriptomic")){
        metatranscriptomic = "-metatranscriptomic"
    } else {
        metatranscriptomic = ""
    }
    threads = Math.max(1, ((task.cpus ?: 1) as int))
    """
    touch ${gsa_mapping_file}
    python ${shared_scripts_dir}/goldstandardfileformat.py -contig -input ${tmp_contig_mapping_file} -genomes ${genome_locations_file} -metadata ${metadata_file} -out ${gsa_mapping_file} -projectDir ${projectDir} ${real_fastq} ${simulator} ${metatranscriptomic} -read_positions ${read_start_positions}
    mkdir --parents ${params.outdir}/hybrid_pooled_gsa
    pigz -p ${threads} -k ${gsa_mapping_file}
    cp ${gsa_mapping_file}.gz ${params.outdir}/hybrid_pooled_gsa/
    """
}
