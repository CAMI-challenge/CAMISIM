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
    take: pooled_gsa_ch
    take: merged_bam_ch
    take: genome_location_file_ch
    take: metadata_ch

    main:

        def types_list = (params.type instanceof List) ? params.type : [params.type]
        def has_se_type = types_list.any { it == "nanosim3" || it == "pbsim3" }
        def has_pe_type = types_list.any { it == "art" || it == "art_modern" || it == "wgsim" }

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
            // row: [sim_type, sample_id, reads_files..., seed]  – single end: 3 elements + seed
            tuple(row[1], (row[2] instanceof List ? row[2].flatten() : row[2]), row[-1])  // (sample_id, read_files, seed)
        }

        pe_reads_seed_ch = reads_seed_ch.filter { row ->
            def sim_type = row[0]
            sim_type == "art" || sim_type == "art_modern" || sim_type == "wgsim"
        }.map { row ->
            // row: [sim_type, sample_id, r1_files, r2_files, seed]
            // r1_files / r2_files may be nested lists after groupTuple → flatten them
            tuple(
                row[1],
                (row[2] instanceof List ? row[2].flatten() : row[2]),
                (row[3] instanceof List ? row[3].flatten() : row[3]),
                row[-1]
            )  // (sample_id, r1, r2, seed)
        }

        if (has_se_type) {
            out_shuffle_se = shuffle(se_reads_seed_ch)
            gs_read_mapping(out_shuffle_se[1].combine(genome_location_file_ch).combine(metadata_ch))
        }
        if (has_pe_type) {
            out_shuffle_pe = shuffle_paired_end(pe_reads_seed_ch)
            gs_read_mapping(out_shuffle_pe[1].combine(genome_location_file_ch).combine(metadata_ch))
        }

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
            .combine(genome_location_file_ch)
            .combine(metadata_ch)
        gs_contig_mapping(gs_contig_ch)

        // anonymize pooled gold standard assembly
        seed_pooled_gsa_ch = seed_file_pooled_gsa_ch.splitCsv(sep:'\t', skip:2)
        shuffle_pooled_gsa(pooled_gsa_ch, seed_pooled_gsa_ch)
        read_start_positions_from_merged_bam(merged_bam_ch)
        pooled_gs_contig_mapping(shuffle_pooled_gsa.out[1], read_start_positions_from_merged_bam.out, genome_location_file_ch, metadata_ch)
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
    tuple val(sample_id), path(read_files), val(seed)

    output:
    tuple val(sample_id), path(anonymous_reads_file)
    tuple val(sample_id), path(tmp_reads_mapping_file)

    script:
    anonymous_reads_file = 'anonymous_reads.fq'
    tmp_reads_mapping_file = 'tmp_reads_mapping.tsv'
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
    mkdir --parents ${params.outdir}/sample_${sample_id}/reads
    pigz -p ${threads} -k ${anonymous_reads_file}
    cp ${anonymous_reads_file}.gz ${params.outdir}/sample_${sample_id}/reads/
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
    tuple val(sample_id), path(first_read_files), path(second_read_files), val(seed)

    output:
    tuple val(sample_id), path(anonymous_reads_file)
    tuple val(sample_id), path(tmp_reads_mapping_file)

    script:
    anonymous_reads_file = 'anonymous_reads.fq'
    tmp_reads_mapping_file = 'tmp_reads_mapping.tsv'
    threads = Math.max(1, ((task.cpus ?: 1) as int))
    """
    set -euo pipefail

    touch ${anonymous_reads_file}
    touch ${tmp_reads_mapping_file}
    cat ${first_read_files} > first_reads.fq
    cat ${second_read_files} > second_reads.fq
    paste -d " " - - - - <first_reads.fq > first_reads_clustered.fq
    paste -d " " - - - - <second_reads.fq > second_reads_clustered.fq
    paste -d ' ' first_reads_clustered.fq second_reads_clustered.fq  > sample${sample_id}_interweaved.fq
    get_seeded_random() { seed="\$1"; openssl enc -aes-256-ctr -pass pass:"\$seed" -nosalt < /dev/zero 2>/dev/null; };
    shuf --random-source=<(get_seeded_random ${seed}) sample${sample_id}_interweaved.fq | tr " " "\\n" | tr -d '\\000' | python3 ${shared_scripts_dir}/anonymizer.py -prefix S${sample_id}R -format fastq -map ${tmp_reads_mapping_file} -out ${anonymous_reads_file}
    mkdir --parents ${params.outdir}/sample_${sample_id}/reads
    pigz -p ${threads} -k ${anonymous_reads_file}
    cp ${anonymous_reads_file}.gz ${params.outdir}/sample_${sample_id}/reads/
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
    tuple val(sample_id), path(tmp_reads_mapping_file), path(genome_locations_file), path(metadata_file)

    output:
    tuple val(sample_id), path(reads_mapping_file)

    script:
    reads_mapping_file = 'reads_mapping.tsv'
    def types_list = (params.type instanceof List) ? params.type : [params.type]
    simulator = ""
    real_fastq = ""
    if (types_list.contains("nanosim3")) {
        if (params.pipeline.equals("metatranscriptomic")) {
            real_fastq = "-nanosim_real_fastq"
        } else if (params.containsKey('simulate_fastq_directly') && params.simulate_fastq_directly) {
            real_fastq = "-nanosim_real_fastq"
        }
    }
    if (types_list.contains("wgsim")) {
        simulator = "-simulator wgsim"
    } else if (types_list.contains("art_modern")) {
        simulator = "-simulator art_modern"
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
    mkdir --parents ${params.outdir}/sample_${sample_id}/reads
    pigz -p ${threads} -k ${reads_mapping_file}
    cp ${reads_mapping_file}.gz ${params.outdir}/sample_${sample_id}/reads/
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
    path gsa_file
    val seed

    output:
    path anonymous_gsa_pooled
    path tmp_reads_mapping_file

    script:
    anonymous_gsa_pooled = 'anonymous_gsa_pooled.fasta'
    tmp_reads_mapping_file = 'tmp_reads_mapping.tsv'
    threads = Math.max(1, ((task.cpus ?: 1) as int))
    """
    touch ${anonymous_gsa_pooled}
    touch ${tmp_reads_mapping_file}
    get_seeded_random() { seed="\$1"; openssl enc -aes-256-ctr -pass pass:"\$seed" -nosalt < /dev/zero 2>/dev/null; };
    cat ${gsa_file} \\
      | paste -d "\\t" - - \\
      | shuf --random-source=<(get_seeded_random ${seed[0]}) \\
      | tr "\\t" "\\n" \\
      | tr -d '\\000' \\
      | python3 ${shared_scripts_dir}/anonymizer.py -prefix PC -format fasta -map ${tmp_reads_mapping_file} -out ${anonymous_gsa_pooled} -s
    mkdir --parents ${params.outdir}
    pigz -p ${threads} -k ${anonymous_gsa_pooled}
    cp ${anonymous_gsa_pooled}.gz ${params.outdir}
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
    path(merged_bam_files)

    output:
    path(filename)

    script:
    filename = "read_start_positions"
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
    def types_list = (params.type instanceof List) ? params.type : [params.type]
    simulator = ""
    real_fastq = ""
    if (types_list.contains("nanosim3")) {
        if (params.containsKey('simulate_fastq_directly') && params.simulate_fastq_directly) {
            real_fastq = "-nanosim_real_fastq"
        }
    }
    if (types_list.contains("wgsim")) {
        simulator = "-simulator wgsim"
    } else if (types_list.contains("art_modern")) {
        simulator = "-simulator art_modern"
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
    path(tmp_contig_mapping_file)
    path(read_start_positions)
    path(genome_locations_file)
    path(metadata_file)

    output:
    path(gsa_mapping_file)

    script:
    gsa_mapping_file = 'gsa_pooled_mapping.tsv'
    def types_list = (params.type instanceof List) ? params.type : [params.type]
    simulator = ""
    real_fastq = ""
    if (types_list.contains("nanosim3")) {
        if (params.containsKey('simulate_fastq_directly') && params.simulate_fastq_directly) {
            real_fastq = "-nanosim_real_fastq"
        }
    }
    if (types_list.contains("wgsim")) {
        simulator = "-simulator wgsim"
    } else if (types_list.contains("art_modern")) {
        simulator = "-simulator art_modern"
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
    mkdir --parents ${params.outdir}
    pigz -p ${threads} -k ${gsa_mapping_file}
    cp ${gsa_mapping_file}.gz ${params.outdir}
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

        // parse seeds: header is 2 lines, then combination_id\tseed
        seed_merged_ch = seed_file_merged_gsa_ch.splitCsv(sep:'\t', skip:2)
            .map { combination_id, seed -> tuple(combination_id.toString(), seed) }

        // index merged gsa channel by stringified combination_id for joining
        merged_gsa_keyed_ch = merged_gsa_ch
            .map { combination_id, sample_ids, gsa -> tuple(combination_id.toString(), sample_ids, gsa) }

        merged_bam_keyed_ch = merged_bam_per_combination
            .map { combination_id, sample_ids, bam -> tuple(combination_id.toString(), bam) }

        // shuffle/anonymize the contigs of the merged GSA per combination
        shuffle_input_ch = merged_gsa_keyed_ch
            .map { combination_id, sample_ids, gsa -> tuple(combination_id, sample_ids, gsa) }
            .join(seed_merged_ch)
        shuffle_merged_gsa(shuffle_input_ch)

        // start positions from the corresponding merged BAM
        read_start_positions_from_merged_combination_bam(merged_bam_keyed_ch)

        // build the gsa_mapping.tsv per combination
        mapping_in_ch = shuffle_merged_gsa.out[1]
            .join(read_start_positions_from_merged_combination_bam.out)
            .combine(genome_location_file_ch)
            .combine(metadata_ch)
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
    tuple val(combination_id), val(sample_ids), path(gsa_file), val(seed)

    output:
    tuple val(combination_id), val(sample_ids), path(anonymous_gsa_file)
    tuple val(combination_id), val(sample_ids), path(tmp_reads_mapping_file)

    script:
    samples_str = sample_ids.join('_')
    anonymous_gsa_file = "anonymous_gsa_merged_samples_${samples_str}.fasta"
    tmp_reads_mapping_file = "tmp_contig_mapping_merged_samples_${samples_str}.tsv"
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
    tuple val(combination_id), path(merged_bam_file)

    output:
    tuple val(combination_id), path(filename)

    script:
    filename = "read_start_positions_combination_${combination_id}"
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
    tuple val(combination_id), val(sample_ids), path(tmp_contig_mapping_file), path(read_start_positions), path(genome_locations_file), path(metadata_file)

    output:
    tuple val(combination_id), val(sample_ids), path(gsa_mapping_file)

    script:
    samples_str = sample_ids.join('_')
    gsa_mapping_file = "gsa_mapping.tsv"
    def types_list = (params.type instanceof List) ? params.type : [params.type]
    simulator = ""
    real_fastq = ""
    if (types_list.contains("nanosim3")) {
        if (params.containsKey('simulate_fastq_directly') && params.simulate_fastq_directly) {
            real_fastq = "-nanosim_real_fastq"
        }
    }
    if (types_list.contains("wgsim")) {
        simulator = "-simulator wgsim"
    } else if (types_list.contains("art_modern")) {
        simulator = "-simulator art_modern"
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
    cp ${gsa_mapping_file}.gz ${params.outdir}/merged_gsa/gsa_mapping_merged_samples_${samples_str}.tsv.gz
    """
}
