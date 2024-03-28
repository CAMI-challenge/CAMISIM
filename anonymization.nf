#!/usr/bin/env nextflow

nextflow.enable.dsl=2

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
    take: samplewise_gsa_ch // tuple val(sample_id), path(gsa)
    take: bam_file_list_per_sample_ch
    take: pooled_gsa_ch
    take: merged_bam_ch
    take: genome_location_file_ch
    take: metadata_ch

    main:

        // get the seed for every sample
        seed_ch = seed_file_read_simulation_ch.splitCsv(sep:'\t', skip:2)

        reads_seed_ch = reads_ch.join(seed_ch)

        if(params.type=="nanosim3") {
            out_shuffle = shuffle(reads_seed_ch)
        } else if(params.type=="art" || params.type=="wgsim") {
            out_shuffle = shuffle_paired_end(reads_seed_ch)
        }

        gs_read_ch = out_shuffle[1].combine(genome_location_file_ch).combine(metadata_ch)
        gs_read_mapping(gs_read_ch)

        // anonymize assembly of every sample
        seed_gsa_ch = seed_file_gsa_ch.splitCsv(sep:'\t', skip:2)
        shuffle_gsa(samplewise_gsa_ch.join(seed_gsa_ch))
        read_start_positions_from_dir_of_bam(bam_file_list_per_sample_ch)

        gs_contig_ch = shuffle_gsa.out[1].join(read_start_positions_from_dir_of_bam.out).combine(genome_location_file_ch).combine(metadata_ch)
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
    container 'biocontainers/biopython:1.81'
    conda "bioconda::biopython"

    input:
    tuple val(sample_id), path(read_files), val(seed)

    output:
    tuple val(sample_id), path(anonymous_reads_file)
    tuple val(sample_id), path(tmp_reads_mapping_file)

    script:
    anonymous_reads_file = 'anonymous_reads.fq'
    tmp_reads_mapping_file = 'tmp_reads_mapping.tsv'
    """
    touch ${anonymous_reads_file}
    touch ${tmp_reads_mapping_file}
    get_seeded_random() { seed="\$1"; openssl enc -aes-256-ctr -pass pass:"\$seed" -nosalt < /dev/zero 2>/dev/null; };
    cat ${read_files} |  sed 'N;N;N;s/\\n/ /g'  | shuf --random-source=<(get_seeded_random ${seed}) | tr " " "\n" | tr -d '\\000' | python3 ${projectDir}/anonymizer.py  -prefix S${sample_id}R -format fastq -map ${tmp_reads_mapping_file} -out ${anonymous_reads_file} -s
    mkdir --parents ${params.outdir}/sample_${sample_id}/reads
    gzip -k ${anonymous_reads_file}
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
    container 'biocontainers/biopython:1.81'
    conda "bioconda::biopython"

    input:
    tuple val(sample_id), path(first_read_files), path(second_read_files), val(seed)

    output:
    tuple val(sample_id), path(anonymous_reads_file)
    tuple val(sample_id), path(tmp_reads_mapping_file)

    script:
    anonymous_reads_file = 'anonymous_reads.fq'
    tmp_reads_mapping_file = 'tmp_reads_mapping.tsv'
    """
    touch ${anonymous_reads_file}
    touch ${tmp_reads_mapping_file}
    cat sample${sample_id}_Genome*1.fq > first_reads.fq
    cat sample${sample_id}_Genome*2.fq > second_reads.fq
    paste -d " " - - - - <first_reads.fq > first_reads_clustered.fq
    paste -d " " - - - - <second_reads.fq > second_reads_clustered.fq
    paste -d ' ' first_reads_clustered.fq second_reads_clustered.fq  > sample${sample_id}_interweaved.fq
    get_seeded_random() { seed="\$1"; openssl enc -aes-256-ctr -pass pass:"\$seed" -nosalt < /dev/zero 2>/dev/null; };
    shuf --random-source=<(get_seeded_random ${seed}) sample${sample_id}_interweaved.fq | tr " " "\n" | tr -d '\\000' | python3 ${projectDir}/anonymizer.py -prefix S${sample_id}R -format fastq -map ${tmp_reads_mapping_file} -out ${anonymous_reads_file}
    mkdir --parents ${params.outdir}/sample_${sample_id}/reads
    gzip -k ${anonymous_reads_file}
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
    container 'biocontainers/biopython:1.81'
    conda "bioconda::biopython"

    input:
    tuple val(sample_id), path(tmp_reads_mapping_file), path(genome_locations_file), path(metadata_file)


    output:
    tuple val(sample_id), path(reads_mapping_file)

    script:
    reads_mapping_file = 'reads_mapping.tsv'
    wgsim = ""
    real_fastq = ""
    if(params.type.equals("nanosim3")) {
        if(params.simulate_fastq_directly){
            real_fastq = "-nanosim_real_fastq"
        }
    } else if(params.type.equals("wgsim")){
            wgsim = "-wgsim"
    }
    """
    touch ${reads_mapping_file}
    python ${projectDir}/scripts/goldstandardfileformat.py -input ${tmp_reads_mapping_file} -genomes ${genome_locations_file} -metadata ${metadata_file} -out ${reads_mapping_file} -projectDir ${projectDir} ${real_fastq} ${wgsim}
    mkdir --parents ${params.outdir}/sample_${sample_id}/reads
    gzip -k ${reads_mapping_file}
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
    container 'biocontainers/biopython:1.81'
    conda "bioconda::biopython"

    input:
    tuple val(sample_id), path(read_files), val(seed)

    output:
    tuple val(sample_id), path(anonymous_gsa_file)
    tuple val(sample_id), path(tmp_reads_mapping_file)

    script:
    anonymous_gsa_file = 'anonymous_gsa.fasta'
    tmp_reads_mapping_file = 'tmp_reads_mapping.tsv'
    """
    touch ${anonymous_gsa_file}
    touch ${tmp_reads_mapping_file}
    get_seeded_random() { seed="\$1"; openssl enc -aes-256-ctr -pass pass:"\$seed" -nosalt < /dev/zero 2>/dev/null; };
    cat ${read_files} |  sed 'N;N;N;s/\\n/ /g'  | shuf --random-source=<(get_seeded_random ${seed}) | tr " " "\n" | tr -d '\\000' | python3 ${projectDir}/anonymizer.py  -prefix S${sample_id}C -format fasta -map ${tmp_reads_mapping_file} -out ${anonymous_gsa_file} -s
    mkdir --parents ${params.outdir}/sample_${sample_id}/contigs
    gzip -k ${anonymous_gsa_file}
    cp ${anonymous_gsa_file}.gz ${params.outdir}/sample_${sample_id}/contigs/
    """
}

/*
* This process shuffles and anonymizes the pooled gsa.
* Takes:
*    A list with the paths to all read files grouped by sample id and the generated seed.
* Output:
*    The anonymous read file for the given sample.
*    The temp reads mapping file for the given sample, containing the read id and the anonymous read id.
 */
process shuffle_pooled_gsa {
    container 'biocontainers/biopython:1.81'
    conda "bioconda::biopython"

    input:
    path(read_files)
    val(seed)

    output:
    tuple path(anonymous_gsa_pooled)
    tuple path(tmp_reads_mapping_file)

    script:
    anonymous_gsa_pooled = 'anonymous_gsa_pooled.fasta'
    tmp_reads_mapping_file = 'tmp_reads_mapping.tsv'
    """
    touch ${anonymous_gsa_pooled}
    touch ${tmp_reads_mapping_file}
    get_seeded_random() { seed="\$1"; openssl enc -aes-256-ctr -pass pass:"\$seed" -nosalt < /dev/zero 2>/dev/null; };
    cat ${read_files} |  sed 'N;N;N;s/\\n/ /g'  | shuf --random-source=<(get_seeded_random ${seed[0]}) | tr " " "\n" | tr -d '\\000' | python3 ${projectDir}/anonymizer.py -prefix PC -format fasta -map ${tmp_reads_mapping_file} -out ${anonymous_gsa_pooled} -s
    mkdir --parents ${params.outdir}
    gzip -k ${anonymous_gsa_pooled}
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
    container 'biocontainers/samtools:1.19.2--h50ea8bc_1'
    conda 'bioconda::samtools'
    
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
        samtools view "\$bamfile" | awk '{print \$1 "\\t" \$4}' >> ${filename}
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
    container 'biocontainers/samtools:1.19.2--h50ea8bc_1'
    conda 'bioconda::samtools'
    
    input:
    path(merged_bam_files)

    output:
    path(filename)

    script:
    filename = "read_start_positions"
    """
    set -o pipefail
    samtools view ${merged_bam_files} | awk '{print \$1 "\\t" \$4}' >> ${filename}
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
process gs_contig_mapping {
    container 'biocontainers/biopython:1.81'
    conda "bioconda::biopython"

    input:
    tuple val(sample_id), path(tmp_contig_mapping_file), path(read_start_positions), path(genome_locations_file), path(metadata_file)


    output:
    tuple val(sample_id), path(gsa_mapping_file)

    script:
    gsa_mapping_file = 'gsa_mapping.tsv'
    reads_mapping_file = 'reads_mapping.tsv'
    wgsim = ""
    real_fastq = ""
    if(params.type.equals("nanosim3")) {
        if(params.simulate_fastq_directly){
            real_fastq = "-nanosim_real_fastq"
        }
    } else if(params.type.equals("wgsim")){
            wgsim = "-wgsim"
    }
    """
    touch ${gsa_mapping_file}
    python ${projectDir}/scripts/goldstandardfileformat.py -contig -input ${tmp_contig_mapping_file} -genomes ${genome_locations_file} -metadata ${metadata_file} -out ${gsa_mapping_file} -projectDir ${projectDir} ${real_fastq} ${wgsim} -read_positions ${read_start_positions}
    mkdir --parents ${params.outdir}/sample_${sample_id}/contigs
    gzip -k ${gsa_mapping_file}
    cp ${gsa_mapping_file}.gz ${params.outdir}/sample_${sample_id}/contigs/
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
    container 'biocontainers/biopython:1.81'
    conda "bioconda::biopython"

    input:
    path(tmp_contig_mapping_file)
    path(read_start_positions)
    path(genome_locations_file)
    path(metadata_file)


    output:
    path(gsa_mapping_file)

    script:
    gsa_mapping_file = 'gsa_pooled_mapping.tsv'
    reads_mapping_file = 'reads_mapping.tsv'
    params.simulate_fastq_directly = false
    wgsim = ""
    real_fastq = ""
    if(params.type.equals("nanosim3")) {
        if(params.simulate_fastq_directly){
            real_fastq = "-nanosim_real_fastq"
        }
    } else if(params.type.equals("wgsim")){
            wgsim = "-wgsim"
    }
    """
    touch ${gsa_mapping_file}
    python ${projectDir}/scripts/goldstandardfileformat.py -contig -input ${tmp_contig_mapping_file} -genomes ${genome_locations_file} -metadata ${metadata_file} -out ${gsa_mapping_file} -projectDir ${projectDir} ${real_fastq} ${wgsim} -read_positions ${read_start_positions}
    mkdir --parents ${params.outdir}
    gzip -k ${gsa_mapping_file}
    cp ${gsa_mapping_file}.gz ${params.outdir}
    """
}
