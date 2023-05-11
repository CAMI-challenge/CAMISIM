#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/** 
* This workflow simulates reads for every sample.
* Takes:
*     genome_location_ch: Path to the file containing the genome locations.
      genome_distribution_ch: Paths to files containing the distributions of every genome for every sample.
* Emits: 
*     A channel containing the merged fastq file and the merged bam file over all genomes for every sample.
**/
workflow anonymization {

    take: reads_ch

    main:
        if(params.type=="nanosim3") {
            shuffle(reads_ch)
        } else if(params.type=="art") {
            shuffle_paired_end(reads_ch)
        }
}

/*
* This process merges all given bam files specified in the pooled_gsa parameter.
* Takes:
*     A list with the paths to all bam files, that should be merged, if the condition is fullfilled.
* Output:
*     The path to the merged bam file.
 */
process shuffle {

    input:
    tuple val(sample_id), path(read_files)

    output:
    path file_name

    script:
    file_name = 'shuffled_reads_sample_'.concat(sample_id).concat('.').concat(read_files[0].getExtension())
    """
    get_seeded_random() { seed="\$1"; openssl enc -aes-256-ctr -pass pass:"\$seed" -nosalt < /dev/zero 2>/dev/null; };
    cat ${read_files} |  sed 'N;N;N;s/\\n/ /g'  | shuf --random-source=<(get_seeded_random ${params.seed}) | awk '{print \$1 "\\n" \$2 "\\n" \$3 "\\n" \$4}' > ${file_name}
    """
}

/*
* This process merges all given bam files specified in the pooled_gsa parameter.
* Takes:
*     A list with the paths to all bam files, that should be merged, if the condition is fullfilled.
* Output:
*     The path to the merged bam file.
 */
process shuffle_paired_end {

    input:
    tuple val(sample_id), path(first_read_files), path(second_read_files)

    output:
    tuple val(sample_id), path(file_name)

    script:
    file_name = 'shuffled_reads_sample_'.concat(sample_id).concat('.').concat(first_read_files[0].getExtension())
    """
    cat ${first_read_files} > first_reads.fq
    cat ${second_read_files} > second_reads.fq
    paste -d " " - - - - <first_reads.fq > first_reads_clustered.fq
    paste -d " " - - - - <second_reads.fq > second_reads_clustered.fq
    paste -d ' ' first_reads_clustered.fq second_reads_clustered.fq  > sample${sample_id}_interweaved.fq
    get_seeded_random() { seed="\$1"; openssl enc -aes-256-ctr -pass pass:"\$seed" -nosalt < /dev/zero 2>/dev/null; };
    shuf --random-source=<(get_seeded_random ${params.seed}) sample${sample_id}_interweaved.fq | tr " " "\n" > ${file_name}
    """
}
