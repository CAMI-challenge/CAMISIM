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
    take: genome_location_ch
    take: seed_file_read_simulation_ch

    main:

        // get the seed for every sample
        seed_ch = seed_file_read_simulation_ch.splitCsv(sep:'\t', skip:2)

        reads_seed_ch = reads_ch.join(seed_ch)

        if(params.type=="nanosim3") {
            out_shuffle = shuffle(reads_ch)
        } else if(params.type=="art" || params.type=="wgsim") {
            out_shuffle = shuffle_paired_end(reads_ch)
        }

        gsa_read_mapping(out_shuffle[1], params.genome_locations_file, params.metadata_file)
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

    conda "bioconda::biopython"

    input:
    tuple val(sample_id), path(first_read_files), path(second_read_files)

    output:
    tuple val(sample_id), path(anonymous_reads_file)
    tuple val(sample_id), path(tmp_reads_mapping_file)

    script:
    anonymous_reads_file = 'anonymous_reads.fq'
    tmp_reads_mapping_file = 'tmp_reads_mapping.tsv'
    """
    touch ${anonymous_reads_file}
    touch ${tmp_reads_mapping_file}
    cat ${first_read_files} > first_reads.fq
    cat ${second_read_files} > second_reads.fq
    paste -d " " - - - - <first_reads.fq > first_reads_clustered.fq
    paste -d " " - - - - <second_reads.fq > second_reads_clustered.fq
    paste -d ' ' first_reads_clustered.fq second_reads_clustered.fq  > sample${sample_id}_interweaved.fq
    get_seeded_random() { seed="\$1"; openssl enc -aes-256-ctr -pass pass:"\$seed" -nosalt < /dev/zero 2>/dev/null; };
    shuf --random-source=<(get_seeded_random ${params.seed}) sample${sample_id}_interweaved.fq | tr " " "\n" | tr -d '\\000' | python3 ${projectDir}/anonymizer.py -prefix S${sample_id}R -format fastq -map ${tmp_reads_mapping_file} -out ${anonymous_reads_file}
    mkdir --parents ${projectDir}/nextflow_out/sample_${sample_id}/reads
    gzip -k ${anonymous_reads_file}
    cp ${anonymous_reads_file}.gz ${projectDir}/nextflow_out/sample_${sample_id}/reads/
    """
}

process gsa_read_mapping {

    conda "bioconda::biopython"

    input:
    tuple val(sample_id), path(tmp_reads_mapping_file)
    path(genome_locations_file)
    path(metadata_file)


    output:
    tuple val(sample_id), path(reads_mapping_file)

    script:
    reads_mapping_file = 'reads_mapping.tsv'
    """
    touch ${reads_mapping_file}
    python ${projectDir}/scripts/goldstandardfileformat.py -input tmp_reads_mapping.tsv -genomes ${genome_locations_file} -metadata ${metadata_file} -out ${reads_mapping_file} -projectDir ${projectDir}
    mkdir --parents ${projectDir}/nextflow_out/sample_${sample_id}/reads
    gzip -k ${reads_mapping_file}
    cp ${reads_mapping_file}.gz ${projectDir}/nextflow_out/sample_${sample_id}/reads/
    """
}
