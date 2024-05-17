#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/** 
* This workflow does the binning of the gold standard assembly.
**/
workflow binning {

    take: samplewise_gsa_ch // tuple val(sample_id), path(gsa)
    take: bam_file_list_per_sample_ch
    take: pooled_gsa_ch
    take: merged_bam_ch
    take: genome_location_file_ch
    take: metadata_ch

    main:

        read_start_positions_from_dir_of_bam(bam_file_list_per_sample_ch)
        binning_per_sample(samplewise_gsa_ch.join(read_start_positions_from_dir_of_bam.out), genome_location_file_ch, metadata_ch)

        read_start_positions_from_merged_bam(merged_bam_ch)
        binning_pooled_gsa(pooled_gsa_ch, read_start_positions_from_merged_bam.out, genome_location_file_ch, metadata_ch)
}

/*
* This process parses 'read' start positions from bam files in a directory.
* Takes:
*   The list of bam files per sample id.
* Output:
*    A file containing the read start posotions for the given sample.
 */
process read_start_positions_from_dir_of_bam {
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
*
* This process parses 'read' start positions from a merged bam file.
* Takes:
*   The bam file.
* Output:
*    A file containing the read start posotions.
 */
process read_start_positions_from_merged_bam {
    conda 'bioconda::samtools'
    
    input:
    path(merged_bam_file)

    output:
    path(filename)

    script:
    filename = "read_start_positions"
    """
    set -o pipefail
    samtools view ${merged_bam_file} | awk '{print \$1 "\\t" \$4}' >> ${filename}
    """
}

/*
* This process does the binning for per sample.
* Takes:
*   The gsa for the given sample and the read positions per sample.
*   A file containing all reference genome locations.
*   The metadata file.
 */
process binning_per_sample {
    conda "bioconda::biopython"

    publishDir "${params.outdir}/sample_${sample_id}/contigs/", mode : 'copy'

    input:
    tuple val(sample_id), path(gsa), path(read_positions)
    path(genome_locations_file)
    path(metadata_file)

    script:
    gsa_mapping_file = 'gsa_mapping.tsv'
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
    goldstandardfileformat.py -binning -read_positions ${read_positions} -genomes ${genome_locations_file} -metadata ${metadata_file} -out ${gsa_mapping_file} -projectDir ${projectDir} -input ${gsa} ${real_fastq} ${wgsim}
    gzip -k ${gsa_mapping_file}
    """
}

/*
* This process does the binning for the pooled gsa.
* Takes:
*   The pooled gsa for the given sample and the read positions of the merged bam file.
*   A file containing all reference genome locations.
*   The metadata file.
 */
process binning_pooled_gsa {
    conda "bioconda::biopython"

    publishDir "${params.outdir}/pooled_gsa/", mode : 'copy'

    input:
    path(gsa)
    path(read_positions)
    path(genome_locations_file)
    path(metadata_file)

    script:
    gsa_mapping_file = 'pooled_gsa_mapping.tsv'
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
    goldstandardfileformat.py -binning -read_positions ${read_positions} -genomes ${genome_locations_file} -metadata ${metadata_file} -out ${gsa_mapping_file} -projectDir ${projectDir} -input ${gsa} ${real_fastq} ${wgsim}
    gzip -k ${gsa_mapping_file}
    """
}
