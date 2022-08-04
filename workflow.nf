#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Defining the module / subworkflow path, and include the elements
 */
read_simulator_folder = "./read_simulators/"
include { read_simulator_nansoim3 } from "${read_simulator_folder}/read_simulator_nansoim3"

// this channel holds the files with the specified distributions of sample
genome_distribution_ch = Channel.fromPath( "./nextflow_defaults/distribution_0.txt" )

// this channel holds the file with the specified locations of the genomes
genome_location_ch = Channel.from( "./nextflow_defaults/genome_locations.tsv" )

/*
 * This is the main workflow and starting point of this nextflow pipeline.
 */
workflow {
   
    // this channel holds the genome location map (key = genome_id, value = absolute path to genome)
    genome_location_ch = Channel.fromPath( "./nextflow_defaults/genome_locations.tsv" ).splitCsv(sep:'\t')
    
    // this channel holds the genome distribution map (key = genome_id, value = distribution)
    genome_distribution_split_ch = genome_distribution_ch.splitCsv(sep:'\t')
    
    // joining of the channels results in new map: key = genome_id, first value = absolute path to genome, second value = distribution
    genome_location_distribution_ch = genome_location_ch.join(genome_distribution_split_ch)
    
    // simulate the reads according to the given simulator type
    if(params.type.equals("nanosim3")) {
        // simulate the reads with nanosim3
        sam_files_channel = read_simulator_nansoim3(genome_location_distribution_ch)
    }

    // convert sam files to sorted bam files
    bam_files_channel = sam_to_bam(sam_files_channel)

    generate_gold_standart_assembly(bam_files_channel)

    generate_gold_standart_assembly.out.view()
}

/* 
* This process converts a sam to a bam file.
* Takes:
*     A tuple with key = genome_id, first value = the sam file to convert, second value = the reference genome (fasta).
* Output:
*     A tuple with key = genome_id, first value = a sorted bam file, second value = the reference genome (fasta).
 */
process sam_to_bam {

    conda 'bioconda::samtools'

    input:
    tuple val(genome_id), path(sam_file), path(fasta_file)

    output:
    tuple val(genome_id), path('*.bam'), path(fasta_file)

    script:
    """
    samtools view -bS ${sam_file} -o alignment_to_sort.bam
    samtools sort -o ${genome_id}.bam alignment_to_sort.bam
    """

}

/*
* This process generates a gold standard assembly for one genome.
* Takes:
*     A tuple with key = genome_id, first value = a sorted bam file, second value = the reference genome (fasta).
* Output:
*     A fasta file with the gold standard assembly of the given genome.
 */
process generate_gold_standart_assembly {

    conda 'bioconda::samtools'

    input:
    tuple val(genome_id), path(bam_file), path(reference_fasta_file)
    
    output:
    path('*.fasta')
    
    script:
    """
    perl -- ${projectDir}/scripts/bamToGold.pl -st samtools -r ${reference_fasta_file} -b ${bam_file} -l 1 -c 1 >> ${genome_id}_gsa.fasta
    """
}
