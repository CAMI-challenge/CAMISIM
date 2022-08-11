#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Defining the module / subworkflow path, and include the elements
 */

 // include read simulator here:
read_simulator_folder = "./read_simulators/"
// include read simulator nanosim3
include { read_simulator_nansoim3 } from "${read_simulator_folder}/read_simulator_nansoim3"

// include workflow for generating gold standard assemblies
include { gold_standard_assembly } from "${projectDir}/gold_standard_assembly"


// this channel holds the files with the specified distributions for every sample
genome_distribution_ch = Channel.fromPath( "./nextflow_defaults/distribution_*.txt" )

// this channel holds the file with the specified locations of the genomes
genome_location_ch = Channel.from( "./nextflow_defaults/genome_locations.tsv" )

/*
 * This is the main workflow and starting point of this nextflow pipeline.
 */
workflow {
   
    // this channel holds the genome location map (key = genome_id, value = absolute path to genome)
    genome_location_ch = Channel.fromPath( "./nextflow_defaults/genome_locations.tsv" ).splitCsv(sep:'\t')

    // create a channel, that holds the path to the genome distribution file by the sample id (key = sample id, first value = genome_id, second value = distribution)
    genome_distribution_by_sample_ch = get_distribution_file_by_sample(genome_distribution_ch).splitCsv(sep:'\t').flatten().collate(3)
    
    // crossing and filtering of the channels results in new map: key = sample_id, first value = genome_id, second value = path to genome, third value = distribution
    genome_location_distribution_ch = unique_genome_id(genome_location_ch.cross(genome_distribution_by_sample_ch).flatten().collate(5))
    
    // simulate the reads according to the given simulator type
    if(params.type.equals("nanosim3")) {
        // simulate the reads with nanosim3
        sam_files_channel = read_simulator_nansoim3(genome_location_distribution_ch)
    }

    // convert sam files to sorted bam files
    bam_files_channel = sam_to_bam(sam_files_channel)

    // generate a gold standard assemblies
    gold_standard_assembly(bam_files_channel)
}

/* 
* This process converts a sam to a bam file.
* Takes:
*     A tuple with key = sample_id, first value = genome id, second value = sam file to convert, third value = the reference fasta file.
* Output:
*     A tuple with key = sample_id, first value = genome_id, second value = a sorted bam file, third value = the reference genome (fasta).
 */
process sam_to_bam {

    conda 'bioconda::samtools'

    input:
    tuple val(sample_id), val(genome_id), path(sam_file), path(fasta_file)

    output:
    tuple val(sample_id), val(genome_id), path('sample*.bam'), path(fasta_file)

    script:
    """
    samtools view -bS ${sam_file} -o alignment_to_sort.bam
    samtools sort -o sample${sample_id}_${genome_id}.bam alignment_to_sort.bam
    """

}

/* 
* This process takes a path to a distribution file and returns a tuple holding the sample Id and the path to the given file.
* Takes:
*     Path to a distribution file.
* Output:
*     A tuple with key = sample_id, value = path to the given distribution file.
 */
process get_distribution_file_by_sample {

    input:
    path distribution_file

    output:
    tuple path(distribution_file), val(sample_id)

    script:
    file_name = distribution_file.toString()
    sample_id = file_name.split('_')[1].split('.txt')[0].toInteger()
    """
    """
}

/* 
* This process takes a tuple and deletes the second entry of the genome_id in it.
* Takes:
*     A tuple with key = genome_id, first value = path to genome, second value = genome_id, third value = distribution, fourth value = sample_id.
* Output:
*     A tuple with key = sample_id, first value = genome_id, second value = path to genome, third value = distribution.
 */
process unique_genome_id {

    input:
    tuple val(genome_id), path(genome_location), val(genome_id), val(distribution), val(sample_id)

    output:
    tuple val(sample_id), val(genome_id), path(genome_location), val(distribution)

    script:
    """
    """
}
