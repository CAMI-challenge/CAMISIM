#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Defining the module / subworkflow path, and include the elements
 */

// include sample wise simulation
include { sample_wise_simulation } from "${projectDir}/sample_wise_simulation"

// this channel holds the files with the specified distributions for every sample
genome_distribution_file_ch = Channel.fromPath( "./nextflow_defaults/distribution_*.txt" )

// this channel holds the file with the specified locations of the genomes
genome_location_file_ch = Channel.fromPath( "./nextflow_defaults/genome_locations.tsv" )

/*
 * This is the main workflow and starting point of this nextflow pipeline.
 */
workflow {

    // simulate reads sample wise
    sample_wise_simulation(genome_location_file_ch, genome_distribution_file_ch)
    // this workflow has two output channels: one bam file per sample and one fasta file per sample
    merged_bam_per_sample = sample_wise_simulation.out[0]
    gsa_for_all_reads_of_one_sample_ch = sample_wise_simulation.out[1]
}
