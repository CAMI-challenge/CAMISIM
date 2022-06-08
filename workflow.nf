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
        read_simulator_nansoim3(genome_location_distribution_ch)
    }
}
