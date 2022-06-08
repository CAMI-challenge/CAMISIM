#!/usr/bin/env nextflow

nextflow.enable.dsl=2

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
    
    // simulate reads via nanosim3
    simulated_reads_ch = simulate_reads_nanosim3(genome_location_distribution_ch, 7408, 651524512, 0.05)
}


process simulate_reads_nanosim3 {
	
    input:
    // map: key = genome_id, first value = absolute path to genome, second value = distribution
    tuple val(genome_id), path(fasta_file), val(abundance)
    // this value has been calculated using the script tools/nanosim_profile/get_mean from the values in nanosim_profile and is used for getting the correct number of reads
    val fragment_size_mean
    // the specified seed
    val seed
    // the specified total size in Gbp
    val total_size
    
    output:
    0
    
    script:
    
    number_of_reads = (total_size*1000000000) * abundance.toFloat() / fragment_size_mean
    number_of_reads = number_of_reads.round(0).toInteger()
    """
    /home/jfunk/CAMISIM/NanoSim/NanoSim/src/simulator.py genome -n ${number_of_reads} -rg ${fasta_file} -o reads -c /home/jfunk/CAMISIM/code/CAMISIM_2/CAMISIM/read_simulators/nanosim_profile/training --seed ${seed} -dna_type linear
    """
}
