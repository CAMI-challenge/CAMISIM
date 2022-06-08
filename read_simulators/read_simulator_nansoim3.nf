workflow read_simulator_nansoim3 {

    take: genome_location_distribution_ch
    main:
        // simulate reads via nanosim3
        simulated_reads_ch = simulate_reads_nanosim3(genome_location_distribution_ch)
}

process simulate_reads_nanosim3 {
	
    input:
    // map: key = genome_id, first value = absolute path to genome, second value = distribution
    tuple val(genome_id), path(fasta_file), val(abundance)
    
    output:
    0
    
    script:
    seed = params.seed
    total_size = params.size
    
    number_of_reads = (total_size*1000000000) * abundance.toFloat() / params.fragment_size_mean_nanosim
    number_of_reads = number_of_reads.round(0).toInteger()
    """
    /home/jfunk/CAMISIM/NanoSim/NanoSim/src/simulator.py genome -n ${number_of_reads} -rg ${fasta_file} -o reads -c /home/jfunk/CAMISIM/code/CAMISIM_2/CAMISIM/read_simulators/nanosim_profile/training --seed ${seed} -dna_type linear
    """
}
