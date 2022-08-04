/** 
* This workflow simulates reads via nanosim3 and converts the resulting sam files into bam files.
* Takes:
*     genome_location_distribution_ch: A channel containing tuples with key = genome_id, first value = absolute path to genome, second value = distribution.
* Emits: bam_file_channel: A channel containing tuples containing the genome id, the simulated bam file and the reference fasta file.
**/
workflow read_simulator_nansoim3 {

    take: genome_location_distribution_ch
    main:
        // simulate reads via nanosim3
        simulated_reads_ch = simulate_reads_nanosim3(genome_location_distribution_ch)

        bam_file_channel = sam_from_reads(simulated_reads_ch)
    emit:
        bam_file_channel
}

/**
* This process simulates reads with nanosim3.
* Input:
*     Tuple containing key = genome_id, first value = path to genome, second value = distribution.
* Output:
*     Tuple containing key = genome_id, first value = path to error profile, second value = path to fasta file with the aligned reads, 
*         third value = path to fasta file with the aligned reads, fourth value = path to reference genome
**/
process simulate_reads_nanosim3 {

    conda 'anaconda::scikit-learn=0.21.3=py37hd81dba3_0 bioconda::nanosim=3.0'
	
    input:
    tuple val(genome_id), path(fasta_file), val(abundance)
    
    output:
    tuple val(genome_id), path('*_error_profile'), path("*_aligned_reads.fasta"), path("*_unaligned_reads.fasta"), path(fasta_file)
    
    script:
    seed = params.seed
    total_size = params.size
    
    number_of_reads = (total_size*1000000000) * abundance.toFloat() / params.fragment_size_mean_nanosim
    number_of_reads = number_of_reads.round(0).toInteger()
    """
    simulator.py genome -n ${number_of_reads} -rg ${fasta_file} -o ${genome_id} -c ${projectDir}/tools/nanosim_profile/training --seed ${seed} -dna_type linear
    """
}

/**
* This process converts a the simulated reads into a sam file.
* Input:
*     Tuple containing key = genome_id, first value = path to error profile, second value = path to fasta file with the aligned reads, 
*         third value = path to fasta file with the aligned reads, fourth value = path to reference genome
* Output:
*    Tuple containing key = genome_id, first value = path to sam file, second value = path to reference genome
**/
process sam_from_reads {

    input: 
    // the error profile files that need to be converted into a sam format
    tuple val(genome_id), val(error_profile), path(aligned_reads), path(unaligned_reads), path(fasta_file)

    output:
    tuple val(genome_id), path('*.sam'), path(fasta_file)

    script:
    """
    ${projectDir}/read_simulators/sam_from_reads.py ${error_profile} ${aligned_reads} ${unaligned_reads} ${fasta_file}
    """
}
