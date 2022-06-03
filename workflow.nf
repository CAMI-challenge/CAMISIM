#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Queue channel -> channel in which data is consumed (used up) - dsl2??
//genome_fasta_ch = Channel.fromPath( "/home/jfunk/CAMISIM/Nextflow/CAMISIM/input/GCA_000092825.1_ASM9282v1.fa" )
genome_distribution_ch = Channel.fromPath( "/home/jfunk/CAMISIM/Nextflow/CAMISIM/input/distribution_0.txt" )
//genome_fasta_distribution_ch = Channel.of(["/home/jfunk/CAMISIM/Nextflow/CAMISIM/input/GCA_000092825.1_ASM9282v1.fa", 0.01842308570184231])

// muss ein value channel werden
genome_location_ch = Channel.from( "/home/jfunk/CAMISIM/Nextflow/CAMISIM/input/genome_locations.tsv" )


workflow {

    /// .map ({ line -> [line.split('\t')[0], line.split('\t')[1]] })
   
    genome_location_ch = Channel.fromPath( "/home/jfunk/CAMISIM/Nextflow/CAMISIM/input/genome_locations.tsv" ).splitCsv(sep:'\t')
    
    //.map ({ line -> [line.split('\t')[0], line.split('\t')[1]] })
    genome_distribution_split_ch = genome_distribution_ch.splitCsv(sep:'\t')
    
    genome_location_distribution_ch = genome_location_ch.join(genome_distribution_split_ch)
    
    //genome_location_distribution_ch_1 = genome_location_ch.mix(genome_distribution_split_ch)
    
    //genome_location_distribution_ch = genome_location_distribution_ch_1.groupTuple()

    //ToDo: Konstanten in einer config Datei festhalten
    
    genome_location_distribution_ch_test = Channel.of(["Genome6.0", "/home/jfunk/CAMISIM/Nextflow/CAMISIM/input/GCA_000006785.2_ASM678v2.fa", 0.01842308570184231])
    
    simulated_reads_ch = simulate_reads_nanosim3(genome_location_distribution_ch, 7408, 651524512, 0.05)
}


process simulate_reads_nanosim3 {
	
    input:
    // Der Pfad zur Fasta File und die abundance muessen als tuple gegeben sein, es koennen nicht verschiedene channel sein, da sie sonst unterschiedlich kombiniert werden wuerden.
    tuple val(genome_id), path(fasta_file), val(abundance)
    val fragment_size_mean
    val seed
    val total_size
    
    output:
    0
    
    script:
    print genome_id
    print fasta_file
    print abundance
    float x = Float.valueOf(abundance)
    number_of_reads = x * (total_size*1000000000) / fragment_size_mean
    number_of_reads = number_of_reads.round(0)
    """
    /home/jfunk/CAMISIM/NanoSim/NanoSim/src/simulator.py genome -n ${number_of_reads} -rg ${fasta_file} -o reads -c /home/jfunk/CAMISIM/code/CAMISIM_2/CAMISIM/read_simulators/nanosim_profile/training --seed ${seed} -dna_type linear
    """
}

// Doku Nextflow???
