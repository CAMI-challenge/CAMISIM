/** 
* This workflow simulates reads via nanosim3 and converts the resulting sam files into bam files.
* Takes:
*     A channel containing tuples with key = genome_id, first value = path to genome, second value = distribution, third value = sample_id, fourth value = seed.
* Emits: 
*     A channel containing tuples with key = sample_id, first value = genome id, second value = simulated bam file, third value = the reference fasta file.
**/
workflow read_simulator_nanosim3 {

    take: genome_location_distribution_ch
    take: read_length_ch
    main:
        // simulate reads via nanosim3
        simulate_reads_nanosim3(genome_location_distribution_ch, read_length_ch)
        simulated_reads_ch = simulate_reads_nanosim3.out[0]

        sam_file_channel = sam_from_reads(simulated_reads_ch)
        bam_file_channel = sam_to_bam(sam_file_channel)
    emit:
        bam_file_channel
        simulate_reads_nanosim3.out[1].groupTuple()
}

/**
* This process simulates reads with nanosim3.
* Input:
*     Tuple containing key = genome_id, first value = path to genome, second value = distribution, third value = sample_id, fourth value = seed.
* Output:
*     Tuple containing key = sample_id, first value = genome_id, second value = path to error profile, third value = path to fasta file with the aligned reads, 
*         fourth value = path to fasta file with the aligned reads, fifth value = path to reference genome.
**/
process simulate_reads_nanosim3 {

    conda 'anaconda::scikit-learn=0.21.3=py37hd81dba3_0 bioconda::nanosim=3.0'
	
    input:
    tuple val(genome_id), val(sample_id), path(fasta_file), val(abundance), val (seed)
    val(read_length_ch)
    
    output:
    tuple val(sample_id), val(genome_id), path('*_error_profile'), path("*_aligned_reads.fasta"), path("*_unaligned_reads.fasta"), path(fasta_file)
    tuple val(sample_id), path("*_aligned_reads.fasta")
    
    script:
    total_size = params.size
    profile = params.base_profile_name
    number_of_reads = (total_size*(10**9)) * abundance.toFloat() / read_length_ch.toFloat()
    number_of_reads = number_of_reads.round(0).toInteger()
    // nanosim seed cannot be > 2**32 -1
    Long used_seed = (seed as Long) % 2**32 - 1

    /**
    String log = "---- sample id: ".concat(sample_id)
    log = log.concat("  genome id: ").concat(genome_id)
    log = log.concat("   fasta file: ").concat(fasta_file.baseName)
    log = log.concat("  abundance: ").concat(abundance)
    log = log.concat("    seed: ").concat(seed)
    log = log.concat("    used_seed: ").concat(Long.toString(used_seed))
    log = log.concat("    number_of_reads: ").concat(Integer.toString(number_of_reads))
    log = log.concat("    profile: ").concat(profile)
    print(log)
    **/

    """
    simulator.py genome -n ${number_of_reads} -rg ${fasta_file} -o sample${sample_id}_${genome_id} -c ${profile} --seed ${used_seed} -dna_type linear

    mkdir --parents ${projectDir}/nextflow_out/sample_${sample_id}/reads/fastq/
    cp *_aligned_reads.fasta ${projectDir}/nextflow_out/sample_${sample_id}/reads/fastq/
    """
}

/**
* This process converts a the simulated reads into a sam file.
* Input:
*     Tuple containing key = sample_id, first value = genome_id, second value = path to error profile, third value = path to fasta file with the aligned reads, 
*         fourth value = path to fasta file with the aligned reads, fifth value = path to reference genome.
* Output:
*    Tuple containing key = sample_id, first value = genome_id, second value = path to sam file, third value = path to reference genome.
**/
process sam_from_reads {

    conda "bioconda::biopython"

    input:
    tuple val(sample_id), val(genome_id), val(error_profile), path(aligned_reads), path(unaligned_reads), path(fasta_file)

    output:
    tuple val(sample_id), val(genome_id), path('*.sam'), path(fasta_file)

    script:
    """
    ${projectDir}/read_simulators/sam_from_reads.py ${error_profile} ${aligned_reads} ${unaligned_reads} ${fasta_file}
    """
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
