/**
* This workflow simulates reads using the wgsim read simuator
* Takes:
*     A channel containing tuples with key = genome_id, first value = path to genome, second value = distribution, third value = sample_id, fourth value = factor.
*     A channel containing the read length.
* Emits: 
*     A channel containing tuples with key = sample_id, first value = genome id, second value = simulated bam file, third value = the reference fasta file.
**/

workflow read_simulator_wgsim {

    take: genome_location_distribution_ch
    take: read_length_ch
    main:
        simulate_reads_wgsim(genome_location_distribution_ch, read_length_ch)
        sam_file_channel = simulate_reads_wgsim.out[0]

        bam_file_channel = sam_to_bam(sam_file_channel)
    emit:
        bam_file_channel
        simulate_reads_wgsim.out[1].groupTuple()
}

/**
* This process simulates short reads with wgsim
* Input:
*     Tuple containing key = genome_id, first value = path to genome, second value = distribution, third value = sample_id, fourth value = factor.
* Output:
*     Tuple containing key = sample_id, first value = genome_id, second value = path to error profile, third value = path to fasta file with the aligned reads, 
*         fourth value = path to fasta file with the aligned reads, fifth value = path to reference genome.
**/

process simulate_reads_wgsim {

    conda 'bioconda::wgsim'

    input:
    tuple val(genome_id), val(sample_id), path(fasta_file), val(abundance), val (seed)
    val(read_length_ch)
    
    output:
    tuple val(sample_id), val(genome_id), path('*.sam'), path(fasta_file)
    tuple val(sample_id), path('*.01.fq'), path('*.02.fq')
    
    script:
    total_size = params.size
    fragment_size = params.fragment_size_mean
    fragment_size_sd = params.fragment_size_sd
    error_rate = params.base_error_rate
    read_length = read_length_ch
    number_of_reads = (total_size*(10**9)) * abundance.toFloat() / read_length_ch.toFloat()
    number_of_reads = number_of_reads.round(0).toInteger()
    """
    wgsim -d ${fragment_size} -s ${fragment_size_sd} -N ${number_of_reads} -1 ${read_length} -2 ${read_length} -S ${seed} -e ${error_rate} -r 0 -R 0 ${fasta_file} sample${sample_id}_${genome_id}.01.fq sample${sample_id}_${genome_id}.02.fq 
    mkdir --parents ${projectDir}/nextflow_out/sample_${sample_id}/reads/
    mkdir --parents ${projectDir}/nextflow_out/sample_${sample_id}/reads/fastq/
    cp *.sam ${projectDir}/nextflow_out/sample_${sample_id}/reads/
    cp sample${sample_id}_${genome_id}*.fq ${projectDir}/nextflow_out/sample_${sample_id}/reads/fastq/
    """
    /**
    @TODO: Maybe add the option to add ALL options of wgsim in the profile
    **/
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
* TODO: Copy from nanosim simulator - make externally evailable
*/
