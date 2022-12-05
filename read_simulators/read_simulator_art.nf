/**
* This workflow simulates reads using the ART read simuator
* Takes:
*     A channel containing tuples with key = sample_id, first value = genome_id, second value = path to genome, third value = distribution.
* Emits: 
*     A channel containing tuples with key = sample_id, first value = genome id, second value = simulated bam file, third value = the reference fasta file.
**/

workflow read_simulator_art {

    take: genome_location_distribution_ch
    take: read_length_ch
    main:
        sam_file_channel = simulate_reads_art(genome_location_distribution_ch, read_length_ch)
        bam_file_channel = sam_to_bam(sam_file_channel)
    emit:
        bam_file_channel
}

/**
* This process simulates short reads with ART
* Input:
*     Tuple containing key = sample_id, first value = genome_id, second value = path to genome, third value = distribution.
* Output:
*     Tuple containing key = sample_id, first value = genome_id, second value = path to error profile, third value = path to fasta file with the aligned reads, 
*         fourth value = path to fasta file with the aligned reads, fifth value = path to reference genome.
**/

process simulate_reads_art {
    
    conda 'bioconda::art=2016.06.05' // TODO: check version and dependencies (gsl, libcblas, libgcc-ng, libstdcxx-ng)
    
    input:
    tuple val(genome_id), path(fasta_file), val(abundance), val(sample_id)
    val(read_length_ch)
    
    output:
    tuple val(sample_id), val(genome_id), path('*.sam'), path(fasta_file) 
   
    script:
    seed = params.seed
    fragment_size_mean = params.fragment_size_mean
    fragment_size_sd = params.fragment_size_sd
    profile = params.base_profile_name
    fold_coverage = abundance // TODO is the abundance already normalised?
    """
    art_illumina -sam -na -i ${fasta_file} -l ${read_length} -m ${fragment_size_mean} -s ${fragment_size_sd} -f ${fold_coverage} -o sample${sample_id}_${genome_id} -1 ${profile}1.txt -2 ${profile}2.txt -rs ${seed}
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
    echo ${sam_file} >> debug
    samtools view -bS ${sam_file} -o alignment_to_sort.bam
    samtools sort -o sample${sample_id}_${genome_id}.bam alignment_to_sort.bam
    """

}
