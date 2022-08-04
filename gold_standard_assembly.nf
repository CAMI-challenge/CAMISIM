/** 
* This workflow generates a gold standard assembly for every genome of a sample. 
* Takes:
*     A channel containing tuples with key = genome_id, first value = a sorted bam file, second value = the reference genome (fasta).
**/
workflow gold_standard_assembly {

    take: bam_files_channel
    main:
        // simulate reads via nanosim3
        gsa_for_every_genome_channel = generate_gold_standard_assembly(bam_files_channel)
}

/*
* This process generates a gold standard assembly for one genome.
* Takes:
*     A tuple with key = genome_id, first value = a sorted bam file, second value = the reference genome (fasta).
* Output:
*     A fasta file with the gold standard assembly of the given genome.
 */
process generate_gold_standard_assembly {

    conda 'bioconda::samtools'

    input:
    tuple val(genome_id), path(bam_file), path(reference_fasta_file)
    
    output:
    path('*.fasta')
    
    script:
    """
    perl -- ${projectDir}/scripts/bamToGold.pl -st samtools -r ${reference_fasta_file} -b ${bam_file} -l 1 -c 1 >> ${genome_id}_gsa.fasta
    """
}