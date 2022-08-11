/** 
* This workflow generates a gold standard assembly for every genome of a sample. 
* Takes:
*     A channel containing tuples with key = genome_id, first value = a sorted bam file, second value = the reference genome (fasta).
**/
workflow gold_standard_assembly {

    take: bam_files_channel
    main:
        // generate gold standard assembly for all genomes of one sample
        // the buffer size needs to be set to the total number of genomes, so that nextflow waits with the execution of the next process until 
        // all genomes of a sample are processed
        gsa_for_every_genome_channel = generate_gold_standard_assembly(bam_files_channel).buffer(size:params.genomes_total)

        // create a fasta file holding all gold standard assemblies of all reads of one sample
        gsa_for_all_reads_of_one_sample_ch = get_fasta_for_sample(gsa_for_every_genome_channel)
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

/*
* This process creates a fasta file holding the content of all given fasta files. The order will be determined by the order they are processed by
* nextflow.
* Takes:
*     The paths to all fasta files, that need to be combined.
* Output:
*     One fasta file holding the content of all given fasta files.
 */
process get_fasta_for_sample {

    input:
    path fasta_files

    output:
    path 'gsa.fasta'
        
    script:
    """
    cat ${fasta_files} > gsa.fasta
    """
}

