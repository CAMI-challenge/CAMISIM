/** 
* This workflow generates a gold standard assemblies.
* Takes:
*     A channel containing tuples with key = sample_id, first value = genome_id, second value = a sorted bam file, third value = the reference genome (fasta).
**/
workflow gold_standard_assembly {

    take: bam_files_channel
    main:
        // generate gold standard assembly for every genome
        gsa_for_every_genome_ch = generate_gold_standard_assembly(bam_files_channel)
        
        // grouping the gold standard assemblies by sample id results in new tuple: key = sample_id, values = path to all gsa of this samples reads
        grouped_gsa_for_every_genome_ch = gsa_for_every_genome_ch.groupTuple()

        // create fasta files holding all gsa of one samples reads
        gsa_for_all_reads_of_one_sample_ch = get_fasta_for_sample(grouped_gsa_for_every_genome_ch)
}

/*
* This process generates a gold standard assembly for one genome.
* Takes:
*     A tuple with key = sample_id, first_value = genome_id, second value = a sorted bam file, third value = the reference genome (fasta).
* Output:
*     A uple with key = sample_id, value = path to fasta file with the gold standard assembly of the given genome.
 */
process generate_gold_standard_assembly {

    conda 'bioconda::samtools'

    input:
    tuple val(sample_id),val(genome_id), path(bam_file), path(reference_fasta_file)
    
    output:
    tuple val(sample_id), path(file_name)
    
    script:
    file_name = 'sample'.concat(sample_id.toString()).concat('_').concat(genome_id).concat('_gsa.fasta')
    """
    perl -- ${projectDir}/scripts/bamToGold.pl -st samtools -r ${reference_fasta_file} -b ${bam_file} -l 1 -c 1 >> ${file_name}
    """
}

/*
* This process creates a fasta file holding the content of all given fasta file. The order will be determined by the order they are processed by
* nextflow.
* Takes:
*     A tuple with key = sample_id, value = the paths to all fasta files, that need to be combined.
* Output:
*     One fasta file holding the content of all given fasta files.
 */
process get_fasta_for_sample {

    input:
    tuple val(sample_id), path(fasta_files)

    output:
    path file_name
        
    script:
    file_name = 'sample'.concat(sample_id.toString()).concat('_gsa.fasta')
    """
    cat ${fasta_files} > ${file_name}
    """
}

