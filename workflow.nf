#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Defining the module / subworkflow path, and include the elements
 */

// include sample wise simulation
include { sample_wise_simulation } from "${projectDir}/sample_wise_simulation"

// this channel holds the files with the specified distributions for every sample
genome_distribution_file_ch = Channel.fromPath( "./nextflow_defaults/distribution_*.txt" )

// this channel holds the file with the specified locations of the genomes
genome_location_file_ch = Channel.fromPath( "./nextflow_defaults/genome_locations.tsv" )

// this channel holds the ncbi tax dump
ncbi_taxdump_file_ch = Channel.fromPath( "./tools/ncbi-taxonomy_20170222.tar.gz" )

/*
 * This is the main workflow and starting point of this nextflow pipeline.
 */
workflow {

    // build ncbi taxonomy from given tax dump
    number_of_samples = genome_distribution_file_ch.count()
    buildTaxonomy(ncbi_taxdump_file_ch.combine(number_of_samples))
    
    // simulate reads sample wise
    sample_wise_simulation(genome_location_file_ch, genome_distribution_file_ch)
    // this workflow has two output channels: one bam file per sample and one fasta file per sample
    merged_bam_per_sample = sample_wise_simulation.out[0].collect()
    gsa_for_all_reads_of_one_sample_ch = sample_wise_simulation.out[1]   

    // merge the bam files of the required samples
    merged_bam_file = merge_bam_files(merged_bam_per_sample)

    reference_fasta_files_ch = genome_location_file_ch.splitCsv(sep:'\t').map { a -> a[1] }

    generate_pooled_gold_standard_assembly(merged_bam_file.combine(reference_fasta_files_ch).groupTuple())

}

/*
* This process merges all given bam files specified in the pooled_gsa parameter.
* Takes:
*     A list with the paths to all bam files, that should be merged, if the condition is fullfilled.
* Output:
*     The path to the merged bam file.
 */
process merge_bam_files {

    conda 'bioconda::samtools'

    input:
    path bam_files

    output:
    path file_name

    script:
    file_name = 'merged.bam'
    compression = 5
    memory = 1

    bam_to_merge = ''

    bam_files.each { 

        bam_file_name = (String) it
        sample_id = bam_file_name.split('_')[1][0].toInteger()
        
        if(sample_id in params.pooled_gsa){
            bam_to_merge = bam_to_merge.concat(' ').concat(bam_file_name)
        }
    }
    """
    samtools merge -u - ${bam_to_merge} | samtools sort -l ${compression} -m ${memory}G -o ${file_name} -O bam
    samtools index ${file_name}
    """
}

/*
* This process generates the pooled gold standard assembly for serveral samples.
* Takes:
*     A tuple with first_value = a sorted bam file and second value = the reference genome (fasta).
* Output:
*     The path to fasta file with the pooled gold standard assembly.
 */
process generate_pooled_gold_standard_assembly {

    conda 'bioconda::samtools'

    input:
    tuple path(bam_file), path(reference_fasta_files)

    output:
    path file_name
    
    script:
    file_name = 'gsa_pooled.fasta'
    """
    cat ${reference_fasta_files} > reference.fasta
    perl -- ${projectDir}/scripts/bamToGold.pl -st samtools -r reference.fasta -b ${bam_file} -l 1 -c 1 >> ${file_name}
    """
}

/*
* This process builds the taxonomy profile for every sample with the given distribution and the ncbi tax dump. The generated profiles will be
* copied to the out directory.
* Takes:
*     The tuple with first_value = zipped ncbi tax dump to build the profile from  and second value = the number of samples.
* Output:
*     The paths to all taxonomic profiles.
 */
process buildTaxonomy {

    input:
    tuple path(dmp), val(number_of_samples)

    output:
    path 'taxonomic_profile_*.txt'

    script:
    index_number_of_samples = number_of_samples - 1
    """
    tar -xf ${dmp}
    ${projectDir}/build_ncbi_taxonomy.py **/names.dmp **/merged.dmp **/nodes.dmp ${number_of_samples} ${projectDir}/nextflow_defaults/distribution_{0..${index_number_of_samples}}.txt
    mkdir --parents ${projectDir}/nextflow_out/
    cp taxonomic_profile_*.txt ${projectDir}/nextflow_out/
    """
}