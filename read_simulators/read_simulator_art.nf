/**
* This workflow simulates reads using the ART read simuator
* Takes:
*     A channel containing tuples with key = genome_id, first value = path to genome, second value = distribution, third value = sample_id, fourth value = factor.
*     A channel containing the read length.
* Emits: 
*     A channel containing tuples with key = sample_id, first value = genome id, second value = simulated bam file, third value = the reference fasta file.
**/

workflow read_simulator_art {

    take: genome_location_distribution_ch
    take: read_length_ch
    main:
        sam_ch = simulate_reads_art(genome_location_distribution_ch, read_length_ch)
        create_bam(sam_ch)
    emit:
        create_bam.out[0]
        simulate_reads_art.out[1].groupTuple()
}

/**
* This process simulates short reads with ART
* Input:
*     Tuple containing key = genome_id, first value = path to genome, second value = distribution, third value = sample_id, fourth value = factor.
* Output:
*     Tuple containing key = sample_id, first value = genome_id, second value = path to error profile, third value = path to fasta file with the aligned reads, 
*         fourth value = path to fasta file with the aligned reads, fifth value = path to reference genome.
**/
process simulate_reads_art {

    scratch true
    container 'quay.io/biocontainers/art:2016.06.05--h589041f_9'
    conda 'bioconda::art=2016.06.05 conda-forge::gsl=2.7 bioconda::samtools' // TODO: check version and dependencies (gsl, libcblas, libgcc-ng, libstdcxx-ng)
    
    input:
    tuple val(genome_id), val(sample_id), path(fasta_file), val(abundance), val(seed), val(factor), path(profile1), path(profile2)
    val(read_length_ch)
    
    output:
    tuple val(sample_id), val(genome_id), path("sample${sample_id}_${genome_id}.sam"), path(fasta_file)
    tuple val(sample_id), path('*1.fq'), path('*2.fq')
   
    script:
    fragment_size_mean = params.fragment_size_mean
    fragment_size_sd = params.fragment_size_sd
    profile = params.base_profile_name
    factor_float_value = Double.valueOf(factor)
    fold_coverage = Double.valueOf(abundance) * factor_float_value // TODO is the abundance already normalised?

    /**
    String log = "---- sample id: ".concat(sample_id)
    log = log.concat("  genome id: ").concat(genome_id)
    log = log.concat("   fasta file: ").concat(fasta_file.baseName)
    log = log.concat("  fragment_size_mean: ").concat(Integer.toString(fragment_size_mean))
    log = log.concat("    fragment_size_sd: ").concat(Integer.toString(fragment_size_sd))
    log = log.concat("    profile: ").concat(profile)
    log = log.concat("    factor_float_value: ").concat(Double.toString(factor_float_value))
    log = log.concat("    fold_coverage: ").concat(Double.toString(fold_coverage))
    log = log.concat("    seed: ").concat(Double.toString(seed))
    print(log)
    **/

    """
    art_illumina -sam -na -i ${fasta_file} -l ${read_length_ch} -m ${fragment_size_mean} -s ${fragment_size_sd} -f ${fold_coverage} -o sample${sample_id}_${genome_id} -1 ${profile2} -2 ${profile2} -rs ${seed}
    for file in sample${sample_id}_${genome_id}*.fq; do gzip -k "\$file"; done
    mkdir --parents ${projectDir}/nextflow_out/sample_${sample_id}/reads/fastq/
    cp sample${sample_id}_${genome_id}*.fq.gz ${projectDir}/nextflow_out/sample_${sample_id}/reads/fastq/
    """
}

process create_bam {
    container 'quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1'
    input:
    tuple val(sample_id), val(genome_id), path(sam_path), path(fasta_file)
    tuple val(sample_id), path(read1_path), path(read2_path)
    
    output:
    tuple val(sample_id), val(genome_id), path("sample${sample_id}_${genome_id}.bam"), path(fasta_file)

    script:
    """
    samtools view -bS sample${sample_id}_${genome_id}.sam | samtools sort -o sample${sample_id}_${genome_id}.bam
    mkdir --parents ${projectDir}/nextflow_out/sample_${sample_id}/bam/
    cp sample${sample_id}_${genome_id}.bam ${projectDir}/nextflow_out/sample_${sample_id}/bam/
    """
}
