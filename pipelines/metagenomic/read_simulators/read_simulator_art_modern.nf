/**
* This workflow simulates reads using the ART_modern read simuator
* Takes:
*     A channel containing tuples with key = genome_id, first value = path to genome, second value = distribution, third value = sample_id, fourth value = factor.
*     A channel containing the read length.
* Emits: 
*     A channel containing tuples with key = sample_id, first value = genome id, second value = simulated bam file, third value = the reference fasta file.
**/

workflow read_simulator_art_modern {

    take: genome_location_distribution_ch
    main:
        simulate_reads_art_modern(genome_location_distribution_ch)
    emit:
        simulate_reads_art_modern.out[0]
        simulate_reads_art_modern.out[1]
}

/**
* This process simulates short reads with ART_modern
* Input:
*     Tuple containing key = genome_id, first value = path to genome, second value = distribution, third value = sample_id, fourth value = factor.
* Output:
*     Tuple containing key = sample_id, first value = genome_id, second value = path to error profile, third value = path to fasta file with the aligned reads, 
*         fourth value = path to fasta file with the aligned reads, fifth value = path to reference genome.
**/
process simulate_reads_art_modern {

    scratch true
    
    conda 'bioconda::art_modern bioconda::seqtk bioconda::samtools bioconda::seqkit conda-forge::pigz'
    
    input:
    tuple val(genome_id), val(sample_id), path(fasta_file), val(abundance), val(seed), val(factor)
    
    output:
    tuple val(sample_id), val(genome_id), path("sample${sample_id}_${genome_id}_art_modern.bam"), path(fasta_file)
    tuple val(sample_id), path('*1.fq'), path('*2.fq')
   
    script:
    fragment_size_mean = params.art_modern.fragment_size_mean
    fragment_size_sd = params.art_modern.fragment_size_sd
    profile = params.art_modern.base_profile_name
    read_length = params.art_modern.profile_read_length
    factor_float_value = Double.valueOf(factor)
    fold_coverage = Double.valueOf(abundance) * factor_float_value // TODO is the abundance already normalised?
    threads = Math.max(1, ((task.cpus ?: 1) as int))

    /**
    String log = "---- sample id: ".concat(sample_id)
    log = log.concat("  genome id: ").concat(genome_id)
    log = log.concat("   fasta file: ").concat(fasta_file.baseName)
    log = log.concat("  fragment_size_mean: ").concat(Integer.toString(fragment_size_mean))
    log = log.concat("    fragment_size_sd: ").concat(Integer.toString(fragment_size_sd))
    log = log.concat("    profile: ").concat(profile)
    log = log.concat("    factor_float_value: ").concat(Double.toString(factor_float_value))
    log = log.concat("    fold_coverage: ").concat(Double.toString(fold_coverage))
    log = log.concat("    seed: ").concat(seed)
    print(log)
    **/

    """
    art_modern --mode wgs --i-file ${fasta_file} --read_len ${read_length} \
        --pe_frag_dist_mean ${fragment_size_mean} --pe_frag_dist_std_dev ${fragment_size_sd} \
        --i-fcov ${fold_coverage} --min_qual 15 \
        --qual_file_1 ${profile}1.txt --qual_file_2 ${profile}2.txt --i-seed ${seed} \
        --o-fastq sample${sample_id}_${genome_id}_art_modern.fq \
        --o-sam sample${sample_id}_${genome_id}_art_modern.bam --o-sam-write_bam
    samtools sort -o tmp.bam sample${sample_id}_${genome_id}_art_modern.bam
    mv tmp.bam sample${sample_id}_${genome_id}_art_modern.bam
    mkdir --parents ${params.outdir}/sample_${sample_id}/bam/art_modern/
    cp sample${sample_id}_${genome_id}_art_modern.bam ${params.outdir}/sample_${sample_id}/bam/art_modern/
    for file in sample${sample_id}_${genome_id}_art_modern*.fq;
    do
    seqtk seq "\$file" -1 > "\${file%.fq}1.fq"
    seqtk seq "\$file" -2 > "\${file%.fq}2.fq"
    pigz -p ${threads} -c "\${file%.fq}1.fq" > "\${file%.fq}1.fq.gz";
    pigz -p ${threads} -c "\${file%.fq}2.fq" > "\${file%.fq}2.fq.gz";
    done
    mkdir --parents ${params.outdir}/sample_${sample_id}/reads/fastq/art_modern/
    cp sample${sample_id}_${genome_id}_art_modern*1.fq.gz ${params.outdir}/sample_${sample_id}/reads/fastq/art_modern/
    cp sample${sample_id}_${genome_id}_art_modern*2.fq.gz ${params.outdir}/sample_${sample_id}/reads/fastq/art_modern/
    """
}
