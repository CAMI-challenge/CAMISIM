scripts_dir = "${projectDir}/pipelines/metatranscriptomic/scripts"
shared_scripts_dir = "${projectDir}/pipelines/shared/scripts"

workflow read_simulator_pbsim3 {

    take: genome_location_distribution_ch
    take: read_length_ch
    main:

        // simulate reads via pbsim3
        simulate_reads_pbsim3(genome_location_distribution_ch, read_length_ch)

    emit:
        simulate_reads_pbsim3.out[0]
        simulate_reads_pbsim3.out[1].groupTuple()
}

process simulate_reads_pbsim3 {

    conda 'bioconda::pbsim3=3.0.0 bioconda::gffread=0.12.7 bioconda::gffutils conda-forge::python==3.7 bioconda::pyfaidx bioconda::samtools=1.9'

    input:
    tuple val(genome_id), val(sample_id), path(fasta_distribution_file), val(abundance), path(fasta_file), path(gff_file), val(seed), path(db)
    val(read_length_ch)

    output:
    tuple val(sample_id), val(genome_id), path("sample${sample_id}_${genome_id}.bam"), path(fasta_file)
    tuple val(sample_id), path("sample${sample_id}_${genome_id}.fq")

    
    //output:
    //tuple val(sample_id), val(genome_id), path('*_error_profile'), path("*_aligned_reads.fastq"), path("*_unaligned_reads.fastq"), path(fasta_file), path("transcript_id_to_seq_id.tsv")
    //tuple val(sample_id), path("*_aligned_reads.fastq")
    
    script:
    model = params.model
    size = params.size
    read_length = read_length_ch
    length_mean = params.fragments_size_mean
    length_sd = params.fragment_size_standard_deviation
    method = params.method
    difference_ratio = params.difference_ratio
    // when the  seed gets to big, the simulation fails
    Long used_seed = (seed as Long) % 2**32 - 1

    child_feature_type = params.child_feature_type

    /**
    String log = "---- sample id: ".concat(sample_id)
    log = log.concat("  read length: ").concat(Integer.toString(read_length_ch))
    log = log.concat("  genome id: ").concat(genome_id)
    log = log.concat("   fasta file: ").concat(fasta_file.baseName)
    log = log.concat("  abundance: ").concat(abundance)
    log = log.concat("    seed: ").concat(seed)
    log = log.concat("    length_mean: ").concat(Integer.toString(length_mean))
    log = log.concat("    length_sd: ").concat(Integer.toString(length_sd))
    log = log.concat("    profile: ").concat(profile)
    print(log)
    **/
    """
    python ${scripts_dir}/create_expression_profile_pbsim3.py \
        --seed ${seed} \
        --db ${db} \
        --fasta_file ${fasta_file} \
        --sample_id ${sample_id} \
        --genome_id ${genome_id} \
        --size ${size} \
        --read_length ${read_length} \
        --child_feature_type ${child_feature_type} \
        --fasta_distribution_file ${fasta_distribution_file}

    #pbsim --strategy trans --method ${method} --${method} ${model} --accuracy-mean 0.85 --difference-ratio ${difference_ratio} --transcript ${sample_id}_${genome_id}_expression_profile.tsv --seed ${seed} --prefix sample${sample_id}_${genome_id}_pbsim3 --length-mean ${length_mean} --length-sd ${length_sd}
    pbsim --strategy trans --method ${method} --${method} ${model} --accuracy-mean 0.85 --difference-ratio ${difference_ratio} --transcript sample_${sample_id}_${genome_id}_expression_profile.tsv --seed ${used_seed} --prefix sample${sample_id}_${genome_id}_pbsim3 --length-mean ${length_mean} --length-sd ${length_sd}

    python ${scripts_dir}/maf_converter.py --transcriptome 
    samtools view -bS sample${sample_id}_${genome_id}.sam | samtools sort -o sample${sample_id}_${genome_id}.bam

    gzip -k sample${sample_id}_${genome_id}.fq

    mkdir --parents ${params.outdir}/sample_${sample_id}/bam/
    cp sample${sample_id}_${genome_id}.bam ${params.outdir}/sample_${sample_id}/bam/
    mkdir --parents ${params.outdir}/sample_${sample_id}/reads/fastq/
    cp sample${sample_id}_${genome_id}.fq.gz ${params.outdir}/sample_${sample_id}/reads/fastq/
    """
}