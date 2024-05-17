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
        simulate_reads_art(genome_location_distribution_ch, read_length_ch)
    emit:
        simulate_reads_art.out[0]
        simulate_reads_art.out[1].groupTuple()
}

process simulate_reads_art {

    conda 'bioconda::art=2016.06.05 conda-forge::gsl=2.7 bioconda::samtools=1.13 anaconda::python=3.6.13 bioconda::gffutils=0.9 bioconda::bedtools=2.31.1 bioconda::pyfaidx=0.5.3'

    //publishDir "${params.outdir}/distributions/final_distributions/", pattern: "${sample_id}_${genome_id}_read_counts.tsv", mode: 'copy'

    publishDir "${params.outdir}/sample_${sample_id}/bam/", pattern: "sample${sample_id}_${genome_id}.bam", mode: 'copy'
    // For some reason this does not work.
    //publishDir "${params.outdir}/sample_${sample_id}/reads/fastq/", pattern: "sample${sample_id}_${genome_id}_1.fq.gz", mode: 'copy'
    //publishDir "${params.outdir}/sample_${sample_id}/reads/fastq/", pattern: "sample${sample_id}_${genome_id}_2.fq.gz", mode: 'copy'

    input:
    tuple val(genome_id), val(sample_id), path(fasta_distribution_file), val(abundance), path(fasta_file), path(gff_file), val(seed), path(db)
    val(read_length_ch)

    output:
    tuple val(sample_id), val(genome_id), path("sample${sample_id}_${genome_id}.bam"), path(fasta_file)
    tuple val(sample_id), path("sample${sample_id}_${genome_id}_1.fq"), path("sample${sample_id}_${genome_id}_2.fq")

    script:
    fragment_size_mean = params.fragment_size_mean
    fragment_size_sd = params.fragment_size_sd
    profile = params.base_profile_name
    size = params.size
    read_length = read_length_ch.toFloat()
    child_feature_type = params.child_feature_type
    """
    python ${projectDir}/pipelines/metatranscriptomic/read_simulators/generate_reads_and_modify_sam.py \
        --seed ${seed} \
        --fasta_distribution_file ${fasta_distribution_file} \
        --gff_file ${gff_file} \
        --child_feature_type ${child_feature_type} \
        --fasta_file ${fasta_file} \
        --sample_id ${sample_id} \
        --genome_id ${genome_id} \
        --size ${size} \
        --read_length ${read_length} \
        --fragment_size_mean ${fragment_size_mean} \
        --fragment_size_sd ${fragment_size_sd} \
        --profile ${profile} \
        --db ${db}

    chmod +x ${sample_id}_${genome_id}_commands.sh
    ./${sample_id}_${genome_id}_commands.sh

    cat sample${sample_id}_${genome_id}_*1.fq > sample${sample_id}_${genome_id}_1.fq
    cat sample${sample_id}_${genome_id}_*2.fq > sample${sample_id}_${genome_id}_2.fq

    gzip -k sample${sample_id}_${genome_id}_1.fq
    gzip -k sample${sample_id}_${genome_id}_2.fq

    mkdir --parents ${params.outdir}/sample_${sample_id}/reads/fastq/
    cp sample${sample_id}_${genome_id}_1.fq.gz ${params.outdir}/sample_${sample_id}/reads/fastq/
    cp sample${sample_id}_${genome_id}_2.fq.gz ${params.outdir}/sample_${sample_id}/reads/fastq/

    samtools merge sample${sample_id}_${genome_id}.sam sample${sample_id}_${genome_id}_*.sam
    samtools view -bS sample${sample_id}_${genome_id}.sam | samtools sort -o sample${sample_id}_${genome_id}.bam
    """
}
