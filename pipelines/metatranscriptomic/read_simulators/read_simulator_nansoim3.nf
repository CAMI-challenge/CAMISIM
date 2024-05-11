/** 
* This workflow simulates reads via nanosim3 and converts the resulting sam files into bam files.
* Takes:
*     A channel containing tuples with key = genome_id, first value = path to genome, second value = distribution, third value = sample_id, fourth value = seed.
* Emits: 
*     A channel containing tuples with key = sample_id, first value = genome id, second value = simulated bam file, third value = the reference fasta file.
**/
workflow read_simulator_nanosim3 {

    take: genome_location_distribution_ch
    take: read_length_ch
    main:
        // simulate reads via nanosim3
        out = simulate_reads_nanosim3(genome_location_distribution_ch, read_length_ch)

        read_ch = out[1].groupTuple()

        bam_ch = bam_from_reads(out[0])

    emit:
        bam_ch
        read_ch
}

/**
* This process simulates reads in fastq format with nanosim3.
* Input:
*     Tuple containing key = genome_id, first value = path to genome, second value = distribution, third value = sample_id, fourth value = seed.
* Output:
*     Tuple containing key = sample_id, first value = genome_id, second value = path to error profile, third value = path to fasta file with the aligned reads, 
*         fourth value = path to fasta file with the aligned reads, fifth value = path to reference genome.
**/
process simulate_reads_nanosim3 {

    conda 'anaconda::scikit-learn=0.21.3=py37hd81dba3_0 bioconda::nanosim=3.1.0 bioconda::gffread=0.12.7 bioconda::gffutils=0.9 conda-forge::python==3.7'

    // ToDo For some reason this does not work
    //publishDir "${params.outdir}/sample_${sample_id}/reads/fastq/", pattern: "*.gz", mode: 'copy'
    //publishDir "${params.outdir}/sample_${sample_id}/reads/fastq/", pattern: "*_aligned_reads.fastq.gz", mode: 'copy'

    input:
    tuple val(genome_id), val(sample_id), path(fasta_distribution_file), val(abundance), path(fasta_file), path(gff_file), val(seed)
    val(read_length_ch)

    
    output:
    tuple val(sample_id), val(genome_id), path('*_error_profile'), path("*_aligned_reads.fastq"), path("*_unaligned_reads.fastq"), path(fasta_file), path("transcript_id_to_seq_id.tsv")
    tuple val(sample_id), path("*_aligned_reads.fastq")
    //tuple val(sample_id), val(genome_id), path('*_error_profile'), path("*_aligned_reads.fasta"), path("*_unaligned_reads.fasta"), path(fasta_file), path("transcript_id_to_seq_id.tsv")
    //tuple val(sample_id), path("*_aligned_reads.fasta")
    
    script:
    profile = params.base_profile_name
    size = params.size
    read_length = read_length_ch
    basecaller = params.basecaller
    number_of_reads = (size*(10**9)) * abundance.toFloat() / read_length_ch.toFloat()
    number_of_reads = number_of_reads.round(0).toInteger()
    // nanosim seed cannot be > 2**32 -1
    //Long seed = 632741178
    //Long used_seed = seed % 2**32 - 1
    Long used_seed = (seed as Long) % 2**32 - 1

    /**
    String log = "---- sample id: ".concat(sample_id)
    log = log.concat("  read length: ").concat(Integer.toString(read_length_ch))
    log = log.concat("  genome id: ").concat(genome_id)
    log = log.concat("   fasta file: ").concat(fasta_file.baseName)
    log = log.concat("  abundance: ").concat(abundance)
    log = log.concat("    seed: ").concat(Long.toString(seed))
    log = log.concat("    used_seed: ").concat(Long.toString(used_seed))
    log = log.concat("    number_of_reads: ").concat(Integer.toString(number_of_reads))
    log = log.concat("    profile: ").concat(profile)
    print(log)
    **/
    """
    python <<CODE
    import csv
    import random
    import gffutils

    random.seed(${seed})

    db_file = 'db_file.db'
    db = gffutils.create_db('${gff_file}', dbfn=db_file)

    total_read_count = 0

    # Specify the path for your output file
    output_file = "sample${sample_id}_${genome_id}.gff3"
    with open("${sample_id}_${genome_id}_expression_profile.tsv", 'w') as exp_f:

        exp_f.write("target_id\\test_counts\\ttpm\\n")

        with open("transcript_id_to_seq_id.tsv", 'w') as transcript_id_to_seq_id_file:
            with open(output_file, 'w') as out_f:
                # Iterate through each file
                with open('${fasta_distribution_file}', 'r') as file:
                    reader = csv.reader(file, delimiter='\\t')
                    # Iterate through each row in the TSV file
                    for row in reader:
                        gene_identifier, abundance = row

                        read_count = (${size}*(10**9)) * float(abundance) / ${read_length}

                        if read_count < 1:
                            read_count = 1 if random.random() < read_count else 0

                        # gene_identifier_modified = gene_identifier.split("transcript:")[1]

                        read_count = round(read_count)

                        # exp_f.write(gene_identifier_modified+"\\t0\\t"+str(read_count)+"\\n")
                        exp_f.write(gene_identifier+"\\t0\\t"+str(read_count)+"\\n")

                        total_read_count = total_read_count + read_count

                        gene = db[gene_identifier]

                        out_f.write(str(gene)+"\\n")

                        # transcript_id_to_seq_id_file.write(gene_identifier_modified+"\\t"+str(gene.seqid)+"\\n")
                        # transcript_id_to_seq_id_file.write(gene_identifier+"\\t"+str(gene.seqid)+"\\n")
                        transcript_id_to_seq_id_file.write(gene_identifier+"\\t"+str(gene.seqid)+"\\t"+str(gene.start)+"\\n")

    with open('total_read_count.txt', 'w') as count_file:
        count_file.write(str(total_read_count))

    CODE

    total_read_count=\$(cat total_read_count.txt)

    gffread -F -w transcriptome.fa -g ${fasta_file} sample${sample_id}_${genome_id}.gff3

    simulator.py transcriptome -rt transcriptome.fa -c ${profile} -e ${sample_id}_${genome_id}_expression_profile.tsv -n \$total_read_count --no_model_ir -b ${basecaller} --fastq -r cDNA_1D --seed ${used_seed} -o sample${sample_id}_${genome_id}

    # gzip -k sample${sample_id}_${genome_id}_aligned_reads.fastq
    gzip -k *_aligned_reads.fastq
    mkdir --parents ${params.outdir}/sample_${sample_id}/reads/fastq/
    cp *_aligned_reads.fastq.gz ${params.outdir}/sample_${sample_id}/reads/fastq/
    """
}

/**
* This process creates a bam file from the simulated reads in fastq format.
* Input:
*     Tuple containing key = sample_id, first value = genome_id, second value = path to error profile, third value = path to fasta file with the aligned reads, 
*         fourth value = path to fasta file with the aligned reads, fifth value = path to reference genome.
* Output:
*    Tuple containing key = sample_id, first value = genome_id, second value = path to bam file, third value = path to reference genome.
**/
process bam_from_reads {

    conda "bioconda::biopython=1.70 anaconda::python=3.5.6 bioconda::samtools=1.13 bioconda::gffutils=0.9"

    publishDir "${params.outdir}/sample_${sample_id}/bam/", pattern: "sample*.bam", mode: 'copy'

    input:
    tuple val(sample_id), val(genome_id), val(error_profile), path(aligned_reads), path(unaligned_reads), path(fasta_file), path(map_transcript_seq_file)

    output:
    tuple val(sample_id), val(genome_id), path('sample*.bam'), path(fasta_file)

    script:
    """
    ${projectDir}/read_simulators/sam_from_reads.py ${error_profile} ${aligned_reads} ${unaligned_reads} ${fasta_file} --transcript_seq_id_map ${map_transcript_seq_file} --transcriptome
    samtools view -bS *.sam | samtools sort -o sample${sample_id}_${genome_id}.bam
    """
}