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

    conda 'bioconda::pbsim3=3.0.0 bioconda::gffread=0.12.7 bioconda::gffutils conda-forge::python==3.7 bioconda::pyfaidx bioconda::samtools'

    input:
    tuple val(genome_id), val(sample_id), path(fasta_distribution_file), val(abundance), path(fasta_file), path(gff_file), val(seed)
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
    python <<CODE
    import csv
    import random
    import gffutils
    import pyfaidx

    random.seed(${seed})

    db_file = 'db_file.db'
    db = gffutils.create_db('${gff_file}', dbfn=db_file)
    fasta = pyfaidx.Fasta('${fasta_file}')

    total_read_count = 0

    # Specify the path for your output file
    with open("${sample_id}_${genome_id}_expression_profile.tsv", 'w') as exp_f:

        with open("transcript_id_to_seq_id.tsv", 'w') as transcript_id_to_seq_id_file:
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

                    total_read_count = total_read_count + read_count

                    gene = db[gene_identifier]

                    if read_count > 0:
                        # exp_f.write(gene_identifier_modified+"\\t0\\t"+str(read_count)+"\\n")
                        exp_f.write(gene_identifier+"\\t"+str(read_count)+"\\t0\\t"+str(gene.sequence(fasta))+"\\n")

                        # transcript_id_to_seq_id_file.write(gene_identifier_modified+"\\t"+str(gene.seqid)+"\\n")
                        # transcript_id_to_seq_id_file.write(gene_identifier+"\\t"+str(gene.seqid)+"\\n")
                        transcript_id_to_seq_id_file.write(gene_identifier+"\\t"+str(gene.seqid)+"\\t"+str(gene.start)+"\\n")

    with open('total_read_count.txt', 'w') as count_file:
        count_file.write(str(total_read_count))

    CODE

    #pbsim --strategy trans --method ${method} --${method} ${model} --accuracy-mean 0.85 --difference-ratio ${difference_ratio} --transcript ${sample_id}_${genome_id}_expression_profile.tsv --seed ${seed} --prefix sample${sample_id}_${genome_id}_pbsim3 --length-mean ${length_mean} --length-sd ${length_sd}
    pbsim --strategy trans --method ${method} --${method} ${model} --accuracy-mean 0.85 --difference-ratio ${difference_ratio} --transcript ${sample_id}_${genome_id}_expression_profile.tsv --seed ${used_seed} --prefix sample${sample_id}_${genome_id}_pbsim3 --length-mean ${length_mean} --length-sd ${length_sd}

    python ${projectDir}/pipelines/metatranscriptomic/read_simulators/maf_converter.py --transcriptome 
    samtools view -bS sample${sample_id}_${genome_id}.sam | samtools sort -o sample${sample_id}_${genome_id}.bam

    gzip -k sample${sample_id}_${genome_id}.fq

    mkdir --parents ${projectDir}/out/sample_${sample_id}/bam/
    cp sample${sample_id}_${genome_id}.bam ${projectDir}/out/sample_${sample_id}/bam/
    mkdir --parents ${projectDir}/out/sample_${sample_id}/reads/fastq/
    cp sample${sample_id}_${genome_id}.fq.gz ${projectDir}/out/sample_${sample_id}/reads/fastq/
    """
}