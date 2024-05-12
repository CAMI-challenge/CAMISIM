/*
 * Defining the module / subworkflow path, and include the elements
 */

// include read simulator here:
read_simulator_folder = "${projectDir}/pipelines/metatranscriptomic/read_simulators/"
include { read_simulator_art } from "${read_simulator_folder}/read_simulator_art"
include { read_simulator_nanosim3 } from "${read_simulator_folder}/read_simulator_nansoim3"
include { read_simulator_pbsim3 } from "${read_simulator_folder}/read_simulator_pbsim3"

include { normalise_abundance_meta_t; normalise_abundance_to_size; count_bases} from "${projectDir}/distribution"

/** 
* This workflow simulates reads for every sample.
* Takes:
*     genome_location_ch: Path to the file containing the genome locations.
      genome_distribution_ch: Paths to files containing the distributions of every genome for every sample.
* Emits: 
*     A channel containing the merged fastq file and the merged bam file over all genomes for every sample.
**/
workflow sample_wise_simulation {

    take: genome_location_file_ch
    take: distribution_file_ch
    take: gene_distribution_file_ch
    take: read_length_ch
    take: seed_file_ch
    take: annotation_file_ch
    
    main:

        // this channel holds the genome location map (key = genome_id, value = absolute path to genome)
        genome_location_ch = genome_location_file_ch
            .splitCsv(sep:'\t') // get genome id and relatvie path from genome location file
            .map { genome_id, path ->
                def abs_path
                if (new File(path).isAbsolute()) { // if the path is an absolute path return it as is
                    abs_path = path
                } else { // else expand relative paths to absolute paths and send to genome_location_ch
                    abs_path = file("${projectDir}/${path}").toAbsolutePath().toString()
                }
                return [genome_id, abs_path]
            }

        // The simulated ART reads (version 016.06.05) doesn't contain the whole header of the reference genome, if there is a space in the header. They then just contain the 
        // substring before the first occurance of the space. In that case the gold standard assembly doesn't work because there are no matching IDs.
        // As a workaround we change the headers of the reference genomes by just selecting the part before the first occurance of a space, if there is a space in the header.
        if(params.type.equals("art")) {
            genome_location_ch = remove_spaces_from_reference_genome(genome_location_ch)
        }    

        // get the seed for every genome
        seed_ch = seed_file_ch.splitCsv(sep:'\t', skip:2)

        // normalize the distributions and read the results from the generated file
        normalised_distribution_ch = normalise_abundance_meta_t(distribution_file_ch.map { a -> tuple(a[1], tuple (a[0], a[2])) }.groupTuple())
            .map { file -> tuple(file.baseName.split('_')[2], file) }.splitCsv(sep:'\t').map { a -> tuple(a[1][0], a[1][1],a[0]) }

        // calculate gene expression
        final_gene_distr_ch = get_final_gene_distr(gene_distribution_file_ch.join(normalised_distribution_ch, by: [0,2]))

        // get_read_count(final_gene_distr_ch)

        location_distribution_seed_ch = final_gene_distr_ch.combine(genome_location_ch, by: 0).combine(annotation_file_ch, by: 0).combine(seed_ch, by: [0,1])

        // read simulation
        if(params.type.equals("art")) {

            read_simulator_art(location_distribution_seed_ch, read_length_ch)

            bam_files_channel = read_simulator_art.out[0]
            reads_ch = read_simulator_art.out[1]

            get_fastq_for_sample_paired_end(reads_ch)

        } else if(params.type.equals("nanosim3")) {

            // simulate the reads with nanosim3
            read_simulator_nanosim3(location_distribution_seed_ch, read_length_ch)

            bam_files_channel = read_simulator_nanosim3.out[0]
            reads_ch = read_simulator_nanosim3.out[1]

            get_fastq_for_sample_single_end(reads_ch)

        } else if (params.type.equals("pbsim3")) {

            // simulate the reads with pbsim3
            read_simulator_pbsim3(location_distribution_seed_ch, read_length_ch)

            bam_files_channel = read_simulator_pbsim3.out[0]
            reads_ch = read_simulator_pbsim3.out[1]

            get_fastq_for_sample_single_end(reads_ch)

        }

        // generate gold standard assembly for every genome and copy it into output folder
        gsa_for_every_genome_ch = generate_gold_standard_assembly(bam_files_channel)

        // grouping the gold standard assemblies by sample id results in new tuple: key = sample_id, values = path to all gsa of this samples reads
        grouped_gsa_for_every_genome_ch = gsa_for_every_genome_ch.groupTuple()

        // create fasta files holding all gsa of one samples reads
        get_fasta_for_sample(grouped_gsa_for_every_genome_ch)

        // group bam files by sample id
        bam_files_by_sample_ch = bam_files_channel.groupTuple().map { a -> tuple(a[0], a[2]) }
        // create bam files holding all bam files of one samples reads 
        merge_bam_files(bam_files_by_sample_ch)

    emit: merge_bam_files.out
    emit: get_fasta_for_sample.out
    emit: reads_ch
    emit: bam_files_by_sample_ch
}

process get_final_gene_distr {

    conda 'anaconda::python=3.6'

    publishDir "${params.outdir}/distributions/final_distributions/", pattern: "${genome_id}_${sample_id}_final_distribution.tsv", mode: 'copy'

    input:
    tuple val(genome_id), val(sample_id), path(gene_distribution_file), val(genome_distribution)

    output:
    tuple val(genome_id), val(sample_id), path("${genome_id}_${sample_id}_final_distribution.tsv"), val(genome_distribution)

    script:
    """
    python <<CODE
    import csv

    genome_distribution = float("${genome_distribution}")
    input_file = "${gene_distribution_file}"
    output_file = "${genome_id}_${sample_id}_final_distribution.tsv"

    # Process the gene distribution file
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\\t')
        writer = csv.writer(outfile, delimiter='\\t')

        for row in reader:
            gene_id, abundance = row
            adjusted_abundance = float(abundance) * genome_distribution
            writer.writerow([gene_id, f'{adjusted_abundance:.10f}'])

    CODE
    """
}

/*
* This process splits all strings in a fasta header line at space character.\
* It then changes the header to the first substring.
*     
 */
process remove_spaces_from_reference_genome {

    input:
    tuple val(genome_id), path(fasta_file)

    output:
    tuple val(genome_id), path(fasta_file)

    script:
    """
    if [ -e ${fasta_file} ]
    then
    mv ${fasta_file} ${fasta_file}_to_rename
    awk '{if(\$0 ~ /^>/) {split(\$0,a," "); print a[1]} else {print \$0}}' ${fasta_file}_to_rename > ${fasta_file}
    else
    echo "File not found: ${fasta_file}"
    exit 1
    fi
    """
}

/*
* This process writes all fastq files into a single one.
 */
process get_fastq_for_sample_paired_end {

    // For some reason this does not work.
    //publishDir "${params.outdir}/sample_${sample_id}/reads/fastq", pattern: "sample_${sample_id}_01.fq.gz", mode: 'copy'
    //publishDir "${params.outdir}/sample_${sample_id}/reads/fastq", pattern: "sample_${sample_id}_02.fq.gz", mode: 'copy'

    input:
    tuple val(sample_id), path(first_read_files), path(second_read_files)

    script:
    """
    # cat ${first_read_files} > sample_${sample_id}_01.fq
    # cat ${second_read_files} > sample_${sample_id}_02.fq

    # Sort files before concatenation to ensure reproducibility
    ls -1 ${first_read_files} | sort | xargs cat > sample_${sample_id}_01.fq
    ls -1 ${second_read_files} | sort | xargs cat > sample_${sample_id}_02.fq

    # Compress the concatenated files
    gzip sample_${sample_id}_01.fq
    gzip sample_${sample_id}_02.fq

    mkdir --parents ${params.outdir}/sample_${sample_id}/reads/fastq
    cp sample_${sample_id}_01.fq.gz ${params.outdir}/sample_${sample_id}/reads/fastq/
    cp sample_${sample_id}_02.fq.gz ${params.outdir}/sample_${sample_id}/reads/fastq/
    """
}

/*
* This process writes all fastq files into a single one.
 */
process get_fastq_for_sample_single_end {

    // For some reason this does not work.
    // publishDir "${params.outdir}/sample_${sample_id}/reads/fastq", pattern: "*.gz", mode: 'copy'

    input:
    tuple val(sample_id), path(read_files)

    script:
    """
    # cat ${read_files} > sample_${sample_id}.fq

    # Sort files before concatenation to ensure reproducibility
    ls -1 ${read_files} | sort | xargs cat > sample_${sample_id}.fq

    # Compress the concatenated files
    gzip sample_${sample_id}.fq

    mkdir --parents ${params.outdir}/sample_${sample_id}/reads/fastq
    cp sample_${sample_id}.fq.gz ${params.outdir}/sample_${sample_id}/reads/fastq/
    """
}

/*
* This process generates a gold standard assembly for one genome.
* Takes:
*     A tuple with key = sample_id, first_value = genome_id, second value = a sorted bam file, third value = the reference genome (fasta).
* Output:
*     A Tuple with key = sample_id, value = path to fasta file with the gold standard assembly of the given genome.
 */
process generate_gold_standard_assembly {

    conda 'bioconda::samtools=1.20 conda-forge::perl=5.32.1'

    input:
    tuple val(sample_id),val(genome_id), path(bam_file), path(reference_fasta_file)

    output:
    tuple val(sample_id), path(file_name)

    when: params.gsa  // Only execute this process when gsa is set to true
    
    script:
    file_name = 'sample'.concat(sample_id.toString()).concat('_').concat(genome_id).concat('_gsa.fasta')
    """
    perl -- ${projectDir}/scripts/bamToGold.pl -st samtools -r ${reference_fasta_file} -b ${bam_file} -l 1 -c 1 >> ${file_name}
    mkdir --parents ${params.outdir}/sample_${sample_id}/gsa
    gzip -k ${file_name}
    cp ${file_name}.gz ${params.outdir}/sample_${sample_id}/gsa/
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
    tuple val(sample_id), path(file_name)
        
    script:
    file_name = 'sample'.concat(sample_id.toString()).concat('_gsa.fasta')
    """
    # cat ${fasta_files} > ${file_name}

    # Sort files before concatenation to ensure reproducibility
    ls -1 ${fasta_files} | sort | xargs cat > ${file_name}

    mkdir --parents ${params.outdir}/sample_${sample_id}/contigs
    gzip -k ${file_name}
    cp ${file_name}.gz ${params.outdir}/sample_${sample_id}/contigs/gsa.fasta.gz
    """
}

/*
* This process merges all given bam files with samtools.
* Takes:
*     A tuple with key = sample_id, value = the paths to all bam files, that need to be combined.
* Output:
*     The path to the merged bam file.
 */
process merge_bam_files {

    conda 'bioconda::samtools=1.13'

    input:
    tuple val(sample_id), path(bam_files)

    output:
    tuple val(sample_id), path(file_name)

    script:
    file_name = 'sample_'.concat(sample_id.toString()).concat('.bam')
    compression = 5
    memory = 1
    """
    samtools merge -u - ${bam_files} | samtools sort -l ${compression} -m ${memory}G -o ${file_name} -O bam
    samtools index ${file_name}
    """
}