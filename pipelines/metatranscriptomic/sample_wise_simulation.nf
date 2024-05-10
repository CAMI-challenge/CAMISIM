/*
 * Defining the module / subworkflow path, and include the elements
 */

// include read simulator here:
read_simulator_folder = "${projectDir}/pipelines/metatranscriptomic/read_simulators/"
include { read_simulator_art } from "${read_simulator_folder}/read_simulator_art"

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

        }
}

process get_final_gene_distr {

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

    mkdir --parents ${params.outdir}/distributions/final_distributions/
    cp ${genome_id}_${sample_id}_final_distribution.tsv ${params.outdir}/distributions/final_distributions/
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

    publishDir "${params.outdir}/sample_${sample_id}/reads/fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(first_read_files), path(second_read_files)

    script:
    """
    cat ${first_read_files} > sample_${sample_id}_01.fq
    cat ${second_read_files} > sample_${sample_id}_02.fq
    """
}
