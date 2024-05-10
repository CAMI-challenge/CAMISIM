#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Defining the module / subworkflow path, and include the elements
 */

// include sample wise simulation
include { sample_wise_simulation } from "${projectDir}/pipelines/metatranscriptomic/sample_wise_simulation"

/*
 * This is the main workflow and starting point of this nextflow pipeline.
 */
workflow metatranscriptomic {

    if(params.seed != null) {
            seed = params.seed
        } else {
            seed = get_random_seed()
        }

    annotation_file_ch = Channel.fromPath(params.gene_annotations_file)


    // this channel holds the genome location map (key = genome_id, value = absolute path to genome)
    annotation_file_ch = annotation_file_ch
        .splitCsv(sep:'\t') // get genome id and relatvie path from genome location file
        .map { genome_id, gff_path ->
            def abs_gff_path
            if (new File(gff_path).isAbsolute()) { // if the path is an absolute path return it as is
                abs_gff_path = gff_path
            } else { // else expand relative paths to absolute paths and send to genome_location_ch
                abs_gff_path = file("${projectDir}/${gff_path}").toAbsolutePath().toString()
            }
            return [genome_id, abs_gff_path]
        }

    // this channel holds the file with the specified locations of the genomes
    genome_location_file_ch = Channel.fromPath(params.genome_locations_file)

    // random seed generation
    get_seed(genome_location_file_ch, seed)

    distribute_gene_abundance(annotation_file_ch)
    gene_distribution_file_ch = distribute_gene_abundance.out[0].flatMap { item ->
        def genomeName = item[0]
        def files = item[1]

        // Transform each file into a new channel item
        files.collect { file ->
            def index = file.baseName.split("_")[-1] // Extracts the numerical part after the last underscore and before the '.tsv'
            return [genomeName, file, index] // Creates the new channel item
        }
    }

    // calculate the genome distributions for each sample for one community
    genome_distribution_file_ch = getCommunityDistribution(genome_location_file_ch, seed).flatten().map { file -> tuple(file.baseName.split('_')[1], file) }.splitCsv(sep:'\t').map { a -> tuple(a[1][0], a[0], a[1][1]) }

    distribution_file_ch = genome_distribution_file_ch.combine(distribute_gene_abundance.out[1], by: 0).map{ item -> tuple( item[0], item[1], Float.valueOf(item[2])*Float.valueOf(item[3]))}

    read_length_ch = params.read_length

    // simulate reads sample wise
    sample_wise_simulation(genome_location_file_ch, distribution_file_ch, gene_distribution_file_ch, read_length_ch, get_seed.out[0], annotation_file_ch)
}

/*
* This process distributes the abundance of the genes of the given genome.
*
*/
process distribute_gene_abundance {

    //publishDir "${params.outdir}/distributions/gene_distributions/", mode : 'copy'

    conda 'bioconda::gffutils=0.12 anaconda::python=3.6'

    input:
    tuple val(genome_id), path(gene_annotations_file)

    output:
    tuple val(genome_id), path("distribution_*.tsv")
    tuple val(genome_id), stdout

    script:
    number_of_samples = params.number_of_samples
    seed = params.seed
    mode = params.mode
    mu = params.mu
    sigma = params.sigma

    gauss_mu = ""
    gauss_sigma = ""
    gene_sigma = ""

    if(mode.equals("timeseries")) {
        gauss_mu = "-gauss_mu " + params.gauss_mu
        gauss_sigma = "-gauss_sigma " + params.gauss_sigma
    } else if(mode.equals("replicate")) {
        gene_sigma = "-gene_sigma " + params.gene_sigma
    }
    """
    python ${projectDir}/pipelines/metatranscriptomic/scripts/get_gene_abundance.py -annotation_file ${gene_annotations_file} -mode ${mode} -number_of_samples ${number_of_samples} -seed ${seed} -log_mu ${mu} -log_sigma ${sigma} ${gauss_mu} ${gauss_sigma} ${gene_sigma}
    """
}

process get_random_seed {

    output:
    stdout

    script:
    """
    #!/usr/bin/env python
    import random

    randomnumber = random.randint(0, ((2**32)-1))

    print(randomnumber, end = "")
    """
}

/*
* This process returns a file containing a random seed for every genome generated from the given seed in the config file.
* In case the simulated reads will be anonymized, it also returns a file containing a random seed for every sample generated from the given seed in the config file.
* Output:
*     The file with the given seed per samle in CSV format.
 */
process get_seed {

    publishDir "${params.outdir}/seed/", mode : 'copy'

    input:
    path (genome_locations)
    val(seed)

    output:
    path ('seed.txt')
    path ('seed_read_anonymisation.txt') optional true
    path ('seed_gsa_anonymisation.txt') optional true
    path ('seed_pooled_gsa_anonymisation.txt') optional true

    script:
    count_samples = params.number_of_samples
    if (params.anonymization){
        param_anonym = "-anonym_seed"
    } else {
        param_anonym = ""
    }
    """
    ${projectDir}/get_seed.py -seed ${seed} -count_samples ${count_samples} -file_genome_locations ${genome_locations} ${param_anonym}
    """
}

/*
* This process calculates the distribution of the genomes for one community.
* Takes: The file with the location to the drawn genomes.
*     
* Output: A file for each sample with the calculcated distributions.
*     
 */
process getCommunityDistribution {

    //publishDir "${params.outdir}/distributions/genome_distributions/", mode : 'copy'

    input:
    path(file_path_of_drawn_genome_location)
    val(seed)

    output:
    path 'distribution_*.txt'

    script:
    number_of_samples = params.number_of_samples
    mode = params.genome_mode
    log_mu = params.genome_log_mu
    log_sigma = params.genome_log_sigma
    gauss_mu = params.genome_gauss_mu
    gauss_sigma = params.genome_gauss_sigma
    """
    python ${projectDir}/get_community_distribution.py ${number_of_samples} ${file_path_of_drawn_genome_location} ${mode} ${log_mu} ${log_sigma} ${gauss_mu} ${gauss_sigma} False ${seed}
    """
}