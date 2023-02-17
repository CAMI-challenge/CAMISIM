/*
 * Defining the module / subworkflow path, and include the elements
 */

// include read simulator here:
read_simulator_folder = "./read_simulators/"
// include read simulator nanosim3
include { read_simulator_nansoim3 } from "${read_simulator_folder}/read_simulator_nansoim3"
include { read_simulator_art } from "${read_simulator_folder}/read_simulator_art"

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
    take: genome_distribution_file_ch
    take: read_length_ch
    main:

        // this channel holds the genome location map (key = genome_id, value = absolute path to genome)
        genome_location_ch = genome_location_file_ch.splitCsv(sep:'\t')

        // create a channel, that holds the path to the genome distribution file by the sample id (key = genome_id, first value = distribution, second value = sample id)
        genome_distribution_ch = genome_distribution_file_ch.flatten().map { file -> tuple(file.baseName.split('_')[1], file) }.splitCsv(sep:'\t').map { a -> tuple(a[0], tuple(a[1][0], a[1][1])) }.groupTuple()
        
        // normalise the abundances of all genomes to 1 for every sample.
        normalised_distribution_ch = normalise_abundance(genome_distribution_ch)
        
        // create a channel, that holds the genome distribution by genome id and sample id (key = genome_id, first value = distribution, second value = sample id)
        genome_distribution_ch = normalised_distribution_ch.map { file -> tuple(file.baseName.split('_')[2], file) }.splitCsv(sep:'\t').map { a -> tuple(a[1][0], a[1][1],a[0]) }

        // combining of the channels results in new map: key = genome_id, first value = path to genome, second value = distribution, third value = sample_id
        genome_location_distribution_ch = genome_location_ch.combine(genome_distribution_ch, by: 0)

        if(params.type.equals("art")) {
            // simulate the reads with art 

            // create a channel that holds: key = sample_id, first value = distribution file, second value = file with all genome locations
            genome_distribution_location_ch = genome_distribution_file_ch.flatten().map { file -> tuple(file.baseName.split('_')[1], file) }.combine(genome_location_file_ch)

            // get the multiplication factor to calculate the fold coverage later (key = sample id, value = factor)
            factor_for_sample_id_ch = get_multiplication_factor(genome_distribution_location_ch)
            
            // join the two channel: key = genome_id, first value = path to genome, second value = distribution, third value = sample_id, fourth value = factor
            genome_location_distribution_factor_ch = genome_location_distribution_ch.map { tuple( it[3], *it ) }.combine(factor_for_sample_id_ch, by: 0 ).map { it[1..-1] }

            bam_files_channel = read_simulator_art(genome_location_distribution_factor_ch, read_length_ch)
        }
        if(params.type.equals("nanosim3")) {
            // simulate the reads with nanosim3
            bam_files_channel = read_simulator_nansoim3(genome_location_distribution_ch, read_length_ch)
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
}

/*
* This process generates a gold standard assembly for one genome.
* Takes:
*     A tuple with key = sample_id, first_value = genome_id, second value = a sorted bam file, third value = the reference genome (fasta).
* Output:
*     A Tuple with key = sample_id, value = path to fasta file with the gold standard assembly of the given genome.
 */
process generate_gold_standard_assembly {

    conda 'bioconda::samtools'

    input:
    tuple val(sample_id),val(genome_id), path(bam_file), path(reference_fasta_file)

    output:
    tuple val(sample_id), path(file_name)
    
    script:
    file_name = 'sample'.concat(sample_id.toString()).concat('_').concat(genome_id).concat('_gsa.fasta')
    """
    perl -- ${projectDir}/scripts/bamToGold.pl -st samtools -r ${reference_fasta_file} -b ${bam_file} -l 1 -c 1 >> ${file_name}
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
    path file_name
        
    script:
    file_name = 'sample'.concat(sample_id.toString()).concat('_gsa.fasta')
    """
    cat ${fasta_files} > ${file_name}
    mkdir --parents ${projectDir}/nextflow_out/sample_${sample_id}/contigs
    cp ${file_name} ${projectDir}/nextflow_out/sample_${sample_id}/contigs/anonymous_gsa.fasta
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

    conda 'bioconda::samtools'

    input:
    tuple val(sample_id), path(bam_files)

    output:
    path file_name

    script:
    file_name = 'sample_'.concat(sample_id.toString()).concat('.bam')
    compression = 5
    memory = 1
    """
    samtools merge -u - ${bam_files} | samtools sort -l ${compression} -m ${memory}G -o ${file_name} -O bam
    samtools index ${file_name}
    """
}

/*
* This process normalises the abundance for the given abundance map for one sample to 1.
* Takes:
*     A tuple with key = sample_id, value = a map with key = genome id and value = abundance.
* Output:
*     The path to file with the normalised abundances.
 */
process normalise_abundance {

    input:
    tuple val(sample_id), val(abundance_map)

    output:
    path file_name

    script:
    file_name = 'normalised_distributions_'.concat(sample_id).concat('.txt')

    double abundance_sum = 0.0

    abundance_map.each { 

        double abundance = Double.parseDouble((String) it[1])
        abundance_sum = abundance_sum + abundance
    }

    String output = ''

    abundance_map.eachWithIndex { item, index ->

        double abundance = Double.parseDouble((String) item[1])
        normalised_abundance = abundance / abundance_sum

        if(index!=0){
            output = output.concat('\n')
        }

        output = output.concat((String) item[0]).concat('\t').concat(Double.toString(normalised_abundance))
    }
    """
    echo "${output}" > ${file_name}
    """
}

/*
* This process calculates the multiplication factor for every every sample.
* This factor is needed in some read simulators for the calculation of the fold coverage. The factor has the same value for every genome of one sample.
* Takes:
*     A tuple with key = sample_id, first value = the file with all genome locations, second value = the file with the distributions of this sample.
* Output:
*     A tuple with key = sample_id, value = the calculated multiplication factor.
 */
process get_multiplication_factor {

    input:
    tuple val(sample_id), path(file_path_distribution), path(genome_locations)

    output:
    //tuple val(sample_id), path('factor.txt')
    tuple val(sample_id), stdout

    script:
    factor = 0
    fragment_size_mean = params.fragment_size_mean
    fragment_size_standard_deviation = params.fragment_size_sd
    total_size = (params.size*1000000000)

    """
    ${projectDir}/calculate_multiplication_factor.py ${fragment_size_mean} ${fragment_size_standard_deviation} ${total_size} ${genome_locations} ${file_path_distribution}
    """
}
