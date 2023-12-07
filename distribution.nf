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
    final_file_name = 'distribution_'.concat(sample_id).concat('.txt')

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
    mkdir --parents ${params.outdir}/distributions/
    cp ${file_name} ${params.outdir}/distributions/${final_file_name}
    """
}

/*
* This process normalises the given abundance with the size of the genome in number of bases.
* Takes:
*     A tuple with key = genome_id, first value = distribution, second value = sample_id, third value = size.
* Output:
*     A tuple with key = genome_id, first value = normalised distribution, second value = sample_id, third value = size.
 */
process normalise_abundance_to_size {

    input:
    tuple val(genome_id), val(distribution), val(sample_id), val(size)

    output:
    tuple val(genome_id), val(sample_id), val(normalised_distribution)

    script:
    normalised_distribution = size.toFloat() * distribution.toFloat()
    """
    """
}

/*
* This process calculates the size of the given genome in number of bases.
* Takes:
*     A tuple with key = genome_id, first value = genome location, second value = distribution, third value = sample_id.
* Output:
*     A tuple with key = genome_id, first value = distribution, second value = sample_id, third value = the calculated size.
 */
process count_bases {

    input:
    tuple val(genome_id), path(genome_location), val(sample_id), val(distribution)

    output:
    tuple val(genome_id), val(distribution), val(sample_id), stdout

    script:
    // bash variables need to be escaped
    // if not, nextflow considers it a nextflow variable
    """
    grep -v ">" ${genome_location} | wc | awk '{print \$3-\$1}'
    """
}