/** 
* This workflow simulates reads via nanosim3 and converts the resulting sam files into bam files.
* Takes:
*     A channel containing tuples with key = genome_id, first value = path to genome, second value = distribution, third value = sample_id, fourth value = seed.
* Emits: 
*     A channel containing tuples with key = sample_id, first value = genome id, second value = simulated bam file, third value = the reference fasta file.
**/
workflow metagenomesimulation_from_profile {

    main:
        
        get_genomes(params.biom_profile, params.number_of_samples, params.reference_genomes, params.seed, params.gauss_mu, params.gauss_sigma, 
            params.max_strains_per_otu, params.no_replace, params.fill_up)

        loc_ch = get_genomes.out[0] //.splitCsv(sep:'\t').map { a -> tuple(a[1].split("/")[-1].split(".fa")[0], a[0]) }
        //fa_ch = get_genomes.out[1].flatten().map { file -> tuple(file.baseName, file) }
        //abundance_ch = get_genomes.out[2].map { file -> tuple(file.splitCsv(sep:'\t')) }
        //abundance_ch = get_genomes.out[2].flatten().map { file -> tuple(file.baseName.split('_')[1], file) }.splitCsv(sep:'\t').map { a -> tuple(a[0], tuple(a[1][0], a[1][1])) }.groupTuple()
        abundance_ch = get_genomes.out[1].flatten()
        dump_ch = get_genomes.out[2]
        meta_data_ch = get_genomes.out[3]

        // combining of the channels results in new map: key = genome_id, value = path to genome
        //genome_location_ch = fa_ch.combine(loc_ch, by: 0).map { a -> tuple(a[2], a[1]) }

    //emit: abundance_ch
    //emit: genome_location_ch
    emit: abundance_ch
    emit: loc_ch
    emit: dump_ch
    emit: meta_data_ch
}

/*
* This process calculates the distribution of the genomes for one community.
* Takes: The file with the location to the drawn genomes.
*     
* Output: A file for each sample with the calculcated distributions.
*     
 */
process get_genomes {

    input:
    path(biom_profile)
    val(number_of_samples)
    path(reference_genomes)
    val(seed)
    val(mu)
    val(sigma)
    val(max_strains)
    val(no_replace)
    val(fill_up)

    output:
    path "genome_to_id.tsv"
    path "abundance_*.tsv"
    path "*.tar.gz"
    path "metadata.tsv"

    script:

    additional_references = "None"

    if(!params.additional_references.isEmpty()) {
        additional_references = params.additional_references
    }

    """
    mkdir --parents ${projectDir}/nextflow_out/internal/genomes/
    python ${projectDir}/get_genomes.py ${biom_profile} ${number_of_samples} ${reference_genomes} ${seed} ${mu} ${sigma} ${max_strains} False ${no_replace} ${fill_up} ${projectDir}/scripts/split_fasta.pl ${projectDir}/nextflow_out/internal/ ${additional_references}
    cp metadata.tsv ${projectDir}/nextflow_out/internal/metadata.tsv
    cp genome_to_id.tsv ${projectDir}/nextflow_out/internal/genome_to_id.tsv
    cp abundance*.tsv ${projectDir}/nextflow_out/internal/
    """
}