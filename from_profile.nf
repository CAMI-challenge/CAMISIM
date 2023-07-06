/** 
* This workflow designs a community based on the given biom profile. It updates the ncbi dump, downloads the genomes and calculates the abundances.
* Emits: 
*     A channel containing the abundance files for every sample.
*     A channel containing the genome location file channel.
*     A channel containing the downloaded ncbi dump.
*     A channel containing the metadata file.
**/
workflow metagenomesimulation_from_profile {

    main:
        
        get_genomes(params.biom_profile, params.number_of_samples, params.reference_genomes, params.seed, params.gauss_mu, params.gauss_sigma, 
            params.max_strains_per_otu, params.no_replace, params.fill_up)

        loc_ch = get_genomes.out[0]
        abundance_ch = get_genomes.out[1].flatten()
        dump_ch = get_genomes.out[2]
        meta_data_ch = get_genomes.out[3]

    emit: abundance_ch
    emit: loc_ch
    emit: dump_ch
    emit: meta_data_ch
}

/*
* This process designs a community based on the given biom profile.
*     
* Output: A file holding the genome id to path to genome file.
*         An abundance file for every sample.
*         The downloaded zipped ncbi dump.
*         The metadata file.
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