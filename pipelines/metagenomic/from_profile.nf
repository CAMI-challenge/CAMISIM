scripts_dir = "${projectDir}/pipelines/metagenomic/scripts"
shared_scripts_dir = "${projectDir}/pipelines/shared/scripts"

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

        reference_genomes = params.containsKey('reference_genomes') ? params.reference_genomes : "${projectDir}/tools/assembly_summary_complete_genomes.txt"
        min_strains_per_otu = params.containsKey('min_strains_per_otu') ? params.min_strains_per_otu : 1
        max_strains_per_otu = params.containsKey('max_strains_per_otu') ? params.max_strains_per_otu : 2
        no_replace = params.containsKey('no_replace') ? params.no_replace : true
        fill_up = params.containsKey('fill_up') ? params.fill_up : false
        max_rank = params.containsKey('max_rank') ? params.max_rank : "family"
        prioritize_additional_genomes = params.containsKey('prioritize_additional_genomes') ? params.prioritize_additional_genomes : false
        additional_genomes_quality_file = params.containsKey('additional_genomes_quality_file') ? params.additional_genomes_quality_file : ""
        additional_genomes_max_contamination = params.containsKey('additional_genomes_max_contamination') ? params.additional_genomes_max_contamination : 5
        additional_genomes_min_completeness = params.containsKey('additional_genomes_min_completeness') ? params.additional_genomes_min_completeness : 95
        additional_genomes_max_num_contigs = params.containsKey('additional_genomes_max_num_contigs') ? params.additional_genomes_max_num_contigs : 200
        
        get_genomes(params.biom_profile, params.number_of_samples, reference_genomes, params.seed, params.gauss_mu, params.gauss_sigma,
            min_strains_per_otu, max_strains_per_otu, no_replace, fill_up,
            params.ncbi_taxdump_file, max_rank, prioritize_additional_genomes,
            additional_genomes_quality_file, additional_genomes_max_contamination,
            additional_genomes_min_completeness, additional_genomes_max_num_contigs)

        loc_ch = get_genomes.out[0]
        abundance_ch = get_genomes.out[1].flatten()
        meta_data_ch = get_genomes.out[2]

    emit: abundance_ch
    emit: loc_ch
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

    cpus 1

    conda "biom-format=2.1.17 ete4=4.3.0 python=3.11 h5py"

    input:
    path(biom_profile)
    val(number_of_samples)
    path(reference_genomes)
    val(seed)
    val(mu)
    val(sigma)
    val(min_strains)
    val(max_strains)
    val(no_replace)
    val(fill_up)
    val(ncbi_taxdump_file)
    val(max_rank)
    val(prioritize_additional_genomes)
    val(additional_genomes_quality_file)
    val(additional_genomes_max_contamination)
    val(additional_genomes_min_completeness)
    val(additional_genomes_max_num_contigs)

    output:
    path "genome_to_id.tsv"
    path "abundance_*.tsv"
    path "metadata.tsv"
    path "genome_filename_mapping.tsv"

    script:

    additional_references = "None"
    additional_genomes_quality = "None"

    if(params.containsKey('additional_references') && !params.additional_references.isEmpty()) {
        additional_references = params.additional_references
    }

    if(!additional_genomes_quality_file.isEmpty()) {
        additional_genomes_quality = additional_genomes_quality_file
    }

    """
    mkdir --parents ${params.outdir}/source_genomes/
    mkdir --parents ${params.outdir}/internal/
    mkdir --parents ${params.outdir}/distributions/

    python ${scripts_dir}/get_genomes.py "${biom_profile}" ${number_of_samples} "${reference_genomes}" ${seed} ${mu} ${sigma} ${min_strains} ${max_strains} False ${no_replace} ${fill_up} "${scripts_dir}/split_fasta.pl" "${params.outdir}/source_genomes/" "${additional_references}" "${additional_genomes_quality}" "${ncbi_taxdump_file}" "${max_rank}" ${prioritize_additional_genomes} ${additional_genomes_max_contamination} ${additional_genomes_min_completeness} ${additional_genomes_max_num_contigs}
    cp metadata.tsv ${params.outdir}/internal/metadata.tsv
    cp genome_to_id.tsv ${params.outdir}/internal/genome_to_id.tsv
    cp genome_to_id.tsv ${params.outdir}/internal/genome_locations.tsv
    cp genome_filename_mapping.tsv ${params.outdir}/internal/genome_filename_mapping.tsv
    cp abundance_*.tsv ${params.outdir}/distributions/
    """
}
