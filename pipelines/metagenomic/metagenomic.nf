#!/usr/bin/env nextflow

nextflow.enable.dsl=2

scripts_dir = "${projectDir}/pipelines/metagenomic/scripts"
shared_scripts_dir = "${projectDir}/pipelines/shared/scripts"

/*
 * Defining the module / subworkflow path, and include the elements
 */

// include sample wise simulation
include { sample_wise_simulation } from "${projectDir}/pipelines/metagenomic/sample_wise_simulation"

// include from profile metagenome simulation
include { metagenomesimulation_from_profile } from "${projectDir}/pipelines/metagenomic/from_profile"

// include anonymization
include { anonymization } from "${projectDir}/pipelines/shared/anonymization"
include { anonymize_merged_gsa } from "${projectDir}/pipelines/shared/anonymization"
include { anonymize_hybrid_gsa } from "${projectDir}/pipelines/shared/anonymization"
include { anonymize_hybrid_pooled_gsa } from "${projectDir}/pipelines/shared/anonymization"

// include binning
include { binning } from "${projectDir}/pipelines/shared/binning"

/*
 * Resolves which pipeline step(s) to run.
 *
 * "step" is the primary control and must be one of:
 *   "community_design" - only design the community, then stop.
 *   "reads_simulate"   - skip the community design and simulate reads from pre-generated files.
 *   "all"              - community design followed by read simulation.
 *
 * When "step" is empty/unset the value is derived from the legacy
 * biom_profile / distribution_files / just_community_design parameters, so existing
 * configs keep behaving exactly as before.
 */
def resolve_step() {

    def valid_steps = ['community_design', 'reads_simulate', 'all']

    def requested = params.containsKey('step') ? params.step?.toString()?.trim() : ''

    def step
    if (requested) {
        step = requested
    } else if (params.biom_profile.isEmpty() && !params.distribution_files.isEmpty()) {
        // legacy: no profile + provided distribution files => resume into read simulation
        step = 'reads_simulate'
    } else if (params.containsKey('just_community_design') && params.just_community_design) {
        // legacy: stop after the community design
        step = 'community_design'
    } else {
        step = 'all'
    }

    if (!valid_steps.contains(step)) {
        throw new IllegalArgumentException("Invalid 'step' value: '${step}'. Must be one of ${valid_steps}.")
    }

    // "reads_simulate" consumes the files produced by a previous community design run.
    if (step == 'reads_simulate') {
        ['distribution_files', 'genome_locations_file', 'metadata_file'].each { key ->
            def value = params.containsKey(key) ? params[key]?.toString()?.trim() : ''
            if (!value) {
                throw new IllegalArgumentException("step='reads_simulate' requires '${key}' to point to the pre-generated file(s) from the community design step.")
            }
        }
    }

    return step
}

/*
 * This is the main workflow and starting point of this nextflow pipeline.
 */
workflow metagenomic {

    // Resolve which step(s) to run (see resolve_step for the mapping of the
    // legacy biom_profile / distribution_files / just_community_design parameters).
    step = resolve_step()

    // contamination_asvs: reads are still simulated for these ASVs, but they are excluded
    // from all ground-truth files. Matching is on the genome_id ('<asv>' or '<asv>.<strain>'),
    // so it applies to biom-designed communities AND to step='reads_simulate' runs whose
    // pre-generated genome_locations/metadata/abundance files use that ASV naming.
    if (params.containsKey('contamination_asvs') && params.contamination_asvs) {
        if (!(params.contamination_asvs instanceof List)) {
            error "params.contamination_asvs must be a list of ASV ids, e.g. ['hASV1'], but got '${params.contamination_asvs}' (${params.contamination_asvs.getClass().name})."
        }
        log.info "contamination_asvs=${params.contamination_asvs}: reads for these ASVs (and their '<asv>.<strain>' genomes) will be simulated but excluded from all gold standards, mappings, per-genome coverage, and the (renormalized) taxonomic profile."
    }

    if(params.seed != null) {
            seed = params.seed
        } else {
            seed = get_random_seed()
        }

    // If no NCBI taxonomy database is given it will be downloaded.
    if(params.ncbi_taxdump_file.isEmpty()) {
        ncbi_taxdump_file_ch = download_NCBI_taxdump()
    } else {
        // this channel holds the ncbi tax dump
        ncbi_taxdump_file_ch = Channel.fromPath(params.ncbi_taxdump_file)
    }

    // ============================ COMMUNITY DESIGN ============================
    // Community design runs for steps "community_design" and "all".

    // Community design from a given biom profile.
    if((step == 'community_design' || step == 'all') && !params.biom_profile.isEmpty()) {

        // start the simulation
        metagenomesimulation_from_profile()

        // get the output channel
        genome_distribution_file_ch = metagenomesimulation_from_profile.out[0]
        genome_location_file_ch = metagenomesimulation_from_profile.out[1]
        metadata_ch = metagenomesimulation_from_profile.out[2]
    }

    // Community design from a given genome list (optionally with strain simulation).
    if((step == 'community_design' || step == 'all') && params.biom_profile.isEmpty()) {

        // this channel holds the file with the specified locations of the genomes
        genome_location_file_ch = Channel.fromPath(params.genome_locations_file)

        metadata_ch = Channel.fromPath(params.metadata_file)

            // if there are more genomes requested than inputted, simulate strains
            if (! (params.genomes_total==params.genomes_real)){

                prepare_strain_simulation(params.genomes_total, params.genomes_real, seed, metadata_ch, params.max_strains_per_otu, genome_location_file_ch)

                // strain simulation is either with gff files or without possible
                if(params.id_to_gff_file.isEmpty()) {

                    // this channel holds the genome location map (key = genome_id, value = absolute path to genome)
                    strain_simulation_ch = prepare_strain_simulation.out
                        .splitCsv(sep:'\t') // get genome id and relatvie path from genome location file
                        .map { genome_id, path, amount, seed ->
                            def abs_path
                            if (new File(path).isAbsolute()) { // if the path is an absolute path return it as is
                                abs_path = path
                            } else { // else expand relative paths to absolute paths and send to genome_location_ch
                                abs_path = file("${projectDir}/${path}").toAbsolutePath().toString()
                            }
                            return [genome_id, abs_path, amount, seed]
                            }.filter { genome_id, abs_path, amount, seed -> amount.toInteger() > 1 }
                            .combine(metadata_ch.splitCsv(sep:'\t', header: false), by:0)

                    strain_simulation_without_gff(strain_simulation_ch)

                    added_genome_location_ch = strain_simulation_without_gff.out[0]
                    added_metadata_ch = strain_simulation_without_gff.out[1]

                } else {
                    // this channel holds the genome location map (key = genome_id, value = absolute path to genome)
                    strain_simulation_ch = prepare_strain_simulation.out
                        .splitCsv(sep:'\t') // get genome id and relatvie path from genome location file
                        .map { genome_id, path, amount, seed, gff ->
                            def abs_path
                            def abs_path_2
                            if (new File(path).isAbsolute()) { // if the path is an absolute path return it as is
                                abs_path = path
                            } else { // else expand relative paths to absolute paths and send to genome_location_ch
                                abs_path = file("${projectDir}/${path}").toAbsolutePath().toString()
                            }

                            if (new File(gff).isAbsolute()) { // if the path is an absolute path return it as is
                                abs_path_2 = gff
                            } else { // else expand relative paths to absolute paths and send to genome_location_ch
                                abs_path_2 = file("${projectDir}/${gff}").toAbsolutePath().toString()
                            }

                        
                            return [genome_id, abs_path, amount, seed, abs_path_2]
                            }.filter { genome_id, abs_path, amount, seed, abs_path_2 -> amount.toInteger() > 1 }
                            .combine(metadata_ch.splitCsv(sep:'\t', header: false), by:0)

                    strain_simulation_with_gff(strain_simulation_ch)

                    added_genome_location_ch = strain_simulation_with_gff.out[0]
                    added_metadata_ch = strain_simulation_with_gff.out[1]
                }
                
                // merge the metadata files together
                merge_metadata_files(added_genome_location_ch.collect(), added_metadata_ch.collect(), genome_location_file_ch, metadata_ch)

                // this channel holds the file with the specified locations of the genomes
                genome_location_file_ch = merge_metadata_files.out[0]

                metadata_ch = merge_metadata_files.out[1]
            }

            // calculate the genome distributions for each sample for one community
            genome_distribution_file_ch = getCommunityDistribution(genome_location_file_ch, seed).flatten()
    }

    // Stop the pipeline if only the community design step is requested.
    if(step == 'community_design'){
        println "Simulation stopping after community design steps."
        return
    }

    // ============================ READ SIMULATION INPUTS ============================
    // For "reads_simulate" the community design is skipped and the files produced by a previous
    // community design run are used directly (their presence is validated in resolve_step).
    if(step == 'reads_simulate') {

        // this channel holds the files with the specified distributions for every sample
        genome_distribution_file_ch = Channel.fromPath(params.distribution_files)

        // this channel holds the file with the specified locations of the genomes
        genome_location_file_ch = Channel.fromPath(params.genome_locations_file)

        metadata_ch = Channel.fromPath(params.metadata_file)
    }

    // build ncbi taxonomy from given tax dump
    number_of_samples_ch = Channel.from(params.number_of_samples)
    buildTaxonomy(number_of_samples_ch.concat(ncbi_taxdump_file_ch.concat(genome_distribution_file_ch)).toList().map { it -> [ it[0], it[1], it[2..-1] ] }, metadata_ch)

    def types_list = (params.type instanceof List) ? params.type : [params.type]

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

    // make all sequences from input genomes (also strain simulated ones) unique and move them to an output location
    genome_location_file_ch = cleanup_and_filter_sequences(genome_location_file_ch, genome_location_ch.collect(flat: false))
    
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

    // random seed generation
    get_seed(genome_location_file_ch, seed)
    // get the text file with the seeds needed for the read simulation
    seed_file_read_simulation_ch = get_seed.out[0]

    // simulate reads sample wise
    sample_wise_simulation(genome_location_ch, genome_location_file_ch, genome_distribution_file_ch, seed_file_read_simulation_ch)

    // out[0]: tuple(sim_type, sample_id, bam)      – per-type merged BAM per sample
    // out[1]: tuple(sim_type, sample_id, gsa_fasta) – per-type GSA per sample
    // out[2]: tuple(sim_type, sample_id, …reads…)   – reads per type per sample (grouped)
    // out[3]: tuple(sim_type, sample_id, [bam_paths]) – per-type BAM list per sample
    // out[4]: tuple(sample_id, bam)                 – all-types combined BAM per sample (for hybrid GSA)
    merged_bam_per_type_per_sample_ch = sample_wise_simulation.out[0]
    gsa_per_type_per_sample_ch        = sample_wise_simulation.out[1]
    reads_per_type_ch                 = sample_wise_simulation.out[2]
    bam_files_by_sample_per_type_ch   = sample_wise_simulation.out[3]
    combined_bam_per_sample_ch        = sample_wise_simulation.out[4]

    // ---- Resolve hybrid types and BAMs ----
    // hybrid = true  → use all simulated types
    // hybrid = [..] → intersect with simulated types; only types present in both lists are hybridised
    // hybrid = false (or absent) → no hybrid
    def hybrid_types = []
    if (params.containsKey('hybrid') && params.hybrid) {
        if (params.hybrid instanceof Boolean) {
            hybrid_types = types_list   // true: all types
        } else if (params.hybrid instanceof List) {
            hybrid_types = params.hybrid.intersect(types_list)  // list: only requested+simulated
        }
    }
    hybrid_bam_ch = Channel.empty()
    if (hybrid_types.size() > 1) {
        if (hybrid_types == types_list) {
            // All types selected: reuse the already-merged combined BAM directly
            hybrid_bam_ch = combined_bam_per_sample_ch
        } else {
            // Subset of types: filter per-type BAMs to the hybrid set, then merge per sample
            filtered_bam_ch = merged_bam_per_type_per_sample_ch.filter { sim_type, sample_id, bam -> hybrid_types.contains(sim_type) }.map { sim_type, sample_id, bam -> tuple(sample_id, bam) }.groupTuple()
            hybrid_bam_ch = merge_bam_files_hybrid(filtered_bam_ch)
        }
    }

    // extract file paths from the tuples to create the reference_fasta_files_ch
    reference_fasta_files_ch = genome_location_ch.map { a -> a[1] }

    // ---- pooled GSA: collect all per-type merged BAMs across samples and types ----
    if (params.pooled_gsa instanceof Boolean && params.pooled_gsa) {
        merged_bam_file = merge_bam_files(merged_bam_per_type_per_sample_ch.map { sim_type, sid, bam -> bam }.collect())
    } else if (params.pooled_gsa instanceof List) {
        merged_bam_file = merge_bam_files(merged_bam_per_type_per_sample_ch
            .filter { sim_type, sid, bam -> params.pooled_gsa*.toString().contains(sid.toString()) }.map { sim_type, sid, bam -> bam }.collect())
    }

    // Initialize a single channel to accumulate all BAMs requiring coverage calculations
    coverages_input_ch = Channel.empty()

    if (params.gsa) {
        coverages_input_ch = coverages_input_ch.mix(merged_bam_per_type_per_sample_ch.map { sim_type, sample_id, bam ->
            tuple(bam, "sample_${sample_id}_${sim_type}", "${params.outdir}/sample_${sample_id}/contigs/${sim_type}") }.combine(genome_location_file_ch.first()))
    }

    // Generate merged GSAs for custom sample combinations
    if (params.containsKey('merged_gsa_combinations') && params.merged_gsa_combinations instanceof List && params.merged_gsa_combinations.size() > 0) {
        // Fail fast on malformed combinations. Without this, an unknown sample id
        // is silently dropped during the BAM merge: if some ids in an entry are
        // valid the merged GSA is built from only those (yet still named after the
        // full requested list), and only an entry whose ids are *all* invalid would
        // incidentally fail later in samtools merge. Validate up front instead.
        def valid_sample_ids = (0..<params.number_of_samples).collect { it.toString() } as Set
        params.merged_gsa_combinations.eachWithIndex { combination, idx ->
            if (!(combination instanceof List) || combination.isEmpty()) {
                error "merged_gsa_combinations entry ${idx} (${combination}) must be a non-empty list of sample ids."
            }
            def invalid_sample_ids = combination.collect { it.toString() }.findAll { !valid_sample_ids.contains(it) } as Set
            if (invalid_sample_ids) {
                error "merged_gsa_combinations entry ${idx} (${combination}) references invalid sample id(s) ${invalid_sample_ids}. Valid sample ids are 0..${params.number_of_samples - 1}."
            }
        }

        combinations_ch = Channel.fromList(params.merged_gsa_combinations.withIndex()).map { combination, idx -> tuple(idx, combination*.toString()) }

        // Mix per-type and hybrid-type BAMs into a single typed channel
        typed_bams_per_sample_ch = merged_bam_per_type_per_sample_ch.map { sim_type, sample_id, bam -> tuple(sim_type, sample_id, bam) }
        if (hybrid_types.size() > 1) {
            def hybrid_label = "hybrid_" + hybrid_types.join('_')
            typed_bams_per_sample_ch = typed_bams_per_sample_ch.mix(
                hybrid_bam_ch.map { sample_id, bam -> tuple(hybrid_label, sample_id, bam) }
            )
        }

        // Group bams to merge by combination and simulation label
        bams_to_merge_ch = combinations_ch
            .combine(typed_bams_per_sample_ch)
            .filter { combination_id, sample_ids, sim_label, sample_id, bam -> sample_ids.contains(sample_id.toString()) }
            .map { combination_id, sample_ids, sim_label, sample_id, bam -> tuple(tuple(sim_label, combination_id, sample_ids), bam) }
            .groupTuple()
            .map { key, bams -> tuple(key[0], key[1], key[2], bams) }

        // Use combined typed BAMs for the cross-sample merged GSAs
        merged_bam_per_combination = merge_bam_files_by_combination_typed(bams_to_merge_ch)

        merged_gsa_ch = generate_merged_gold_standard_assembly(merged_bam_per_combination.combine(reference_fasta_files_ch.collect().map { [it] }))

        // ---- Coverage for merged-combo GSAs ----
        coverages_input_ch = coverages_input_ch.mix(merged_bam_per_combination.map { sim_label, combination_id, sample_ids, bam ->
            def samples_str = sample_ids.join('_')
            tuple(bam, "merged_samples_${samples_str}_${sim_label}", "${params.outdir}/merged_gsa") }.combine(genome_location_file_ch.first()))
    }

    if (params.pooled_gsa) {

        // ---- Per-type pooled GSA ----
        // For each simulated type, collect its per-sample BAMs and generate a per-type pooled GSA.
        // Output files: pooled_gsa/gsa_pooled_<type>.fasta.gz
        if (params.pooled_gsa instanceof Boolean && params.pooled_gsa) {
            per_type_pooled_bam_ch = merged_bam_per_type_per_sample_ch.map { sim_type, sid, bam -> tuple(sim_type, bam) }.groupTuple()
        } else if (params.pooled_gsa instanceof List) {
            per_type_pooled_bam_ch = merged_bam_per_type_per_sample_ch.filter { sim_type, sid, bam -> params.pooled_gsa*.toString().contains(sid.toString()) }.map { sim_type, sid, bam -> tuple(sim_type, bam) }.groupTuple()
        }
        per_type_pooled_merged_bam_ch = merge_bam_files_per_type_pooled(per_type_pooled_bam_ch)
        generate_pooled_gold_standard_assembly_per_type(per_type_pooled_merged_bam_ch.combine(reference_fasta_files_ch.collect().map { [it] }))

        // ---- Coverage for per-type pooled GSAs ----
        coverages_input_ch = coverages_input_ch.mix(per_type_pooled_merged_bam_ch.map { sim_type, bam ->
            tuple(bam, "pooled_${sim_type}", "${params.outdir}/pooled_gsa") }.combine(genome_location_file_ch.first()))

        // ---- Hybrid pooled GSA ----
        // Only combine the types listed in `hybrid`; save to hybrid_pooled_gsa/ directory.
        if (hybrid_types.size() > 1) {
            if (hybrid_types == types_list) {
                // All types: use the already-merged combined BAM across all samples
                hybrid_pooled_bam_ch = merged_bam_file
            } else {
                // Subset of types: filter pooled per-type BAMs to only the hybrid types, then merge
                hybrid_type_bams_ch = per_type_pooled_merged_bam_ch.filter { sim_type, bam -> hybrid_types.contains(sim_type) }.map { sim_type, bam -> bam }.collect()
                hybrid_pooled_bam_ch = merge_bam_files_hybrid_pooled(hybrid_type_bams_ch)
            }
            generate_hybrid_pooled_gold_standard_assembly(hybrid_pooled_bam_ch.combine(reference_fasta_files_ch.collect().map { [it] }), hybrid_types)

            // ---- Coverage for hybrid pooled GSA ----
            coverages_input_ch = coverages_input_ch.mix(hybrid_pooled_bam_ch.map { bam ->
                def types_suffix = (hybrid_types instanceof List ? hybrid_types : [hybrid_types]).join('_')
                tuple(bam, "pooled_hybrid_${types_suffix}", "${params.outdir}/hybrid_pooled_gsa") }.combine(genome_location_file_ch.first()))
        }

        // ---- Hybrid GSA per sample ----
        if (hybrid_types.size() > 1) {
            generate_hybrid_gold_standard_assembly(hybrid_bam_ch.combine(reference_fasta_files_ch.collect().map { [it] }), hybrid_types)

            // ---- Coverage for hybrid per-sample GSAs ----
            coverages_input_ch = coverages_input_ch.mix(hybrid_bam_ch.map { sample_id, bam ->
                def types_suffix = (hybrid_types instanceof List ? hybrid_types : [hybrid_types]).join('_')
                tuple(bam, "sample_${sample_id}_hybrid_${types_suffix}", "${params.outdir}/sample_${sample_id}/hybrid_gsa") }.combine(genome_location_file_ch.first()))
        }

        // if requested, anonymize reads, gsa and pooled gsa
        if(params.anonymization) {
            // Pass the first per-type pooled GSA output as the "pooled_gsa" channel for anonymization.
            // (Anonymization of all per-type pooled GSAs and the hybrid pooled GSA would require
            //  further refactoring of the anonymization workflow.)
            anonymization(reads_per_type_ch, get_seed.out[1], get_seed.out[2], get_seed.out[3], gsa_per_type_per_sample_ch, bam_files_by_sample_per_type_ch, generate_pooled_gold_standard_assembly_per_type.out, per_type_pooled_merged_bam_ch, genome_location_file_ch, metadata_ch)
            if (params.containsKey('merged_gsa_combinations') && params.merged_gsa_combinations instanceof List && params.merged_gsa_combinations.size() > 0) {
                anonymize_merged_gsa(merged_gsa_ch, merged_bam_per_combination, get_seed.out[4], genome_location_file_ch, metadata_ch)
            }
            if (hybrid_types.size() > 1) {
                anonymize_hybrid_gsa(generate_hybrid_gold_standard_assembly.out, hybrid_bam_ch, get_seed.out[5], genome_location_file_ch, metadata_ch, hybrid_types)
                anonymize_hybrid_pooled_gsa(generate_hybrid_pooled_gold_standard_assembly.out, hybrid_pooled_bam_ch, get_seed.out[3], genome_location_file_ch, metadata_ch, hybrid_types)
            }
        } else { // if no anonymization is requested, create binning gold standard
            binning(gsa_per_type_per_sample_ch, bam_files_by_sample_per_type_ch, generate_pooled_gold_standard_assembly_per_type.out.map { sim_type, f -> f }.first(), per_type_pooled_merged_bam_ch.map { sim_type, bam -> bam }.first(), genome_location_file_ch, metadata_ch)
        }
    }

    // Trigger all coverage computations at once through a single channel
    compute_coverage(coverages_input_ch)
}

/*
* This process downloads the NCBI taxonomy database.
*
*/
process download_NCBI_taxdump {

    conda 'conda-forge::ete3'

    output:
    path "taxdump.tar.gz"

    script:
    """
    #!/usr/bin/env python
    import os
    from ete3 import NCBITaxa

    # Initialize NCBITaxa
    ncbi = NCBITaxa()

    # Update taxonomy database
    ncbi.update_taxonomy_database()

    # Create the output directory if it does not exist
    output_dir = "${params.outdir}/internal/genomes/"
    os.makedirs(output_dir, exist_ok=True)

    # Copy the downloaded taxdump to the output directory
    taxdump_file = "./*.tar.gz"
    os.system(f"cp {taxdump_file} {output_dir}")
    """
}

/*
* This process merges all given bam files specified in the pooled_gsa parameter.
* Takes:
*     A list with the paths to all bam files that should be merged if the condition is fulfilled.
* Output:
*     The path to the merged bam file.
 */
process merge_bam_files {

    conda 'bioconda::samtools'

    input:
    path bam_files

    output:
    path file_name

    script:
    file_name = 'merged.bam'
    compression = 5
    memory = 1
    threads_for_sort = Math.max(1, ((task.cpus ?: 1) as int))

    bam_to_merge = ''

    bam_files.each {

        bam_file_name = (String) it
        sample_id = bam_file_name.split('_')[1][0].toInteger()

        //if(sample_id in params.pooled_gsa){
        bam_to_merge = bam_to_merge.concat(' ').concat(bam_file_name)
        //}
    }
    """
    samtools merge -u - ${bam_to_merge} | samtools sort -@ ${threads_for_sort} -l ${compression} -m ${memory}G -o ${file_name} -O bam
    """
}

/*
* This process merges BAM files from a selected subset of types for one sample,
* used as input for the hybrid gold standard assembly.
*/
process merge_bam_files_hybrid {

    conda 'bioconda::samtools'

    input:
    tuple val(sample_id), path(bam_files)

    output:
    tuple val(sample_id), path(file_name)

    script:
    file_name = "sample_${sample_id}_hybrid.bam"
    compression = 5
    memory = 1
    threads_for_sort = Math.max(1, ((task.cpus ?: 1) as int))
    """
    samtools merge -u - ${bam_files} | samtools sort -@ ${threads_for_sort} -l ${compression} -m ${memory}G -o ${file_name} -O bam
    """
}

/*
* Merges per-sample BAMs of a single sequencing type across all pooled samples,
* producing one merged BAM per type for the per-type pooled GSA.
*/
process merge_bam_files_per_type_pooled {

    conda 'bioconda::samtools'

    input:
    tuple val(sim_type), path(bam_files)

    output:
    tuple val(sim_type), path(file_name)

    script:
    file_name = "pooled_${sim_type}.bam"
    compression = 5
    memory = 1
    threads_for_sort = Math.max(1, ((task.cpus ?: 1) as int))
    """
    samtools merge -u - ${bam_files} | samtools sort -@ ${threads_for_sort} -l ${compression} -m ${memory}G -o ${file_name} -O bam
    """
}

/*
* Merges the per-type pooled BAMs of the hybrid types into a single BAM,
* used as input for the hybrid pooled gold standard assembly.
*/
process merge_bam_files_hybrid_pooled {

    conda 'bioconda::samtools'

    input:
    path bam_files

    output:
    path file_name

    script:
    file_name = "pooled_hybrid.bam"
    compression = 5
    memory = 1
    threads_for_sort = Math.max(1, ((task.cpus ?: 1) as int))
    """
    samtools merge -u - ${bam_files} | samtools sort -@ ${threads_for_sort} -l ${compression} -m ${memory}G -o ${file_name} -O bam
    """
}

/*
* This process merges BAM files for a specific combination of samples.
* Takes:
*     combination_id: An identifier for the combination (index)
*     sample_ids: List of sample IDs to include in this combination
*     all_bam_tuples: All BAM file tuples (sample_id, bam_path)
* Output:
*     A tuple with combination_id and the path to the merged BAM file.
*/
process merge_bam_files_by_combination {

    conda 'bioconda::samtools'

    input:
    tuple val(combination_id), val(sample_ids)
    val all_bam_tuples

    output:
    tuple val(combination_id), val(sample_ids), path(file_name)

    script:
    file_name = "merged_combination_${combination_id}.bam"
    compression = 5
    memory = 1
    threads_for_sort = Math.max(1, ((task.cpus ?: 1) as int))

    // Filter BAM files that belong to the specified sample IDs
    bam_to_merge = all_bam_tuples
        .findAll { sample_id, bam_path -> sample_ids.contains(sample_id.toString()) }
        .sort { a, b -> sample_ids.indexOf(a[0].toString()) <=> sample_ids.indexOf(b[0].toString()) }
        .collect { sample_id, bam_path -> bam_path }
        .join(' ')
    """
    samtools merge -u - ${bam_to_merge} | samtools sort -@ ${threads_for_sort} -l ${compression} -m ${memory}G -o ${file_name} -O bam
    """
}

/*
* This process merges BAM files for a specific combination of samples, grouped by simulation type/label.
*/
process merge_bam_files_by_combination_typed {

    conda 'bioconda::samtools'

    input:
    tuple val(sim_label), val(combination_id), val(sample_ids), path(bam_files)

    output:
    tuple val(sim_label), val(combination_id), val(sample_ids), path(file_name)

    script:
    file_name = "merged_combination_${combination_id}_${sim_label}.bam"
    compression = 5
    memory = 1
    threads_for_sort = Math.max(1, ((task.cpus ?: 1) as int))
    """
    samtools merge -u - ${bam_files} | samtools sort -@ ${threads_for_sort} -l ${compression} -m ${memory}G -o ${file_name} -O bam
    """
}

/*
* This process generates a merged gold standard assembly for a custom sample combination.
* Takes:
*     A tuple with combination_id, sample_ids list, merged BAM file, and reference FASTA files.
* Output:
*     The path to the FASTA file with the merged gold standard assembly.
*/
process generate_merged_gold_standard_assembly {

    conda 'bioconda::samtools conda-forge::pigz'

    input:
    tuple val(sim_label), val(combination_id), val(sample_ids), path(bam_file), path(reference_fasta_files)

    output:
    tuple val(sim_label), val(combination_id), val(sample_ids), path(file_name)

    script:
    // Create a descriptive filename with the sample IDs and labels
    samples_str = sample_ids.join('_')
    file_name = "gsa_merged_samples_${samples_str}_${sim_label}.fasta"
    threads = Math.max(1, ((task.cpus ?: 1) as int))
    """
    cat ${reference_fasta_files} > reference.fasta
    samtools faidx reference.fasta
    python ${shared_scripts_dir}/bamToGold.py -st samtools -r reference.fasta -b ${bam_file} -l 1 -c 1 >> ${file_name}
    mkdir --parents ${params.outdir}/merged_gsa
    pigz -p ${threads} -k ${file_name}
    cp ${file_name}.gz ${params.outdir}/merged_gsa/
    """
}

/*
* This process generates the pooled gold standard assembly for one specific
* sequencing type, across all pooled samples.
* Output file: pooled_gsa/gsa_pooled_<type>.fasta.gz
*/
process generate_pooled_gold_standard_assembly_per_type {

    conda 'bioconda::samtools conda-forge::pigz'

    input:
    tuple val(sim_type), path(bam_file), path(reference_fasta_files)

    output:
    tuple val(sim_type), path(file_name)

    script:
    file_name = "gsa_pooled_${sim_type}.fasta"
    threads = Math.max(1, ((task.cpus ?: 1) as int))
    """
    cat ${reference_fasta_files} > reference.fasta
    samtools faidx reference.fasta
    python ${shared_scripts_dir}/bamToGold.py -st samtools -r reference.fasta -b ${bam_file} -l 1 -c 1 >> ${file_name}
    mkdir --parents ${params.outdir}/pooled_gsa
    pigz -p ${threads} -k ${file_name}
    cp ${file_name}.gz ${params.outdir}/pooled_gsa/
    """
}

/*
* This process generates the hybrid pooled gold standard assembly combining
* only the types listed in `hybrid`, across all pooled samples.
* Output file: hybrid_pooled_gsa/gsa_pooled_hybrid_<types>.fasta.gz
*/
process generate_hybrid_pooled_gold_standard_assembly {

    conda 'bioconda::samtools conda-forge::pigz'

    input:
    tuple path(bam_file), path(reference_fasta_files)
    val hybrid_types

    output:
    path file_name

    script:
    types_suffix = (hybrid_types instanceof List ? hybrid_types : [hybrid_types]).join('_')
    file_name = "gsa_pooled_hybrid_${types_suffix}.fasta"
    threads = Math.max(1, ((task.cpus ?: 1) as int))
    """
    cat ${reference_fasta_files} > reference.fasta
    samtools faidx reference.fasta
    python ${shared_scripts_dir}/bamToGold.py -st samtools -r reference.fasta -b ${bam_file} -l 1 -c 1 >> ${file_name}
    mkdir --parents ${params.outdir}/hybrid_pooled_gsa
    pigz -p ${threads} -k ${file_name}
    cp ${file_name}.gz ${params.outdir}/hybrid_pooled_gsa/
    """
}

/*
* This process generates a hybrid gold standard assembly for one sample,
* combining reads from the selected hybrid sequencing types via their merged BAM.
* Input:
*     tuple(sample_id, combined_bam, [reference_fasta_files])
*     val hybrid_types_str  – underscore-joined sorted type names, e.g. "iss_nanosim3"
* Output:
*     tuple(sample_id, hybrid_gsa_fasta)
*/
process generate_hybrid_gold_standard_assembly {

    conda 'bioconda::samtools conda-forge::pigz'

    input:
    tuple val(sample_id), path(bam_file), path(reference_fasta_files)
    val hybrid_types

    output:
    tuple val(sample_id), path(file_name)

    script:
    types_suffix = (hybrid_types instanceof List ? hybrid_types : [hybrid_types]).join('_')
    file_name = "sample${sample_id}_hybrid_gsa_${types_suffix}.fasta"
    threads = Math.max(1, ((task.cpus ?: 1) as int))
    """
    cat ${reference_fasta_files} > reference.fasta
    samtools faidx reference.fasta
    python ${shared_scripts_dir}/bamToGold.py -st samtools -r reference.fasta -b ${bam_file} -l 1 -c 1 >> ${file_name}
    mkdir --parents ${params.outdir}/sample_${sample_id}/hybrid_gsa
    pigz -p ${threads} -k ${file_name}
    cp ${file_name}.gz ${params.outdir}/sample_${sample_id}/hybrid_gsa/
    """
}

/*
* This process builds the taxonomy profile for every sample with the given distribution and the ncbi tax dump. The generated profiles will be
* copied to the out directory.
* Takes:
*     The tuple with first_value = zipped ncbi tax dump to build the profile from  and second value = the number of samples.
* Output:
*     The paths to all taxonomic profiles.
 */
process buildTaxonomy {
    
    conda 'python'

    input:
    tuple val(number_of_samples), path(dmp), path(distribution_files)
    path(metadata_ch)

    output:
    path 'taxonomic_profile_*.txt'

    script:
    index_number_of_samples = number_of_samples - 1
    // Contamination ASVs are excluded from the taxonomic profile (and the remaining
    // taxa are renormalized). Passed as a comma-separated string; "" disables it.
    contamination_csv = (params.containsKey('contamination_asvs') && params.contamination_asvs) ? params.contamination_asvs.join(',') : ''
    """
    tar -xf ${dmp}
    [ -f **/names.dmp ] && mv **/names.dmp ./names.dmp
    [ -f **/merged.dmp ] && mv **/merged.dmp ./merged.dmp
    [ -f **/nodes.dmp ] && mv **/nodes.dmp ./nodes.dmp
    ${scripts_dir}/build_ncbi_taxonomy.py names.dmp merged.dmp nodes.dmp ${number_of_samples} ${metadata_ch} "${contamination_csv}" ${distribution_files}
    mkdir --parents ${params.outdir}
    cp taxonomic_profile_*.txt ${params.outdir}
    """
}

process get_random_seed {

    conda 'python'

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
* This process prepares the strain simulation by calculating the genome amounts.
*/
process prepare_strain_simulation {

    conda "conda-forge::python=3.11.5 conda-forge::numpy"

    input:
    val genomes_total
    val genomes_real
    val seed
    path metadata
    val max_strains_per_otu
    path id_to_genome_file

    output:
    path "genome_id_to_file_amount_gff.tsv"

    script:
    if(params.id_to_gff_file.isEmpty()){
        gff = ""
    } else {
        gff += "--id_to_gff_file ${id_to_gff_file}"
    }
    """
    python ${scripts_dir}/prepare_strain_simulation.py -genomes_total ${genomes_total} -genomes_real ${genomes_real} -seed ${seed} -metadata ${metadata} -max_strains_per_otu ${max_strains_per_otu} -id_to_genome_file ${id_to_genome_file} ${gff}
    """
}

/*
* This process simulates strains using sgEvolver by using an empty gff file.
*/
process strain_simulation_without_gff {

    conda "bioconda::perl-bioperl conda-forge::biopython=1.83 conda-forge::python=3.11.5"

    input:
    tuple val(genome_id), path(fasta), val(amount), val(seed), val(OTU), val(NCBI_ID), val(novelty_category)

    output:
    path "genome_id_to_file_path_${genome_id}.tsv"
    path "meta_table_${genome_id}.tsv"

    script:
    strain_simulation_template = params.strain_simulation_template
    """
    # Run the Perl script
    touch empty_gff.gff
    ${projectDir}/scripts/sgEvolver/simujobrun.pl ${fasta} empty_gff.gff ${seed} ${strain_simulation_template}

    # Run the Python script
    python ${scripts_dir}/pick_random_strains.py ${amount} ${genome_id} ${NCBI_ID} ${novelty_category} ${OTU} ${strain_simulation_template} ${params.outdir}

    mkdir --parents ${params.outdir}/source_genomes/

    # Read the TSV file and copy each file to its destination
    while IFS=\$'\t' read -r genome_id dest_path; do
        base_name=\$(basename "\$dest_path")
        cp "\$base_name" "\$dest_path"
    done < "genome_id_to_file_path_${genome_id}.tsv"
    """
}

/*
* This process simulates strains using sgEvolver by using the specified gff files.
*/
process strain_simulation_with_gff {

    conda 'perl python'

    input:
    tuple val(genome_id), path(fasta), val(amount), val(seed), path(gff), val(OTU), val(NCBI_ID), val(novelty_category)

    output:
    path "genome_id_to_file_path_${genome_id}.tsv"
    path "meta_table_${genome_id}.tsv"
    path "sequence_id_map_genome_${genome_id}.txt", optional: true

    script:
    strain_simulation_template = params.strain_simulation_template
    """
    # Run the Perl script
    ${projectDir}/scripts/sgEvolver/simujobrun.pl ${fasta} ${gff} ${seed} ${strain_simulation_template}

    # Run the Python script
    python ${scripts_dir}/pick_random_strains.py ${amount} ${genome_id} ${NCBI_ID} ${novelty_category} ${OTU} ${strain_simulation_template} ${params.outdir}

    # Read the TSV file and copy each file to its destination
    while IFS=\$'\t' read -r genome_id dest_path; do
        base_name=\$(basename "\$dest_path")
        cp "\$base_name" "\$dest_path"
    done < "genome_id_to_file_path_${genome_id}.tsv"
    """
}

/*
* This process merges the metadata files created during the strain simulation with the user specified ones.
*/
process merge_metadata_files {

    input:
    path genome_id_to_file_paths
    path meta_tables
    path genome_location_file
    path metadata_table

    output:
    path 'merged_genome_location.tsv'
    path 'merged_meta_data.tsv'

    script:
    """
    # Write the content of genome_location to a new file
    cp ${genome_location_file} merged_genome_location.tsv

    # Append the content of each genome_id_to_file_path to the new file
    cat ${genome_id_to_file_paths} >> merged_genome_location.tsv

    # Write the content of genome_location to a new file
    cp ${metadata_table} merged_meta_data.tsv

    # Append the content of each genome_id_to_file_path to the new file
    cat ${meta_tables} >> merged_meta_data.tsv

    mkdir --parents ${params.outdir}/internal/
    cp merged_genome_location.tsv ${params.outdir}/internal/genome_locations.tsv
    cp merged_meta_data.tsv ${params.outdir}/internal/meta_data.tsv
    """
}

process cleanup_and_filter_sequences {

    conda "conda-forge::biopython=1.83 conda-forge::python=3.11.5"

    input:
    path genome_id_to_file_path
    val genome_entries

    output:
    path "internal_*"

    script:
    normalized_genome_entries = genome_entries.collect { entry ->
        if (entry instanceof List || entry instanceof Object[]) {
            if (entry.size() >= 2) {
                return [entry[0], entry[1]]
            }
        }
        throw new IllegalArgumentException("cleanup_and_filter_sequences expected [genome_id, path] entries but received: ${entry}")
    }
    absolute_genome_id_to_file_path = normalized_genome_entries.collect { "${it[0]}\t${it[1]}" }.join('\n')
    """
    mkdir --parents ${params.outdir}/source_genomes/
    mkdir --parents ${params.outdir}/internal/

cat > absolute_${genome_id_to_file_path} <<'EOF'
${absolute_genome_id_to_file_path}
EOF

    python ${scripts_dir}/clean_up_sequences.py absolute_${genome_id_to_file_path} ${params.outdir}/source_genomes/ internal_${genome_id_to_file_path}

    cp ./out_genomes/* ${params.outdir}/source_genomes/
    cp internal_${genome_id_to_file_path} ${params.outdir}/internal/genome_locations.tsv
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

    conda 'python'

    input:
    path(file_path_of_drawn_genome_location)
    val(seed)

    output:
    path 'distribution_*.txt'

    script:
    number_of_samples = params.number_of_samples
    mode = params.mode
    log_mu = params.log_mu
    log_sigma = params.log_sigma
    gauss_mu = params.gauss_mu
    gauss_sigma = params.gauss_sigma
    verbose = params.verbose
    """
    python ${shared_scripts_dir}/get_community_distribution.py ${number_of_samples} ${file_path_of_drawn_genome_location} ${mode} ${log_mu} ${log_sigma} ${gauss_mu} ${gauss_sigma} ${verbose} ${seed}
    mkdir --parents ${params.outdir}/distributions/
    cp distribution_*.txt ${params.outdir}/distributions/
    """
}

/*
* This process returns a file containing a random seed for every genome generated from the given seed in the config file.
* In case the simulated reads will be anonymized, it also returns a file containing a random seed for every sample generated from the given seed in the config file.
* Output:
*     The file with the given seed per samle in CSV format.
 */
process get_seed {

    conda 'python'

    input:
    path (genome_locations)
    val(seed)

    output:
    path ('seed.txt')
    path ('seed_read_anonymisation.txt'), optional: true
    path ('seed_gsa_anonymisation.txt'), optional: true
    path ('seed_pooled_gsa_anonymisation.txt'), optional: true
    path ('seed_merged_gsa_anonymisation.txt'), optional: true
    path ('seed_hybrid_gsa_anonymisation.txt'), optional: true

    script:
    count_samples = params.number_of_samples
    if (params.anonymization){
        param_anonym = "-anonym_seed"
    } else {
        param_anonym = ""
    }
    merged_count = 0
    if (params.containsKey('merged_gsa_combinations') && params.merged_gsa_combinations instanceof List) {
        merged_count = params.merged_gsa_combinations.size()
    }
    hybrid_count = 0
    if (params.containsKey('hybrid') && params.hybrid) {
        hybrid_count = params.number_of_samples
    }
    """
    ${shared_scripts_dir}/get_seed.py -seed ${seed} -count_samples ${count_samples} -file_genome_locations ${genome_locations} ${param_anonym} -merged_count ${merged_count} -hybrid_count ${hybrid_count}
    mkdir --parents ${params.outdir}/seed/
    cp seed*.txt ${params.outdir}/seed/
    """
}

/*
* This process computes the mean coverage per genome from a BAM file and writes
* a TSV file <label>_coverage.tsv with two columns:
*   genome_fasta_basename  mean_coverage
*
* It mirrors the logic of get_coverage.sh but runs inside the Nextflow workflow
* so it is automatically triggered for every gold standard assembly BAM.
*
* Input:
*   bam              – path to the (sorted, indexed) BAM file
*   label            – string used to name the output file, e.g. "sample_0_art"
*   out_subdir       – output directory path (e.g. ${params.outdir}/sample_0/contigs/art)
*   genome_locations – the genome_locations.tsv file (genome_id <TAB> fasta_path)
* Output:
*   The coverage TSV file.
*/
process compute_coverage {

    conda 'bioconda::samtools conda-forge::python'

    input:
    tuple path(bam_file), val(label), val(out_subdir), path(genome_locations)

    output:
    path(cov_file)

    script:
    cov_file = "${label}_coverage.tsv"
    """
    set -euo pipefail

    samtools coverage ${bam_file} \
      | python3 ${scripts_dir}/aggregate_coverage.py "${genome_locations}" "${cov_file}"

    mkdir --parents "${out_subdir}"
    cp "${cov_file}" "${out_subdir}/"
    """
}
