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

// include binning
include { binning } from "${projectDir}/pipelines/shared/binning"

/*
 * This is the main workflow and starting point of this nextflow pipeline.
 */
workflow metagenomic {

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

    // if this parameter is set, the metagenome simulation hast to be from the given profile
    if(!params.biom_profile.isEmpty()) {

        // start the simulation
        metagenomesimulation_from_profile()

        // get the output channel
        genome_distribution_file_ch = metagenomesimulation_from_profile.out[0]
        genome_location_file_ch = metagenomesimulation_from_profile.out[1]
        metadata_ch = metagenomesimulation_from_profile.out[2]

        // Stop the pipeline if just community design steps are needed
        if(params.just_community_design){
            println "Simulation stopping after community design steps."
            return
        }

    } else { // not from profile

        // this channel holds the file with the specified locations of the genomes
        genome_location_file_ch = Channel.fromPath(params.genome_locations_file)

        metadata_ch = Channel.fromPath(params.metadata_file) 

        // if there are distribution files given for each sample use those
        if(params.distribution_files.isEmpty()) {

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

            // Stop the pipeline if just community design steps are needed
            if(params.just_community_design){
                println "Simulation stopping after community design steps."
                return
            }

        // otherwise, create distribution files for each sample
        } else {
        
            // this channel holds the files with the specified distributions for every sample
            genome_distribution_file_ch = Channel.fromPath(params.distribution_files)
        }

    }    
    // build ncbi taxonomy from given tax dump
    number_of_samples_ch = Channel.from(params.number_of_samples)
    buildTaxonomy(number_of_samples_ch.concat(ncbi_taxdump_file_ch.concat(genome_distribution_file_ch)).toList().map { it -> [ it[0], it[1], it[2..-1] ] }, metadata_ch)
    

    if(params.type.equals("nanosim3")) {

        if(params.read_length != null) {
            read_length_ch = params.read_length
        } else {
            read_length_ch = calculate_Nanosim_read_length(params.base_profile_name) // this takes very long
        }
    } else {
        read_length_ch = params.profile_read_length
    }

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
    genome_location_file_ch = cleanup_and_filter_sequences(genome_location_file_ch, genome_location_ch.map { it[1] }.collect())
    
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
    sample_wise_simulation(genome_location_ch, genome_location_file_ch, genome_distribution_file_ch, read_length_ch, seed_file_read_simulation_ch)
    // this workflow has two output channels: one bam file per sample and one fasta file per sample
    merged_bam_per_sample = sample_wise_simulation.out[0]
    gsa_for_all_reads_of_one_sample_ch = sample_wise_simulation.out[1]    

    // extract file paths from the tuples to create the reference_fasta_files_ch
    reference_fasta_files_ch = genome_location_ch.map { a -> a[1] }

    // merge the bam files required for the pooled gsa
    if (params.pooled_gsa instanceof Boolean && params.pooled_gsa) {
        merged_bam_file = merge_bam_files(merged_bam_per_sample.map { it[1] }.collect())
    } else if (params.pooled_gsa instanceof List) {
        merged_bam_file = merge_bam_files(merged_bam_per_sample.filter { params.pooled_gsa*.toString().contains(it[0]) }.map { it[1] }.collect())
    }

    generate_pooled_gold_standard_assembly(merged_bam_file.combine(reference_fasta_files_ch).groupTuple())    

    // if requested, anonymize reads, gsa and pooled gsa
    if(params.anonymization) {
        anonymization(sample_wise_simulation.out[2], get_seed.out[1], get_seed.out[2], get_seed.out[3], gsa_for_all_reads_of_one_sample_ch, sample_wise_simulation.out[3], generate_pooled_gold_standard_assembly.out, merged_bam_file, genome_location_file_ch, metadata_ch)
    } else { // if no anonymization is requested, create binning gold standard
        binning(gsa_for_all_reads_of_one_sample_ch, sample_wise_simulation.out[3], generate_pooled_gold_standard_assembly.out, merged_bam_file, genome_location_file_ch, metadata_ch)
    }
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
* This process calculates the average read length of Nanosim reads from the pickle of the predefined profile
*
*/
process calculate_Nanosim_read_length {
    // TODO: Packages which are needed multiple times should be loaded only once
    conda 'conda-forge::scikit-learn=0.21.3=py37* conda-forge::joblib=1.2.0'

    input:
    val profile

    output:
    stdout

    script:
    """
    #!/usr/bin/env python
    import joblib
    import sys
    import numpy as np
    from scipy.integrate import quad
    from sklearn.neighbors import KernelDensity

    read_length_file = "${profile}_aligned_reads.pkl"
    #default is {prefix}_aligned_reads.pkl

    kde = joblib.load(read_length_file) # length is stored as joblib pkl

    # the kd has a probability density function from which we can get mean and variance via integration
    # it is the log density function though, need to np.exp
    pdf = lambda x : np.exp(kde.score_samples([[x]]))[0]
    mean = quad(lambda x: x * pdf(x), a=-np.inf, b=np.inf)[0]
    print(mean)
    """
}

/*
* This process merges all given bam files specified in the pooled_gsa parameter.
* Takes:
*     A list with the paths to all bam files, that should be merged, if the condition is fullfilled.
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

    bam_to_merge = ''

    bam_files.each {

        bam_file_name = (String) it
        sample_id = bam_file_name.split('_')[1][0].toInteger()
        
        //if(sample_id in params.pooled_gsa){
        bam_to_merge = bam_to_merge.concat(' ').concat(bam_file_name)
        //}
    }
    """
    samtools merge -u - ${bam_to_merge} | samtools sort -l ${compression} -m ${memory}G -o ${file_name} -O bam
    """
}

/*
* This process generates the pooled gold standard assembly for serveral samples.
* Takes:
*     A tuple with first_value = a sorted bam file and second value = the reference genome (fasta).
* Output:
*     The path to fasta file with the pooled gold standard assembly.
 */
process generate_pooled_gold_standard_assembly {

    conda 'bioconda::samtools'

    input:
    tuple path(bam_file), path(reference_fasta_files)

    output:
    path file_name
    
    script:
    file_name = 'gsa_pooled.fasta'
    """
    cat ${reference_fasta_files} > reference.fasta
    perl -- ${shared_scripts_dir}/bamToGold.pl -st samtools -r reference.fasta -b ${bam_file} -l 1 -c 1 >> ${file_name}
    mkdir --parents ${params.outdir}/pooled_gsa
    gzip -k ${file_name}
    cp ${file_name}.gz ${params.outdir}/pooled_gsa/
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

    input:
    tuple val(number_of_samples), path(dmp), path(distribution_files)
    path(metadata_ch)

    output:
    path 'taxonomic_profile_*.txt'

    script:
    index_number_of_samples = number_of_samples - 1
    """
    tar -xf ${dmp}
    [ -f **/names.dmp ] && mv **/names.dmp ./names.dmp
    [ -f **/merged.dmp ] && mv **/merged.dmp ./merged.dmp
    [ -f **/nodes.dmp ] && mv **/nodes.dmp ./nodes.dmp
    ${scripts_dir}/build_ncbi_taxonomy.py names.dmp merged.dmp nodes.dmp ${number_of_samples} ${metadata_ch} ${distribution_files}
    mkdir --parents ${params.outdir}
    cp taxonomic_profile_*.txt ${params.outdir}
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

    input:
    tuple val(genome_id), path(fasta), val(amount), val(seed), path(gff)
    val seed

    output:
    tuple val(genome_id), path("genome_id_to_file_path_${genome_id}.tsv")
    tuple val(genome_id), path("meta_table_${genome_id}.tsv")
    tuple val(genome_id), path("sequence_id_map_genome_${genome_id}.txt")

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
    path genomes

    output:
    path genome_id_to_file_path

    script:
    """
    mkdir --parents ${params.outdir}/source_genomes/
    mkdir --parents ${params.outdir}/internal/

    touch internal_${genome_id_to_file_path}

    python ${scripts_dir}/clean_up_sequences.py ${genome_id_to_file_path} ${params.outdir}/source_genomes/ internal_${genome_id_to_file_path}

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
    ${shared_scripts_dir}/get_seed.py -seed ${seed} -count_samples ${count_samples} -file_genome_locations ${genome_locations} ${param_anonym}
    mkdir --parents ${params.outdir}/seed/
    cp seed*.txt ${params.outdir}/seed/
    """
}
