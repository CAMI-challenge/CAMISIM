/*
 * Defining the module / subworkflow path, and include the elements
 */

scripts_dir = "${projectDir}/pipelines/metagenomic/scripts"
shared_scripts_dir = "${projectDir}/pipelines/shared/scripts"

// include read simulator here:
read_simulator_folder = "${projectDir}/pipelines/metagenomic/read_simulators/"
// include read simulator nanosim3
include { read_simulator_nanosim3 } from "${read_simulator_folder}/read_simulator_nansoim3"
include { read_simulator_art } from "${read_simulator_folder}/read_simulator_art"
include { read_simulator_art_modern } from "${read_simulator_folder}/read_simulator_art_modern"
include { read_simulator_wgsim } from "${read_simulator_folder}/read_simulator_wgsim"
include { read_simulator_iss } from "${read_simulator_folder}/read_simulator_iss"
include {
    normalise_abundance as normalise_abundance_standard
    normalise_abundance as normalise_abundance_size_norm
    normalise_abundance_to_size
    count_bases
} from "${projectDir}/pipelines/shared/distribution"

/** 
* This workflow simulates reads for every sample.
* Takes:
*     genome_location_ch: Path to the file containing the genome locations.
      genome_distribution_ch: Paths to files containing the distributions of every genome for every sample.
* Emits: 
*     A channel containing the merged fastq file and the merged bam file over all genomes for every sample.
**/
/*
* A genome_id belongs to a contamination ASV when it equals the ASV id or is one of
* its strain genomes '<asv>.<n>'. The trailing '.' prevents 'hASV1' matching 'hASV10'.
* Mirrors _is_contamination_genome in build_ncbi_taxonomy.py / goldstandardfileformat.py.
*/
def is_contamination_genome(genome_id, contamination_asvs) {
    def gid = genome_id.toString()
    return contamination_asvs.any { asv -> def a = asv.toString(); gid == a || gid.startsWith(a + '.') }
}

workflow sample_wise_simulation {

    take: genome_location_ch
    take: genome_location_file_ch
    take: genome_distribution_file_ch
    take: seed_file_ch
    main:

        def types_list = (params.type instanceof List) ? params.type : [params.type].collect { it.toString() }.unique()

        // get the seed for every genome
        seed_ch = seed_file_ch.splitCsv(sep:'\t', skip:2)

        // get tuple with key = sample id, first value = genome_id, second value = distribution from distribution files
        distribution_file_ch = genome_distribution_file_ch.flatten().map { file -> tuple(file.baseName.split('_')[1], file) }.splitCsv(sep:'\t').map { a -> tuple(a[1][0], a[0], a[1][1]) }

        // for read simulators that need size-normalised abundance (nanosim3, wgsim, iss):
        def needs_size_norm = types_list.any { it == "nanosim3" || it == "wgsim" || it == "iss" }
        if (needs_size_norm) {
            genome_location_distribution_ch_tmp = genome_location_ch.combine(distribution_file_ch, by: 0)
            counted_bases_ch = count_bases(genome_location_distribution_ch_tmp)
            size_norm_distribution_file_ch = normalise_abundance_to_size(counted_bases_ch)
            genome_size_ch = counted_bases_ch.map { tuple(it[0], it[2], it[3].toString().trim()) }

            size_norm_normalised_distribution_ch = normalise_abundance_size_norm(size_norm_distribution_file_ch.map { a -> tuple(a[1], tuple(a[0], a[2])) }.groupTuple())
                .map { file -> tuple(file.baseName.split('_')[2], file) }.splitCsv(sep:'\t').map { a -> tuple(a[1][0], a[1][1], a[0]) }.filter { it[1] != '0.0' }
            location_distribution_seed_size_norm_ch = genome_location_ch.combine(size_norm_normalised_distribution_ch, by: 0).map { a -> tuple(a[0], a[3], a[1], a[2]) }.combine(seed_ch, by:[0,1])
        }

        def needs_coverage_factor = types_list.any { it == "art" || it == "art_modern" }
        // standard non-size-normalised distributions are only needed by art / art_modern
        if (needs_coverage_factor) {
            normalised_distribution_ch = normalise_abundance_standard(distribution_file_ch.map { a -> tuple(a[1], tuple(a[0], a[2])) }.groupTuple())
                .map { file -> tuple(file.baseName.split('_')[2], file) }.splitCsv(sep:'\t').map { a -> tuple(a[1][0], a[1][1], a[0]) }.filter { it[1] != '0.0' }

            genome_distribution_location_ch = genome_distribution_file_ch.flatten().map { file -> tuple(file.baseName.split('_')[1], file) }.combine(genome_location_file_ch)
        }

        // ---- accumulators across all types ----
        // Each type contributes to these channels; they are mixed together.
        all_type_bam_ch              = Channel.empty()  // tuple(sim_type, sample_id, genome_id, bam, ref_fasta)  – per genome, pre-merge
        all_type_reads_ch            = Channel.empty()  // tuple(sim_type, sample_id, …reads…)
        all_type_bam_by_sample_ch    = Channel.empty()  // tuple(sim_type, sample_id, [bam_paths])  – after groupTuple per type
        paired_end_reads_to_write_ch = Channel.empty()  // tuple(sim_type, sample_id, r1_files, r2_files)
        single_end_reads_to_write_ch = Channel.empty()  // tuple(sim_type, sample_id, read_files)

        // ======================================================================
        // Run each type
        // ======================================================================
        for (sim_type in types_list) {
            if (sim_type == "art") {
                // The simulated ART reads (version 016.06.05) doesn't contain the whole header of the reference genome, if there is a space in the header. They then just contain the
                // substring before the first occurance of the space. In that case the gold standard assembly doesn't work because there are no matching IDs.
                // As a workaround we change the headers of the reference genomes by just selecting the part before the first occurance of a space, if there is a space in the header.
                // TODO does this happen for ART_modern?
                art_genome_location_ch = remove_spaces_from_reference_genome(genome_location_ch)

                art_location_distribution_seed_ch = art_genome_location_ch.combine(normalised_distribution_ch, by: 0).map { a -> tuple(a[0], a[3], a[1], a[2]) }.combine(seed_ch, by:[0,1])
                art_factor_for_sample_id_ch = get_multiplication_factor_art(genome_distribution_location_ch)
                art_genome_location_distribution_factor_ch = art_location_distribution_seed_ch.map { tuple(it[1], *it) }.combine(art_factor_for_sample_id_ch, by: 0).map { it[1..-1] }

                read_simulator_art(art_genome_location_distribution_factor_ch)
                type_bam_ch   = read_simulator_art.out[0].map { sid, gid, bam, ref -> tuple("art", sid, gid, bam, ref) }
                type_reads_ch = read_simulator_art.out[1].map { sid, r1, r2 -> tuple("art", sid, r1, r2) }

                // get_fastq_for_sample_paired_end_typed(read_simulator_art.out[1], sim_type)
                paired_end_reads_to_write_ch = paired_end_reads_to_write_ch.mix(type_reads_ch)
                all_type_bam_ch   = all_type_bam_ch.mix(type_bam_ch)
                all_type_reads_ch = all_type_reads_ch.mix(type_reads_ch)

            } else if (sim_type == "art_modern") {
                // combine genome location + normalised distribution + seed
                location_distribution_seed_ch = genome_location_ch.combine(normalised_distribution_ch, by: 0).map { a -> tuple(a[0], a[3], a[1], a[2]) }.combine(seed_ch, by:[0,1])
                factor_for_sample_id_ch = get_multiplication_factor_art_modern(genome_distribution_location_ch)
                genome_location_distribution_factor_ch = location_distribution_seed_ch.map { tuple(it[1], *it) }.combine(factor_for_sample_id_ch, by: 0).map { it[1..-1] }

                read_simulator_art_modern(genome_location_distribution_factor_ch)
                type_bam_ch   = read_simulator_art_modern.out[0].map { sid, gid, bam, ref -> tuple("art_modern", sid, gid, bam, ref) }
                type_reads_ch = read_simulator_art_modern.out[1].map { sid, r1, r2 -> tuple("art_modern", sid, r1, r2) }

                // get_fastq_for_sample_paired_end_typed(read_simulator_art_modern.out[1], sim_type)
                paired_end_reads_to_write_ch = paired_end_reads_to_write_ch.mix(type_reads_ch)
                all_type_bam_ch   = all_type_bam_ch.mix(type_bam_ch)
                all_type_reads_ch = all_type_reads_ch.mix(type_reads_ch)

            } else if (sim_type == "nanosim3") {
                location_distribution_seed_size_ch = location_distribution_seed_size_norm_ch.combine(genome_size_ch, by:[0,1])

                read_simulator_nanosim3(location_distribution_seed_size_ch)
                type_bam_ch   = read_simulator_nanosim3.out[0].map { sid, gid, bam, ref -> tuple("nanosim3", sid, gid, bam, ref) }
                type_reads_ch = read_simulator_nanosim3.out[1].map { sid, reads -> tuple("nanosim3", sid, reads) }

                // get_fastq_for_sample_single_end_typed(read_simulator_nanosim3.out[1], sim_type)
                single_end_reads_to_write_ch = single_end_reads_to_write_ch.mix(type_reads_ch)
                all_type_bam_ch   = all_type_bam_ch.mix(type_bam_ch)
                all_type_reads_ch = all_type_reads_ch.mix(type_reads_ch)

            } else if (sim_type == "wgsim") {
                read_simulator_wgsim(location_distribution_seed_size_norm_ch)
                type_bam_ch   = read_simulator_wgsim.out[0].map { sid, gid, bam, ref -> tuple("wgsim", sid, gid, bam, ref) }
                type_reads_ch = read_simulator_wgsim.out[1].map { sid, r1, r2 -> tuple("wgsim", sid, r1, r2) }

                // get_fastq_for_sample_paired_end_typed(read_simulator_wgsim.out[1], sim_type)
                paired_end_reads_to_write_ch = paired_end_reads_to_write_ch.mix(type_reads_ch)
                all_type_bam_ch   = all_type_bam_ch.mix(type_bam_ch)
                all_type_reads_ch = all_type_reads_ch.mix(type_reads_ch)

            } else if (sim_type == "iss") {
                read_simulator_iss(location_distribution_seed_size_norm_ch)
                type_bam_ch   = read_simulator_iss.out[0].map { sid, gid, bam, ref -> tuple("iss", sid, gid, bam, ref) }
                type_reads_ch = read_simulator_iss.out[1].map { sid, r1, r2 -> tuple("iss", sid, r1, r2) }

                paired_end_reads_to_write_ch = paired_end_reads_to_write_ch.mix(type_reads_ch)
                all_type_bam_ch   = all_type_bam_ch.mix(type_bam_ch)
                all_type_reads_ch = all_type_reads_ch.mix(type_reads_ch)
            }
        }

        // Contamination ASVs: their reads remain in all_type_reads_ch and the fastq
        // channels above, but drop their per-genome BAMs here so they are excluded from
        // every downstream gold standard assembly, gsa_mapping, and coverage file (all of
        // which derive from all_type_bam_ch). Reads output is intentionally left untouched.
        if (params.containsKey('contamination_asvs') && params.contamination_asvs) {
            all_type_bam_ch = all_type_bam_ch.filter { read_type, sid, gid, bam, ref -> !is_contamination_genome(gid, params.contamination_asvs) }
        }

        get_fastq_for_sample_paired_end_typed(paired_end_reads_to_write_ch.groupTuple(by: [0, 1]))
        get_fastq_for_sample_single_end_typed(single_end_reads_to_write_ch.groupTuple(by: [0, 1]))

        // ---- per-type: generate GSA per genome, then merge to per-type per-sample GSA ----
        // strip sim_type from the BAM tuples for generate_gold_standard_assembly (unchanged process signature)
        typed_bam_for_gsa_ch = all_type_bam_ch.map { read_type, sid, gid, bam, ref -> tuple(read_type, sid, gid, bam, ref) }
        gsa_per_type_per_genome_ch = generate_gold_standard_assembly_typed(typed_bam_for_gsa_ch)

        // group by (sim_type, sample_id) -> all per-genome GSAs for that type+sample
        grouped_gsa_per_type_ch = gsa_per_type_per_genome_ch.groupTuple(by: [0, 1])
        gsa_per_type_per_sample_ch = get_fasta_for_sample_typed(grouped_gsa_per_type_ch)

        // ---- per-type: merge BAMs per sample ----
        // group: (sim_type, sample_id) -> [bam_paths]
        bam_files_by_sample_per_type_ch = all_type_bam_ch.map { read_type, sid, gid, bam, ref -> tuple(read_type, sid, bam) }.groupTuple(by: [0, 1])
        merged_bam_per_type_per_sample_ch = merge_bam_files_typed(bam_files_by_sample_per_type_ch)

        // ---- combined (all-types): merge BAMs per sample for hybrid GSA ----
        all_bam_flat_ch = all_type_bam_ch.map { read_type, sid, gid, bam, ref -> tuple(sid, bam) }
        all_bam_by_sample_ch = all_bam_flat_ch.groupTuple()
        combined_bam_per_sample_ch = merge_bam_files_combined(all_bam_by_sample_ch)

        // ---- reads grouped per-type per-sample (used by anonymization) ----
        reads_per_type_grouped_ch = all_type_reads_ch.groupTuple(by: [0, 1])

    emit: merged_bam_per_type_per_sample_ch        // [0]: tuple(read_type, sample_id, bam)
    emit: gsa_per_type_per_sample_ch               // [1]: tuple(read_type, sample_id, gsa_fasta)
    emit: reads_per_type_grouped_ch                // [2]: tuple(read_type, sample_id, …reads…)
    emit: bam_files_by_sample_per_type_ch          // [3]: tuple(read_type, sample_id, [bam_paths])
    emit: combined_bam_per_sample_ch               // [4]: tuple(sample_id, bam) – all types combined
}

/*
* This process generates a gold standard assembly for one genome.
* Takes:
*     A tuple with key = sample_id, first_value = genome_id, second value = a sorted bam file, third value = the reference genome (fasta).
* Output:
*     A Tuple with key = sample_id, value = path to fasta file with the gold standard assembly of the given genome.
 */
process generate_gold_standard_assembly {

    conda 'bioconda::samtools conda-forge::pigz'

    input:
    tuple val(sample_id),val(genome_id), path(bam_file), path(reference_fasta_file)

    output:
    tuple val(sample_id), path(file_name)

    when: params.gsa  // Only execute this process when gsa is set to true

    script:
    file_name = 'sample'.concat(sample_id.toString()).concat('_').concat(genome_id).concat('_gsa.fasta')
    threads = Math.max(1, ((task.cpus ?: 1) as int))
    """
    samtools faidx ${reference_fasta_file}
    python ${shared_scripts_dir}/bamToGold.py -st samtools -r ${reference_fasta_file} -b ${bam_file} -l 1 -c 1 >> ${file_name}
    mkdir --parents ${params.outdir}/sample_${sample_id}/gsa
    pigz -p ${threads} -k ${file_name}
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

    conda 'conda-forge::pigz'

    input:
    tuple val(sample_id), path(fasta_files)

    output:
    tuple val(sample_id), path(file_name)

    script:
    file_name = 'sample'.concat(sample_id.toString()).concat('_gsa.fasta')
    threads = Math.max(1, ((task.cpus ?: 1) as int))
    """
    set -euo pipefail

    # Sort files before concatenation to ensure reproducibility
    printf '%s\\n' ${fasta_files} | sort | xargs -r cat -- > ${file_name}

    mkdir --parents ${params.outdir}/sample_${sample_id}/contigs
    pigz -p ${threads} -k ${file_name}
    cp ${file_name}.gz ${params.outdir}/sample_${sample_id}/contigs/gsa.fasta.gz
    """
}

/*
* Generates a gold standard assembly for one genome for one specific sequencing type.
* The sim_type tag is carried through so per-type GSAs stay separate.
*/
process generate_gold_standard_assembly_typed {

    conda 'bioconda::samtools conda-forge::pigz'

    input:
    tuple val(sim_type), val(sample_id), val(genome_id), path(bam_file), path(reference_fasta_file)

    output:
    tuple val(sim_type), val(sample_id), path(file_name)

    when: params.gsa

    script:
    file_name = "sample${sample_id}_${genome_id}_${sim_type}_gsa.fasta"
    threads = Math.max(1, ((task.cpus ?: 1) as int))
    """
    samtools faidx ${reference_fasta_file}
    python ${shared_scripts_dir}/bamToGold.py -st samtools -r ${reference_fasta_file} -b ${bam_file} -l 1 -c 1 >> ${file_name}
    mkdir --parents ${params.outdir}/sample_${sample_id}/gsa/${sim_type}
    pigz -p ${threads} -k ${file_name}
    cp ${file_name}.gz ${params.outdir}/sample_${sample_id}/gsa/${sim_type}/
    """
}

/*
* Creates a single FASTA holding all per-genome GSAs for one type+sample.
*/
process get_fasta_for_sample_typed {

    conda 'conda-forge::pigz'

    input:
    tuple val(sim_type), val(sample_id), path(fasta_files)

    output:
    tuple val(sim_type), val(sample_id), path(file_name)

    script:
    file_name = "sample${sample_id}_${sim_type}_gsa.fasta"
    threads = Math.max(1, ((task.cpus ?: 1) as int))
    """
    set -euo pipefail

    printf '%s\\n' ${fasta_files} | sort | xargs -r cat -- > ${file_name}

    mkdir --parents ${params.outdir}/sample_${sample_id}/contigs/${sim_type}
    pigz -p ${threads} -k ${file_name}
    cp ${file_name}.gz ${params.outdir}/sample_${sample_id}/contigs/${sim_type}/gsa_${sim_type}.fasta.gz
    """
}

/*
* Merges all BAM files for one type+sample into a single sorted BAM.
*/
process merge_bam_files_typed {

    conda 'bioconda::samtools'

    input:
    tuple val(sim_type), val(sample_id), path(bam_files)

    output:
    tuple val(sim_type), val(sample_id), path(file_name)

    script:
    file_name = "sample_${sample_id}_${sim_type}.bam"
    compression = 5
    memory = 1
    threads_for_sort = Math.max(1, ((task.cpus ?: 1) as int))
    """
    samtools merge -u - ${bam_files} | samtools sort -@ ${threads_for_sort} -l ${compression} -m ${memory}G -o ${file_name} -O bam
    samtools index ${file_name}
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
    tuple val(sample_id), path(file_name)

    script:
    file_name = 'sample_'.concat(sample_id.toString()).concat('.bam')
    compression = 5
    memory = 1
    threads_for_sort = Math.max(1, ((task.cpus ?: 1) as int))
    """
    samtools merge -u - ${bam_files} | samtools sort -@ ${threads_for_sort} -l ${compression} -m ${memory}G -o ${file_name} -O bam
    samtools index ${file_name}
    """
}

/*
* Merges BAM files from ALL types for one sample into a single sorted BAM.
* Used as input for the hybrid gold standard assembly.
*/
process merge_bam_files_combined {

    conda 'bioconda::samtools'

    input:
    tuple val(sample_id), path(bam_files)

    output:
    tuple val(sample_id), path(file_name)

    script:
    file_name = "sample_${sample_id}_combined.bam"
    compression = 5
    memory = 1
    threads_for_sort = Math.max(1, ((task.cpus ?: 1) as int))
    """
    samtools merge -u - ${bam_files} | samtools sort -@ ${threads_for_sort} -l ${compression} -m ${memory}G -o ${file_name} -O bam
    samtools index ${file_name}
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
process get_multiplication_factor_art {

    conda "conda-forge::biopython"

    input:
    tuple val(sample_id), path(file_path_distribution), path(genome_locations)

    output:
    tuple val(sample_id), stdout

    script:
    factor = 0
    total_size = new BigDecimal(params.size).multiply(new BigDecimal(10**9))
    """
    python ${scripts_dir}/get_multiplication_factor.py ${total_size} "${genome_locations}" "${file_path_distribution}" "${projectDir}"
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
process get_multiplication_factor_art_modern {

    conda "conda-forge::biopython"

    input:
    tuple val(sample_id), path(file_path_distribution), path(genome_locations)

    output:
    tuple val(sample_id), stdout

    script:
    factor = 0
    total_size = new BigDecimal(params.size).multiply(new BigDecimal(10**9))
    """
    python ${scripts_dir}/get_multiplication_factor.py ${total_size} "${genome_locations}" "${file_path_distribution}" "${projectDir}"
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
* Writes all fastq files (single-end) for one type+sample into one file.
* sim_type is embedded in the output filename to avoid collisions.
*/
process get_fastq_for_sample_single_end_typed {

    conda 'conda-forge::pigz'

    input:
    tuple val(sim_type), val(sample_id), path(read_files)

    script:
    """
    set -euo pipefail

    # Sort files before concatenation to ensure reproducibility
    printf '%s\\n' ${read_files} | sort | xargs -r cat -- | pigz -p ${task.cpus} -c > sample_${sample_id}_${sim_type}.fq.gz

    mkdir --parents ${params.outdir}/sample_${sample_id}/reads/fastq/${sim_type}/
    cp sample_${sample_id}_${sim_type}.fq.gz ${params.outdir}/sample_${sample_id}/reads/fastq/${sim_type}/
    """
}

/*
* Writes all fastq files (paired-end) for one type+sample into one pair of files.
* sim_type is embedded in the output filenames to avoid collisions.
*/
process get_fastq_for_sample_paired_end_typed {

    conda 'conda-forge::pigz'

    input:
    tuple val(sim_type), val(sample_id), path(first_read_files), path(second_read_files)

    script:
    """
    set -euo pipefail

    # Sort files before concatenation to ensure reproducibility
    printf '%s\\n' ${first_read_files}  | sort | xargs -r cat -- | pigz -p ${task.cpus} -c > sample_${sample_id}_${sim_type}_01.fq.gz
    printf '%s\\n' ${second_read_files} | sort | xargs -r cat -- | pigz -p ${task.cpus} -c > sample_${sample_id}_${sim_type}_02.fq.gz

    mkdir --parents ${params.outdir}/sample_${sample_id}/reads/fastq/${sim_type}/
    cp sample_${sample_id}_${sim_type}_01.fq.gz ${params.outdir}/sample_${sample_id}/reads/fastq/${sim_type}/
    cp sample_${sample_id}_${sim_type}_02.fq.gz ${params.outdir}/sample_${sample_id}/reads/fastq/${sim_type}/
    """
}
