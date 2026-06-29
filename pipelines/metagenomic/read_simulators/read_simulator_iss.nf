scripts_dir = "${projectDir}/pipelines/metagenomic/scripts"
shared_scripts_dir = "${projectDir}/pipelines/shared/scripts"

/**
* This workflow simulates reads using the InSilicoSeq (iss) read simulator
* Takes:
*     A channel containing tuples with key = genome_id, first value = path to genome, second value = distribution, third value = sample_id, fourth value = seed.
* Emits:
*     A channel containing tuples with key = sample_id, first value = genome id, second value = simulated bam file, third value = the reference fasta file.
**/

workflow read_simulator_iss {

    take: genome_location_distribution_ch
    main:
        simulate_reads_iss(genome_location_distribution_ch)
    emit:
        simulate_reads_iss.out[0]
        simulate_reads_iss.out[1]
}

/**
* This process simulates short reads with iss
* Input:
*     Tuple containing key = genome_id, first value = path to genome, second value = distribution, third value = sample_id, fourth value = seed.
* Output:
*     Tuple containing key = sample_id, first value = genome_id, second value = path to bam file, third value = path to reference genome.
**/

process simulate_reads_iss {

    conda 'conda-forge::python=3.10 conda-forge::biopython=1.83 bioconda::insilicoseq=2.0.1 bioconda::minimap2 bioconda::samtools conda-forge::pigz'

    input:
    tuple val(genome_id), val(sample_id), path(fasta_file), val(abundance), val(seed)

    output:
    tuple val(sample_id), val(genome_id), path("sample${sample_id}_${genome_id}_iss.bam"), path(fasta_file), optional: true
    tuple val(sample_id), path("sample${sample_id}_${genome_id}_iss1.fq"), path("sample${sample_id}_${genome_id}_iss2.fq"), optional: true

    script:
    total_size = params.size
    model = params.iss.base_profile_name
    Long used_seed = ((seed as Long) % (2**32 - 1))
    threads = Math.max(1, ((task.cpus ?: 1) as int))

    """
    # Calculate read length and number of reads dynamically from the profile
    read_length=\$(python -c "import numpy as np; m = np.load('${model}', allow_pickle=True); print(int(m['read_length']))")
    number_of_reads=\$(python -c "import math; total_size = ${total_size}; abundance = ${abundance}; read_length = \${read_length}; print(round((total_size * 10**9) * abundance / read_length))")

    echo "DEBUG: read_length = \$read_length"
    echo "DEBUG: number_of_reads = \$number_of_reads"

    if [ "\$number_of_reads" -eq 0 ]; then
        echo "INFO: Skipping ISS for sample${sample_id}_${genome_id} because number_of_reads is 0."
        exit 0
    fi

    # iss requires number_of_reads >= 2
    if [ "\$number_of_reads" -eq 1 ]; then
        echo "INFO: Increasing ISS number_of_reads from 1 to 2 for sample${sample_id}_${genome_id}."
        number_of_reads=2
    fi

    # Run InSilicoSeq simulator
    # Use one cpu because parallel iss can fail for read counts < cpus
    iss generate --genomes ${fasta_file} --model ${model} --n_reads \${number_of_reads} --seed ${used_seed} --cpus 1 --output sample${sample_id}_${genome_id}_iss_temp

    if [ ! -s sample${sample_id}_${genome_id}_iss_temp_R1.fastq ] || [ ! -s sample${sample_id}_${genome_id}_iss_temp_R2.fastq ]; then
        echo "ERROR: ISS did not create non-empty paired FASTQ files for sample${sample_id}_${genome_id}." >&2
        exit 1
    fi

    # Rename outputs to follow standard naming conventions
    mv sample${sample_id}_${genome_id}_iss_temp_R1.fastq sample${sample_id}_${genome_id}_iss1.fq
    mv sample${sample_id}_${genome_id}_iss_temp_R2.fastq sample${sample_id}_${genome_id}_iss2.fq

    # Align generated reads back to reference using minimap2 to produce coordinate-sorted BAM file
    minimap2 -t ${threads} -ax sr ${fasta_file} sample${sample_id}_${genome_id}_iss1.fq sample${sample_id}_${genome_id}_iss2.fq | samtools view -bS | samtools sort -o sample${sample_id}_${genome_id}_iss.bam

    # Organize outputs
    mkdir --parents ${params.outdir}/sample_${sample_id}/bam/iss/
    cp sample${sample_id}_${genome_id}_iss.bam ${params.outdir}/sample_${sample_id}/bam/iss/

    pigz -p ${threads} -k sample${sample_id}_${genome_id}_iss1.fq
    pigz -p ${threads} -k sample${sample_id}_${genome_id}_iss2.fq

    mkdir --parents ${params.outdir}/sample_${sample_id}/reads/fastq/iss/
    cp sample${sample_id}_${genome_id}_iss1.fq.gz ${params.outdir}/sample_${sample_id}/reads/fastq/iss/
    cp sample${sample_id}_${genome_id}_iss2.fq.gz ${params.outdir}/sample_${sample_id}/reads/fastq/iss/
    """
}
