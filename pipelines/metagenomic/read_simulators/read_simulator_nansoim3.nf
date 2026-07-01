scripts_dir = "${projectDir}/pipelines/metagenomic/scripts"
shared_scripts_dir = "${projectDir}/pipelines/shared/scripts"

/** 
* This workflow simulates reads via nanosim3 and converts the resulting sam files into bam files.
* Takes:
*     A channel containing tuples with key = genome_id, first value = path to genome, second value = distribution, third value = sample_id, fourth value = seed.
* Emits: 
*     A channel containing tuples with key = sample_id, first value = genome id, second value = simulated bam file, third value = the reference fasta file.
**/
workflow read_simulator_nanosim3 {

    take: genome_location_distribution_ch
    main:
        // simulate reads via nanosim3

        // Precompute safe_max and join it back to the original channel by IDs
        precompute_limit(genome_location_distribution_ch)

        // This joins [genome_id, sample_id, path, abundance, seed, size]
        // with [genome_id, sample_id, safe_max]
        ch_with_limit = genome_location_distribution_ch.join(precompute_limit.out, by: [0, 1])

        // NanoSim derives read length from the trained model (-c) and the -max cap;
        // the previously-computed mean read length was never used by simulator.py, so it
        // (and the sklearn-pickle load it required) has been removed.

        // simulate reads in fastq format with nanosim directly
        if(params.nanosim3.simulate_fastq_directly) {

            simulate_reads_fastq_nanosim3(ch_with_limit)
            read_ch = simulate_reads_fastq_nanosim3.out[1]

            bam_ch = bam_from_reads_fastq(simulate_reads_fastq_nanosim3.out[0])

        } else { // simulate reads in fasta format with nanosim and convert them later
            simulate_reads_fasta_nanosim3(ch_with_limit)

            bam_ch = bam_from_reads_fasta(simulate_reads_fasta_nanosim3.out[0])[0]
            read_ch = bam_from_reads_fasta.out[1]
        }
        
    emit:
        bam_ch
        read_ch
}

/**
* This process calculates the safe -max parameter for NanoSim.
**/
process precompute_limit {

    conda "conda-forge::python=3.10"

    input:
    tuple val(genome_id), val(sample_id), path(fasta_file), val(abundance), val (seed), val(genome_size)

    output:
    tuple val(genome_id), val(sample_id), stdout

    script:
    """
    #!/usr/bin/env python
    import sys

    def read_fasta_max_len(fasta_path):
        max_len = 0
        current_len = 0
        try:
            with open(fasta_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line: continue
                    if line.startswith('>'):
                        if current_len > max_len:
                            max_len = current_len
                        current_len = 0
                    else:
                        current_len += len(line)
                if current_len > max_len:
                    max_len = current_len
        except FileNotFoundError:
            sys.exit(1)
        return max_len

    max_contig = read_fasta_max_len("${fasta_file}")
    safe_max = max_contig - 200
    if safe_max < 50:
        safe_max = max(50, max_contig - 10)

    print(safe_max, end='')
    """
}

/**
* This process simulates reads in fasta format with nanosim3.
* Input:
*     Tuple containing key = genome_id, first value = path to genome, second value = distribution, third value = sample_id, fourth value = seed.
* Output:
*     Tuple containing key = sample_id, first value = genome_id, second value = path to error profile, third value = path to fasta file with the aligned reads, 
*         fourth value = path to fasta file with the aligned reads, fifth value = path to reference genome.
**/
process simulate_reads_fasta_nanosim3 {

    conda 'conda-forge::scikit-learn=0.23.2 conda-forge::numpy=1.23.5 bioconda::nanosim=3.2'

    input:
    tuple val(genome_id), val(sample_id), path(fasta_file), val(abundance), val (seed), val(genome_size), val(safe_max)

    output:
    tuple val(sample_id), val(genome_id), path('*_error_profile'), path("*_aligned_reads.fasta"), path("*_unaligned_reads.fasta"), path(fasta_file)
    
    script:
    total_size = new BigDecimal(params.size).multiply(new BigDecimal(10**9))
    profile = params.nanosim3.base_profile_name
    coverage = total_size.multiply(new BigDecimal(abundance.toString())).divide(new BigDecimal(genome_size.toString()), java.math.MathContext.DECIMAL64)
    coverage = coverage.stripTrailingZeros()
    // nanosim seed cannot be > 2**32 -1
    Long used_seed = ((seed as Long) % (2**32 - 1))

    /**
    String log = "---- sample id: ".concat(sample_id)
    log = log.concat("  genome id: ").concat(genome_id)
    log = log.concat("   fasta file: ").concat(fasta_file.baseName)
    log = log.concat("  abundance: ").concat(abundance)
    log = log.concat("    seed: ").concat(seed)
    log = log.concat("    used_seed: ").concat(Long.toString(used_seed))
    log = log.concat("    coverage: ").concat(coverage.toPlainString())
    log = log.concat("    profile: ").concat(profile)
    print(log)
    **/

    """
    simulator.py genome -x ${coverage} -rg ${fasta_file} -o sample${sample_id}_${genome_id}_nanosim3 -c ${profile} --seed ${used_seed} -dna_type linear -min ${params.nanosim3.min_read_length} -max ${safe_max}
    """
}

/**
* This process simulates reads in fastq format with nanosim3.
* Input:
*     Tuple containing key = genome_id, first value = path to genome, second value = distribution, third value = sample_id, fourth value = seed.
* Output:
*     Tuple containing key = sample_id, first value = genome_id, second value = path to error profile, third value = path to fasta file with the aligned reads, 
*         fourth value = path to fasta file with the aligned reads, fifth value = path to reference genome.
**/
process simulate_reads_fastq_nanosim3 {

    conda 'conda-forge::scikit-learn=0.23.2 conda-forge::numpy=1.23.5 bioconda::nanosim=3.2 conda-forge::pigz'

    input:
    tuple val(genome_id), val(sample_id), path(fasta_file), val(abundance), val (seed), val(genome_size), val(safe_max)

    output:
    tuple val(sample_id), val(genome_id), path('*_error_profile'), path("*_aligned_reads.fastq"), path("*_unaligned_reads.fastq"), path(fasta_file)
    tuple val(sample_id), path("*_aligned_reads.fastq")
    
    script:
    total_size = new BigDecimal(params.size).multiply(new BigDecimal(10**9))
    profile = params.nanosim3.base_profile_name
    coverage = total_size.multiply(new BigDecimal(abundance.toString())).divide(new BigDecimal(genome_size.toString()), java.math.MathContext.DECIMAL64)
    coverage = coverage.stripTrailingZeros()
    // nanosim seed cannot be > 2**32 -1
    Long used_seed = ((seed as Long) % (2**32 - 1))
    threads = Math.max(1, ((task.cpus ?: 1) as int))

    /**
    echo "sample_id: ${sample_id}" > debug_inputs.txt
    echo "genome_id: ${genome_id}" >> debug_inputs.txt
    echo "fasta_file: ${fasta_file}" >> debug_inputs.txt
    echo "abundance: ${abundance}" >> debug_inputs.txt
    echo "genome_size: ${genome_size}" >> debug_inputs.txt
    echo "read_length_ch: ${read_length_ch}" >> debug_inputs.txt
    echo "seed: ${seed}" >> debug_inputs.txt
    echo "used_seed: ${used_seed}" >> debug_inputs.txt
    echo "total_size: ${total_size}" >> debug_inputs.txt
    echo "coverage: ${coverage}" >> debug_inputs.txt
    echo "profile: ${profile}" >> debug_inputs.txt
    **/

    """
    simulator.py genome -x ${coverage} -rg ${fasta_file} -o sample${sample_id}_${genome_id}_nanosim3 -c ${profile} --seed ${used_seed} -dna_type linear --fastq -t ${threads} -min ${params.nanosim3.min_read_length} -max ${safe_max}
    mkdir --parents ${params.outdir}/sample_${sample_id}/reads/fastq/nanosim3/
    for file in *_aligned_reads.fastq; do pigz -p ${threads} -k "\$file"; done
    cp *_aligned_reads.fastq.gz ${params.outdir}/sample_${sample_id}/reads/fastq/nanosim3/
    """
}

/**
* Tis process creates a bam file from the simulated reads in fasta format.
* In a second step it converts the fasta read file into a fq file.
* Input:
*     Tuple containing key = sample_id, first value = genome_id, second value = path to error profile, third value = path to fasta file with the aligned reads, 
*         fourth value = path to fasta file with the aligned reads, fifth value = path to reference genome.
* Output:
*    Tuple containing key = sample_id, first value = genome_id, second value = path to bam file, third value = path to reference genome.
**/
process bam_from_reads_fasta {

    conda "conda-forge::biopython bioconda::samtools conda-forge::pigz"

    input:
    tuple val(sample_id), val(genome_id), val(error_profile), path(aligned_reads), path(unaligned_reads), path(fasta_file)

    output:
    tuple val(sample_id), val(genome_id), path('sample*.bam'), path(fasta_file)
    tuple val(sample_id), path("sample${sample_id}_${genome_id}_nanosim3.fq")

    script:
    threads = Math.max(1, ((task.cpus ?: 1) as int))

    """
    ${shared_scripts_dir}/sam_from_reads.py ${error_profile} ${aligned_reads} ${unaligned_reads} ${fasta_file} --stdout | samtools view -bS | samtools sort -o sample${sample_id}_${genome_id}_nanosim3.bam
    mkdir --parents ${params.outdir}/sample_${sample_id}/bam/nanosim3/
    cp sample*.bam ${params.outdir}/sample_${sample_id}/bam/nanosim3/
    mkdir --parents ${params.outdir}/sample_${sample_id}/reads/fastq/nanosim3/
    for file in sample${sample_id}_${genome_id}_nanosim3.fq; do pigz -p ${threads} -k "\$file"; done
    cp sample${sample_id}_${genome_id}_nanosim3.fq.gz ${params.outdir}/sample_${sample_id}/reads/fastq/nanosim3/
    """
}

/**
* This process creates a bam file from the simulated reads in fastq format.
* Input:
*     Tuple containing key = sample_id, first value = genome_id, second value = path to error profile, third value = path to fasta file with the aligned reads, 
*         fourth value = path to fasta file with the aligned reads, fifth value = path to reference genome.
* Output:
*    Tuple containing key = sample_id, first value = genome_id, second value = path to bam file, third value = path to reference genome.
**/
process bam_from_reads_fastq {

    conda "conda-forge::biopython bioconda::samtools"

    input:
    tuple val(sample_id), val(genome_id), val(error_profile), path(aligned_reads), path(unaligned_reads), path(fasta_file)

    output:
    tuple val(sample_id), val(genome_id), path('sample*.bam'), path(fasta_file)

    script:
    """
    ${shared_scripts_dir}/sam_from_reads.py ${error_profile} ${aligned_reads} ${unaligned_reads} ${fasta_file} --fastq --stdout | samtools view -bS -F 4 | samtools sort -o sample${sample_id}_${genome_id}_nanosim3.bam
    mkdir --parents ${params.outdir}/sample_${sample_id}/bam/nanosim3/
    cp sample*.bam ${params.outdir}/sample_${sample_id}/bam/nanosim3/
    """
}

/*
* This process calculates the average read length of Nanosim reads from the pickle of the predefined profile
*
*/
process calculate_Nanosim_read_length {
    // TODO: Packages which are needed multiple times should be loaded only once
    // conda 'conda-forge::scikit-learn=0.23 conda-forge::numpy=1.23 conda-forge::joblib=1.2.0'
    conda 'conda-forge::scikit-learn=0.23.2 conda-forge::numpy=1.23.5 bioconda::nanosim=3.2'


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
