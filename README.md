[![Smoke Test](https://github.com/CAMI-challenge/CAMISIM/actions/workflows/ci.yml/badge.svg)](https://github.com/CAMI-challenge/CAMISIM/actions/workflows/ci.yml)
[![Unit Tests](https://github.com/CAMI-challenge/CAMISIM/actions/workflows/tests.yml/badge.svg)](https://github.com/CAMI-challenge/CAMISIM/actions/workflows/tests.yml)

# CAMISIM

CAMISIM simulates realistic shotgun metagenome and metatranscriptome datasets from microbial communities. It models abundance distributions, generates reads using established simulators, and produces gold standard assemblies and binning references for benchmarking.

Originally developed for the [Critical Assessment of Metagenome Interpretation (CAMI)](http://microbiome-cosi.org/cami) challenge, CAMISIM is suitable for general use. Please [open an issue](https://github.com/CAMI-challenge/CAMISIM/issues) if you run into problems.

## What's new in v2.0

CAMISIM 2.0 is a complete rewrite using [Nextflow](https://www.nextflow.io/) (DSL2), replacing the previous Python-based pipeline. Key changes:

- Nextflow-managed workflow with support for local, SLURM, and cloud executors
- Conda/Mamba-based dependency management
- Metatranscriptomic simulation pipeline (new)
- Configurable read simulators: ART (Illumina), NanoSim3 (Nanopore), WGSIM, pbsim3 (PacBio)

The previous standalone Python version is available via the git tag `1.31-final` but will not receive updates. Use `convert_config.py` to migrate v1 config files to v2 format.

## Quick start

### Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2)
- [Conda](https://docs.conda.io/) or [Mamba](https://mamba.readthedocs.io/) (for automatic environment management)

### Run a metagenomic simulation

```bash
nextflow run main.nf --pipeline metagenomic
```

### Run a metatranscriptomic simulation

```bash
nextflow run main.nf --pipeline metatranscriptomic
```

### Use a custom config

```bash
nextflow run main.nf --pipeline metagenomic --config path/to/your.config
```

Both pipelines ship with default example data under `nextflow_defaults/` so you can do a test run out of the box.

## Key parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `size` | `0.05` | Sample size in Gbp |
| `type` | `art` | Read simulator (`art`, `nanosim3`, `wgsim`, `pbsim3`) |
| `number_of_samples` | `2` | Number of samples to simulate |
| `gsa` | `true` | Generate gold standard assembly |
| `anonymization` | `true` | Anonymize read and contig IDs |
| `seed` | `632741178` | Random seed for reproducibility |
| `biom_profile` | `""` | BIOM profile for community design (optional) |

See [Configuration File Options](https://github.com/CAMI-challenge/CAMISIM/wiki/Configuration-File-Options) for the full parameter reference.

## Input files

| File | Format | Description |
|------|--------|-------------|
| `genome_locations.tsv` | TSV | `genome_ID<tab>path/to/genome.fa` |
| `meta_data.tsv` | TSV with header | `genome_ID<tab>OTU<tab>NCBI_ID<tab>novelty_category` |
| `distribution_N.txt` | TSV | `genome_ID<tab>abundance_fraction` (one per sample) |
| `gene_annotation_locations.tsv` | TSV | `genome_ID<tab>path/to/annotations.gff3` (metatranscriptomic only) |

## Project structure

```
CAMISIM/
  main.nf                          # Entry point
  nextflow.config                   # Global config (selects pipeline)
  pipelines/
    metagenomic/                    # Metagenomic simulation pipeline
      metagenomic.nf                #   Main workflow
      sample_wise_simulation.nf     #   Per-sample read simulation
      read_simulators/              #   ART, NanoSim3, WGSIM modules
      scripts/                      #   Python utilities
      config/                       #   Default parameters
    metatranscriptomic/             # Metatranscriptomic simulation pipeline
      metatranscriptomic.nf         #   Main workflow (+ gene expression)
      read_simulators/              #   ART, NanoSim3, pbsim3 modules
      scripts/                      #   Python utilities
      config/                       #   Default parameters
    shared/                         # Shared modules (both pipelines)
      anonymization.nf              #   Read/contig anonymization
      binning.nf                    #   Read-to-genome mapping
      distribution.nf               #   Abundance normalization
      scripts/                      #   Shared Python utilities
  nextflow_defaults/                # Example input data for test runs
  tools/                            # Bundled simulators and reference data
  tests/                            # Unit tests (pytest)
```

## Documentation

- [User manual](https://github.com/CAMI-challenge/CAMISIM/wiki/User-manual)
- [Configuration File Options](https://github.com/CAMI-challenge/CAMISIM/wiki/Configuration-File-Options)
- [File Formats](https://github.com/CAMI-challenge/CAMISIM/wiki/File-Formats)

## Citation

If you use CAMISIM, please cite:

> Fritz\*, Hofmann\*, et al. (2019). **CAMISIM: Simulating metagenomes and microbial communities.** *Microbiome*, 7:17. doi:[10.1186/s40168-019-0633-6](https://doi.org/10.1186/s40168-019-0633-6)

You may also cite the CAMI benchmark paper:

> Sczyrba\*, Hofmann\*, Belmann\*, et al. (2017). **Critical Assessment of Metagenome Interpretation — a benchmark of metagenomics software.** *Nature Methods*, 14(11):1063–1071. doi:[10.1038/nmeth.4458](https://doi.org/10.1038/nmeth.4458)

## License

Apache License 2.0. See [LICENSE.txt](LICENSE.txt).
