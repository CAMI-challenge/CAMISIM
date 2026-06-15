# CAMISIM 2.0 — User Manual

CAMISIM models the abundance distributions of microbial communities and simulates the
corresponding shotgun **metagenome** and **metatranscriptome** datasets (reads, alignments,
gold‑standard assemblies and taxonomic/binning truth). This version is a
[Nextflow](https://www.nextflow.io/) (DSL2) reimplementation.

This manual is a step‑by‑step operational guide written from a careful review of the
current code. It covers:

1. [Architecture & how it runs](#1-architecture--how-it-runs)
2. [Prerequisites & installation](#2-prerequisites--installation)
3. [Quick start](#3-quick-start)
4. [How configuration works](#4-how-configuration-works)
5. [Metagenomic pipeline — step by step](#5-metagenomic-pipeline--step-by-step)
6. [Metatranscriptomic pipeline — step by step](#6-metatranscriptomic-pipeline--step-by-step)
7. [Configuration parameter reference](#7-configuration-parameter-reference)
8. [Input file formats](#8-input-file-formats)
9. [Output files & directory layout](#9-output-files--directory-layout)
10. [FAQ](#10-faq)

> Throughout, `${projectDir}` is the repository root (the directory containing `main.nf`),
> regardless of where you launch Nextflow from. Code references are given as
> `file:line` so you can jump to the source.

---

## 1. Architecture & how it runs

```
nextflow run main.nf  --pipeline <metagenomic|metatranscriptomic>  [params...]
        │
        ├─ nextflow.config        # picks which *.config to load (or your --config)
        │
        └─ main.nf                # includes & calls the chosen sub‑workflow
              ├─ pipelines/metagenomic/metagenomic.nf
              │     ├─ community design  (from genome list OR from a BIOM profile)
              │     ├─ read simulation   (art_modern | art | nanosim3 | wgsim)
              │     ├─ gold standard assemblies (per‑sample, pooled, merged)
              │     └─ anonymization  OR  binning
              └─ pipelines/metatranscriptomic/metatranscriptomic.nf
                    ├─ genome‑level + gene‑level expression design
                    ├─ read simulation  (art_modern | art | nanosim3 | pbsim3)
                    └─ gold standard assemblies + anonymization/binning
```

* The single entry point is `main.nf`. It branches **only** on `params.pipeline`
  (`main.nf:6-17`). The value must be exactly `metagenomic` or `metatranscriptomic`;
  any other value runs an empty workflow that produces nothing (there is no error branch).
* Each process declares its own conda environment, so tools (samtools, art, nanosim,
  wgsim, …) are installed automatically on first run.

---

## 2. Prerequisites & installation

| Requirement | Notes |
|---|---|
| **Nextflow** ≥ 22 (DSL2) | Tested with Nextflow 25.x. Install per the [Nextflow docs](https://www.nextflow.io/docs/latest/getstarted.html). |
| **conda / mamba** | `conda` is **enabled by default** (`pipelines/shared/config/conda.config`). Channels are `conda-forge`, `bioconda`. Each process pulls its own env on first run. |
| **git clone** of this repo | All scripts, profiles, toy data and the sgEvolver binary are bundled. |

**The repo does NOT bundle simulator binaries.** `tools/` ships only *profiles / models /
reference tables*; the actual binaries (`art_illumina`, `art_modern`, NanoSim `simulator.py`,
`pbsim`, `wgsim`, `samtools`, …) are resolved by conda at runtime. What `tools/` contains:

* `tools/art_illumina-2.3.6/profiles/` — ART Illumina quality profiles (default
  `ART_MBARC-26_HiSeq_R1.txt` / `R2.txt`).
* `tools/nanosim_profile/` — NanoSim 3 pre‑trained models (default
  `nanosim323_r10_hc_lomanLab`).
* `tools/assembly_summary_complete_genomes.txt` — default reference‑genome catalog for
  BIOM‑profile mode.
* `tools/ncbi-taxonomy_20170222.tar.gz` (and `…_20180226.tar.gz`) — default NCBI taxdump.
* `scripts/sgEvolver/` — the sgEvolver strain‑simulation tool + template (`simulation_dir/`).

**Conda env caching (recommended on clusters).** Uncomment and set `conda.cacheDir` in
`pipelines/shared/config/conda.config` (or in your external config) to an existing
directory so envs are reused across runs.

---

## 3. Quick start

> ⚠️ **Read this first — the default `step` is not what you may expect.** Out of the box the
> metagenomic pipeline resolves `step` to **`reads_simulate`**, because the shipped
> `distribution_files` glob is non‑empty and `biom_profile` is empty
> (`pipelines/metagenomic/metagenomic.nf:46-47`). To run a *full* design + simulation with the
> defaults, pass `--step all` explicitly. See [§5.1](#51-execution-control-the-step-parameter).

```bash
# Full metagenomic run (community design + read simulation) with bundled toy data
nextflow run main.nf --pipeline metagenomic --step all

# Community design only (design the community, write the tables, then stop)
nextflow run main.nf --pipeline metagenomic --step community_design

# Read simulation only — reuse pre‑generated community‑design files
nextflow run main.nf --pipeline metagenomic --step reads_simulate \
    --distribution_files "/abs/path/distribution_*.txt" \
    --genome_locations_file /abs/path/genome_locations.tsv \
    --metadata_file /abs/path/meta_data.tsv

# Use a full external config (replaces the built‑in config — see §4)
nextflow run main.nf --pipeline metagenomic --config configs/metax_soil_SRP261862.config

# Resume a previous run (reuses cached work/) — deterministic when seed is fixed
nextflow run main.nf --pipeline metagenomic --step all --seed 632741178 -resume

# Metatranscriptomic run — see §6 for important caveats about the default simulator
nextflow run main.nf --pipeline metatranscriptomic
```

---

## 4. How configuration works

### 4.1 Config selection (`nextflow.config`)

`nextflow.config` defines two top‑level params: `pipeline` (default `metagenomic`) and
`config` (default empty). The include logic (`nextflow.config:37-49`) is:

* **If `--config <path>` is set** → only that file is included (`nextflow.config:37-39`).
  The path is resolved by `resolveExternalConfigPath` in this order:
  **absolute path (if it exists) → launch directory → `${projectDir}` → otherwise throw**
  (`nextflow.config:8-35`).
  **Important:** when you use `--config`, the built‑in config chain
  (`metagenomic.config` → `art_modern.config`/`distribution.config`/`conda.config`) is
  **bypassed**. Your external config must be self‑contained (define everything) or itself
  `includeConfig` the building blocks. The shipped example configs
  (`metax_bench_data.config`, `configs/metax_soil_SRP261862.config`) are full standalone
  configs.
* **Else if `pipeline == metagenomic`** → includes
  `pipelines/metagenomic/config/metagenomic.config`.
* **Else if `pipeline == metatranscriptomic`** → includes
  `pipelines/metatranscriptomic/config/metatranscriptomic.config`.

### 4.2 The built‑in metagenomic config chain

`pipelines/metagenomic/config/metagenomic.config` declares the core params and then
conditionally includes (`metagenomic.config:77-101`):

* one simulator config based on `params.type` — `art.config` / `art_modern.config` /
  `nanosim.config` / `wgsim.config`;
* `profile.config` **only if `biom_profile` is non‑empty**;
* always `distribution.config` and `conda.config`.

A consequence: `fragment_size_mean`, `fragment_size_sd`, `profile_read_length` and
`base_profile_name` live in the *simulator* configs, not the core config — only the matching
one is loaded.

### 4.3 Overriding parameters

Three ways, increasing precedence/scope:

1. **Edit the config files** in `pipelines/<pipeline>/config/`.
2. **Pass `-c extra.config`** — Nextflow *merges* this on top of the loaded config. Good for
   cluster resources / `workDir` (see `camisim.cfg`, `configs/camisim_soil_SRP261862.cfg`).
3. **`--config full.config`** — *replaces* the built‑in config entirely (§4.1).
4. **CLI `--param value`** or `-params-file params.yaml` — overrides individual params.

> `-c` **merges**; `--config` **replaces**. Don't confuse the two. A small resource‑only
> `.cfg` (just `workDir` + `withName` overrides) must be supplied with `-c`, not `--config`.

### 4.4 Reproducibility, seed and `-resume`

* `-resume` is the standard Nextflow flag — it reuses cached `work/` outputs whenever a
  task's inputs and script are unchanged.
* The master seed is `params.seed` if set, otherwise the `get_random_seed` process draws a
  random 32‑bit integer (`metagenomic.nf:82-86`, `metatranscriptomic.nf:26-30`).
* Both pipelines **default to a fixed seed** (`632741178`), so default runs are
  deterministic and fully resumable. To get a fresh random community each run, set
  `seed = null`.
* All downstream random seeds (per‑genome read seeds, per‑sample anonymization seeds, …)
  are derived deterministically from this single master seed by `get_seed.py`, so the
  *order* of derivation matters — see [§5.7](#57-seeds).

### 4.5 Executors / clusters

There are **no Nextflow `profiles {}` blocks** in this repo (so `-profile` does nothing).
The default executor is local. To run on SLURM, either uncomment the executor block at the
top of `metagenomic.config` / `metatranscriptomic.config`, or supply an external config with
a `process { executor='slurm'; memory=…; time=… }` block (as the example soil configs do).
The default per‑process wall time is `process.time = '2h'`.

---

## 5. Metagenomic pipeline — step by step

The metagenomic workflow (`pipelines/metagenomic/metagenomic.nf`) runs in two phases —
**community design** and **read simulation** — gated by `step`.

### 5.1 Execution control: the `step` parameter

`resolve_step()` (`metagenomic.nf:37-71`) determines what runs. Valid values:

| `step` | What runs |
|---|---|
| `community_design` | Design the community, write the tables, then **stop** (prints "Simulation stopping after community design steps."). |
| `reads_simulate` | **Skip** community design; simulate reads from pre‑generated files. **Requires** `distribution_files`, `genome_locations_file`, `metadata_file` (else the run aborts in `resolve_step`). |
| `all` | Community design **then** read simulation. |
| `""` (default) | **Auto‑detect** (legacy mapping). |

Auto‑detection when `step` is empty (`metagenomic.nf:44-54`):

* if `biom_profile` is empty **and** `distribution_files` is non‑empty → `reads_simulate`;
* else if `just_community_design == true` → `community_design`;
* else → `all`.

> Because the shipped default `distribution_files` is a non‑empty glob
> (`distribution.config:6`) and `biom_profile` is empty, the **out‑of‑the‑box default is
> `reads_simulate`**. Pass `--step all` for a full design+simulate, or clear
> `distribution_files`.

`just_community_design` is **deprecated** and honored only when `step` is empty.
The metatranscriptomic pipeline has **no** `step` mechanism — it always designs then simulates.

### 5.2 Community design — Path A: from a genome list (default)

Triggered when `biom_profile` is empty (`metagenomic.nf:112`). You supply:

* `genome_locations_file` — TSV `genome_id <TAB> path_to_fasta` (no header). Relative paths
  are resolved against `${projectDir}`.
* `metadata_file` — TSV **with header** `genome_ID <TAB> OTU <TAB> NCBI_ID <TAB> novelty_category`.

Flow:

1. Load both files into channels (`metagenomic.nf:115-117`).
2. **Optional strain simulation** when `genomes_total != genomes_real` (`metagenomic.nf:120`).
   Defaults are `2 == 2`, so **no strain simulation by default**. When enabled,
   `prepare_strain_simulation.py` decides how many strains per source genome (respecting
   `max_strains_per_otu`), then sgEvolver (`scripts/sgEvolver/simujobrun.pl` +
   `pick_random_strains.py`) synthesizes evolved strains using the
   `strain_simulation_template`. New strains are merged back into the genome‑location and
   metadata tables by `merge_metadata_files`.
3. **Abundances**: `getCommunityDistribution` runs `get_community_distribution.py`, which
   draws per‑sample relative abundances according to `mode` (see below) and writes
   `distribution_*.txt`.

**Abundance models (`mode`)** — `get_community_distribution.py`:

| `mode` | Behaviour |
|---|---|
| `differential` (default) | Each sample drawn **independently** from `lognormal(log_mu, log_sigma)`. |
| `replicates` | Sample 0 lognormal; later samples = sample 0 + `gauss(gauss_mu, gauss_sigma)` noise. |
| `timeseries_normal` | Sequential Gaussian random walk; a genome that hits 0 goes extinct. |
| `timeseries_lognormal` | Sample n = mean of (sample n‑1, a fresh lognormal draw). |

Notes: `log_*`/`gauss_*` are **cast to int**, so use integer values. An invalid `mode`
string silently yields sample‑0‑only abundances (no error). Code review indicates an
off‑by‑one in how per‑sample files are indexed in `write_distribution_file`
(`distribution_0.txt` is populated from the last sample column); this does not affect the
*set* of abundances, only the file‑index labelling.

### 5.3 Community design — Path B: from a BIOM profile

Triggered when `biom_profile` is set (`metagenomic.nf:100-109`, subworkflow
`metagenomesimulation_from_profile` in `from_profile.nf`). `get_genomes.py`:

* Reads the BIOM table (`defaults/mini.biom` is a toy example). Observation (row) ids are
  OTUs; each must carry `metadata.taxonomy` as a lineage string `k__…;p__…;…;s__SciName`.
  Columns are samples; cell values are abundances.
* Matches each OTU lineage to NCBI taxids, then draws genomes for it from the
  `reference_genomes` catalog (and optional `additional_references`), respecting
  `min_strains_per_otu`/`max_strains_per_otu`, `no_replace`, `max_rank`, and quality
  thresholds for additional genomes.
* Downloads/copies the selected genome FASTAs, splits the per‑OTU abundance across drawn
  strains (lognormal split using **`gauss_mu`/`gauss_sigma`**), and writes
  `genome_to_id.tsv`, `metadata.tsv`, `genome_filename_mapping.tsv` and per‑sample
  `abundance_*.tsv` (the distribution files for this mode).

`profile.config` (loaded only in this mode) holds all the selection parameters — see §7.

> In BIOM mode `min_strains_per_otu` is used; in genome‑list mode only `max_strains_per_otu`
> is used. `gauss_mu`/`gauss_sigma` double as the within‑OTU lognormal split parameters here.

### 5.4 Taxonomy & sequence cleanup (both paths)

* `download_NCBI_taxdump` runs only if `ncbi_taxdump_file` is empty (default points to the
  bundled `tools/ncbi-taxonomy_20170222.tar.gz`).
* `buildTaxonomy` → `build_ncbi_taxonomy.py` produces one CAMI‑format
  `taxonomic_profile_<i>.txt` per sample (`metagenomic.nf:540-562`).
* `cleanup_and_filter_sequences` (`metagenomic.nf:713-747`) uniquifies/filters every input
  sequence id, copies the cleaned FASTAs to `source_genomes/`, and writes the **final**
  `internal/genome_locations.tsv`. It runs for `all` and `reads_simulate` (not for
  `community_design`, which returns earlier).

### 5.5 Read simulation

`sample_wise_simulation` (`sample_wise_simulation.nf`) dispatches on `params.type`:

| `type` | Tool | Reads | Coverage / read‑count basis |
|---|---|---|---|
| `art_modern` (default) | `art_modern` | paired‑end short | `fold_coverage = abundance × factor` |
| `art` | `art_illumina` (2016.06.05) | paired‑end short | `fold_coverage = abundance × factor` |
| `nanosim3` | NanoSim 3 `simulator.py genome` | single‑end long (ONT) | `coverage = size·1e9 · abundance / genome_size` |
| `wgsim` | `wgsim` | paired‑end short | `N_reads = round(size·1e9 · abundance / read_length)` |

Key mechanics:

* `size` is the **target yield per sample in Gbp** (default `0.05`).
* For `art`/`art_modern`, `get_multiplication_factor` computes a per‑sample coverage scaling
  factor from `size` and genome sizes; for `nanosim3`/`wgsim`, abundance is first
  **weighted by genome length** (`count_bases` + `normalise_abundance_to_size`).
* All abundances are renormalized to sum to 1 per sample; genomes whose normalized
  abundance is exactly `0.0` are dropped (`sample_wise_simulation.nf:66`).
* **Read length**: `profile_read_length` (default 150) for art/art_modern/wgsim. For
  nanosim3, read length is derived from the model by `calculate_Nanosim_read_length`
  ("this takes very long") and capped per genome by `safe_max = max_contig − 200` (with a
  floor of 50, raised to `max_contig − 10` when `max_contig − 200 < 50`);
  nanosim ignores `profile_read_length`.
* **ART quirk**: classic `art` truncates FASTA headers at the first space, so for
  `type=art` references are pre‑processed by `remove_spaces_from_reference_genome`. (An open
  TODO questions whether art_modern needs the same.)
* **wgsim**: error‑free by default (`base_error_rate=0`); mutation/indel rates hardcoded to
  0; `create_cigar=false` writes a trivial `<len>M` CIGAR.
* **nanosim3**: `simulate_fastq_directly=true` (default) emits FASTQ with real quality
  scores; `false` emits FASTA and synthesizes FASTQ. The seed is reduced `mod 2^32−1`.

Example commands (per genome): ART
`art_illumina -sam -na -i <fa> -l <len> -m <mean> -s <sd> -f <fcov> -qL 15 -p -o … -1 <prof>1.txt -2 <prof>2.txt -rs <seed>`;
wgsim `wgsim -d <mean> -s <sd> -N <N> -1 <len> -2 <len> -S <seed> -e <err> -r 0 -R 0 …`.

### 5.6 Gold standard assemblies (GSA)

A GSA is built by `bamToGold.py`, which runs `samtools mpileup -B -Q 0 -f <ref> <bam>` and
emits the **reference bases covered by ≥ 1 read** (`-l 1 -c 1`), broken into contigs at
coverage gaps. FASTA headers are `>{seqname}_from_{start}_to_{stop}_total_{len}`. There are
three scopes:

| Scope | Controlled by | Output |
|---|---|---|
| **Per‑sample** (per genome, then concatenated) | `gsa` (default true) | `sample_<id>/gsa/…_gsa.fasta.gz` and `sample_<id>/contigs/gsa.fasta.gz` |
| **Pooled** (across all/selected samples) | `pooled_gsa` (`true` / `[ids]` / `[]`) | `pooled_gsa/gsa_pooled.fasta.gz` |
| **Merged** (custom groups, e.g. one per patient) | `merged_gsa_combinations` (list of lists) | `merged_gsa/gsa_merged_samples_<ids>.fasta.gz` |

* `pooled_gsa = true` pools all samples; a list `[0,1]` pools only those; `[]` skips. (Note:
  an empty list is falsy in the `if (params.pooled_gsa)` guard, so `[]` truly disables it.)
* `merged_gsa_combinations = [[0,1],[2,3]]` produces one merged GSA per inner list
  (each inner list = one "patient"). Sample ids are validated up front against
  `0..number_of_samples-1`; an out‑of‑range or malformed entry aborts the run with a clear
  error (`metagenomic.nf:271-286`). The merged GSA's covered positions equal the union of
  the member samples' coverage.

### 5.7 Anonymization vs binning

Both run **only inside `if (params.pooled_gsa)`** and are **mutually exclusive**
(`metagenomic.nf:311-322`):

* `anonymization = true` (default) → the `anonymization` subworkflow shuffles + renames
  reads and GSA contigs (deterministic, seeded with openssl AES‑256‑CTR keystream) and emits
  gold‑standard mapping files. If `merged_gsa_combinations` is set, merged GSAs are also
  anonymized (`anonymize_merged_gsa`).
* `anonymization = false` → the `binning` subworkflow produces CAMI binning ground‑truth
  mappings (`@@SEQUENCEID  BINID  TAXID  _LENGTH`) instead.

Anonymous prefixes: reads `S<sample>R`, per‑sample contigs `S<sample>C`, pooled `PC`, merged
`M<combination>C`. Note: the *anonymized* pooled outputs are written at the **top level** of
`outdir` (`anonymous_gsa_pooled.fasta.gz`, `gsa_pooled_mapping.tsv.gz`), whereas the
non‑anonymized pooled outputs go under `pooled_gsa/`.

### 5.8 Seeds

`get_seed.py` (process `get_seed`, `metagenomic.nf:788-819`) derives all seeds from the
master seed and writes them to `outdir/seed/`:

| File | When | Content |
|---|---|---|
| `seed.txt` | always | per‑`(genome, sample)` read‑simulation seeds (2 header lines + rows). |
| `seed_read_anonymisation.txt` | anonymization on | one seed per sample. |
| `seed_gsa_anonymisation.txt` | anonymization on | one seed per sample. |
| `seed_pooled_gsa_anonymisation.txt` | anonymization on | one seed. |
| `seed_merged_gsa_anonymisation.txt` | anonymization on **and** `merged_gsa_combinations` non‑empty | one seed per combination. |

The merged‑GSA seed block is generated **last**, so adding `merged_gsa_combinations` to an
existing run does **not** change any earlier seed — it only appends new ones (relevant for
`-resume`). Because seeds feed downstream processes as *parsed values* (not as file inputs),
re‑running `get_seed` with identical content does not invalidate cached read‑simulation tasks.

---

## 6. Metatranscriptomic pipeline — step by step

The metatranscriptomic workflow (`pipelines/metatranscriptomic/metatranscriptomic.nf`) reuses
the metagenomic scaffolding (seeds, GSA, pooled GSA, anonymization/binning) but adds a
**gene/transcript expression layer**: reads are simulated proportional to *per‑gene
expression* over *spliced transcript sequences* extracted from a GFF3, then re‑mapped back to
**genome coordinates** so the BAM/GSA remain in genome space.

> ⚠️ **The metatranscriptomic pipeline is not fully production‑ready.** The default
> short‑read path is currently broken (verified below). The long‑read paths (`nanosim3`,
> `pbsim3`) are the most complete but require externally installed models. See
> [§6.5 Known limitations](#65-known-limitations--needs-optimization). Treat this part as
> experimental and verify outputs.

### 6.1 Inputs (what's different from metagenomic)

* `genome_locations_file` — TSV `genome_id <TAB> fasta_path` (as in metagenomic).
* **`gene_annotations_file`** — TSV `genome_id <TAB> gff3_path` (the defining extra input).
  genome ids must match `genome_locations_file`.
* `metadata_file` — TSV with header `genome_ID OTU NCBI_ID novelty_category`.
* The GFF3 must contain features of type `feature_type` (default `mRNA`) and child features
  of type `child_feature_type` (default `exon`; `CDS` is also supported, with phase
  handling). There is **no** `step` parameter and **no** BIOM‑profile mode.

Toy data lives under `nextflow_defaults/metatranscriptomic/data/` (6 organisms with paired
`.fa` + `.gff3`), but the default `genome_locations.tsv` / `gene_annotation_locations.tsv`
reference only **3** of them, so a default run simulates 3 genomes.

### 6.2 Design: genome‑level × gene‑level abundance

1. **Genome‑level distribution** — same machinery as metagenomic, using the `genome_*`
   parameters (`genome_mode`, `genome_log_mu`, …). Writes
   `distributions/genome_distributions/distribution_<i>.txt`. Pre‑supplied
   `genome_distribution_files` are used if non‑empty.
2. **Gene‑level distribution** — `distribute_gene_abundance` → `get_gene_abundance.py` builds
   a gffutils DB per genome, assigns each gene a per‑sample abundance (`mode` =
   `differential`/`timeseries`/`replicate`, using `mu`/`sigma`/`gauss_*`/`gene_sigma`),
   weights by gene length and normalizes. Writes
   `distributions/gene_distributions/distribution_<genome>_sample_<i>.tsv`.
3. **Final** — `get_final_gene_distr` multiplies gene abundance × genome abundance →
   `distributions/final_distributions/<genome>_<sample>_final_distribution.tsv`.

### 6.3 Read simulation

`type` selects `art_modern` (default), `art`, `nanosim3`, or `pbsim3` (note: **no wgsim**;
`pbsim3` is metatranscriptomic‑only). For each simulator the spliced transcript is extracted
(`sequence_extractor.py`), reads are simulated per gene, and the alignment coordinates are
shifted from transcript space to genome space:

* **art / art_modern** — `generate_reads_and_modify_sam.py` emits per‑gene `art`/`art_modern`
  commands plus an `awk` that shifts SAM POS/PNEXT by the gene's genomic start.
* **nanosim3** — Trans‑NanoSim transcriptome mode; coordinates remapped via
  `transcript_id_to_seq_id.tsv` in `sam_from_reads.py --transcriptome`. Requires an
  externally installed Trans‑NanoSim model.
* **pbsim3** — `pbsim --strategy trans`; `maf_converter.py --transcriptome` remaps MAF
  coordinates to the genome. Requires an externally installed PBSIM3 model.

### 6.4 Outputs

Same shape as metagenomic per‑sample/pooled outputs, **plus** the three distribution
sub‑folders above. There is **no** `merged_gsa` and **no** `internal/`/`source_genomes/`
(genomes are user‑supplied). Default `outdir` is `${projectDir}/out_mtx`.

### 6.5 Known limitations / needs optimization

These are concrete issues found by code review (file:line given) — the metatranscriptomic
path needs work before it can be relied on:

1. **Default short‑read type is broken (verified).** `read_simulator_art_modern.nf:13`
   declares `workflow read_simulator_art` (copy‑pasted from the art module) but
   `sample_wise_simulation.nf:11` imports it as `read_simulator_art_modern`. The named
   workflow does not exist in that file, so a default run (`type=art_modern`) fails at
   include/execution time.
2. **`--simulator` argument never passed (verified).**
   `generate_reads_and_modify_sam.py:24` makes `--simulator` **required**, but neither
   `read_simulator_art.nf` nor `read_simulator_art_modern.nf` passes it (their command blocks
   end at `--db`). Both short‑read paths therefore error with
   `the following arguments are required: --simulator`.
3. **art_modern conda env lacks samtools.** `read_simulator_art_modern.nf:26` declares
   `bioconda::art_modern bioconda::gffutils bioconda::bedtools bioconda::pyfaidx`, but the
   process body calls `samtools merge/view/sort` — the BAM step would fail even if (1) and
   (2) were fixed.
4. **External, home‑relative model paths.** `nanosim.config` points
   `base_profile_name` at `~/NanoSim/pre-trained_models/…` and `pbsim3.config` points `model`
   at `~/pbsim3_data/data/QSHMM-RSII.model`. These are **not bundled**, and the `~` is passed
   verbatim into the tool (not reliably expanded). You must install the models and set
   absolute paths.
5. **NanoSim reproducibility / approximate counts.** The config comments state nanosim3 reads
   are not yet reproducible and read counts are approximate. The seed reduction here is
   `(seed as Long) % 2**32 - 1` (`read_simulator_nansoim3.nf:68`, `read_simulator_pbsim3.nf:44`),
   which by Groovy precedence is `(seed mod 2^32) − 1` — *not* `mod (2^32−1)` as the parenthesized
   metagenomic version does — so the effective seed is silently shifted (e.g. 632741178 → 632741177).
6. **Fragile genome‑coordinate remapping.** The transcript→genome SAM shift is an `awk`
   one‑liner appended to a generated shell script (`generate_reads_and_modify_sam.py:79`);
   there's no validation, and gene‑id sanitization for filenames only handles `/`, `(`, `)`.
7. **gffutils version mismatch.** The DB is built with `gffutils=0.12`
   (`metatranscriptomic.nf:132`) but read with `gffutils=0.9` in the art/nanosim modules —
   cross‑version sqlite schema reads may misbehave.
8. **`sequence_extractor.py` assumes well‑formed GFF3.** Only `+`/`-` strands are handled —
   a `.`/unset strand leaves the result undefined and raises. (CDS phase parsing is safe: a
   non‑digit phase falls back to 0.)
9. **Minor**: a typo'd `scripts_dir` (`…/metatransciptomic/…`) in
   `sample_wise_simulation.nf:5` (currently harmless), several commented‑out `publishDir`
   blocks replaced by manual `cp` (bypasses Nextflow's resume guarantees), and an unused
   `normalise_abundance_to_size` import.

**Practical guidance today:** if you must run the metatranscriptomic pipeline, the
long‑read paths (`--type nanosim3` or `--type pbsim3`) are the most complete, *after* you
install the required models and set absolute model paths. The short‑read (`art`/`art_modern`)
paths need code fixes (items 1–3) first.

---

## 7. Configuration parameter reference

Defaults shown are the shipped values. "Config file" is where the parameter is declared.

### 7.1 Top‑level (`nextflow.config`)

| Param | Default | Meaning |
|---|---|---|
| `pipeline` | `metagenomic` | `metagenomic` or `metatranscriptomic`. Selects the workflow and the config. |
| `config` | `""` | Optional external config path. When set, **replaces** the built‑in config (resolved absolute → launch dir → projectDir). |

### 7.2 Metagenomic core (`pipelines/metagenomic/config/metagenomic.config`)

| Param | Default | Meaning |
|---|---|---|
| `outdir` | `${projectDir}/nextflow_out` | Output directory (created if missing). |
| `size` | `0.05` | Target yield **per sample in Gbp**. Scales coverage / read counts. |
| `type` | `art_modern` | Read simulator: `art_modern` / `art` / `nanosim3` / `wgsim`. |
| `number_of_samples` | `2` | Number of samples; valid sample ids are `0..N-1`. |
| `gsa` | `true` | Build per‑genome/per‑sample gold‑standard assemblies. |
| `pooled_gsa` | `true` | `true` (all) / list of ids (subset) / `[]` (skip). Also gates anonymization/binning. |
| `merged_gsa_combinations` | `[]` | List of lists of sample ids; one merged GSA per inner list. |
| `anonymization` | `true` | `true` → anonymize; `false` → binning ground truth. |
| `seed` | `632741178` | Master RNG seed; `null` → random 32‑bit seed. |
| `biom_profile` | `""` | If set, design from a BIOM profile (and load `profile.config`). |
| `genome_locations_file` | `…/nextflow_defaults/metagenomic/genome_locations.tsv` | genome_id→FASTA map (genome‑list mode). |
| `metadata_file` | `…/nextflow_defaults/metagenomic/meta_data.tsv` | per‑genome metadata (with header). |
| `ncbi_taxdump_file` | `${projectDir}/tools/ncbi-taxonomy_20170222.tar.gz` | NCBI taxdump; empty → auto‑download. |
| `max_strains_per_otu` | `2` | Max strains drawn per OTU. |
| `min_strains_per_otu` | `1` | Min strains per OTU (**BIOM mode only**). |

### 7.3 Distribution / steps (`pipelines/metagenomic/config/distribution.config`)

| Param | Default | Meaning |
|---|---|---|
| `step` | `""` | `community_design` / `reads_simulate` / `all` / `""` (auto). See §5.1. |
| `just_community_design` | `false` | **Deprecated**; honored only when `step==""`. |
| `distribution_files` | `…/distribution_*.txt` | Pre‑generated per‑sample abundance files (glob/list). Non‑empty default ⇒ default step = `reads_simulate`. |
| `mode` | `differential` | Abundance model: `differential` / `replicates` / `timeseries_normal` / `timeseries_lognormal`. |
| `log_mu` | `1` | Lognormal mean (int). |
| `log_sigma` | `2` | Lognormal sd (int). |
| `gauss_mu` | `1` | Gaussian noise mean (int). Also lognormal `mu` for within‑OTU split in BIOM mode. |
| `gauss_sigma` | `1` | Gaussian noise sd (int). Also lognormal `sigma` in BIOM mode. |
| `verbose` | `false` | Keep `false` — `true` puts abundance generation in an infinite loop (`get_community_distribution.py:275-288`). |
| `genomes_total` | `2` | Total genomes after strain simulation; `!= genomes_real` triggers sgEvolver (must be `>= genomes_real`; equality = no strain simulation). |
| `genomes_real` | `2` | Number of real (un‑evolved) input genomes. |
| `id_to_gff_file` | `""` | Optional genome_id→GFF map; selects the with‑GFF strain path. |
| `strain_simulation_template` | `${projectDir}/scripts/sgEvolver/simulation_dir` | sgEvolver template dir. |

### 7.4 Simulator configs (loaded per `type`)

`art.config` / `art_modern.config` (`metagenomic`):

| Param | Default | Meaning |
|---|---|---|
| `profile_read_length` | `150` | Read length. |
| `fragment_size_mean` | `270` | Mean paired‑end fragment size. |
| `fragment_size_sd` | `27` | Fragment size sd. |
| `base_profile_name` | `${projectDir}/tools/art_illumina-2.3.6/profiles/ART_MBARC-26_HiSeq_R` | Quality‑profile prefix; `1.txt`/`2.txt` are appended. |

`nanosim.config` (`metagenomic`): `base_profile_name = …/tools/nanosim_profile/nanosim323_r10_hc_lomanLab`,
`simulate_fastq_directly = true`.

`wgsim.config`: `profile_read_length=150`, `fragment_size_mean=270`, `fragment_size_sd=27`,
`base_error_rate=0`, `create_cigar=false` (no `base_profile_name`).

### 7.5 BIOM‑profile selection (`profile.config`, loaded only when `biom_profile` set)

| Param | Default | Meaning |
|---|---|---|
| `reference_genomes` | `${projectDir}/tools/assembly_summary_complete_genomes.txt` | Candidate genome catalog (`taxid  name  ftp/path`). |
| `no_replace` | `true` | Sample genomes without replacement. |
| `fill_up` | `false` | Assign leftover genomes to unmatched OTUs. |
| `additional_references` | `""` | Extra genomes (requires a quality file). |
| `prioritize_additional_genomes` | `false` | Prefer high‑quality additional genomes. |
| `additional_genomes_quality_file` | `""` | Quality table for additional genomes. |
| `additional_genomes_max_contamination` | `5` | Quality threshold. |
| `additional_genomes_min_completeness` | `90` | Quality threshold. |
| `additional_genomes_max_num_contigs` | `200` | Quality threshold. |
| `max_rank` | `family` | Highest taxonomic rank to draw genomes from. |

### 7.6 Metatranscriptomic (`pipelines/metatranscriptomic/config/…`)

Core (`metatranscriptomic.config`): `outdir=${projectDir}/out_mtx`, `size=0.05`,
`type=art_modern` (`art_modern`/`art`/`pbsim3`/`nanosim3`), `feature_type=mRNA`,
`child_feature_type=exon`, `number_of_samples=2`, `seed=632741178`,
`genome_locations_file`, `gene_annotations_file`, `metadata_file`, `gsa`, `pooled_gsa`,
`anonymization`. (No `step`, no `merged_gsa_combinations`, no `biom_profile`.)

Gene/genome distribution (`distribution.config`): gene‑level `mode=differential`, `mu=1`,
`sigma=2`, `gauss_mu=0`, `gauss_sigma=1`, `gene_sigma=1`; genome‑level
`genome_distribution_files`, `genome_mode=differential`, `genome_log_mu=1`,
`genome_log_sigma=2`, `genome_gauss_mu=1`, `genome_gauss_sigma=1`.

Simulators: `art.config`/`art_modern.config` (`read_length=150`, fragment params, ART
profile); `nanosim.config` (external model prefix, `read_length=1125`, `basecaller=guppy`);
`pbsim3.config` (`method=qshmm`, external `model`, `difference_ratio=6:55:39`,
`fragments_size_mean=1500`, `read_length=1000`, `fragment_size_standard_deviation=500`).

---

## 8. Input file formats

| Input | Format |
|---|---|
| `genome_locations_file` | TSV, **no header**, 2 columns: `genome_id <TAB> fasta_path`. Relative paths resolved against `${projectDir}`. |
| `metadata_file` | TSV, **with header**: `genome_ID <TAB> OTU <TAB> NCBI_ID <TAB> novelty_category`. `novelty_category=plasmid` routes to a plasmid lineage. |
| `distribution_files` | One TSV per sample, **no header**: `genome_id <TAB> relative_abundance`. |
| `biom_profile` | BIOM table (JSON BIOM 1.0 or HDF5). Rows = OTUs with `metadata.taxonomy` lineage `k__…;…;s__SciName`; columns = samples. |
| `ncbi_taxdump_file` | gzipped tar of NCBI taxdump (`names.dmp`, `nodes.dmp`, `merged.dmp`). |
| `reference_genomes` (BIOM mode) | TSV: `NCBI_taxid <TAB> scientific_name <TAB> genome_path_or_ftp`. Header optional (alias‑detected). |
| `additional_references` (BIOM mode) | TSV: `taxid <TAB> name <TAB> path [<TAB> novelty]`; needs a quality file with `genome_path, completeness, contamination, num_contigs`. |
| `id_to_gff_file` (strain sim) | TSV: `genome_id <TAB> gff_path`. |
| `gene_annotations_file` (metatranscriptomic) | TSV: `genome_id <TAB> gff3_path`. GFF3 must contain `feature_type` and `child_feature_type` features. |

---

## 9. Output files & directory layout

### 9.1 Metagenomic (`outdir` default `${projectDir}/nextflow_out`)

```
<outdir>/
├── internal/
│   ├── genome_locations.tsv          # final genome_id→FASTA map (after cleanup)
│   ├── meta_data.tsv                 # ONLY when strain simulation runs (genomes_total != genomes_real)
│   ├── metadata.tsv                  # BIOM mode
│   ├── genome_to_id.tsv              # BIOM mode
│   ├── genome_filename_mapping.tsv   # BIOM mode
│   └── genomes/                      # downloaded taxdump (only if ncbi_taxdump_file empty)
├── source_genomes/                   # cleaned / downloaded reference FASTAs
├── distributions/
│   ├── distribution_<i>.txt          # per-sample normalized abundances (genome-list mode)
│   └── abundance_<i>.tsv             # per-sample abundances (BIOM mode)
├── seed/
│   ├── seed.txt
│   ├── seed_read_anonymisation.txt          # anonymization=true
│   ├── seed_gsa_anonymisation.txt           # anonymization=true
│   ├── seed_pooled_gsa_anonymisation.txt    # anonymization=true
│   └── seed_merged_gsa_anonymisation.txt    # anonymization=true AND merged_gsa_combinations set
├── taxonomic_profile_<i>.txt         # one CAMI taxonomic profile per sample
├── sample_<i>/
│   ├── bam/   sample<i>_<genome>.bam            # per-genome alignments
│   ├── reads/
│   │   ├── fastq/   sample<i>_<genome>*.fq.gz   # per-genome reads
│   │   │            sample_<i>_01.fq.gz / _02.fq.gz   # merged paired (art/art_modern/wgsim)
│   │   │            sample_<i>.fq.gz                  # merged single-end (nanosim3)
│   │   ├── anonymous_reads.fq.gz       # anonymization=true
│   │   └── reads_mapping.tsv.gz        # anonymization=true
│   ├── gsa/      sample<i>_<genome>_gsa.fasta.gz   # per-genome GSA (gsa=true)
│   └── contigs/
│       ├── gsa.fasta.gz                # per-sample combined GSA
│       ├── anonymous_gsa.fasta.gz      # anonymization=true
│       └── gsa_mapping.tsv.gz          # anonymization (gs_contig_mapping) OR binning
├── pooled_gsa/
│   ├── gsa_pooled.fasta.gz
│   └── pooled_gsa_mapping.tsv.gz       # binning path (anonymization=false)
├── merged_gsa/                         # only if merged_gsa_combinations set
│   ├── gsa_merged_samples_<ids>.fasta.gz
│   ├── anonymous_gsa_merged_samples_<ids>.fasta.gz   # anonymization=true
│   └── gsa_mapping_merged_samples_<ids>.tsv.gz       # anonymization=true
├── anonymous_gsa_pooled.fasta.gz       # TOP LEVEL, anonymization=true
└── gsa_pooled_mapping.tsv.gz           # TOP LEVEL, anonymization=true
```

Mapping‑file columns: `reads_mapping.tsv` =
`#anonymous_read_id, genome_id, tax_id, read_id`; contig `gsa_mapping.tsv` (anonymization
path) = `#anonymous_contig_id, genome_id, tax_id, contig_id, number_reads, start_position,
end_position`; binning `gsa_mapping.tsv` = `@@SEQUENCEID  BINID  TAXID  _LENGTH`.

### 9.2 Metatranscriptomic (`outdir` default `${projectDir}/out_mtx`)

Same `seed/`, `sample_<i>/…`, `pooled_gsa/`, and top‑level anonymized pooled outputs as
metagenomic, **plus**:

```
<outdir>/distributions/
├── genome_distributions/distribution_<i>.txt
├── gene_distributions/distribution_<genome>_sample_<i>.tsv
└── final_distributions/<genome>_<sample>_final_distribution.tsv
```

No `merged_gsa/`, `internal/` or `source_genomes/`.

### 9.3 Files to reuse for a `reads_simulate` re‑run

To re‑simulate reads from a previous design, point `step=reads_simulate` at:
`internal/genome_locations.tsv`, `distributions/distribution_*.txt` (or
`distributions/abundance_*.tsv` in BIOM mode), and a metadata file. For metadata use
`internal/metadata.tsv` (BIOM mode) or `internal/meta_data.tsv` **only if strain simulation
ran** (`genomes_total != genomes_real`); in the default genome‑list path that file is not
written, so reuse your **original** `metadata_file` instead.

---

## 10. FAQ

**Q1. I ran `--pipeline metagenomic` with the defaults and it skipped community design — why?**
Because `step` auto‑resolves to `reads_simulate` out of the box: `biom_profile` is empty and
the default `distribution_files` glob is non‑empty (`metagenomic.nf:46-47`). Use
`--step all` to design + simulate, or clear `distribution_files`.

**Q2. How do I run only the community design, or only the read simulation?**
`--step community_design` stops after writing `internal/`, `source_genomes/`,
`distributions/`, `seed/` and the taxonomic profiles. `--step reads_simulate` skips design and
**requires** `distribution_files`, `genome_locations_file`, `metadata_file` (else it aborts).

**Q3. Nothing happened / no outputs and no error.**
`params.pipeline` was probably not exactly `metagenomic` or `metatranscriptomic` — `main.nf`
has no error branch and silently runs an empty workflow for any other value.

**Q4. How do I make runs reproducible? Does `-resume` work?**
Keep `seed` fixed (the default `632741178` is fixed). All derived seeds come from it
deterministically, so `-resume` reuses cached tasks. Set `seed = null` for a fresh random
community each run (which breaks resume‑equality across runs).

**Q5. I added `merged_gsa_combinations` and re‑ran with `-resume`; `get_seed` re‑ran but
everything else was cached. Is that a problem?**
No. `get_seed`'s command embeds the combination count, so changing it forces `get_seed` to
re‑run — but the merged‑seed block is appended *last*, so `seed.txt` and the other seed files
are byte‑identical. Seeds reach downstream processes as parsed values, so cached
read‑simulation/GSA tasks stay valid. Only the new merged GSAs are computed.

**Q6. What does `pooled_gsa` accept, and how is it different from `merged_gsa_combinations`?**
`pooled_gsa` makes **one** GSA over all samples (`true`), a subset (`[0,1]`), or none (`[]`).
`merged_gsa_combinations` is a **list of lists** that makes one GSA *per inner list* (e.g.
one per patient): `[[0,1],[2,3]]`. Sample ids must be in `0..number_of_samples-1`, else the
run errors immediately.

**Q7. A merged/pooled GSA looks empty or smaller than expected.**
A GSA contains only reference bases **covered by ≥ 1 read** (`bamToGold.py`, `-l 1 -c 1`).
Low coverage (small `size`, low abundance) yields little or no GSA. The merged GSA's covered
positions equal the *union* of the member samples' coverage.

**Q8. I set `pooled_gsa = []` to disable pooled GSA — did anonymization/binning still run?**
No. Anonymization and binning run only inside `if (params.pooled_gsa)`. With `[]` (or
`false`) that block is skipped, so no anonymized/binning mapping is produced — although
per‑sample and per‑genome GSA FASTAs are still written.

**Q9. Where did the anonymized pooled files go? They're not under `pooled_gsa/`.**
The *anonymized* pooled outputs are written at the **top level** of `outdir`
(`anonymous_gsa_pooled.fasta.gz`, `gsa_pooled_mapping.tsv.gz`). Only the non‑anonymized
`gsa_pooled.fasta.gz` (and binning `pooled_gsa_mapping.tsv.gz`) live under `pooled_gsa/`.

**Q10. How do I choose / change the read simulator and its error model?**
Set `type` (`art_modern` default, `art`, `nanosim3`, `wgsim` for metagenomic). The matching
simulator config is loaded automatically. Change `base_profile_name` to point at a different
ART profile prefix (the `<prefix>1.txt`/`2.txt` files must exist) or NanoSim model prefix
(all model files must exist). `wgsim` has no external profile — only `base_error_rate`.

**Q11. How is the amount of data controlled?**
`size` is the **target yield per sample in Gbp** (default `0.05`). It scales fold coverage
(art/art_modern), nanosim coverage, and wgsim read counts. Increase it for larger datasets.

**Q12. How do I add more genomes than I have inputs for (strain simulation)?**
Set `genomes_total > genomes_real`. sgEvolver synthesizes the extra strains
(`max_strains_per_otu` caps strains per OTU). With the defaults `2 == 2`, no strain
simulation runs. Note: the **with‑GFF** strain path (`id_to_gff_file` set) appears
under‑maintained relative to the no‑GFF path — prefer leaving `id_to_gff_file` empty unless
you've verified that branch.

**Q13. How do I supply my own community as a taxonomic profile instead of a genome list?**
Set `biom_profile` to a BIOM file. CAMISIM then selects genomes from `reference_genomes`
matching each OTU lineage. This loads `profile.config` (selection params). The example soil
configs (`metax_bench_data.config`, `configs/metax_soil_SRP261862.config`) show full
profile‑mode setups.

**Q14. Do I need internet access? What about the NCBI taxonomy?**
Conda needs network on first run to build envs. If `ncbi_taxdump_file` is empty, CAMISIM
downloads a fresh taxdump via ete3 (network needed); otherwise it uses the bundled
`tools/ncbi-taxonomy_20170222.tar.gz` (offline). BIOM‑mode genome downloads also need network
unless your `reference_genomes` paths are local.

**Q15. How do I run on a SLURM cluster?**
Add a `process { executor='slurm'; memory=…; time=… }` block (uncomment it in the config, or
use an external config), and tune per‑process resources with a `-c resources.cfg`
(`withName:` overrides, like `configs/camisim_soil_SRP261862.cfg`). Set `conda.cacheDir` to a
persistent path so envs are reused. The default `process.time` is `2h` — raise it for large
or long‑read runs.

**Q16. What's the difference between `-c` and `--config`?**
`-c file` *merges* on top of the loaded config (use it for resources/`workDir`). `--config
file` (i.e. `params.config`) *replaces* the entire built‑in config chain, so that file must
define everything. A small resource‑only `.cfg` must be passed with `-c`.

**Q17. Why are some sample/genome ids dropped from simulation?**
Genomes whose normalized abundance is exactly `0.0` are filtered out before simulation
(`sample_wise_simulation.nf:66`). This is expected for very rare genomes under some
distribution settings.

**Q18. The NanoSim step "takes very long" before reads are produced.**
For `nanosim3`, the read length is computed from the trained model by
`calculate_Nanosim_read_length` (integration over the model's KDE), which is slow — this is
expected. The actual read length comes from the model, not `profile_read_length`.

**Q19. The metatranscriptomic run fails immediately — is it broken?**
The **default short‑read path (`type=art_modern`) is currently broken** (a workflow‑name
mismatch and a missing required `--simulator` argument; see §6.5). Until those are fixed, use
`--type nanosim3` or `--type pbsim3` *after* installing their models and setting absolute
model paths (the configs ship `~/...` placeholders that are not bundled).

**Q20. Which output do I feed to an assembler / binner as the "ground truth"?**
Per‑sample reads: `sample_<i>/reads/fastq/sample_<i>_*.fq.gz` (or anonymized
`anonymous_reads.fq.gz`). Gold‑standard assembly: per‑sample `contigs/gsa.fasta.gz`,
`pooled_gsa/gsa_pooled.fasta.gz`, or `merged_gsa/gsa_merged_samples_*.fasta.gz`. Truth
mappings: `reads_mapping.tsv.gz`, `gsa_mapping.tsv.gz`; taxonomic truth:
`taxonomic_profile_<i>.txt`.

**Q21. Can I convert an old CAMISIM 1 config?**
`convert_config.py` (repo root) converts a v1 `.ini` to Nextflow params (no correctness
guarantee). It hardcodes some values (including a `conda.cacheDir` path) you must edit
afterwards. See `defaults/CAMISIM1.3_config.ini` and `CAMI1_documentation/` for the legacy
format.

**Q22. `step=reads_simulate` aborted with "requires '…' to point to the pre‑generated
file(s)".**
That step needs all three of `distribution_files`, `genome_locations_file`, `metadata_file`
set to existing files (`metagenomic.nf:60-68`). Point them at a previous run's
`distributions/distribution_*.txt` and `internal/genome_locations.tsv`. For the metadata file,
reuse your original `metadata_file` (the default genome‑list path does **not** write
`internal/meta_data.tsv` unless strain simulation ran — see §9.3).

---

*Generated from a code review of the CAMISIM 2.0 Nextflow pipeline. For the upstream project
see the [CAMISIM GitHub wiki](https://github.com/CAMI-challenge/CAMISIM/wiki) and the
citation in `README.md`.*
