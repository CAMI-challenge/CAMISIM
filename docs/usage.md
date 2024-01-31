Contents
===========

* [Introduction](#introduction)  
* [Installation](#installation) 
	* [Dependencies](#dependencies)
* [Usage](#usage)

Introduction
===========
With CAMISIM you can create simulated metagenome data sets out of taxonomic profiles or *de novo* from a list of genomes. If a taxonomic profile is used as input, the output data set is created from the NCBI complete genomes, reflecting the input profile as closely as possible and will contain the same number of samples as the input profile, if not specified otherwise. If the community is designed *de novo*, a user-defined number of the complete genomes is used for creating a community which maximizes genome novelty as well as phylogenetic spread.

CAMISIM can generate four types of simulated metagenome samples *de novo*:
- Individual simulate metagenome sample.  
  Uses a taxonomic profile sampled from a log-normal distribution
- A time series of simulated metagenome samples.  
  Uses a taxonomic profile sampled from a log-normal distribution with Gaussian noise, a normal distribution added to successively generated samples.
- A set of replicate simulated metagenome samples.  
  Uses a taxonomic profile sampled from a log-normal distribution, with Gaussian noise repeatedly added to the original log-normal values.
- Differential abundance metagenome samples.  
  Uses a taxonomic profile sampled from a log-normal distribution.

Installation
===========
### From Source
Download CAMISIM from the [github repository](https://github.com/CAMI-challenge/CAMISIM) into any directory and make sure the following [[Dependencies|Dependencies]] are fullfilled.
The metagenome simulation pipeline can be downloaded using git:

    git clone https://github.com/CAMI-challenge/CAMISIM

# Dependencies
## Operation system

The pipeline has only been tested on UNIX systems and might fail to run on any other systems (and no support for these will be provided).

## Software
The following software is required for the simulation pipeline to run:

#### [nextflow](https://github.com/nextflow-io/nextflow/releases)

Nextflow is the workflow framework used for this software. Nextflow can also be installed via conda (https://anaconda.org/bioconda/nextflow).

#### [mamba](https://github.com/mamba-org/mamba)

CAMISIM uses mamba to install and configure all needed software packages. Mamba can also be installed via conda (https://anaconda.org/conda-forge/mamba).

#### [conda](https://docs.conda.io/en/latest/) optional

CAMISIM can also use conda instead of mamba to install and configure all needed software packages. This is not recommended for performance reasons. If it is still needed to use conda, adjust it in the config files ([[configuration file|Configuration-File-Options]]).

## Hardware

The simulation will be conducted primarily in the work directory created by the pipeline.
The results will then be copied to the specified output directory. Be sure to have enough space at both locations.
Required hard drive space and RAM can vary a lot. Very small simulations can be run on a laptop with 4GB RAM, while realistic simulations of several hundred genomes, given a realistic metagenome size, can easily require several hundreds of gigabyte of both RAM and HD space, the chosen metagenome size being the relevant factor.

## Resources

A database dump of the [NCBI](http://www.ncbi.nlm.nih.gov/guide/taxonomy/) taxonomy is included, current versions can be downloaded from the NCBI [FTP-Server](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/).

### Genomes
If the community design should be performed *de novo*, genomes in [fasta](https://en.wikipedia.org/wiki/FASTA) format to be sampled from are needed. Otherwise they will be downloaded from the NCBI complete genomes.

The *de novo* community design needs two files to run: 

- A file containing, tab separated, a genome identifier and that path to the file of the genome.

- A [[meta data file|meta-data-file-format] that contains, tab separated and with header, genome identifier, novelty categorization, otu assignment and a taxonomic classification.

### [NCBI taxdump](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)
At minimum the following files are required: "nodes.dmp", "merged.dmp", "names.dmp"

USAGE:
===========

_de novo_:

    >> nextflow run workflow.nf

or _from\_profile_:

    >> nextflow run workflow.nf --params.biom_profile "${projectDir}/defaults/mini.biom"
   
This uses the `nextflow.config` found in CAMISIM's main directory.

# Metagenome simulation 
Multiple kinds of simulated metagenome data sets can be created, either using an existing 16S rRNA profile or _de novo_: Single sample, differential abundance or time series data sets with 
different insert sizes and complexities. 
All data sets can be simulated as paired-end sequencing from Illumina HiSeq or other machines as well as other technologies (like ONT or PacBio) and error rates.
For generation of data sets and their samples, the following steps are undertaken:  
After genome data validation (step 1), the community composition is designed according to specified criteria or the given taxonomic profile (step 2) and the metagenome data sets simulated (step 3). 
Creation of the gold standards (step 4) represents the last part of the pipeline. 

## Step 1. Data preprocessing and validation
For the _de novo_ metagenome simulation, all sequences shorter than 1000 base pairs are removed from the provided genome assemblies and sequences validated to contain only valid characters.
The input genomes can be draft genomes in fasta format and as such contain beside the bases ACGT also ambiguous DNA characters such as 'RYWSMKHBVDN'.

## Step 2. Community Design

### Genome selection

If the input is a taxonomic profile in [BIOM](http://biom-format.org/) format, the taxa appearing in it, are matched to closely related, complete genomes from the NCBI RefSeq or additional reference files with e.g. new genomes, such that the simulated metagenome (genomes and their abundances) reflects the input profile as closely as possible. This is done via the matching of scientific names in the input profile. Given the taxonomic names of the OTUs in the BIOM profile (this is required), a NCBI taxonomy ID is inferred. For all reference genomes, a taxonomy ID is also required. In the next step, the lineages between the OTUs and the genomes are compared and an OTU is mapped to a genome which has the smallest path between OTU taxonomy ID and genome taxonomy ID. Since the shortest path might be really long, a threshold is set to family level: An OTU will not get mapped to a genome, which is not in the same family, i.e. for every mapped genome it is guaranteed that they are at least of the same family.
Usually this leaves some OTUs unassigned, either because their scientific name did not yield a NCBI taxonomy ID or because for the taxonomy ID no complete genomes are available. It is possible to assign "random" genomes from the reference collection(s) to the unmapped OTUs to increase complexity via an extra option.
 
For the _de novo_ community design, a specified number of sequences (either circular elements or genomes) is selected according to certain criteria, see below, from all input sequences (e.g. circular elements or genomes) to generate a specific data set. 
The random selection of genomes for a particular data set is guided based on their membership in novelty category (Taxonomic annotation step 4.) and OTUs: 
To increase taxonomic spread, for a specified overall number of genomes to be included in a particular data set, the selection algorithm attempts to draw the same number of genomes from each novelty category (new strain, new species, new genus, new family, new order). 
If a category does not have enough genomes to thus enable coming up with the specified genome number, more genomes are drawn from the other categories, in equal amounts, if possible. 
This is repeated until the specified number of genomes has been sampled. 
To model strain level diversity within a species, for every OTU, all available genomes are included in one particular data set up to five genomes at most. 
More than five genomes per OTU are only included, if the specified data set size requires more genomes to be added, that are not available otherwise.

### Creating the community abundances  
For data set simulation, genome and circular element abundance distributions are required for every data set. 
For a single sample data set, the abundance distribution is created by sampling from a log-normal distribution with mu set to 1 and sigma to 2 by default. 
To reflect a differential abundance experiment, two or more abundance samples can be drawn independently from a log-normal distribution with the same parameter settings.
For consecutive samples of a time series, abundances for the next sample are calculated by sampling new values from the log-normal distribution, adding these to the previous value and dividing by two.
The absolute abundance of the circular elements can be set to be like 15 times as high as the abundance of the microbial genomes, to replicate high natural plasmid abundances. 
For each sample, a taxonomic profile for higher ranking taxa is calculated based on the genome abundance values and their taxonomic IDs.

## Step 3. Metagenome sample simulation
Based on the genome abundance distribution, the genome sizes and the specified metagenome sample output size, for every genome, the coverage in the simulated samples is calculated. 
The chosen read simulator then generates reads, using the technology and properties (like read length, insert size or error rates if applicable) as specified by the user, sampling each genome proportionally to the specified abundance in the community. 
A read simulation of a dataset can be run a second time using a different mean insert size using the same genome abundances as before to replicate the use of a mix of technologies on the same samples. 
The size of every sample of each dataset can be set to any size.

Outputs of this step are for every sample of a data set, a FASTQ and a BAM file, which specifies the location of every read in the input sequences.

## Step 4. Creation of the gold standards and anonymization
A gold standard assembly is made with SAMtools for each sample individually and for all samples of each dataset together (pooled gold standard), which included every sequence position with a coverage of at least one, respectively. 
Specifically, the SAM file obtained from the read simulation specifying the location of each read on the reference genomes is used with SAMtools for calculation of the per-site coverage. 
Sequences are then broken into gold standard contig sequences at zero coverage sites. 
The contig sequences are then labeled with a unique dataset and sample identifying prefix and an increasing index as suffix. 
In the last step, sequences of the simulated data sets are shuffled using the unix command 'shuf' and anonymized, thus creating the challenge data sets. 
