### Output
#### {out}/distributions/distribution_{i}.txt
'{i}' is the index for each sample that is to be generated.
* Column 1: genome_ID
* Column 2: abundance

'genome_ID' is the identifier of the genomes used.  
'abundance' is the relative abundance of a genome to be simulated. 'abundance' does not reflect the amount of genetic data of a genome, but the amount of genomes.  
In a set of two genomes, with both having a abundance of 0.5 but one genome is double the size of the other, the bigger genome will be 66% of the genetic data in the simulated metagenome.

#### {out}/internal/genome_to_id.tsv
If the metagenome is simulated from profile this file is present.
It contains a list of genomes paths to the copies in the output directory in the 'genomes' folder.
* Column 1: genome_ID
* Column 2: file path

#### {out}/internal/meta_data.tsv
If the metagenome is simulated from profile this file is present.
Merged [[meta data|User-manual#meta_datatsv]] of genomes of each community that are actually used for the simulation.

#### {out}/internal/genomes/*.fa
If the metagenome is simulated from profile this files are present.
The used genomes are stored here.

#### {out}/internal/genomes/taxdump.tar.gz
If no Taxonomy dump from the NCBI was specified in the nextflow.config, the downloaded and used one by CAMISIM will be stored here.

#### {out}/sample_{i}/bam/sample{i}_{genome_id}.bam
The bam files generated based on reads generated from the read simulator.

#### {out}/sample_{i}/reads/sample{i}_{genome_id}_*.fq
The simulated reads.

#### {out}/sample_{i}/reads/sample{i}_*.fq
A file containing all reads of this sample.

#### {out}/sample_{i}/reads/anonymous_reads.fq
If anonymization is done, this will be the only fastq file.

#### {out}/sample_{i}/reads/reads_mapping.tsv
Mapping of reads for evaluation

* Column 1: anonymous read id
* Column 2: genome id
* Column 3: taxonomic id
* Column 4: read id

#### {out}/sample_{i}/contigs/anonymous_gsa.fasta
Fasta file with perfect assembly of reads of this sample

#### {out}/sample_{i}/contigs/gsa_mapping.tsv
Mapping of contigs for evaluation

* Column 1: anonymous contig id
* Column 2: genome id
* Column 3: taxonomic id
* Column 4: sequence id of the original genome (in 'source_genomes' folder)
* Column 5: number of reads used in the contig
* Column 6: start position
* Column 7: end position

#### {out}/anonymous_gsa_pooled.fasta
Fasta file with perfect assembly of reads from all samples

#### {out}/gsa_pooled_mapping.tsv
Mapping of contigs from pooled reads fo_validate_raw_genomesr evaluation.

* Column 1: anonymous_contig_id
* Column 2: genome id
* Column 3: taxonomic id
* Column 4: sequence id of the original genome (in 'source_genomes' folder)
* Column 5: number of reads used in the contig
* Column 6: start position
* Column 7: end position

#### {out}/taxonomic_profile_{i}.txt
Taxonomic profile for each sample

#### {out}/seed/seed*.txt
The seeds used for the metagenomesimulation.
If the results of the simulation need to be reproduced, use the used_initial_seed from the {out}/seed/seed.txt in the next run.
The seed can be [[configured|Configuration-File-Options]] in the nextflow.config file.
