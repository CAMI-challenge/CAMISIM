GS(A) File Format
================

##GSA mapping
The header of the gold standard assembly mapping file starts with a '#' to mark it as a line that is to be ignored when read.  

    #anonymous_contig_id	genome_id	tax_id	contig_id	number_reads	start_position	end_position

'anonymous_contig_id' is a column of anonymized contig ids of the original contig ids of the gold standard assembly fasta file.  
'genome_id' is the id of the original genome.  
'taxid' is the taxonomic classification of the genome.  
'contig_id' is the column of the original contig ids of the gold standard assembly fasta file.  
'number_reads' is the column with the number of reads each contig was assembled with.

##GS mapping
The header of the gold standard mapping file starts with a '#' to mark it as a line that is to be ignored when read.  

    #anonymous_read_id	genome_id	tax_id	read_id

'anonymous_read_id' is a column of anonymized read ids of the original read ids of the fastq file with the simulated reads.  
'genome_id' is the id of the original genome.  
'taxid' is the taxonomic classification of the genome.  
'read_id' is the column of the original read ids of the fastq file with the simulated reads.
