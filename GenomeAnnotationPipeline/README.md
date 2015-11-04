GenomeAnnotationPipeline
=================================

Annotation of genomes by clustering of (16S) marker genes

###NAME
	GenomeAnnotationPipeline

###FILE
	./run.py

###FUNCTIONS
	ani_of_genomes_and_novelty_prediction(parser)
		The fourth step is to calculate the average nucleotide identity.
		To lessen the calculation burden, only genomes of sequences within the same cluster as an unpublished genome (marker gene) are compared.
		In case only SILVA sequences are in a cluster, no comparison can be done.
		The tool Mummer is used for the genome comparison, specifically nucmer.
		Novelty predictions are made only for genomes with ani's better than 96%
		ani > 96% -> same species
		ani > 98% -> same strain
		input:
		- meta data table with a list of the genomes and data of the previous step, the output table path will be used as input.
		- working directory which contains the mothur formatted file with the clusters (mothur_otu.txt)
		- file containing a list of fasta file paths, the genomes that someone wants to cluster.
		- file containing a list of reference fasta file paths.
		output:
		- meta data table with a list of the genomes, with columns added that contain ani based tax prediction, rank and novelty prediction

	classification_of_genomes_and_novelty_prediction(parser)
		As the third step, the unpublished genomes are classified based on the clusters they are found in.
		Since clusters were made in 0.01 distance steps, the classification can be done using the smallest clusters first, using bigger ones if a classification can not be made.
		If a marker gene of an unpublished genome is found in a cluster together with references, a common taxon that 90% of sequences agree with will be the predicted taxon.
		The 90% is arbitrary chosen and is required because of taxonomic inconsistencies.
		When a specific rank is checked for agreement, sequences with unknown classification on that rank are ignored.
		TODO: check for taxonomic consitency on higher ranks for those!
		Novelty prediction is based on the predicted taxon's rank. a high rank (phylum, order, class) with low distance can be a strong indicator for taxonomic inconsistencies.
		But it could also be caused by sequences that are not fully classified, yet.
		input:
		- meta data table with a list of the genomes that are to be classified
		- working directory where the results will be saved and which contains the mothur formatted file with the clusters
		output:
		- meta data table with a list of the genomes, with columns added that contain cluster based tax prediction, rank and novelty prediction
		
	gene_alignment_and_clustering(parser)
		The second step is to align 16S sequences and clustering them.
		All is done using "mothur".
		The alignment requires a high quality reference alignment (e.g. from SILVA) and
		is done using the "gotoh" algorithm. (needleman, and blast are possible alternatives, but not implemented here)
		Also using "mothur", empty or uninformative columns are removed.
		When calculating distances (similar to DNADist) multi nucleotide gaps will be counted as one, gaps at the end are ignored.
		To add even more references, the distances of the reference alignment will be merged with those of the working data.
		These were precalculated and only the missing distances between the working data to the reference alignment need to be calculated.
		The clustering will be done based on the final distances using the "Furthest neighbor" algorithm.
		Mothur outputs cluster in steps up to a cutoff. The size of the steps can be chosen, by default 0.01 steps.
		input:
		- fasta formatted reference alignment (e.g. from SILVA)
		- a threshold (between 0 and 1), that will be used for the clustering. 0.03 is default.
		- working directory where the results will be saved and which contains the fasta formatted file with the extracted marker genes of all genomes
		- the number of processors that are available for parallel processing. The program "mothur" can use several cores for the alignments and distance calculations.
		output:
		- a mothur formatted file containing the clusters, from unique up to the given threshold
	
	marker_gene_extraction(parser)
		The first step is to find and extract 16S marker gene sequences. The sequences are found using "hmmsearch" and extracted based on the given positions.
		Two hmmer can currently be used. HMMER3.0 with a profile from 2010 and "rnammer" using HMMER2.0 with a profile from 2006.
		A third one using HMMER3.0 is still to be evaluated.
		So far it seems like rnammer provides better(more) results, but is also very slow.
		input:
		- file containing a list of fasta file paths, for example the genomes that someone wants to cluster.
		- file containing a list of reference fasta file paths or alternatively, a fasta formated file containing the already extracted marker genes of the reference genomes.
		- working directory where the results will be saved (intermediate files will be worked on in the designated /tmp folder)
		- the number of processors that are available for parallel processing. The program "parallel" is used to process several genomes at the same time.
		output:
		- fasta formatted file containing the extracted marker genes of all genomes
	
	my_main()
		16S marker genes extraction, clustering, classification and novelty prediction
		
		NAME run.py
		
		FILE  install_path/run.py
		
		DESCRIPTION
		Script for extracting 16S rRNA sequences from a set of genomes (new ones and reference genomes),
		clustering of 16S sequences by first doing a multiple alignment, and then calculation the distance.
		If the marker genes of the unpublished genomes clustered with one or more of the reference genomes, a taxonomic id is assign to them.
		Finally the ani is calculated for those unpublished genomes, that got clustered with reference genomes.
		
		The three stages in short:
		1. Extraction of 16S sequences
		2. Creation of a multiple alignment of 16S sequences, distance matrix calculation and clustering
		3. Classification of genomes based on clustering and novelty prediction.
		4. ANI calculation and novelty prediction.
		
		input:
		- file containing a list of fasta file paths, for example the genomes that someone wants to cluster.
		- file containing a list of reference fasta file paths. (Step3) alternatively , a fasta formated file containing the already extracted marker genes of the reference genomes.
		- a threshold (between 0 and 1), that will be used for the clustering. 0.03 is default.
		- meta data table with a list of the genomes that are to be classified
		- working directory where results of each stage will be saved
		- the number of processors that are available for parallel processing.
		output:
		- meta data table with a list of the genomes, with columns added that contain tax predictions, average nucleotide identity and novelty predictions
		
		TODO
		- Evaluate third HMMER option
		- confusion matrix of different cutoffs
		- have all logfiles in the working directory
	
###AUTHOR
    hofmann


###USAGE
	run.py [-h] [-irf INPUT_REFERENCE_FNA_FILE] [-ir INPUT_REFERENCE_FILE]
				  [-i INPUT_FILE] [-wd WORKING_DIRECTORY] [-th THRESHOLD]
				  [-p PROCESSORS] [-s STEP] [-im METADATA_TABLE_IN]
				  [-om METADATA_TABLE_OUT]
	
	taxonomic classify genomes by obtaining cluster based on 16S marker genes from genomes
	
	optional arguments:
	
	  -h, --help  
							show this help message and exit  
	  
	  -irf INPUT_REFERENCE_FNA_FILE, --input_reference_fna_file INPUT_REFERENCE_FNA_FILE  
							path to a fasta file containing the 16S marker genes  
							of the reference genomes  
							
	  -ir INPUT_REFERENCE_FILE, --input_reference_file INPUT_REFERENCE_FILE  
							path to a file containing list of reference genomes  
							file paths: <id>\t<path>. NO COLUMN NAMES!  
							
	  -i INPUT_FILE, --input_file INPUT_FILE  
							path to a file containing tab separated list of  
							genomes and their file paths: <id>\t<path>. NO COLUMN  
							NAMES!  
							
	  -wd WORKING_DIRECTORY, --working_directory WORKING_DIRECTORY  
							folder that will contain the log-files, found marker  
							genes and also a file in mothur format containing the  
							clustering  
							
	  -th THRESHOLD, --threshold THRESHOLD  
							threshold for clustering.  Default: 0.03  
							
	  -p PROCESSORS, --processors PROCESSORS  
							number of processors to be used. Default: 2  
							
	  -s STEP, --step STEP  available options: 0-4:  
							0 -> Full run through,  
							1 -> Marker gene extraction,  
							2 -> Gene alignment and clustering,  
							3 -> Classification of Genomes and novelty prediction  
							4 -> Average Nucleotide Identity calculation  
							Default: 0  
							
	  -im METADATA_TABLE_IN, --metadata_table_in METADATA_TABLE_IN  
							path to file containing tab separated list of genomes  
							and their file path  
							
	  -om METADATA_TABLE_OUT, --metadata_table_out METADATA_TABLE_OUT  
							path to file containing tab separated list of genomes  
							and their file path  
