#!/usr/bin/env python

__author__ = 'hofmann'

import sys
import os
from scripts.MetaDataTable.metadatatable import MetadataTable
from scripts.argumenthandler import ArgumentHandler
from scripts.NcbiTaxonomy.ncbitaxonomy import NcbiTaxonomy
from scripts.MGAnnotate.mgannotate import MGAnnotate
import ani_prediction_to_meta_table


class MetagenomeSimulationPipeline(ArgumentHandler):
	"""
	Pipeline for the generation of a simulated metagenome
	"""

	def my_main(self):
		"""16S marker genes extraction, clustering, classification and novelty prediction

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
	"""

		# example:
		# usage = '''Example usage:
		# 	# do it all but ani
		# 	{0} -c config.cfg
		#
		# 	# if marker gene extraction already finished
		# 	{0} -c config.cfg -s 2
		# 	'''.format(sys.argv[0])

		options = ArgumentHandler()
		if options._verbose:
			self._logger.info(options.to_string())
		if not options.is_valid():
			self._logger.info("Abort")
			return

		if options.novelty_only:
			reference_map_table = MetadataTable(separator=self._separator, logfile=self._logfile, verbose=self._verbose)
			reference_map_table.read(options.reference_genome_locations_file, column_names=False)
			ref_genome_ids = set(reference_map_table.get_column(0))
			metadata_table = MetadataTable(separator=self._separator, logfile=self._logfile, verbose=self._verbose)
			metadata_table.read(options.metadata_table_in)
			taxonomy = NcbiTaxonomy(options.ncbi_reference_directory, False)
			refernce_ncbi_id_set = set([gid.split('.')[0] for gid in ref_genome_ids])
			mgcluster = MGAnnotate()
			mgcluster.establish_novelty_categorisation(
				taxonomy, refernce_ncbi_id_set, metadata_table, options.column_name_cluster_prediction,
				options.column_name_cluster_novelty)
			metadata_table.write(options.metadata_table_out)
			return

		if options.stage == 0 or options.stage == 1:
			if not self.marker_gene_extraction():
				sys.exit(1)

		if options.stage == 0 or options.stage == 2:
			if not self.gene_alignment_and_clustering():
				sys.exit(1)

		if options.stage == 0 or options.stage == 3:
			if not self.classification_of_genomes_and_novelty_prediction():
				sys.exit(1)

		if options.stage == 4:
			if not self.ani_of_genomes_and_novelty_prediction():
				sys.exit(1)

	def marker_gene_extraction(self):
		"""The first step is to find and extract 16S marker gene sequences. The sequences are found using "hmmsearch" and extracted based on the given positions.
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
	"""
		assert isinstance(self, ArgumentHandler)
		from scripts.MGExtract.mgextract import MGExtract
		mg_extract = MGExtract(
			mg_analyse_executable=os.path.join(self.pipeline_directory, "tools", "16SDetector", "run.py"),
			filename_query_genome_file_paths=self.query_genomes_location_file,
			filename_reference_genome_file_paths=self.reference_genome_locations_file,
			filename_reference_marker_genes=self.reference_markergene_file_path,
			config_path=self.config_file_path,
			filename_iid_map=os.path.join(self.project_directory, self.filename_internal_id_map),
			max_processors=self.processors,
			temp_directory=self.temp_directory,
			debug=self.debug_mode)

		return mg_extract.gather_markergenes(
			hmmer=self.hmmer,
			mg_type="16S",
			output_file=os.path.join(self.project_directory, self.file_mg_16s))

	def gene_alignment_and_clustering(self):
		"""The second step is to align 16S sequences and clustering them.
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
	"""
		assert isinstance(self, ArgumentHandler)
		project_directory = self.project_directory

		input_file = os.path.join(project_directory, self.file_mg_16s)
		if not os.path.isfile(input_file):
			self._logger.error("'-o': File not found: {}".format(input_file))
			return False

		from scripts.MGCluster.mgcluster import MGCluster
		mg_cluster = MGCluster(
			mothur_executable=self.binary_mothur,
			directory_silva_reference=self.silva_reference_directory,
			max_processors=self.processors,
			temp_directory=self.temp_directory,
			debug=self.debug_mode)

		return mg_cluster.cluster(
			marker_gene_fasta=os.path.join(self.project_directory, self.file_mg_16s),
			output_cluster_file=os.path.join(self.project_directory, self.file_cluster_mg_16s),
			distance_cutoff=self.distance_cutoff,
			precision=self.precision,
			method=self.cluster_method)

	def classification_of_genomes_and_novelty_prediction(self):
		"""As the third step, the unpublished genomes are classified based on the clusters they are found in.
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
	"""

		assert isinstance(self, ArgumentHandler)
		from scripts.MGAnnotate.mgannotate import MGAnnotate
		mgcluster = MGAnnotate()
		mgcluster.run(self)

	def ani_of_genomes_and_novelty_prediction(self):
		"""The fourth step is to calculate the average nucleotide identity.
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
	"""
		assert isinstance(self, ArgumentHandler)
		return ani_prediction_to_meta_table.my_main(self)

	def create_meta_table_from_fasta_path_file(self):
		metadata_table = MetadataTable(separator=self._separator, logfile=self._logfile, verbose=self._verbose)
		metadata_table.read(self.reference_genome_locations_file, column_names=False)
		if metadata_table.get_number_of_rows() == 0:
			return False
		id_column = metadata_table.get_column(0)
		metadata_table.clear()
		metadata_table.insert_column(id_column, self.column_name_unpublished_genomes_id)
		metadata_table.write(self.metadata_table_in)
		return True


if __name__ == "__main__":
	pipeline = MetagenomeSimulationPipeline()
	pipeline.my_main()
