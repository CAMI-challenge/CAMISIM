#!/usr/bin/env python

__author__ = 'Peter Hofmann'
__version__ = "0.0.6"

import os
import traceback
from scripts.MetaDataTable.metadatatable import MetadataTable
from scripts.argumenthandler_ga import ArgumentHandler
from scripts.NcbiTaxonomy.ncbitaxonomy import NcbiTaxonomy
from scripts.MGCluster.mgcluster import MGCluster
from scripts.MGAnnotate.mgannotate import MGAnnotate
from scripts.MGExtract.mgextract import MGExtract
from scripts.MGAnnotate.mothurcluster import MothurCluster
from scripts.MGAnnotate.taxonomiccluster import TaxonomicCluster


class GenomeAnnotation(ArgumentHandler):
	"""
	Pipeline for the annotation of genomes
	"""

	_label = "GenomeAnnotationPipeline"

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

		@rtype: None
		"""

		# example:
		# usage = '''Example usage:
		# 	# do it all but ani
		# 	{0} -c config.cfg
		#
		# 	# if marker gene extraction already finished
		# 	{0} -c config.cfg -s 2
		# 	'''.format(sys.argv[0])

		# self._logger.info(self.to_string())
		if not self._input_valid():
			self._logger.info("Aborted")
			return
		self._logger.info("Starting")
		try:
			if self._validate_genomes:
				self._validate_raw_genomes()

			if self._phase == 0 or self._phase == 1:
				self.marker_gene_extraction()

			if self._phase == 0 or self._phase == 2:
				self.gene_alignment_and_clustering()

			if self._phase == 0 or self._phase == 3:
				self.marker_gene_annotation()

		except (KeyboardInterrupt, SystemExit, Exception, ValueError, OSError) as e:
			self._logger.debug("\n{}\n".format(traceback.format_exc()))
			if len(e.args) > 0:
				self._logger.error(e.args[0])
			self._logger.error("Aborted")
		except AssertionError:
			self._logger.error("Aborted")
		else:
			self._logger.info("Finished")

		if not self._debug:
			self._project_file_folder_handler.remove_directory_temp()
		else:
			self._logger.info("Temporary data stored at:\n{}".format(self._project_file_folder_handler.get_tmp_wd()))

	def _validate_raw_genomes(self):
		"""
		Validate format raw and reference genomes

		@rtype: None
		"""
		self._logger.info("Validating Genomes")
		meta_data_table = MetadataTable(
			separator=self._separator,
			logfile=self._logfile,
			verbose=self._verbose)

		are_valid = True
		meta_data_table.read(self._file_path_query_genomes_location_file, column_names=False)
		list_of_file_paths = meta_data_table.get_column(1)
		if not self._validate_format(list_of_file_paths, file_format="fasta", sequence_type="dna", ambiguous=True):
			are_valid = False

		meta_data_table.read(self._file_path_reference_genome_locations, column_names=False)
		list_of_file_paths = meta_data_table.get_column(1)
		if not self._validate_format(list_of_file_paths, file_format="fasta", sequence_type="dna", ambiguous=True):
			are_valid = False

		if not are_valid:
			msg = "Invalid genomes found!"
			self._logger.error(msg)
			raise RuntimeError(msg)
		self._logger.info("Validating Genomes Done")

	def _validate_format(self, list_of_file_paths, file_format="fasta", sequence_type="dna", ambiguous=True):
		"""
		Validate file format of a list of fasta files

		@param list_of_file_paths: List of fasta file paths
		@type list_of_file_paths: list[str|unicode]
		@param file_format: 'fasta' or 'fastq'
		@type file_format: str | unicode
		@param sequence_type: 'dna' or 'rna' or 'protein'
		@type sequence_type: str | unicode
		@param ambiguous: If true ambiguous characters are valid
		@type ambiguous: bool

		@return: True if all valid
		@rtype: bool
		"""
		result = True
		for file_path in list_of_file_paths:
			if not self.validate_sequence_file(file_path, file_format, sequence_type, ambiguous):
				result = False
		return result

	def marker_gene_extraction(self):
		"""
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

		@rtype: None
		"""
		assert isinstance(self, ArgumentHandler)
		mg_extract = MGExtract(
			mg_analyse_executable=self._get_mg_analyse_executable(),
			file_path_query_genome_file_paths=self._file_path_query_genomes_location_file,
			file_path_reference_genome_file_paths=self._file_path_reference_genome_locations,
			file_path_name_reference_marker_genes=self._file_path_reference_markergene,
			config_path=self._file_path_config,
			file_path_map_reference_genome_id_to_tax_id=self._file_path_map_reference_genome_id_to_tax_id,
			max_processors=self._max_processors,
			temp_directory=self._project_file_folder_handler.get_tmp_wd(),
			separator=self._separator, logfile=self._logfile, verbose=self._verbose, debug=self._debug)

		mg_extract.gather_markergenes(
			hmmer=self._hmmer,
			mg_type="16S",
			file_path_output=self._project_file_folder_handler.get_file_path_mg_16s(),
			file_path_map_uid_sid=self._project_file_folder_handler.get_file_path_internal_id_map())

		# merge silva iid with genome iid
		data_table_iid_mapping_silva = MetadataTable(separator=self._separator, logfile=self._logfile, verbose=self._verbose)
		file_path_silva_map = os.path.join(self._silva_reference_directory, MGCluster.get_file_name_of_map())
		data_table_iid_mapping_silva.read(file_path_silva_map)
		data_table_iid_mapping = MetadataTable(separator=self._separator, logfile=self._logfile, verbose=self._verbose)
		data_table_iid_mapping.read(self._project_file_folder_handler.get_file_path_internal_id_map())
		data_table_iid_mapping.concatenate(data_table_iid_mapping_silva, strict=False)
		data_table_iid_mapping.write(self._project_file_folder_handler.get_file_path_internal_id_map())

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

		@rtype: None
	"""
		assert self.validate_file(self._project_file_folder_handler.get_file_path_mg_16s())

		mg_cluster = MGCluster(
			mothur_executable=self._binary_mothur,
			directory_silva_reference=self._silva_reference_directory,
			max_processors=self._max_processors,
			temp_directory=self._directory_temp,
			logfile=self._logfile, verbose=self._verbose, debug=self._debug)

		mg_cluster.cluster(
			marker_gene_fasta=self._project_file_folder_handler.get_file_path_mg_16s(),
			output_cluster_file=self._project_file_folder_handler.get_file_path_cluster_mg_16s(),
			distance_cutoff=self._distance_cutoff,
			precision=self._precision,
			method=self._cluster_method)

	def marker_gene_annotation(self):
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

		@rtype: None
		"""
		# set of taxonomic ids of well known genomes
		data_table = MetadataTable(separator=self._separator, logfile=self._logfile, verbose=self._verbose)
		data_table.read(self._file_path_map_reference_genome_id_to_tax_id)
		list_of_refernce_ncbi_id = data_table.get_column(1)

		# mapping of all internal ids
		# data_table_iid_mapping_silva = MetadataTable(
		# 	separator=self._separator, logfile=self._logfile, verbose=self._verbose)
		# file_path_silva_map = os.path.join(self._silva_reference_directory, MGCluster.get_file_name_of_map())
		# data_table_iid_mapping_silva.read(file_path_silva_map)
		data_table_iid_mapping = MetadataTable(separator=self._separator, logfile=self._logfile, verbose=self._verbose)
		data_table_iid_mapping.read(self._project_file_folder_handler.get_file_path_internal_id_map())
		# data_table_iid_mapping.concatenate(data_table_iid_mapping_silva, strict=False)

		mg_annotate = MGAnnotate(
			# ncbi_reference_directory=self._ncbi_reference_directory,
			file_path_query_genomes_location=self._file_path_query_genomes_location_file,
			file_path_reference_genomes_location=self._file_path_reference_genome_locations,
			file_path_reference_taxid_map=self._file_path_map_reference_genome_id_to_tax_id,
			file_path_nucmer=self._file_path_nucmer,
			column_name_genome_id=self._column_name_genome_id,
			column_name_otu=self._column_name_otu_id,
			column_name_novelty_category=self._column_name_cluster_novelty,
			column_name_ncbi=self._column_name_ncbi,
			column_name_scientific_name=self._column_name_cluster_scientific_name,
			column_name_ani=self._column_name_ani,
			column_name_ani_novelty=self._column_name_ani_novelty,
			column_name_ani_ncbi=self._column_name_ani_compare,
			column_name_ani_scientific_name=self._column_name_ani_scientific_name,
			temp_directory=self._directory_temp, max_processors=self._max_processors,
			separator=self._separator, logfile=self._logfile, verbose=self._verbose, debug=self._debug
		)

		metadata_table = MetadataTable(separator=self._separator, logfile=self._logfile, verbose=self._verbose)
		metadata_table.read(self._metadata_table_in, column_names=True)
		metadata_table.remove_empty_columns()

		list_query_gid = metadata_table.get_column(self._column_name_genome_id)
		if list_query_gid is None:
			msg = "Meta data file does not contain the required header '{}'".format(self._column_name_genome_id)
			self._logger.error(msg)
			raise IOError(msg)

		taxonomy = NcbiTaxonomy(self._ncbi_reference_directory, verbose=self._verbose, logfile=self._logfile)

		mothur_cluster = MothurCluster(
			self._precision, iid_gid_mapping=data_table_iid_mapping.get_map(0, 1),
			logfile=self._logfile, verbose=self._verbose, debug=self._debug)
		mothur_cluster.read(self._project_file_folder_handler.get_file_path_cluster_mg_16s(), list_query_gid)

		taxonomy_cluster = TaxonomicCluster(
			mothur_cluster, taxonomy, iid_tid_map=data_table_iid_mapping.get_map(0, 2),
			set_reference_genome_ncbi=set(list_of_refernce_ncbi_id),
			logfile=self._logfile, verbose=self._verbose, debug=self._debug)

		if self._annotate_classify:
			self._logger.info("Taxonomic classification")
			# also, novelty based clustering
			mg_annotate.taxonomic_classification(
				metadata_table, mothur_cluster, taxonomy_cluster, taxonomy, self._classification_distance_minimum)
			self._logger.info("Taxonomic classification Done")

		if self._annotate_novelty:
			self._logger.info("Novelty categorisation")
			# novelty by comparing with reference taxonomic ids
			mg_annotate.novelty_categorisation(taxonomy, set(list_of_refernce_ncbi_id), metadata_table)
			self._logger.info("Novelty categorisation Done")

		if self._annotate_otu:
			self._logger.info("OTU")
			mg_annotate.set_otu_id(metadata_table, mothur_cluster, self._otu_distance)
			self._logger.info("OTU Done")

		if self._annotate_ani:
			self._logger.info("Calculating ANI")
			mg_annotate.calculate_ani(
				mothur_cluster, taxonomy, metadata_table, self._distance_cutoff, self._ani_minimum_alignment)
			self._logger.info("Calculating ANI Done")
		metadata_table.write(self._project_file_folder_handler.get_file_path_meta_data_table(), column_names=True)

	def create_meta_table(self, file_path_metadata_table):
		"""
		Generate a input metadata file with genome ids only

		@param file_path_metadata_table:
		@type file_path_metadata_table: str|unicode

		@rtype: None
		"""
		metadata_table = MetadataTable(separator=self._separator, logfile=self._logfile, verbose=self._verbose)
		metadata_table.read(self._file_path_reference_genome_locations, column_names=False)
		if metadata_table.get_number_of_rows() == 0:
			raise ValueError("Invalid file content")
		id_column = metadata_table.get_column(0)
		metadata_table.clear()
		metadata_table.insert_column(id_column, self._column_name_genome_id)
		metadata_table.write(file_path_metadata_table, column_names=True)


if __name__ == "__main__":
	pipeline = GenomeAnnotation(
		args=None, separator="\t",
		column_name_genome_id="genome_ID", column_name_otu="OTU", column_name_novelty_category="novelty_category",
		column_name_ncbi="NCBI_ID")
	pipeline.my_main()
