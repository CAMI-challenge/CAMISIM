__author__ = 'hofmann'

from scripts.MGAnnotate.mothurcluster import MothurCluster
from scripts.MGAnnotate.taxonomiccluster import TaxonomicCluster
from scripts.NcbiTaxonomy.ncbitaxonomy import NcbiTaxonomy
from scripts.MetaDataTable.metadatatable import MetadataTable
from scripts.Validator.validator import Validator
from scripts.MGAnnotate.anim import ANIm
from scripts.MGAnnotate.novelty import Novelty


class MGAnnotate(Validator):

	_label = "MGAnnotate"

	def __init__(
		self,
		file_path_query_genomes_location, file_path_reference_genomes_location, file_path_reference_taxid_map,
		file_path_nucmer=None,
		column_name_genome_id="genome_ID", column_name_otu="OTU", column_name_novelty_category="novelty_category",
		column_name_ncbi="NCBI_ID", column_name_scientific_name="SCIENTIFIC_NAME",
		column_name_threshold="prediction_threshold",
		column_name_ani="ANI", column_name_ani_novelty="ANI_NOVELTY_CATEGORY",
		column_name_ani_ncbi="ANI_NCBI_ID", column_name_ani_scientific_name="ANI_SCIENTIFIC_NAME",
		temp_directory=None, max_processors=1, separator="\t", logfile=None, verbose=False, debug=False):
		"""
		Constructor

		@param separator: Expected column separator in metadata files
		@type separator: str|unicode
		@param logfile: file handler or file path to a log file
		@type logfile: file | FileIO | StringIO | basestring
		@param verbose: Not verbose means that only warnings and errors will be past to stream
		@type verbose: bool
		@param debug: Display debug messages
		@type debug: bool
		"""
		assert file_path_nucmer is None or self.validate_file(file_path_nucmer, executable=True)
		# assert self.validate_dir(ncbi_reference_directory)
		assert isinstance(max_processors, (int, long))
		assert self.validate_number(max_processors, minimum=1)
		assert isinstance(separator, basestring)
		assert isinstance(column_name_genome_id, basestring)
		assert isinstance(column_name_otu, basestring)
		assert isinstance(column_name_novelty_category, basestring)
		assert isinstance(column_name_ncbi, basestring)
		assert isinstance(column_name_scientific_name, basestring)
		assert isinstance(column_name_threshold, basestring)
		assert isinstance(column_name_ani, basestring)
		assert isinstance(column_name_ani_novelty, basestring)
		assert isinstance(column_name_ani_ncbi, basestring)
		assert isinstance(column_name_ani_scientific_name, basestring)
		super(MGAnnotate, self).__init__(logfile=logfile, verbose=verbose, debug=debug)
		self._column_name_genome_id = column_name_genome_id
		self._column_name_otu_id = column_name_otu
		self._column_name_ncbi_prediction = column_name_ncbi
		self._column_name_ncbi_scientific_name = column_name_scientific_name
		self._column_name_novelty = column_name_novelty_category
		self._column_name_threshold = column_name_threshold
		self._column_name_ani = column_name_ani
		self._column_name_ani_novelty = column_name_ani_novelty
		self._column_name_ani_ncbi = column_name_ani_ncbi
		self._column_name_ani_scientific_name = column_name_ani_scientific_name
		self._separator = separator
		# self._ncbi_reference_directory = ncbi_reference_directory
		self._file_path_nucmer = file_path_nucmer
		self._max_processors = max_processors
		self._file_path_query_genomes_location = file_path_query_genomes_location
		self._file_path_reference_genomes_location = file_path_reference_genomes_location
		self._file_path_reference_taxid_map = file_path_reference_taxid_map
		self._tmp_dir = temp_directory

	def taxonomic_classification(self, metadata_table, mothur_cluster, taxonomy_cluster, taxonomy, classification_distance_minimum):
		"""
		Taxonomic classification of genomes

		@param metadata_table: Handler of MetadataTable
		@type metadata_table: MetadataTable
		@param mothur_cluster: Handler of MothurCluster
		@type mothur_cluster: MothurCluster
		@param taxonomy_cluster: Handler of TaxonomicCluster
		@type taxonomy_cluster: TaxonomicCluster
		@param taxonomy: Handler of NcbiTaxonomy
		@type taxonomy: NcbiTaxonomy
		@param classification_distance_minimum: Minimum genetic distance in percent
		@type classification_distance_minimum: int | float

		@rtype: None
		"""
		assert isinstance(metadata_table, MetadataTable)
		assert isinstance(mothur_cluster, MothurCluster)
		assert isinstance(taxonomy_cluster, TaxonomicCluster)
		assert isinstance(taxonomy, NcbiTaxonomy)
		assert isinstance(classification_distance_minimum, (int, float))
		column_minimum_threshold = metadata_table.get_empty_column()
		column_novely_threshold = metadata_table.get_empty_column()
		column_support = metadata_table.get_empty_column()
		column_cutoff = metadata_table.get_empty_column()
		column_ncbi_prediction = metadata_table.get_empty_column()
		column_science_name = metadata_table.get_empty_column()
		column_novelty = metadata_table.get_empty_column()
		list_query_gid = metadata_table.get_column(self._column_name_genome_id)
		if list_query_gid is None:
			msg = "Meta data file does not contain the required header '{}'".format(list_query_gid)
			self._logger.error(msg)
			raise IOError(msg)
		# _____statistic = {}
		lowest_predicted_novelty = {}

		classification_distance = float(classification_distance_minimum)
		max_threshold = mothur_cluster.get_max_threshold()
		if max_threshold == "unique" or float(classification_distance) > float(max_threshold):
			classification_distance = max_threshold
			self._logger.warning("Minimum classification distance unavailable, changed to {}!".format(classification_distance))

		sorted_lists_of_cutoffs = mothur_cluster.get_sorted_lists_of_thresholds()
		prediction_thresholds = mothur_cluster.get_prediction_thresholds(minimum=classification_distance)
		for cluster_cutoff in sorted_lists_of_cutoffs:
			if cluster_cutoff == "unique":
				continue
			self._logger.debug("Threshold {}".format(cluster_cutoff))
			cluster_cutoff = float(cluster_cutoff)

			for row_index, query_gid in enumerate(list_query_gid):
				if row_index not in lowest_predicted_novelty:
					lowest_predicted_novelty[row_index] = {"novelty": '', "support": '', "threshold": '', "minimum": ''}
				if not mothur_cluster.element_exists(cluster_cutoff, query_gid):
					# self._logger.debug("'{}' not found!".format(query_gid))
					# if no marker gene was found it will not be in the clustering
					continue
				if query_gid == "":
					continue
				separator = ""
				list_of_cluster_id, list_of_cluster = mothur_cluster.get_cluster_of_threshold_of_gid(cluster_cutoff, query_gid)
				if len(list_of_cluster) > 1:
					separator = ";"
				predicted__ncbi = []
				predicted_science_name = []
				for cluster in list_of_cluster:
					ncbi_prediction, novelty, support = taxonomy_cluster.predict_tax_id_of(cluster, lowest_predicted_novelty[row_index])
					if ncbi_prediction is None:
						continue
					if cluster_cutoff not in prediction_thresholds or ncbi_prediction in predicted__ncbi or column_ncbi_prediction[row_index] != "":
						continue
					column_cutoff[row_index] = str(cluster_cutoff)
					predicted__ncbi.append(ncbi_prediction)
					predicted_science_name.append(taxonomy.get_scientific_name(ncbi_prediction))

				if column_ncbi_prediction[row_index] != "":
					continue
				column_ncbi_prediction[row_index] = separator.join(predicted__ncbi)
				column_science_name[row_index] = separator.join(predicted_science_name)

		# unknown_novelty = False
		for row_index, query_gid in enumerate(list_query_gid):
			if lowest_predicted_novelty[row_index]["novelty"] is not '':
				column_novelty[row_index] = "new_" + lowest_predicted_novelty[row_index]["novelty"]
				column_novely_threshold[row_index] = lowest_predicted_novelty[row_index]["threshold"]
				column_support[row_index] = lowest_predicted_novelty[row_index]["support"]
				column_minimum_threshold[row_index] = lowest_predicted_novelty[row_index]["minimum"]
			# else:
			# 	unknown_novelty = True

		metadata_table.insert_column(column_cutoff, self._column_name_threshold)
		metadata_table.insert_column(column_ncbi_prediction, self._column_name_ncbi_prediction)
		metadata_table.insert_column(column_science_name, self._column_name_ncbi_scientific_name)
		metadata_table.insert_column(column_novelty, self._column_name_novelty)
		# metadata_table.insert_column(column_novely_threshold, "novelty threshold")
		# metadata_table.insert_column(column_support, "novelty support")
		# metadata_table.insert_column(column_minimum_threshold, "minimum threshold")

	def novelty_categorisation(self, taxonomy, reference_ncbi_id_set, metadata_table):
		"""
		Predict novelty of a genome

		@param taxonomy: Handler of NcbiTaxonomy
		@type taxonomy: NcbiTaxonomy
		@param reference_ncbi_id_set: Reference taxonomic ids of known genomes
		@type reference_ncbi_id_set: set[str|unicode]
		@param metadata_table: Handler of MetadataTable
		@type metadata_table: MetadataTable

		@rtype: None
		"""
		assert isinstance(taxonomy, NcbiTaxonomy)
		assert isinstance(reference_ncbi_id_set, set)
		assert isinstance(metadata_table, MetadataTable)
		novelty = Novelty(
			taxonomy,
			column_name_ncbi_id=self._column_name_ncbi_prediction,
			column_name_novelty=self._column_name_novelty,
			separator="\t", logfile=None, verbose=True, debug=False)
		novelty.read_reference(set(reference_ncbi_id_set))
		novelty.compute_novelty(metadata_table)

	def set_otu_id(self, metadata_table, mothur_cluster, otu_distance):
		"""
		Set OTU id based on clusters at a threshold

		@param metadata_table: Handle of MetadataTable
		@type metadata_table: MetadataTable
		@param mothur_cluster: Handle of MothurCluster
		@type mothur_cluster: MothurCluster
		@param otu_distance: Genetic distance in percent
		@type otu_distance: int|float

		@rtype: None
		"""
		assert isinstance(metadata_table, MetadataTable)
		assert isinstance(mothur_cluster, MothurCluster)
		assert isinstance(otu_distance, (int, float))
		list_of_unclustered_elements = set()
		column_name_unpublished_genomes_id = metadata_table.get_column(self._column_name_genome_id)
		number_of_genomes = len(column_name_unpublished_genomes_id)
		max_threshold = mothur_cluster.get_max_threshold()
		if max_threshold == "unique" or float(otu_distance) > float(max_threshold):
			otu_distance = max_threshold
			self._logger.warning("OTU distance unavailable, changed to {}!".format(otu_distance))

		column_otu_id = metadata_table.get_empty_column()
		sorted_lists_of_cutoffs = mothur_cluster.get_sorted_lists_of_thresholds()
		for cluster_cutoff in sorted_lists_of_cutoffs:
			if cluster_cutoff == "unique" or otu_distance == "unique" or float(cluster_cutoff) != float(otu_distance):
				continue
			cluster_cutoff = float(cluster_cutoff)
			for row_index in range(0, number_of_genomes):
				unpublished_genome_id = column_name_unpublished_genomes_id[row_index]
				if not mothur_cluster.element_exists(cluster_cutoff, unpublished_genome_id):
					# if no marker gene was found it will not be in the clustering
					list_of_unclustered_elements.add(unpublished_genome_id)
					continue
				separator = ""
				list_of_cluster_id, list_of_cluster = mothur_cluster.get_cluster_of_threshold_of_gid(cluster_cutoff, unpublished_genome_id)
				if len(list_of_cluster_id) > 1:
					separator = ";"
				column_otu_id[row_index] = separator.join([str(otu_id) for otu_id in sorted(set(list_of_cluster_id))])
		metadata_table.insert_column(column_otu_id, self._column_name_otu_id)
		if len(list_of_unclustered_elements) > 0:
			self._logger.warning("No cluster found for {} ids!".format(len(list_of_unclustered_elements)))

	def calculate_ani(self, mothur_cluster, taxonomy, metadata_table, ani_distance, ani_minimum_alignment=0.9):
		"""
		Calculate the average nucleotide identity.
		To lessen the calculation burden, only genomes of sequences within the same cluster as an unpublished genome (marker gene) are compared.
		In case only SILVA sequences are in a cluster, no comparison can be done.
		The tool Mummer is used for the genome comparison, specifically nucmer.
		Novelty predictions are made only for genomes with ani's better than 96%
		ani > 96% -> same species
		ani > 98% -> same strain

		@param mothur_cluster: Handler of MothurCluster
		@type mothur_cluster: MothurCluster
		@param taxonomy: Handler of NcbiTaxonomy
		@type taxonomy: NcbiTaxonomy
		@param metadata_table: Handler of MetadataTable
		@type metadata_table: MetadataTable
		@param ani_distance: Looks up if one of the taxids within a cluster at that distance exists
		@type ani_distance: int|long

		@rtype: None
		"""
		assert isinstance(mothur_cluster, MothurCluster)
		assert isinstance(taxonomy, NcbiTaxonomy)
		assert isinstance(metadata_table, MetadataTable)
		assert isinstance(ani_distance, (float, long, int)), "Got '{}'".format(ani_distance)
		assert self.validate_number(ani_distance, minimum=0, maximum=1)
		ani_scientific_name_column = metadata_table.get_empty_column()
		ani_prediction_novelty_column = metadata_table.get_empty_column()
		ani_prediction_column = metadata_table.get_empty_column()
		ani_column = metadata_table.get_empty_column()
		query_genome_ids_column = metadata_table.get_column(self._column_name_genome_id)

		if not mothur_cluster.has_threshold(ani_distance):
			max_threshold = mothur_cluster.get_max_threshold()
			self._logger.warning("Distance '{}' not available using '{}' instead.".format(
				ani_distance, max_threshold))
			if not max_threshold == "unique":
				max_threshold = float(max_threshold)
			ani_distance = max_threshold

		ani_calculator = ANIm(
			file_path_nucmer=self._file_path_nucmer,
			minimum_alignment=ani_minimum_alignment,
			file_path_query_genomes_location=self._file_path_query_genomes_location,
			file_path_reference_genomes_location=self._file_path_reference_genomes_location,
			file_path_reference_taxid_map=self._file_path_reference_taxid_map,
			separator=self._separator, temp_directory=self._tmp_dir, max_processors=self._max_processors,
			logfile=self._logfile, verbose=self._verbose, debug=self._debug)
		list_of_clusters_by_gid = {}
		for gid in query_genome_ids_column:
			list_of_index, clusters = mothur_cluster.get_cluster_of_threshold_of_gid(ani_distance, gid)
			list_of_clusters_by_gid[gid] = clusters
		for query_genomes_id in query_genome_ids_column:
			for clusters in list_of_clusters_by_gid[query_genomes_id]:
				if clusters is None:
					continue
				ani_calculator.add_nucmer_cmd_lines(mothur_cluster.iid_to_gid_list(clusters), [query_genomes_id])

		total_lengths, sim_errors, percent_identity, percent_alignment, ncbi = ani_calculator.calculate_best_anim()

		for row_index, query_genomes_id in enumerate(query_genome_ids_column):
			if query_genomes_id in percent_identity and percent_identity[query_genomes_id] > 0:
				if float(percent_identity[query_genomes_id]) > 0.98:
					ani_prediction_novelty_column[row_index] = "same_strain"
				elif float(percent_identity[query_genomes_id]) > 0.96:
					ani_prediction_novelty_column[row_index] = "same_species"
				ani_prediction_column[row_index] = ncbi[query_genomes_id]
				ani_column[row_index] = str(percent_identity[query_genomes_id])
				science_name = taxonomy.get_scientific_name(ncbi[query_genomes_id])
				if science_name is not None:
					ani_scientific_name_column[row_index] = science_name

		metadata_table.insert_column(ani_column, self._column_name_ani)
		metadata_table.insert_column(ani_prediction_novelty_column, self._column_name_ani_novelty)
		metadata_table.insert_column(ani_prediction_column, self._column_name_ani_ncbi)
		metadata_table.insert_column(ani_scientific_name_column, self._column_name_ani_scientific_name)
