__author__ = 'hofmann'

import os
from scripts.MGAnnotate.mothurcluster import MothurCluster
from scripts.MGAnnotate.taxonomiccluster import TaxonomicCluster
from scripts.NcbiTaxonomy.ncbitaxonomy import NcbiTaxonomy
from scripts.MetaDataTable.metadatatable import MetadataTable
from scripts.Validator.validator import Validator


class MGAnnotate(Validator):

	_label = "MGAnnotate"

	def __init__(self, separator="\t", logfile=None, verbose=False, debug=False):
		super(MGAnnotate, self).__init__(logfile=logfile, verbose=verbose, debug=debug)
		self._separator = separator

	def run(self, options):
		metadata_table = MetadataTable()
		metadata_table.read(options.metadata_table_in)
		metadata_table.remove_empty_columns()

		taxonomy = NcbiTaxonomy(options.ncbi_reference_directory, verbose=self._verbose, logfile=self._logfile)

		sequence_mapping = None
		if options.silva_ref_map_file is not None:
			silva_sequence_map = MetadataTable(separator=self._separator, logfile=self._logfile, verbose=self._verbose)
			silva_sequence_map_filename = os.path.join(options.silva_reference_directory, options.silva_ref_map_file)
			silva_sequence_map.read(silva_sequence_map_filename, column_names=False)
			sequence_mapping = silva_sequence_map.get_map(0, 1)

		cluster_file = os.path.join(options.project_directory, options.file_cluster_mg_16s)
		mothur_cluster = MothurCluster(options.precision, sequence_map=sequence_mapping, logfile=self._logfile)
		mothur_cluster.read(cluster_file)

		taxonomy_cluster = TaxonomicCluster(mothur_cluster, taxonomy)

		column_name_unpublished_genomes_id = metadata_table.get_column(options.column_name_unpublished_genomes_id)
		if column_name_unpublished_genomes_id is None:
			msg = "Meta data file does not contain the required header '{}'".format(column_name_unpublished_genomes_id)
			self._logger.error(msg)
			raise IOError(msg)

		self.taxonomic_prediction(options, metadata_table, mothur_cluster, taxonomy_cluster, taxonomy)
		self.set_otu_id(options, metadata_table, mothur_cluster)
		metadata_table.write(options.metadata_table_out)
		return True

	def taxonomic_prediction(self, options, metadata_table, mothur_cluster, taxonomy_cluster, taxonomy):
		reference_map_table = MetadataTable()
		reference_map_table.read(options.input_reference_file, column_names=False)
		ref_genome_ids = set(reference_map_table.get_column(0))

		column_minimum_threshold = metadata_table.get_empty_column()
		column_novely_threshold = metadata_table.get_empty_column()
		column_support = metadata_table.get_empty_column()
		column_cutoff = metadata_table.get_empty_column()
		column_ncbi_prediction = metadata_table.get_empty_column()
		column_science_name = metadata_table.get_empty_column()
		column_novelty = metadata_table.get_empty_column()
		column_name_unpublished_genomes_id = metadata_table.get_column(options.column_name_unpublished_genomes_id)
		if column_name_unpublished_genomes_id is None:
			msg = "Meta data file does not contain the required header '{}'".format(column_name_unpublished_genomes_id)
			self._logger.error(msg)
			raise IOError(msg)
		# _____statistic = {}
		number_of_genomes = len(column_name_unpublished_genomes_id)
		lowest_predicted_novelty = {}

		classification_distance = float(options.classification_distance_minimum)
		max_threshold = mothur_cluster.get_max_threshold()
		if max_threshold == "unique" or float(classification_distance) > float(max_threshold):
			classification_distance = max_threshold
			self._logger.warning("Minimum classification distance unavailable, changed to {}!".format(classification_distance))

		all_done = False
		sorted_lists_of_cutoffs = mothur_cluster.get_sorted_lists_of_cutoffs()
		prediction_thresholds = mothur_cluster.get_prediction_thresholds(minimum=classification_distance)
		# print sorted_lists_of_cutoffs
		# print prediction_thresholds
		# sys.exit()
		for cluster_cutoff in sorted_lists_of_cutoffs:
			if all_done:
				break
			if cluster_cutoff == "unique":
				continue
			self._logger.info("#threshold {}".format(cluster_cutoff))
			cluster_cutoff = float(cluster_cutoff)

			# if cluster_cutoff not in _____statistic:
			# 	_____statistic[cluster_cutoff] = {"sname": metadata_table.get_empty_column(), "novelty": metadata_table.get_empty_column(), "support": metadata_table.get_empty_column() }

			all_done = True
			for row_index in range(0, number_of_genomes):
				if row_index not in lowest_predicted_novelty:
					lowest_predicted_novelty[row_index] = {"novelty": '', "support": '', "threshold": '', "minimum": ''}
				unpublished_genome_id = column_name_unpublished_genomes_id[row_index]
				if not mothur_cluster.element_exists(cluster_cutoff, unpublished_genome_id):
					# if no marker gene was found it will not be in the clustering
					continue
				if unpublished_genome_id == "":
					continue
				all_done = False
				separator = ""
				list_of_cluster_id, list_of_cluster = mothur_cluster.get_cluster_of_cutoff_of_element(cluster_cutoff, unpublished_genome_id)
				if len(list_of_cluster) > 1:
					separator = ";"
				predicted__ncbi = []
				predicted_science_name = []
				# predicted_novelty = []
				# list_support = []
				for cluster in list_of_cluster:
					# ncbi_prediction, novelty = taxonomy_cluster.get_cluster_ncbi_tax_prediction(cluster, column_name_unpublished_genomes_id, unpublished_genome_ids)
					ncbi_prediction, novelty, support = taxonomy_cluster.predict_tax_id_of(cluster, column_name_unpublished_genomes_id, unpublished_genome_id, lowest_predicted_novelty[row_index], cluster_cutoff, ref_genome_ids)
					if ncbi_prediction is None:
						continue
					# print unpublished_genome_id, novelty
					if cluster_cutoff not in prediction_thresholds or ncbi_prediction in predicted__ncbi or column_ncbi_prediction[row_index] != "":
						continue
					column_cutoff[row_index] = str(cluster_cutoff)
					predicted__ncbi.append(ncbi_prediction)
					predicted_science_name.append(taxonomy.get_scientific_name(ncbi_prediction))
					# predicted_novelty.append("new_" + novelty)

				# 		list_support.append(support)
				# _____statistic[cluster_cutoff]["sname"][row_index] = separator.join(predicted_science_name)
				# _____statistic[cluster_cutoff]["novelty"][row_index] = separator.join(predicted_novelty)
				# _____statistic[cluster_cutoff]["support"][row_index] = separator.join(list_support)
				if column_ncbi_prediction[row_index] != "":
					continue
				column_ncbi_prediction[row_index] = separator.join(predicted__ncbi)
				column_science_name[row_index] = separator.join(predicted_science_name)

		# unknown_novelty = False
		for row_index in range(0, number_of_genomes):
			if lowest_predicted_novelty[row_index]["novelty"] is not '':
				column_novelty[row_index] = "new_" + lowest_predicted_novelty[row_index]["novelty"]
				column_novely_threshold[row_index] = lowest_predicted_novelty[row_index]["threshold"]
				column_support[row_index] = lowest_predicted_novelty[row_index]["support"]
				column_minimum_threshold[row_index] = lowest_predicted_novelty[row_index]["minimum"]
			# else:
			# 	unknown_novelty = True

		metadata_table.insert_column(column_cutoff, options.column_name_cutoff)
		metadata_table.insert_column(column_ncbi_prediction, options.column_name_cluster_prediction)
		metadata_table.insert_column(column_science_name, options.column_name_cluster_scientific_name)
		metadata_table.insert_column(column_novelty, options.column_name_cluster_novelty)
		# metadata_table.insert_column(column_novely_threshold, "novelty threshold")
		# metadata_table.insert_column(column_support, "novelty support")
		# metadata_table.insert_column(column_minimum_threshold, "minimum threshold")

		# if unknown_novelty:
		refernce_ncbi_id_set = set([gid.split('.')[0] for gid in ref_genome_ids])
		self.establish_novelty_categorisation(
			taxonomy, refernce_ncbi_id_set, metadata_table, options.column_name_cluster_prediction, options.column_name_cluster_novelty)

		# for cluster_cutoff in sorted_lists_of_cutoffs:
		# 	if cluster_cutoff == "unique":
		# 		continue
		# 	metadata_table.insert_column(_____statistic[cluster_cutoff]["support"], "{}_support".format(cluster_cutoff))
		# for cluster_cutoff in sorted_lists_of_cutoffs:
		# 	if cluster_cutoff == "unique":
		# 		continue
		# 	metadata_table.insert_column(_____statistic[cluster_cutoff]["sname"], "{}_name".format(cluster_cutoff))
		# for cluster_cutoff in sorted_lists_of_cutoffs:
		# 	if cluster_cutoff == "unique":
		# 		continue
		# 	metadata_table.insert_column(_____statistic[cluster_cutoff]["novelty"], "{}_novelty".format(cluster_cutoff))

		self._logger.info("Taxonomic prediction finished")

	def establish_novelty_categorisation(
		self, taxonomy, refernce_ncbi_id_set, metadata_table, column_name_cluster_prediction, column_name_cluster_novelty):
		from scripts.MGAnnotate.novelty import Novelty
		self._logger.info("Establish novelty categorisation")
		# novelty = Novelty(taxonomy,
		# 				  logger=logger,
		# 				  column_name_ncbi_id=options.column_name_cluster_prediction,
		# 				  column_name_novelty="NOVELTY_CATEGORY_Jessika")
		novelty = Novelty(
			taxonomy,
			column_name_ncbi_id=column_name_cluster_prediction,
			column_name_novelty=column_name_cluster_novelty,
			separator="\t", logfile=None, verbose=True, debug=False)
		novelty.read_reference(refernce_ncbi_id_set)
		novelty.compute_novelty(metadata_table)
		self._logger.info("Done")

	def set_otu_id(self, options, metadata_table, mothur_cluster):
		list_of_unclustered_elements = set()
		column_name_unpublished_genomes_id = metadata_table.get_column(options.column_name_unpublished_genomes_id)
		number_of_genomes = len(column_name_unpublished_genomes_id)
		otu_distance = options.otu_distance
		max_threshold = mothur_cluster.get_max_threshold()
		if max_threshold == "unique" or float(otu_distance) > float(max_threshold):
			otu_distance = max_threshold
			self._logger.warning("OTU distance unavailable, changed to {}!".format(otu_distance))

		column_otu_id = metadata_table.get_empty_column()
		sorted_lists_of_cutoffs = mothur_cluster.get_sorted_lists_of_cutoffs()
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
				list_of_cluster_id, list_of_cluster = mothur_cluster.get_cluster_of_cutoff_of_element(cluster_cutoff, unpublished_genome_id)
				if len(list_of_cluster_id) > 1:
					separator = ";"
				column_otu_id[row_index] = separator.join([str(otu_id) for otu_id in sorted(set(list_of_cluster_id))])
		metadata_table.insert_column(column_otu_id, options.column_name_otu_id)
		if len(list_of_unclustered_elements) > 0:
			self._logger.warning("No cluster found for {} ids!".format(len(list_of_unclustered_elements)))
		self._logger.info("OTU finished")
