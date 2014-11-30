__author__ = 'hofmann'

#from source.ArgumentHandler import ArgumentHandler
from source.MothurCluster import MothurCluster
from source.TaxonomicCluster import TaxonomicCluster
from source.NcbiTaxonomy.NcbiTaxonomy import NcbiTaxonomy
from source.MetaTable import MetaTable
from source.Logger import Logger
import sys
import os


def my_main(options):
	logger = Logger()
	metadata_table = MetaTable(logger=logger)
	metadata_table.read(options.metadata_table_in)
	metadata_table.remove_empty_columns()

	taxonomy = NcbiTaxonomy(options.ncbi_reference_directory, False, logger)

	cluster_file = os.path.join(options.project_directory, options.file_cluster_mg_16s)
	mothur_cluster = MothurCluster(logger=logger)
	mothur_cluster.read(cluster_file)

	taxonomy_cluster = TaxonomicCluster(mothur_cluster, taxonomy, logger)

	column_name_unpublished_genomes_id = metadata_table.get_column(options.column_name_unpublished_genomes_id)
	if column_name_unpublished_genomes_id is None:
		logger.error("Meta data file does not contain the required header '{}'".format(column_name_unpublished_genomes_id))
		sys.exit(1)

	taxonomic_prediction(options, metadata_table, mothur_cluster, taxonomy_cluster, taxonomy, logger)
	set_otu_id(options, metadata_table, mothur_cluster, logger)
	metadata_table.write(options.metadata_table_out)


def taxonomic_prediction(options, metadata_table, mothur_cluster, taxonomy_cluster, taxonomy, logger):
	column_cutoff = metadata_table.get_empty_column()
	column_ncbi_prediction = metadata_table.get_empty_column()
	column_science_name = metadata_table.get_empty_column()
	column_novelty = metadata_table.get_empty_column()
	column_name_unpublished_genomes_id = metadata_table.get_column(options.column_name_unpublished_genomes_id)
	if column_name_unpublished_genomes_id is None:
		logger.error("Meta data file does not contain the required header '{}'".format(column_name_unpublished_genomes_id))
		sys.exit(1)
	#_____statistic = {}
	number_of_genomes = len(column_name_unpublished_genomes_id)
	classification_distance = str(options.classification_distance_minimum)
	all_done = False
	sorted_lists_of_cutoffs = mothur_cluster.get_sorted_lists_of_cutoffs()
	for cluster_cutoff in sorted_lists_of_cutoffs:
		if all_done:
			break
		if cluster_cutoff == "unique":
			continue
		if float(cluster_cutoff) < float(classification_distance):
			continue
		logger.info("#threshold {}".format(cluster_cutoff))

		#if cluster_cutoff not in _____statistic:
		#	_____statistic[cluster_cutoff] = {"sname": metadata_table.get_empty_column(), "novelty": metadata_table.get_empty_column(), "support": metadata_table.get_empty_column() }

		all_done = True
		for row_index in range(0, number_of_genomes):
			unpublished_genome_ids = column_name_unpublished_genomes_id[row_index]
			if not mothur_cluster.element_exists(cluster_cutoff, unpublished_genome_ids):
				# if no marker gene was found it will not be in the clustering
				continue
			if unpublished_genome_ids == "" or column_ncbi_prediction[row_index] != "":
				continue
			all_done = False
			separator = ""
			list_of_cluster_id, list_of_cluster = mothur_cluster.get_cluster_of_cutoff_of_element(cluster_cutoff, unpublished_genome_ids)
			if len(list_of_cluster) > 1:
				separator = ";"
			predicted__ncbi = []
			predicted_science_name = []
			predicted_novelty = []
			#list_support = []
			for cluster in list_of_cluster:
				#ncbi_prediction, novelty = taxonomy_cluster.get_cluster_ncbi_tax_prediction(cluster, column_name_unpublished_genomes_id, unpublished_genome_ids)
				ncbi_prediction, novelty, support = taxonomy_cluster.predict_tax_id_of(cluster, column_name_unpublished_genomes_id, unpublished_genome_ids)
				if ncbi_prediction is not None and ncbi_prediction not in predicted__ncbi:
					column_cutoff[row_index] = str(cluster_cutoff)
					predicted__ncbi.append(ncbi_prediction)
					predicted_science_name.append(taxonomy.get_scientific_name(ncbi_prediction))
					predicted_novelty.append("new_" + novelty)

			#		list_support.append(support)
			#_____statistic[cluster_cutoff]["sname"][row_index] = separator.join(predicted_science_name)
			#_____statistic[cluster_cutoff]["novelty"][row_index] = separator.join(predicted_novelty)
			#_____statistic[cluster_cutoff]["support"][row_index] = separator.join(list_support)

			column_ncbi_prediction[row_index] = separator.join(predicted__ncbi)
			column_science_name[row_index] = separator.join(predicted_science_name)
			column_novelty[row_index] = separator.join(predicted_novelty)

	metadata_table.set_column(column_cutoff, options.column_name_cutoff)
	metadata_table.set_column(column_ncbi_prediction, options.column_name_cluster_prediction)
	metadata_table.set_column(column_science_name, options.column_name_cluster_scientific_name)
	metadata_table.set_column(column_novelty, options.column_name_cluster_novelty)

	#for cluster_cutoff in sorted_lists_of_cutoffs:
	#	if cluster_cutoff == "unique":
	#		continue
	#	metadata_table.set_column(_____statistic[cluster_cutoff]["support"], "{}_support".format(cluster_cutoff))
	#for cluster_cutoff in sorted_lists_of_cutoffs:
	#	if cluster_cutoff == "unique":
	#		continue
	#	metadata_table.set_column(_____statistic[cluster_cutoff]["sname"], "{}_name".format(cluster_cutoff))
	#for cluster_cutoff in sorted_lists_of_cutoffs:
	#	if cluster_cutoff == "unique":
	#		continue
	#	metadata_table.set_column(_____statistic[cluster_cutoff]["novelty"], "{}_novelty".format(cluster_cutoff))

	logger.info("Taxonomic prediction finished")


def set_otu_id(options, metadata_table, mothur_cluster, logger):
	list_of_unclustered_elements = set()
	column_name_unpublished_genomes_id = metadata_table.get_column(options.column_name_unpublished_genomes_id)
	number_of_genomes = len(column_name_unpublished_genomes_id)
	otu_distance = options.otu_distance
	column_otu_id = metadata_table.get_empty_column()
	sorted_lists_of_cutoffs = mothur_cluster.get_sorted_lists_of_cutoffs()
	for cluster_cutoff in sorted_lists_of_cutoffs:
		if cluster_cutoff == "unique" or float(cluster_cutoff) != float(otu_distance):
			continue
		for row_index in range(0, number_of_genomes):
			unpublished_genome_ids = column_name_unpublished_genomes_id[row_index]
			if not mothur_cluster.element_exists(cluster_cutoff, unpublished_genome_ids):
				# if no marker gene was found it will not be in the clustering
				list_of_unclustered_elements.add(unpublished_genome_ids)
				continue
			separator = ""
			list_of_cluster_id, list_of_cluster = mothur_cluster.get_cluster_of_cutoff_of_element(cluster_cutoff, unpublished_genome_ids)
			if len(list_of_cluster_id) > 1:
				separator = ";"
			column_otu_id[row_index] = separator.join([str(otu_id) for otu_id in sorted(set(list_of_cluster_id))])
	metadata_table.set_column(column_otu_id, options.column_name_otu_id)
	if len(list_of_unclustered_elements) > 0:
		logger.warning("No cluster found for {} ids!".format(len(list_of_unclustered_elements)))
	logger.info("OTU finished")