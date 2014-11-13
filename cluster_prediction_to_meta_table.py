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
	logger.info("Reading metadata table file: '{}'".format(options.metadata_table_in))
	metadata_table = MetaTable('\t', logger)
	metadata_table.read(options.metadata_table_in)
	metadata_table.remove_empty_columns()

	logger.info("Loading taxonomic database: '{}'".format(options.ncbi_reference_directory))
	taxonomy = NcbiTaxonomy(options.ncbi_reference_directory, False, logger)

	cluster_file = os.path.join(options.project_directory, options.file_cluster_mg_16s)
	logger.info("Reading mothur cluster file: '{}'".format(cluster_file))
	mothur_cluster = MothurCluster()
	with open(cluster_file, 'r') as file_handle:
		mothur_cluster.read_mothur_clustering_file(file_handle)

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
		logger.info("#cutoff {}".format(cluster_cutoff))
		all_done = True
		for row_index in range(0, number_of_genomes):
			unpublished_genome_ids = column_name_unpublished_genomes_id[row_index]
			if not mothur_cluster.element_exists(cluster_cutoff, unpublished_genome_ids):
				#logger.warning("{}: No cluster found for id '{}'".format(cluster_cutoff, unpublished_genome_ids))
				continue
			if unpublished_genome_ids == "" or column_ncbi_prediction[row_index] != "":
				continue
			all_done = False
			separator = ""
			list_of_cluster_id, list_of_cluster = mothur_cluster.get_cluster_of_cutoff_of_element(cluster_cutoff, unpublished_genome_ids)
			if len(list_of_cluster_id) > 1:
				separator = ";"
			otu_ncbi = []
			otu_name = []
			otu_novelty = []
			for index in range(0, len(list_of_cluster_id)):
				cluster = list_of_cluster[index]
				if cluster is None:
					continue
				ncbi_prediction, novelty = taxonomy_cluster.get_cluster_ncbi_tax_prediction(cluster, column_name_unpublished_genomes_id, unpublished_genome_ids)
				if ncbi_prediction is not None and ncbi_prediction not in otu_ncbi:
					column_cutoff[row_index] = str(cluster_cutoff)
					otu_ncbi.append(ncbi_prediction)
					otu_name.append(taxonomy.get_scientific_name(ncbi_prediction))
					otu_novelty.append("new_" + novelty)
			column_ncbi_prediction[row_index] = separator.join(otu_ncbi)
			column_science_name[row_index] = separator.join(otu_name)
			column_novelty[row_index] = separator.join(otu_novelty)

	metadata_table.set_column(column_cutoff, options.column_name_cutoff)
	metadata_table.set_column(column_ncbi_prediction, options.column_name_cluster_prediction)
	metadata_table.set_column(column_science_name, options.column_name_cluster_scientific_name)
	metadata_table.set_column(column_novelty, options.column_name_cluster_novelty)
	logger.info("Taxonomic prediction finished")


def set_otu_id(options, metadata_table, mothur_cluster, logger):
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
				logger.warning("{}: No cluster found for id '{}'".format(cluster_cutoff, unpublished_genome_ids))
				continue
			separator = ""
			list_of_cluster_id, list_of_cluster = mothur_cluster.get_cluster_of_cutoff_of_element(cluster_cutoff, unpublished_genome_ids)
			if len(list_of_cluster_id) > 1:
				separator = ";"
			column_otu_id[row_index] = separator.join([str(otu_id) for otu_id in list_of_cluster_id])
	metadata_table.set_column(column_otu_id, options.column_name_otu_id)
	logger.info("OTU finished")