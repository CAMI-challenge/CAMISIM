__author__ = 'hofmann'

#from source.ArgumentHandler import ArgumentHandler
from source.MothurCluster import MothurCluster
from source.NcbiTaxonomy.NcbiTaxonomy import NcbiTaxonomy
from source.MetaTable import MetaTable
from source.ANIm import ANIm
from source.Logger import Logger
import sys
import os


def my_main(options, logger=None):
	if logger == None:
		logger = Logger("ANI prediction")
	#ranks = args.ranks.strip().split(',')

	column_name_ani_novelty = options.column_name_ani_novelty
	column_name_ani_taxid = options.column_name_ani_compare
	column_name_ani_scientific_name = options.column_name_ani_scientific_name
	column_name_ani = options.column_name_ani

	logger.info("Reading metadata table file: '{}".format(options.metadata_table_in))
	metadata_table = MetaTable(options.metadata_table_in)

	logger.info("Loading taxonomic database: '{}'".format(options.ncbi_reference_directory))
	taxonomy = NcbiTaxonomy(options.ncbi_reference_directory, False, logger)

	cluster_file = os.path.join(options.project_directory, options.file_cluster_mg_16s)
	logger.info("Reading mothur cluster file: '{}'".format(cluster_file))
	mothur_cluster = MothurCluster()
	with open(cluster_file, 'r') as file_handle:
		mothur_cluster.read_mothur_clustering_file(file_handle)

	ani_scientific_name_column = metadata_table.get_empty_column()
	ani_prediction_novelty_column = metadata_table.get_empty_column()
	ani_prediction_column = metadata_table.get_empty_column()
	ani_column = metadata_table.get_empty_column()
	cutoff_column = metadata_table.get_column(options.column_name_cutoff)
	unpublished_genome_ids_column = metadata_table.get_column(options.column_name_unpublished_genomes_id)

	# TODO: get rid of import package
	nucmer_exe = "importpackage mummer;nucmer"
	with ANIm(options.input_genomes_file, options.input_reference_file, None, nucmer_exe=nucmer_exe, logger=logger, pool_size=options.processors) as ani_calculator:
		list_of_clusters = mothur_cluster.get_clusters_of_elements(cutoff_column, unpublished_genome_ids_column)
		#logger.info("OTUS: {}".format(otus))
		#sys.exit(0)
		for unknown_genomes_id in unpublished_genome_ids_column:
			if list_of_clusters[unknown_genomes_id] is None:
				continue
			ani_calculator.add_nucmer_cmd_lines(list_of_clusters[unknown_genomes_id], [unknown_genomes_id])

		min_lengths, min_sim_errors, min_perc_ids, min_perc_aln, ncbi, ani_ish = ani_calculator.calculate_minimum_anim()

		number_of_unknown = len(unpublished_genome_ids_column)
		for row_index in range(0, number_of_unknown):
			unknown_genomes_id = unpublished_genome_ids_column[row_index]
			if unknown_genomes_id in ani_ish and ani_ish[unknown_genomes_id] > 0:
				if float(ani_ish[unknown_genomes_id]) > 0.98:
					ani_prediction_novelty_column[row_index] = "same_strain"
				elif float(ani_ish[unknown_genomes_id]) > 0.96:
					ani_prediction_novelty_column[row_index] = "same_species"
				ani_prediction_column[row_index] = ncbi[unknown_genomes_id]
				ani_column[row_index] = str(ani_ish[unknown_genomes_id])
				science_name = taxonomy.get_scientific_name(ncbi[unknown_genomes_id])
				if science_name is not None:
					ani_scientific_name_column[row_index] = science_name

	metadata_table.set_column(column_name_ani, ani_column)
	metadata_table.set_column(column_name_ani_novelty, ani_prediction_novelty_column)
	metadata_table.set_column(column_name_ani_taxid, ani_prediction_column)
	metadata_table.set_column(column_name_ani_scientific_name, ani_scientific_name_column)
	metadata_table.write(options.metadata_table_out)
	logger.info("ANI prediction finished")
	sys.exit(0)