__author__ = 'hofmann'

from source.argumenthandler import ArgumentHandler
from source.mothurcluster import MothurCluster
from source.taxonomy.ncbitaxonomy import NcbiTaxonomy
from source.metatable import MetaTable
from source.anim import ANIm
from source.logger import Logger
import os


def my_main(options, logger=None):
	if logger is None:
		logger = Logger("ANI prediction")
	#ranks = args.ranks.strip().split(',')

	assert isinstance(options, ArgumentHandler)

	column_name_ani_novelty = options.column_name_ani_novelty
	column_name_ani_taxid = options.column_name_ani_compare
	column_name_ani_scientific_name = options.column_name_ani_scientific_name
	column_name_ani = options.column_name_ani

	if not os.path.isfile(options.metadata_table_out):
		logger.error("Mothur file with list of clusters not found at: '{}'".format(options.metadata_table_out))
		return False

	metadata_table = MetaTable(logger=logger)
	# metadata_table_out, need data from previous steps
	metadata_table.read(options.metadata_table_out)

	logger.info("Loading taxonomic database: '{}'".format(options.ncbi_reference_directory))
	taxonomy = NcbiTaxonomy(options.ncbi_reference_directory, False, logger)

	cluster_file = os.path.join(options.project_directory, options.file_cluster_mg_16s)
	mothur_cluster = MothurCluster(options.precision, logger=logger)
	mothur_cluster.read(cluster_file)

	ani_scientific_name_column = metadata_table.get_empty_column()
	ani_prediction_novelty_column = metadata_table.get_empty_column()
	ani_prediction_column = metadata_table.get_empty_column()
	ani_column = metadata_table.get_empty_column()
	#cutoff_column = metadata_table.get_column(options.column_name_cutoff)
	unpublished_genome_ids_column = metadata_table.get_column(options.column_name_unpublished_genomes_id)

	# TODO: get rid of import package
	nucmer_exe = "importpackage mummer;nucmer"
	with ANIm(options.input_genomes_file, options.input_reference_file, options.temp_directory, nucmer_exe=nucmer_exe, logger=logger, pool_size=options.processors) as ani_calculator:
		list_of_clusters = mothur_cluster.get_clusters_of_elements(options.distance_cutoff, unpublished_genome_ids_column)
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

	metadata_table.set_column(ani_column, column_name_ani)
	metadata_table.set_column(ani_prediction_novelty_column, column_name_ani_novelty)
	metadata_table.set_column(ani_prediction_column, column_name_ani_taxid)
	metadata_table.set_column(ani_scientific_name_column, column_name_ani_scientific_name)
	metadata_table.write(options.metadata_table_out)
	logger.info("ANI prediction finished")
	return True