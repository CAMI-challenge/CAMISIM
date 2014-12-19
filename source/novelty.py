# -*- coding: utf-8 -*-
import sys
import os
import argparse
from taxonomy import ncbitaxonomy
import metatable
from logger import Logger

__author__ = 'jessika'

"""
	Step 4: Novelty_category determination (Jessika)

	NAME Novelty.py

	FILE /net/metagenomics/projects/cami_2014/reference_data_preparation/Novelty.py

	DESCRIPTION

	Define the taxonomic 'novelty category' for a list of genomes with given NCBI_IDs (taxon IDs, extracted from the metafile.csv) relative to a list of reference strain sequences (extracted from a directory which includes NCBI reference sequences named [taxID].[nr].fna). The taxon ID reflects the lowest rank until which the genome could be placed in the reference taxonomy, meaning that if the rank is above strain level, it is likely new up to the rank below where the taxon ID belongs to. The 'novelty_category' specifies the highest taxonomic rank that a sequenced strain belongs to for which there is no sequenced genome yet in the list of reference strain sequences. This can be even a higher rank than one minus the rank of the taxon ID, as not all taxa which are known have sequenced genomes available. The output is ADD DeTAILS.


	Algorithm:

	if taxon is a strain
		is included  in the reference?  return 'None': return new_strain

	for each rank  from the taxons rank up to rank 'superkingdom':
		 is a taxon from the reference a child of the taxon? return new_(rank -1) // this means that there is a sequenced genomes from this taxon in the reference, and one rank below the rank is new (as apparently not defined in reference taxonomy yet and detectable by us using 16S analysis)

	return new_'superkingdom'

	TODO
		None

	Usage: Novelty.py -ref [directory including ncbi reference sequences] -db [ncbi database] -meta [metafile including column NCBI_ID] -o [output file]

"""


class Novelty():
	def __init__(self, taxonomy=None, logger=None, column_name_ncbi_id="NCBI_ID", column_name_novelty="NOVELTY_CATEGORY"):
		"""
			@param sqlite_db: usually file named "ncbitax_sqlite.db"
		"""
		self._logger = logger
		if logger is None:
			self._logger = Logger("Novelty")

		if taxonomy:
			self._tax = taxonomy
		else:
			self._logger.error("[Novelty] No taxonomy given.")
			sys.exit(1)

		self._column_name_ncbi_id = column_name_ncbi_id
		self._column_name_novelty = column_name_novelty
		#self.ranks = ['strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom', 'root']
		self._ranks = ['strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
		#wrapping:
		self._included_parents_at_rank = dict()
		#for rank in self._ranks:
		#	self._included_parents_at_rank[rank] = set()
		#self._included_parents_at_rank['no rank'] = set()  # includes no ranks and all strain ids
		#self._included_parents_at_rank['root'].add(1)

	def read_reference(self, refernce_ncbi_id_set, excluded=None):
		"""
			extracts all taxonomic IDs from the reference directory and stores all contained IDs for each rank
			@param referenceIdsSet:   including reference IDs
		"""

		self._logger.info("[Novelty] Extracting included parents at each rank for each reference ID.. This may take a while.")

		for ncbi_id in refernce_ncbi_id_set:
			if excluded is not None and ncbi_id in excluded:
				continue

			taxonomic_lineage = self._tax.get_lineage_of_legal_ranks(ncbi_id, ranks=self._ranks)
			if taxonomic_lineage is None:
				continue

			for rank_index in xrange(len(self._ranks)):
				if taxonomic_lineage[rank_index] is None:
					continue
				rank = self._ranks[rank_index]
				if rank not in self._included_parents_at_rank:
					self._included_parents_at_rank[rank] = set()
				self._included_parents_at_rank[rank].add(taxonomic_lineage[rank_index])

		self._logger.info("[Novelty] Reference processing done.")

	def compute_novelty_for_metafile(self, in_meta_file, out_meta_file):
		"""
			computes the novelty_category for each NCBI ID in the metafile and updates it to the output file
			(Note that the metafile must include a header with column name 'NCBI_ID'
								whereas novelty_category is added if it does not exist)
			@param in_meta_file: usually file named 'metadata_table_[version].csv'#
			@param out_meta_file:  file for the output
		"""

		self._logger.info("[Novelty] Processing information from metafile: '{}'".format(in_meta_file))
		meta_table = metatable.MetaTable(logger=self._logger)
		meta_table.read(in_meta_file)
		self.compute_novelty(meta_table)
		meta_table.write(out_meta_file)

	def compute_novelty(self, meta_table):
		"""
			computes the novelty_category for each NCBI ID in the metafile and updates it to the output file
			(Note that the metafile must include a header with column name 'NCBI_ID'
								whereas novelty_category is added if it does not exist)
			@param in_meta_file: usually file named 'metadata_table_[version].csv'#
			@param out_meta_file:  file for the output
		"""

		column_ncbi_id = meta_table.get_column(self._column_name_ncbi_id)
		column_novelty_category = meta_table.get_column(self._column_name_novelty)
		if column_novelty_category is None:
			column_novelty_category = meta_table.get_empty_column()

		number_of_rows = meta_table.get_number_of_rows()
		for row_index in range(number_of_rows):
			list_novelty = []
			if column_ncbi_id[row_index] == '':
				continue
			if column_novelty_category[row_index] != '':
				list_novelty.append(column_novelty_category[row_index][4:])
			list_ncbi_id = column_ncbi_id[row_index].split(';')
			for new_ncbi_id in list_ncbi_id:
				if not new_ncbi_id.isdigit():
					continue
				novelty = self.get_novelty(new_ncbi_id)
				if novelty is None:
					continue
				list_novelty.append(novelty)
				#self._logger.info("[Novelty] {} is not included at rank {}".format(new_ncbi_id, novelty))
			#column_novelty_category[row_index] = ','.join(["new_" + str(novelty) for novelty in list_novelty])
			column_novelty_category[row_index] = "new_" + self.get_lowest_rank(set(list_novelty))

		meta_table.set_column(column_novelty_category, self._column_name_novelty)

	def get_novelty(self, ncbi_id):
		"""
			Compute novelty_category for a new taxon ID according to the following algorithm
			@param ncbi_id: NCBI ID for which the novelty_category should be computed
			@return: novelty category or None
			@rtype: str
		"""

		taxonomic_lineage = self._tax.get_lineage_of_legal_ranks(ncbi_id, ranks=self._ranks)

		if taxonomic_lineage is None:
			self._logger.warning("[Novelty] {} not included in the taxonomy.".format(ncbi_id))
			return None

		for rank_index in xrange(len(self._ranks)):
			if taxonomic_lineage[rank_index] is None:
				continue
			if taxonomic_lineage[rank_index] in self._included_parents_at_rank[self._ranks[rank_index]]:
				return self._ranks[rank_index-1]

		return None
		#return 'superkingdom'

	def get_taxonomic_ids_from_directory(self, directory):
		"""
			search a directory for all files with taxonomic IDS and save them as a set
			@param directory: directory containing sequences named [ID].[nr].fna
			@return: set of the IDs
			@rtype: set
		"""
		directory_list = os.listdir(directory)
		tax_ids = set()
		for item in directory_list:
			if not os.path.isfile(os.path.join(directory, item)):
				continue
			tid = item.split(".")[0]
			if '_' in tid:
				tid = tid.split('_')[0]
			tax_ids.add(tid)

		return tax_ids

	def get_lowest_rank(self, list_of_ranks=set(), ranks=None):
		if ranks is None:
			ranks = self._ranks
		if len(list_of_ranks) == 0:
			return ''
		lowest_rank = list_of_ranks.pop()
		for rank in list_of_ranks:
			if ranks.index(rank) < ranks.index(lowest_rank):
				lowest_rank = rank
		return lowest_rank


#***************** Main ************************
if __name__ =='__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-dref', action='store', help='directory with reference sequence files (named taxID.1.fna)', required=True)
	parser.add_argument('-db', "--ncbi_reference_directory", action='store', help='complete path to ncbi reference directory', required=True)
	parser.add_argument('-mi', "--metadata_file_in", action='store', required=True, help='metadatafile input')
	parser.add_argument('-mo', "--metadata_file_out", action='store', required=True, help='metadatafile output')

	options = parser.parse_args()

	logger = Logger("Novelty")

	taxonomy = ncbitaxonomy.NcbiTaxonomy(options.ncbi_reference_directory, False, logger)

	Nov = Novelty(taxonomy, logger=logger)
	refernceIdsSet = Nov.get_taxonomic_ids_from_directory(options.dref)
	Nov.read_reference(refernceIdsSet)
	Nov.compute_novelty_for_metafile(options.metadata_file_in, options.metadata_file_out)
	logger.info("[Novelty] Done")