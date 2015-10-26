# -*- coding: utf-8 -*-

import os
import argparse
from scripts.NcbiTaxonomy.ncbitaxonomy import NcbiTaxonomy
from scripts.MetaDataTable.metadatatable import MetadataTable
from scripts.Validator.validator import Validator

__author__ = 'jessika'


class Novelty(Validator):
	"""
	Novelty_category determination (Jessika)

	NAME novelty.py

	FILE /net/metagenomics/projects/cami_2014/reference_data_preparation/Novelty.py

	DESCRIPTION

	Define the taxonomic 'novelty category' for a list of genomes with given NCBI_IDs (taxon IDs, extracted from the metafile.csv)
	relative to a list of reference strain sequences (extracted from a directory which includes NCBI reference sequences named [taxID].[nr].fna).
	The taxon ID reflects the lowest rank until which the genome could be placed in the reference taxonomy, meaning that if the rank is above strain level,
	it is likely new up to the rank below where the taxon ID belongs to.
	The 'novelty_category' specifies the highest taxonomic rank that a sequenced strain belongs to
	for which there is no sequenced genome yet in the list of reference strain sequences.
	This can be even a higher rank than one minus the rank of the taxon ID, as not all taxa which are known have sequenced genomes available.
	The output is ADD DeTAILS.


	Algorithm:

	if taxon is a strain
		is included  in the reference?  return 'None': return new_strain

	for each rank  from the taxons rank up to rank 'superkingdom':
		is a taxon from the reference a child of the taxon? return new_(rank -1)
		// this means that there is a sequenced genomes from this taxon in the reference,
		// and one rank below the rank is new (as apparently not defined in reference taxonomy yet and detectable by us using 16S analysis)

	return new_'superkingdom'

	TODO
		None

	Usage: Novelty.py -ref [directory including ncbi reference sequences] -db [ncbi database] -meta [metafile including column NCBI_ID] -o [output file]
	"""

	_label = "Novelty"
	_ranks = ['strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']

	def __init__(
		self, taxonomy, column_name_ncbi_id="NCBI_ID", column_name_novelty="NOVELTY_CATEGORY",
		separator="\t", logfile=None, verbose=True, debug=False):
		"""
		Constructor

		@param taxonomy: Handle to NcbiTaxonomy
		@type taxonomy: NcbiTaxonomy
		@param column_name_ncbi_id: Header name of the ncbi column
		@type column_name_ncbi_id: str | unicode
		@param column_name_novelty: Header name of the novelty prediction column
		@type column_name_novelty: str | unicode
		@param separator: Column separator
		@type separator: str | unicode
		@param logfile: File handler or file path to a log file
		@type logfile: file | FileIO | StringIO | basestring
		@param verbose: Not verbose means that only warnings and errors will be past to stream
		@type verbose: bool
		@param debug: Display debug messages
		@type debug: bool
		"""
		super(Novelty, self).__init__(logfile=logfile, verbose=verbose, debug=debug)
		assert isinstance(taxonomy, NcbiTaxonomy), "No taxonomy given."
		assert isinstance(column_name_ncbi_id, basestring)
		assert isinstance(column_name_novelty, basestring)
		assert isinstance(separator, basestring)
		self._tax = taxonomy
		self._separator = separator
		self._column_name_ncbi_id = column_name_ncbi_id
		self._column_name_novelty = column_name_novelty

		self._included_parents_at_rank = dict()

	def read_reference(self, set_reference_taxonomic_ids, excluded=None):
		"""
		extracts all taxonomic IDs from the reference directory and stores all contained IDs for each rank
		@param set_reference_taxonomic_ids: List of reference ids
		@type set_reference_taxonomic_ids: set[str|unicode]
		@param excluded: List of reference ids
		@type excluded: set[str|unicode]

		@rtype: None
		"""
		assert isinstance(set_reference_taxonomic_ids, set)
		assert excluded is None or isinstance(excluded, set)
		self._logger.info("Extracting included parents at each rank for each reference ID.. This may take a while.")
		for ncbi_id in set_reference_taxonomic_ids:
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
		self._logger.info("Reference processing done.")

	def compute_novelty_for_metafile(self, in_meta_file, out_meta_file):
		"""
		computes the novelty_category for each NCBI ID in the metafile and updates it to the output file
		(Note that the metafile must include a header with column name 'NCBI_ID'
							whereas novelty_category is added if it does not exist)

		@param in_meta_file: filepath to file named 'metadata_table_[version].csv'#
		@type in_meta_file: str | unicode
		@param out_meta_file: file path of the output
		@type out_meta_file: str | unicode

		@rtype: None
		"""
		assert self.validate_file(in_meta_file)
		assert self.validate_file(out_meta_file)
		self._logger.info("Processing information from metafile: '{}'".format(in_meta_file))
		meta_table = MetadataTable(separator=self._separator, logfile=self._logfile, verbose=self._verbose)
		meta_table.read(in_meta_file, column_names=True)
		self.compute_novelty(meta_table)
		meta_table.write(out_meta_file, column_names=True)

	def compute_novelty(self, meta_table):
		"""
		computes the novelty_category for each NCBI ID in the metafile and updates it to the output file
		(Note that the metafile must include a header with column name 'NCBI_ID'
							whereas novelty_category is added if it does not exist)
		@param meta_table: Handle of a MetadataTable
		@type meta_table: MetadataTable
		"""
		assert isinstance(meta_table, MetadataTable)
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
				# self._logger.info("[Novelty] {} is not included at rank {}".format(new_ncbi_id, novelty))
			# column_novelty_category[row_index] = ','.join(["new_" + str(novelty) for novelty in list_novelty])
			column_novelty_category[row_index] = "new_" + self.get_lowest_rank(set(list_novelty))

		meta_table.insert_column(column_novelty_category, self._column_name_novelty)

	def get_novelty(self, ncbi_id):
		"""
		Compute novelty_category for a new taxon ID according to the following algorithm.

		@param ncbi_id: NCBI ID for which the novelty_category should be computed
		@type ncbi_id: str | unicode

		@return: novelty category or None
		@rtype: str
		"""
		assert isinstance(ncbi_id, basestring)
		taxonomic_lineage = self._tax.get_lineage_of_legal_ranks(ncbi_id, ranks=self._ranks)

		if taxonomic_lineage is None:
			self._logger.warning("{} not included in the taxonomy.".format(ncbi_id))
			return None

		for rank_index in xrange(len(self._ranks)):
			if taxonomic_lineage[rank_index] is None:
				continue
			if taxonomic_lineage[rank_index] in self._included_parents_at_rank[self._ranks[rank_index]]:
				return self._ranks[rank_index - 1]
		return None

	def get_taxonomic_ids_from_directory(self, directory):
		"""
		search a directory for all files with taxonomic IDS and return them as a set

		@param directory: directory containing sequences named [ID].[nr].fna
		@type directory: str | unicode

		@return: set of the IDs
		@rtype: set[str | unicode]
		"""
		assert self.validate_dir(directory)
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

	def get_lowest_rank(self, set_of_ranks, ranks=None):
		"""
		Get the lowest rank form a set of ranks

		@param set_of_ranks: List of ranks found in self._ranks
		@type set_of_ranks: set[str|unicode]
		@param ranks: Ordered list of ranks
		@type ranks: list[str|unicode]

		@return: Name of lowest rank
		@rtype: str | unicode
		"""
		assert isinstance(set_of_ranks, set)
		assert ranks is None or isinstance(ranks, list)
		if ranks is None:
			ranks = self._ranks
		if len(set_of_ranks) == 0:
			return ''
		lowest_rank = set_of_ranks.pop()
		for rank in set_of_ranks:
			if ranks.index(rank) < ranks.index(lowest_rank):
				lowest_rank = rank
		return lowest_rank


# ***************** Main ************************
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-dref', action='store', help='directory with reference sequence files (named taxID.1.fna)', required=True)
	parser.add_argument('-db', "--ncbi_reference_directory", action='store', help='complete path to ncbi reference directory', required=True)
	parser.add_argument('-mi', "--metadata_file_in", action='store', required=True, help='metadatafile input')
	parser.add_argument('-mo', "--metadata_file_out", action='store', required=True, help='metadatafile output')

	options = parser.parse_args()

	taxonomy_handle = NcbiTaxonomy(options.ncbi_reference_directory, False)

	Nov = Novelty(taxonomy_handle)
	refernceIdsSet = Nov.get_taxonomic_ids_from_directory(options.dref)
	Nov.read_reference(refernceIdsSet)
	Nov.compute_novelty_for_metafile(options.metadata_file_in, options.metadata_file_out)
