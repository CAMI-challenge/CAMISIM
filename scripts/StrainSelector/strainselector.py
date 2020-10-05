#!/usr/bin/env python3

__author__ = 'hofmann'
__original_author__ = 'eik.dahms@uni-duesseldorf.de'
__version__ = '0.0.2'

import random
from scripts.loggingwrapper import DefaultLogging
from scripts.MetaDataTable.metadatatable import MetadataTable


class StrainSelector(DefaultLogging):
	"""
	algorithm:
	----------

	# of strains per category = C = # strains / # categories

	- read the input table and sort by category and OTU
	- sort categories by increasing size

	- drawing from categories:
		IF (size of category) < C :
		the whole category is taken and C is recalculated
		for the remaining categories

		ELSE:
		go to category:
			- get OTUs random
				FOR EACH:
				- draw random strains until num_per_otu or C is reached
				- non-drawn strains are saved
				- IF end is reached and # drawn strains is not C :
					take (C - #drawn strains) from saved strains
	"""
	def __init__(
		self, column_name_genome_id="genome_ID", column_name_otu="OTU", column_name_novelty_category="novelty_category",
		logfile=None, verbose=True, debug=False, seed=None):
		"""


		@param column_name_genome_id: Column name for genome ids
		@type column_name_genome_id: str
		@param column_name_otu: Column name for
		@type column_name_otu: str
		@param column_name_novelty_category: Column name for
		@type column_name_novelty_category: str
		@param logfile: file handler or file path to a log file
		@type logfile: file | FileIO | StringIO | str
		@param verbose: Not verbose means that only warnings and errors will be past to stream
		@type verbose: bool
		@param debug: Display debug messages
		@type debug: bool
		@param seed: Seed for random module

		@return: Nothing
		@rtype: None
		"""
		super(StrainSelector, self).__init__(logfile=logfile, verbose=verbose, debug=debug)
		assert isinstance(column_name_genome_id, str)
		assert isinstance(column_name_otu, str)
		assert isinstance(column_name_novelty_category, str)
		if seed is not None:
			random.seed(seed)
		self._cats = 0
		self._draw = 0
		self._per_cat = 0
		self._rest = 0
		self._per_otu = 0
		self._column_name_genome_id = column_name_genome_id
		self._column_name_otu = column_name_otu
		self._column_name_novelty_category = column_name_novelty_category

	def get_drawn_genome_id(self, metadata_table, number_of_strains, number_of_strains_per_otu):
		"""
		Get a list of drawn genome ids.

		@param metadata_table: Metadata containing genome id, otu and a novelty category
		@type metadata_table: MetadataTable
		@param number_of_strains: Amount of strains to be drawn
		@type number_of_strains: int
		@param number_of_strains_per_otu:
		@type number_of_strains_per_otu: int

		@return: list of drawn genome ids @raise ValueError:
		@rtype: list[str|unicode]
		"""
		assert isinstance(metadata_table, MetadataTable)
		assert isinstance(number_of_strains, int)
		assert isinstance(number_of_strains_per_otu, int)
		genomes_read_from_metadata = 0
		lost_lines = 0

		categories = dict()

		required_headers = [self._column_name_genome_id, self._column_name_otu, self._column_name_novelty_category]
		if not metadata_table.validate_column_names(required_headers):
			msg = 'A header is missing, needed: {}'.format(', '.join(required_headers))
			self._logger.error(msg)
			raise ValueError(msg)

		column_genome_id = metadata_table.get_column(self._column_name_genome_id)
		column_otu = metadata_table.get_column(self._column_name_otu)
		column_novelty = metadata_table.get_column(self._column_name_novelty_category)
		row_count = metadata_table.get_number_of_rows()

		# read metadata from provided file
		for row_index in range(0, row_count):
			if column_genome_id[row_index] == '' or column_otu[row_index] == '' or column_novelty[row_index] == '':
				lost_lines += 1
				continue
			genomes_read_from_metadata += 1
			novelty = column_novelty[row_index]
			if novelty not in categories:
				categories[novelty] = NoveltyCategory(novelty)
			categories[novelty].add_strain(column_otu[row_index], column_genome_id[row_index])

		if lost_lines > 0:
			msg = "Invalid lines/ unusable genomes: {}".format(lost_lines)
			self._logger.debug(msg)

		if genomes_read_from_metadata < number_of_strains:
			msg = 'Not enough data to draw.'
			self._logger.error(msg)
			raise ValueError(msg)

		# sample OTUs from the pool of available OTUs
		drawn_genome_id = self._draw_strains(categories, number_of_strains, number_of_strains_per_otu)
		return drawn_genome_id

	def _recalc(self, number_of_drawn):
		"""
		Adjust values based on the number of strains drawn from the previous category

		@param number_of_drawn: Number of drawn strains from the previous category
		@rtype: None
		"""
		assert isinstance(number_of_drawn, int)
		self._draw -= number_of_drawn
		self._cats -= 1
		self._per_cat = int(self._draw/self._cats)
		self._rest = self._draw % self._cats

	def _substract(self):
		"""
		Adjust values by subtracting the maximum number of strains that can be drawn from a category.

		@rtype: None
		"""
		self._draw -= self._per_cat
		if self._rest > 0:
			self._rest -= 1

	def _per_category(self):
		"""
		Get maximum number of strains valid to be drawn from a category.

		@rtype: int
		"""
		return self._per_cat + (self._rest > 0)

	def _draw_strains(self, categories, number_of_strains_to_draw, max_amount_of_strains_per_otu):
		"""
		Get list of drawn strains

		@param categories:
		@type: dict[str, NoveltyCategory]
		@param number_of_strains_to_draw: Amount of strains to be drawn.
		@type: int
		@param max_amount_of_strains_per_otu: Maximum number of strains to be drawn from one category.
		@type: int

		@return: @raise ValueError:
		@rtype: list[str|unicode]
		"""
		assert isinstance(categories, dict)
		assert isinstance(number_of_strains_to_draw, int)
		assert isinstance(max_amount_of_strains_per_otu, int)
		number_of_novelty_categories = len(categories)

		self._cats = number_of_novelty_categories
		self._draw = number_of_strains_to_draw
		self._per_cat = int(number_of_strains_to_draw/number_of_novelty_categories)
		self._rest = number_of_strains_to_draw % number_of_novelty_categories
		self._per_otu = max_amount_of_strains_per_otu

		drawn_genome_id = []
		for x in sorted(categories, key=lambda name: categories[name].get_strain_amount()):
			if categories[x].get_strain_amount() < self._per_category():
				# draw all strains
				drawn_genome_id += categories[x].get_all_strains()
				self._recalc(categories[x].get_strain_amount())
			else:
				# draw subset of strains
				drawn_genome_id += categories[x].draw_strains(self._per_category(), self._per_otu)
				self._substract()

		if len(drawn_genome_id) < number_of_strains_to_draw:
			msg = "Could only draw {} samples instead of {}".format(len(drawn_genome_id), number_of_strains_to_draw)
			self._logger.error(msg)
			raise ValueError(msg)
		return drawn_genome_id


class NoveltyCategory(object):
	"""
	Handling list of OTU belonging to one novelty category
	"""

	def __init__(self, name):
		assert isinstance(name, str)
		self._name = name
		self._number_of_strains = 0
		self._otu_list = {}

	def get_name(self):
		"""
		Return name of novelty category

		@return: Name of novelty category
		@rtype: str | unicode
		"""
		return self._name

	def get_strain_amount(self):
		"""
		Return amount of strains within this category

		@return: Amount of strains within this category
		@rtype: int
		"""
		return self._number_of_strains

	def draw_strains(self, total, limit_per_otu):
		"""
		Draw subset of strains from this category

		@param total: Amount to be drawn.
		@type total: int
		@param limit_per_otu: Maximum number of strains to be drawn from an otu
		@type limit_per_otu: int

		@return: List of drawn strain id
		@rtype: list[str|unicode]
		"""
		assert isinstance(total, int)
		assert isinstance(limit_per_otu, int)
		assert total <= self.get_strain_amount()
		assert limit_per_otu > 0
		drawn_strain = []
		overhead = []
		drawn_strain_count_overall = 0
		for otu_id in random.sample(self._otu_list.keys(), len(self._otu_list)):
			drawn_strain_count_otu = 0
			for strain_id in random.sample(self._otu_list[otu_id], len(self._otu_list[otu_id])):
				if drawn_strain_count_otu < limit_per_otu and drawn_strain_count_overall < total:
					drawn_strain.append(strain_id)
					drawn_strain_count_otu += 1
					drawn_strain_count_overall += 1
			overhead += list(set(self._otu_list[otu_id])-set(drawn_strain))

		if drawn_strain_count_overall < total:
			# out += overhead[0:total-drawn_strain_cound_overall]
			#print "\n", total-drawn_strain_count_overall, len(overhead), len(drawn_strain), "\n"
			drawn_strain += random.sample(overhead, total-drawn_strain_count_overall)
		return drawn_strain

	def get_all_strains(self):
		"""
		Return list of all strains from this category

		@return: List of drawn strain id
		@rtype: list[str|unicode]
		"""
		out = []
		for otu_id in self._otu_list:
			for strain_id in self._otu_list[otu_id]:
				out.append(strain_id)
		return out

	def add_strain(self, otu_id, strain_id):
		"""
		Add a strain

		@param otu_id: Identifier of otu.
		@type otu_id: str | unicode
		@param strain_id: Identifier of strain.
		@type strain_id: str | unicode

		@rtype: None
		"""
		assert isinstance(otu_id, str)
		assert isinstance(strain_id, str)
		if otu_id not in self._otu_list:
			self._otu_list[otu_id] = []
		self._otu_list[otu_id].append(strain_id)
		self._number_of_strains += 1
