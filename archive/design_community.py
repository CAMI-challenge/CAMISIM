import sys
import random
import numpy.random as np_random
import os
from collections import Counter

column_name_genome_id = "genome_ID"
column_name_otu = "OTU"
column_name_novelty_category = "novelty_category"
column_name_ncbi = "NCBI_ID"
column_name_source = "source"
filename_distribution_comunity_joint = "distribution_{sample_index}.txt"
cats = 0
draw = 0
per_cat = 0
rest = 0
per_otu = 0
list_of_column_names_from_file=[]
list_of_column_names_new_metadata_table = []
number_of_rows = 0

def get_distribution_file_paths(directory, number_of_samples):
		"""
		Generate directory paths for each sample

		@param directory: Output stream
		@type directory: str | unicode
		@param number_of_samples: Number of samples
		@type number_of_samples: int | long

		@return: list of directories
		@rtype: list[str | unicode]
		"""
		file_path = os.path.join(directory, filename_distribution_comunity_joint)
		return [file_path.format(sample_index=sample_index) for sample_index in range(number_of_samples)]

def get_genome_amounts_geometric_fix(num_real_genomes, max_genome_amount, geometric_probability=0.3):
		"""
		Get amounts of genomes by original genome

		@param num_real_genomes: exact number of real genomes
		@type num_real_genomes: int 
		@param max_genome_amount: Total number of genomes
		@type max_genome_amount: int 

		@return: List of integers representing amount of strains
		@rtype: list[int]
		"""
		assert isinstance(num_real_genomes, int)
		assert isinstance(max_genome_amount, int)

		final_amounts = [1] * num_real_genomes
		index = 0
		while index < len(final_amounts):
			if sum(final_amounts) >= max_genome_amount:
				break
			final_amounts[index] += 1 + np_random.geometric(geometric_probability)
			index += 1

		final_amounts[index-1] -= sum(final_amounts) - max_genome_amount
		return final_amounts

def get_genome_amounts_geometric(probability, max_genome_amount, geometric_probability=0.3):
		"""
		Get amounts of genomes by original genome

		@param probability: Proportion of simulated original genomes
		@type probability: int  | float
		@param max_genome_amount: Total number of genomes
		@type max_genome_amount: int 

		@return: List of integers representing amount of strains
		@rtype: list[int]
		"""
		assert isinstance(probability, (int, float))
		assert 0 <= probability <= 1
		assert isinstance(max_genome_amount, int)

		final_amounts = []
		while sum(final_amounts) < max_genome_amount:
			if random.uniform(0, 1) < probability:
				final_amounts.append(1)
			else:
				amount = 1 + np_random.geometric(geometric_probability)
				final_amounts.append(amount)

		final_amounts[-1] -= sum(final_amounts) - max_genome_amount
		return final_amounts

def get_genome_amounts(probability, max_genome_amount):
		"""
		Get amounts of genomes by original genome

		@param probability: Proportion of simulated original genomes
		@type probability: int  | float
		@param max_genome_amount: Total number of genomes
		@type max_genome_amount: int 

		@return: List of integers representing amount of strains
		@rtype: list[int]
		"""
		assert isinstance(probability, (int, float))
		assert 0 <= probability <= 1
		assert isinstance(max_genome_amount, int)

		genome_amounts = get_genome_amounts_geometric(probability, max_genome_amount)
		diverence = Counter(genome_amounts)[1] / float(len(genome_amounts))
		if max_genome_amount >= 10:
			while abs(diverence - probability) > 0.05:
				# print "need: {}, gotten: {}".format(probability, diverence)
				genome_amounts = get_genome_amounts_geometric(probability, max_genome_amount)
				diverence = Counter(genome_amounts)[1] / float(len(genome_amounts))
		return genome_amounts

def get_genome_amounts(probability, max_genome_amount, num_real_genomes=None, silent=True):
		"""
		Get amounts of genomes by original genome

		@param probability: Proportion of simulated original genomes
		@type probability: int  | float
		@param max_genome_amount: Total number of genomes
		@type max_genome_amount: int 
		@param num_real_genomes: exact number of real genomes
		@type num_real_genomes: int 

		@return:
		@rtype: list[int]
		"""
		assert probability is None or isinstance(probability, (int, float))
		if probability:
			assert 0 <= probability <= 1
		assert isinstance(max_genome_amount, int)
		assert isinstance(num_real_genomes, int)
		assert isinstance(silent, bool)

		if num_real_genomes is not None:
			genome_amounts = get_genome_amounts_geometric_fix(num_real_genomes, max_genome_amount)
		else:
			genome_amounts = get_genome_amounts(probability, max_genome_amount)

		#if not silent:
			#self.print_distribution(genome_amounts)
			#message = "Do you accept this distribution? [y/n]"
			#while not self.get_confirmation(message):
				#if num_real_genomes is not None:
					#genome_amounts = self._get_genome_amounts_geometric_fix(num_real_genomes, max_genome_amount)
				#else:
					#genome_amounts = self._get_genome_amounts(probability, max_genome_amount)
				#self.print_distribution(genome_amounts)
		return genome_amounts

def open_compress(self, file_path, mode='r', compresslevel=5, compression_type=None):
		"""
		Open a file for reading or writing

		@attention: When reading file and compression_type None, type will be guessed.

		@param file_path: Path to file
		@type file_path: str | unicode
		@param mode: mode a file is opened with. 'r' or 'w'
		@type mode: str | unicode
		@param compresslevel: Higher level is slower but likely smaller. 0-9, except zip 0-8.
		@type compresslevel: int
		@param compression_type: "zip", "gz", "bz2",
		@type compression_type: str | unicode

		@return: Return a file object
		@rtype: file
		"""
		assert mode in self._modes, "Unsupported mode '{}'.".format(mode)
		if compression_type is None:
			compression_type = self.get_compression_type(file_path)
		if mode == 'r':
			return self._open[compression_type](file_path, mode=mode)
		elif compression_type == "gz":
			assert self.validate_number(compresslevel, minimum=0, maximum=9)
			return self._open[compression_type](file_path, mode='w', compresslevel=compresslevel)
		elif compression_type == "bz2":
			assert self.validate_number(compresslevel, minimum=0, maximum=9)
			return self._open[compression_type](file_path, mode='w', compresslevel=compresslevel)
		elif compression_type == "zip":
			assert self.validate_number(compresslevel, minimum=0, maximum=8)
			return self._open[compression_type](file_path, mode='w', compression=compresslevel)

def parse_column_names(self, stream_input, separator):
			row = stream_input.readline().rstrip('\n').rstrip('\r')
			list_of_column_names_from_file = row.split(separator)
			assert self._has_unique_columns(list_of_column_names_from_file), "Column names must be unique!"

def read(file_path, is_filepath_metadata_table, separator=None, column_names=False, comment_line=None):
		"""
			Reading comma or tab separated values in a file as table

			@param file_path: path to file to be opened
			@type file_path: str | unicode
			@param separator: default character assumed to separate values in a file
			@type separator: str | unicode
			@param column_names: True if column names available
			@type column_names: bool
			@param comment_line: character or list of character indication comment lines
			@type comment_line: str | unicode | list[str|unicode]

			@return: None
			@rtype: None
		"""

		if comment_line is None:
			comment_line = ['#']
		elif isinstance(comment_line, str):
			comment_line = [comment_line]

		if separator is None:
			separator = "\t"

		assert isinstance(file_path, str)
		assert isinstance(separator, str)
		assert isinstance(comment_line, list)
		assert isinstance(column_names, bool)
		
		meta_table = {}

		with open_compress(file_path) as file_handler:

			# read column names
			if column_names:
				parse_column_names (file_handler, separator)
				for column_name in list_of_column_names_from_file:
					meta_table[column_name] = []

			# read rows
			row_count = 0
			for line in file_handler:
				row_count += 1
				row = line.rstrip('\n').rstrip('\r')
				if line[0] in comment_line or len(row) == 0:
					continue
				number_of_rows += 1
				row_cells = row.split(separator)
				number_of_columns = len(list(list_of_column_names_from_file))
				if number_of_columns != 0 and number_of_columns != len(row_cells):
					msg = "Format error. Bad number of values in line {}".format(row_count)
					raise ValueError(msg)
				for index, value in enumerate(row_cells):
					if column_names:
						column_name = list_of_column_names_from_file[index]
					else:
						column_name = index
						if column_name not in meta_table:
							meta_table[column_name] = []
					meta_table[column_name].append(row_cells[index].rstrip('\n').rstrip('\r'))
				if number_of_columns == 0:
					list_of_column_names_from_file = sorted(meta_table.keys())

		return meta_table

def has_column(column_name, meta_table):
		"""
			Get index of value in a column

			@attention:

			@param column_name: column name
			@type column_name: int | long | str | unicode

			@return: True if column available
			@rtype: bool
		"""
		assert isinstance(column_name, (str, int))

		if column_name in meta_table:
			return True
		else:
			return False

def validate_column_names(list_of_column_names, meta_table):
		"""
			Validate that a list of column names exists in the loaded table

			@attention:

			@param list_of_column_names: column name
			@type list_of_column_names: list[str|unicode]

			@return: True if all column names exist
			@rtype: bool
		"""
		assert isinstance(list_of_column_names, list)

		list_of_invalid_column_names = []
		for column_name in list_of_column_names:
			if not has_column(column_name, meta_table):
				list_of_invalid_column_names.append(column_name)
		if len(list_of_invalid_column_names) > 0:
			return False
		return True

def per_category():
		"""
		Get maximum number of strains valid to be drawn from a category.

		@rtype: int
		"""
		return per_cat + (rest > 0)

def get_all_strains(novelty_category):
		"""
		Return list of all strains from this category

		@return: List of drawn strain id
		@rtype: list[str|unicode]
		"""
		out = []
		for otu_id in novelty_category['otu_list']:
			for strain_id in novelty_category['otu_list'][otu_id]:
				out.append(strain_id)
		return out

def recalc(number_of_drawn):
		"""
		Adjust values based on the number of strains drawn from the previous category

		@param number_of_drawn: Number of drawn strains from the previous category
		@rtype: None
		"""
		assert isinstance(number_of_drawn, int)
		draw -= number_of_drawn
		cats -= 1
		per_cat = int(draw/cats)
		rest = draw % cats

def draw_strains(novelty_category, total, limit_per_otu):
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
		assert total <= novelty_category['number_of_strains']
		assert limit_per_otu > 0
		drawn_strain = []
		overhead = []
		drawn_strain_count_overall = 0
		otu_list = novelty_category['otu_list']
		for otu_id in random.sample(otu_list.keys(), len(otu_list)):
			drawn_strain_count_otu = 0
			for strain_id in random.sample(otu_list[otu_id], len(otu_list[otu_id])):
				if drawn_strain_count_otu < limit_per_otu and drawn_strain_count_overall < total:
					drawn_strain.append(strain_id)
					drawn_strain_count_otu += 1
					drawn_strain_count_overall += 1
			overhead += list(set(otu_list[otu_id])-set(drawn_strain))

		if drawn_strain_count_overall < total:
			# out += overhead[0:total-drawn_strain_cound_overall]
			#print "\n", total-drawn_strain_count_overall, len(overhead), len(drawn_strain), "\n"
			drawn_strain += random.sample(overhead, total-drawn_strain_count_overall)
		return drawn_strain

def substract():
		"""
		Adjust values by subtracting the maximum number of strains that can be drawn from a category.

		@rtype: None
		"""
		draw -= per_cat
		if rest > 0:
			rest -= 1

def draw_strains(categories, number_of_strains_to_draw, max_amount_of_strains_per_otu):
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

		cats = number_of_novelty_categories
		draw = number_of_strains_to_draw
		per_cat = int(number_of_strains_to_draw/number_of_novelty_categories)
		rest = number_of_strains_to_draw % number_of_novelty_categories
		per_otu = max_amount_of_strains_per_otu

		drawn_genome_id = []
		for x in sorted(categories, key=lambda name: categories[name].get_strain_amount()):
			if categories[x]['number_of_strains'] < per_category():
				# draw all strains
				drawn_genome_id += get_all_strains(categories[x])
				recalc(categories[x]['number_of_strains'])
			else:
				# draw subset of strains
				drawn_genome_id += draw_strains(categories[x], per_category(), per_otu)
				substract()

		if len(drawn_genome_id) < number_of_strains_to_draw:
			msg = "Could only draw {} samples instead of {}".format(len(drawn_genome_id), number_of_strains_to_draw)
			raise ValueError(msg)
		return drawn_genome_id

def get_drawn_genome_id(metadata_table, number_of_strains, number_of_strains_per_otu):
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
		assert isinstance(number_of_strains, int)
		assert isinstance(number_of_strains_per_otu, int)
		genomes_read_from_metadata = 0
		lost_lines = 0

		categories = dict()

		required_headers = [column_name_genome_id, column_name_otu, column_name_novelty_category]
		if not validate_column_names(required_headers, metadata_table):
			msg = 'A header is missing, needed: {}'.format(', '.join(required_headers))
			raise ValueError(msg)

		column_genome_id = metadata_table[column_name_genome_id]
		column_otu = metadata_table[column_name_otu]
		column_novelty = metadata_table[column_name_novelty_category]
		row_count = number_of_rows()

		# read metadata from provided file
		for row_index in range(0, row_count):
			if column_genome_id[row_index] == '' or column_otu[row_index] == '' or column_novelty[row_index] == '':
				lost_lines += 1
				continue
			genomes_read_from_metadata += 1
			novelty = column_novelty[row_index]
			if novelty not in categories:
				categories[novelty] = {'name': novelty, 
				 					   'number_of_strains': 0,
				 					   '_otu_list': {}}
			categories[novelty].add_strain(column_otu[row_index], column_genome_id[row_index])

		if lost_lines > 0:
			msg = "Invalid lines/ unusable genomes: {}".format(lost_lines)

		if genomes_read_from_metadata < number_of_strains:
			msg = 'Not enough data to draw.'
			raise ValueError(msg)

		# sample OTUs from the pool of available OTUs
		drawn_genome_id = draw_strains(categories, number_of_strains, number_of_strains_per_otu)
		return drawn_genome_id

def write(meta_table, file_path, separator=None, column_names=False, compression_level=0,
		exclude=None, value_list=None, key_column_name=None):
		"""
			Write tab separated files

			@attention: No comments will be written

			@param file_path: path to file to be opened
			@type file_path: str | unicode
			@param separator: file handler or file path to a log file
			@type separator: str | unicode
			@param column_names: True if column names should be written
			@type column_names: bool
			@param compression_level: any value above 0 will compress files
			@type compression_level: int | long
			@param exclude: If True, rows with a value in the value_list at the key_column_names are removed, False: all others are removed
			@type exclude: None | bool
			@param value_list:
			@type value_list: list[str|unicode]
			@param key_column_name: column name of excluded or included rows
			@type key_column_name: str | unicode

			@return: None
			@rtype: None
		"""

		if separator is None:
			separator = "\t"

		assert isinstance(file_path, str)
		assert isinstance(separator, str)
		assert isinstance(column_names, bool)
		assert isinstance(compression_level, int)
		assert 0 <= compression_level < 10
		assert exclude is None or isinstance(exclude, bool)
		assert value_list is None or isinstance(value_list, list)
		assert key_column_name is None or isinstance(key_column_name, str), "Invalid: {}".format(key_column_name)

		if compression_level > 0:
			file_handler = open_compress(file_path, "w", compression_level)
		else:
			file_handler = open(file_path, "w")

		if column_names:
			if not isinstance(list_of_column_names_from_file[0], str):
				header = separator.join([str(index) for index in list_of_column_names_from_file])
			else:
				header = separator.join(list_of_column_names_from_file)
			file_handler.write(header + '\n')
		for row_number in range(0, number_of_rows):
			if exclude is not None:
				if not exclude and meta_table[key_column_name][row_number] not in value_list:
					continue
				if exclude and meta_table[key_column_name][row_number] in value_list:
					continue

			row = []
			for column_names in list_of_column_names_from_file:
				row.append(str(meta_table[column_names][row_number]))
			file_handler.write(separator.join(row) + '\n')
		file_handler.close()

def get_map(metadata_table, key_column_name, value_column_name, unique_key=True):
		"""
			Keep rows at key values of a column

			@attention:

			@param key_column_name: Column name
			@type key_column_name: str | unicode | int | long
			@param value_column_name: Column name
			@type value_column_name: str | unicode | int | long

			@return: map
			@rtype: dict[str|unicode, str|unicode]

			@raises: KeyError
		"""

		assert isinstance(key_column_name, (str, int))
		assert isinstance(value_column_name, (str, int))
		assert has_column(metadata_table, key_column_name), "Column '{}' not found!".format(key_column_name)
		assert has_column(metadata_table, value_column_name), "Column '{}' not found!".format(value_column_name)

		if key_column_name not in metadata_table:
			# ??? self._logger.error("Column name '{}' not available!".format(key_column_name))
			return None
		if value_column_name not in metadata_table:
			# ??? self._logger.error("Column name '{}' not available!".format(value_column_name))
			return None
		new_map = {}
		if len(metadata_table) < 2:
			return new_map
		row_keys = metadata_table[key_column_name]
		row_values = metadata_table[value_column_name]
		for index, key in enumerate(row_keys):
			if unique_key and key in new_map:
				msg = "Key column is not unique! Key: '{}'".format(key)
				# ??? self._logger.error(msg)
				raise KeyError(msg)
			new_map[key] = row_values[index]
		return new_map

def get_genome_id_to_path_map(file_path_of_file_mapping_genome_id_to_paths, list_of_drawn_genome_id):
		"""
		Get a dictionary mapping genome id to the path of their genome

		@param file_path_of_file_mapping_genome_id_to_paths: File path to file with format 'id \t path'
		@type file_path_of_file_mapping_genome_id_to_paths: str | unicode
		@param list_of_drawn_genome_id: List of genome identifiers
		@type list_of_drawn_genome_id: list[str|unicode]

		@return: genome ids mapped to their gnome file path
		@rtype: dict[str|unicode, str|unicode]
		"""
		genome_id_to_path_map = {}

		mdt = read(file_path_of_file_mapping_genome_id_to_paths)
		# ??? if mdt.get_number_of_rows() > 0:
		if mdt:
			genome_id_to_path_map = get_map(mdt, 0, 1, unique_key=True)
		msg = "'{}' is missing one or more genome id".format(os.path.basename(file_path_of_file_mapping_genome_id_to_paths))
		assert set(genome_id_to_path_map.keys()).issuperset(list_of_drawn_genome_id), msg
		return {genome_id: genome_id_to_path_map[genome_id] for genome_id in list_of_drawn_genome_id}

def get_column(metadata_table, column_name):
		"""
			Get a column

			@attention: use index if no name available

			@param column_name: column name
			@type column_name: int | long | str | unicode

			@return: Cell values of a column
			@rtype: list[str|unicode]
		"""
		assert isinstance(column_name, (str, int))
		assert has_column(metadata_table, column_name), "Column '{}' not found!".format(column_name)
		return list(metadata_table[column_name])


def reduce_rows_to_subset(metadata_table, list_of_values, key_column_name):
		"""
			Keep rows at key values of a column

			@attention:

			@param list_of_values: Cell values of table column
			@type list_of_values: list[str|unicode]
			@param key_column_name: Column name
			@type key_column_name: str | unicode

			@return: Nothing
			@rtype: None
		"""

		assert isinstance(key_column_name, (str, int))
		assert isinstance(list_of_values, list)
		assert has_column(key_column_name, metadata_table), "Column '{}' not found!".format(key_column_name)

		new_meta_table = {}
		for column_name in list_of_column_names_from_file:
			new_meta_table[column_name] = []
		column = get_column(metadata_table, key_column_name)
		for index, value in enumerate(column):
			if value not in list_of_values:
				continue
			for column_name in list_of_column_names_from_file:
				new_meta_table[column_name].append(metadata_table[column_name][index])
		metadata_table = new_meta_table
		metadata_table = len(metadata_table[key_column_name])

		return metadata_table

def insert_column(meta_table, list_of_column_names, list_of_values=None, column_name=None):
		"""
			Insert a new column or overwrite an old one.

			@attention: if column_name exists, it will be overwritten

			@param list_of_values: Cell values of table column
			@type list_of_values: list[str|unicode]
			@param column_name: column name or index
			@type column_name: int | long | str | unicode

			@return: Nothing
			@rtype: None
		"""
		if column_name is None:
			column_name = len(ist_of_column_names)
		assert isinstance(column_name, (str, int))
		# assert len(values) == self._number_of_rows, ""

		if list_of_values is None:
			list_of_values = self.get_empty_column()
		assert isinstance(list_of_values, list)
		assert len(list_of_values) == self._number_of_rows, "Bad amount of values: {}/{}".format(
			len(list_of_values), self._number_of_rows)

		if column_name not in self._list_of_column_names:
			self._list_of_column_names.append(column_name)
		self._meta_table[column_name] = list_of_values

def concatenate(metadata_table, meta_table, strict=True):
		"""
			Concatenate two metadata tables

			@attention:

			@param meta_table: Table to be concatenated. The community metadata table.
			@type meta_table: MetadataTable
			@param strict: if true, both tables must have the same column names, else empty cells will be added where needed
			@type strict: bool

			@return: Nothing
			@rtype: None
		"""
		# metadata_table  self  list_of_column_names_new_metadata_table
		assert isinstance(strict, bool)

		if len(list_of_column_names_new_metadata_table) == 0:
			strict = False
		if strict:
			valid_foreign_column_names = validate_column_names(list_of_column_names_from_file, metadata_table)
			valid_own_column_names = validate_column_names(list_of_column_names_new_metadata_table, meta_table)
			if not valid_foreign_column_names or not valid_own_column_names:
				msg = "Column names are not identical!"
				#self._logger.error(msg)
				raise ValueError(msg)
			for column_name in list_of_column_names_new_metadata_table:
				metadata_table[column_name].extend(meta_table.get_column(column_name))
		else:
			for column_name in list_of_column_names_from_file:
				if column_name not in list_of_column_names_new_metadata_table:
					self.insert_column(self.get_empty_column(), column_name)
				self._meta_table[column_name].extend(meta_table.get_column(column_name))

		self._number_of_rows += meta_table.get_number_of_rows()

		for column_name in self._list_of_column_names:
			if len(self._meta_table[column_name]) < self._number_of_rows:
				self._meta_table[column_name].extend([''] * (self._number_of_rows - len(self._meta_table[column_name])))

def design_community(community, number_of_samples, directory_in_template=None):
		"""
		Design artificial community, of a specific design, with different distributions for each sample

		@param file_path_distributions: File path where distributions will be written to
		@type file_path_distributions: str | unicode
		@param community: Input data for the creation of a community
		@type community: Community
		@param number_of_samples: Amount of samples to be simulated
		@type number_of_samples: int
		@param metadata_table: Will contain metadata of all (simulated) genomes/plasmids drawn
		@type metadata_table: MetadataTable
		@param directory_out_metadata: Metadata tables of separated by chosen and not chosen genomes are written to here
		@type directory_out_metadata: str | unicode
		@param directory_in_template: contains template data for strain simulation
		@type directory_in_template: str | unicode

		@return: Dictionary with drawn genome ids as key and file paths as value
		@rtype: dict[str|unicode, str|unicode]
		"""

		number_of_strains = community['genomes_total']
		genomes_real = community['genomes_real']
		silent = not community['verbose']

		# pick how much a strain will be simulated
		genome_amounts = []
		strain_simulation = None
		if community['simulate_strains']:

			probability = None  # 1-options.communities[community_id]["evolve"]
			genome_amounts = get_genome_amounts(probability ,number_of_strains, genomes_real, silent)
			number_of_strains = len(genome_amounts)

		# draw strains
		meta_table_community = read(community['file_path_metadata_table'], column_names=True)
		
		list_of_drawn_genome_id = get_drawn_genome_id(
			metadata_table=meta_table_community,
			number_of_strains=number_of_strains,
			number_of_strains_per_otu=community.limit_per_otu
			)

		# write unused data to separate file
		old_base_name = os.path.basename(community['file_path_metadata_table'])
		file_prefix, extention = os.path.splitext(old_base_name)
		new_file_name = "unused_c{index}_{prefix}{ext}".format(
			prefix=file_prefix,
			index=community.id,
			ext=extention)
		metadata_new_file_path = os.path.join("./", new_file_name)
		write(meta_table_community,
			metadata_new_file_path,
			exclude=True,
			value_list=list_of_drawn_genome_id,
			key_column_name=column_name_genome_id,
			column_names=True)

		# get path for every genome
		genome_id_to_file_path_gff = None
		if community["file_path_gff_locations"]:
			genome_id_to_file_path_gff = get_genome_id_to_path_map(community["file_path_gff_locations"], list_of_drawn_genome_id)
		genome_id_to_path_map = get_genome_id_to_path_map(community.file_path_genome_locations, list_of_drawn_genome_id)

		# concatenate
		meta_table_community = reduce_rows_to_subset(meta_table_community, list_of_drawn_genome_id, column_name_genome_id)

		# ???
		# in case of de novo community design, there is an empty metadata table at the beginning
		metadata_table = {}
		list_of_column_names_new_metadata_table = []

		concatenate(metadata_table, meta_table_community, strict=False)

		# validate correct format of files
		self._logger.info("Validating raw sequence files!")
		assert self.validate_format(
			list_of_file_paths=genome_id_to_path_map.values(),
			file_format="fasta",
			sequence_type="dna",
			ambiguous=True
			), "Validation of file format failed!"

		# simulate diversity around strains
		if community.simulate_strains:
			genome_id_to_amounts = strain_simulation.get_genome_id_to_amounts(list_of_drawn_genome_id, genome_amounts)
			strain_simulation.simulate_strains(
				meta_table=metadata_table,
				genome_id_to_amounts=genome_id_to_amounts,
				genome_id_to_file_path_genome=genome_id_to_path_map,
				genome_id_to_file_path_gff=genome_id_to_file_path_gff)
			# adopt new list that includes simulated strains
			self._logger.info("Validating simulated sequence files!")
			for genome_id, file_path in genome_id_to_path_map.items():
				if genome_id in list_of_drawn_genome_id:
					continue
				assert self.validate_sequence_file(
					file_path,
					file_format="fasta",
					sequence_type="dna",
					ambiguous=True)
			list_of_drawn_genome_id = genome_id_to_path_map.keys()

		# get community distributions
		population_distribution = PopulationDistribution(
			logfile=self._logfile, verbose=self._verbose, debug=self._debug)
		list_of_distributions = population_distribution.get_lists_of_distributions(
			size_of_population=len(list_of_drawn_genome_id),
			number_of_samples=number_of_samples,
			modus=community.mode,
			log_mu=community.log_mu, log_sigma=community.log_sigma,
			gauss_mu=community.gauss_mu, gauss_sigma=community.gauss_sigma,
			view_distribution=community.verbose
		)

		# move and clean up files (removes sequence description)
		# genome_id_to_total_length = self.move_genome_files(
		#     genome_id_to_path_map,
		#     directory_output=directory_out_genomes,
		#     sequence_min_length=min_sequence_length,
		#     set_of_sequence_names=set_of_sequence_names)

		# write distribution file
		# genome_id_to_distributions = self._get_genome_id_to_distributions(list_of_drawn_genome_id, list_of_distributions)
		assert len(list_of_drawn_genome_id) == len(list_of_distributions)
		genome_id_to_distributions = dict(zip(list_of_drawn_genome_id, list_of_distributions))

		# genome_id_to_file_name = self._get_genome_id_to_file_name(genome_id_to_path_map)
		with open(file_path_distributions, 'w') as stream_out:
			self._write_distribution_file(stream_out=stream_out, genome_id_to_abundance=genome_id_to_distributions)
		return genome_id_to_path_map


def design_samples(list_of_communities, list_of_file_paths_distribution, number_of_samples, directory_in_template=None):
		"""
		Design artificial community, of a specific design, with different distributions for each sample

		@param list_of_file_paths_distribution: Output file path list for distribution files
		@type list_of_file_paths_distribution: list[str | unicode]
		@param list_of_communities: List of input data for the creation of a community
		@type list_of_communities: list[Community]
		@param metadata_table: Will contain metadata of all (simulated) genomes/plasmids drawn
		@type metadata_table: MetadataTable
		@param directory_out_metadata: Metadata tables of separated by chosen and not chosen genomes are written to here
		@type directory_out_metadata: str | unicode
		@param directory_in_template: contains template data for strain simulation
		@type directory_in_template: str | unicode

		@return: Dictionary with drawn genome ids as key and file paths as value
		@rtype: tuple[dict[str|unicode, str|unicode]
		"""
		assert isinstance(list_of_communities, list)
		assert isinstance(list_of_file_paths_distribution, list)

		list_of_comunity_distribution_file_paths = []

		merged_genome_id_to_path_map = {}
		for community in list_of_communities:
			list_of_comunity_distribution_file_paths.append(file_path_output_comunity)
			genome_id_to_path_map = design_community(community, number_of_samples, directory_in_template)
			merged_genome_id_to_path_map.update(genome_id_to_path_map)

		for index_sample, file_path_output in enumerate(list_of_file_paths_distribution):
			self.merge_communities(
				list_of_communities, list_of_comunity_distribution_file_paths, index_sample, file_path_output)

		# delete now obsolete files
		if not self._debug:
			for file_path in list_of_comunity_distribution_file_paths:
				if os.path.isfile(file_path):
					os.remove(file_path)

		return merged_genome_id_to_path_map


def start_design_community(max_processors, number_of_samples, directory_simulation_template, list_of_communities, seed=None):
		"""
		Start designing sample a community

		@return: map genome id to genome file path and list of distribution file paths
		@rtype: tuple[dict[str|unicode, str|unicode], list[str|unicode]]]
		"""

		if seed is not None:
			random.seed(seed)
			np_random.seed(abs(hash(seed)) % 4294967295)  # numpy accepts only 32 bit integers

		list_of_file_paths_distribution = get_distribution_file_paths("./", number_of_samples)
		merged_genome_id_to_path_map = design_samples(list_of_communities, list_of_file_paths_distribution, number_of_samples, directory_simulation_template)
		#     directory_out_distributions=directory_out_distributions,
		self.write_profile_gold_standard(meta_data_table, list_of_file_paths_distribution)

		file_path_metadata = self._project_file_folder_handler.get_genome_metadata_file_path()
		meta_data_table.write(file_path_metadata, column_names=True)
		return merged_genome_id_to_path_map, list_of_file_paths_distribution

# TODO mehrere communities
# main method and entry point of this script
if __name__ == "__main__":

	#TODO  _num_real_genomes? siehe configfilehandler zeile 168 -> ist das optional? wird das ueberhaupt genutzt?
	#TODO id_to_gff_file? siehe configfilehandler zeile 148 -> ist das optional?
	max_processors = sys.argv[1]
	number_of_samples = sys.argv[2]
	strain_simulation_template = sys.argv[3]
	genomes_total = sys.argv[4]
	max_strains_per_otu = sys.argv[5]
	file_path_metadata_table = sys.argv[6]
	file_path_genome_locations = sys.argv[7]
	file_path_gff_locations = sys.argv[8]
	ratio = sys.argv[9]
	mode = sys.argv[10]
	log_mu = sys.argv[11]
	log_sigma = sys.argv[12]
	gauss_mu = sys.argv[13]
	gauss_sigma = sys.argv[14]
	verbose = sys.argv[15]

	num_real_genomes = None

	simulate_strains = False
	if num_real_genomes and num_real_genomes < genomes_total:
		simulate_strains = True

	list_of_communities = [{'genomes_total': genomes_total, 
				 'max_strains_per_otu': max_strains_per_otu,   
				 'file_path_metadata_table': file_path_metadata_table,
				 'file_path_gff_locations': file_path_gff_locations,
				 'ratio': ratio,
				 'mode': mode,
				 'log_mu': log_mu,
				 'log_sigma': log_sigma,
				 'gauss_mu': gauss_mu,
				 'gauss_sigma': gauss_sigma,
				 'simulate_strains' : simulate_strains,
				 'verbose': verbose}]

	start_design_community(max_processors, number_of_samples, strain_simulation_template, list_of_communities, seed=None)