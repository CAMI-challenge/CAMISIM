#!/usr/bin/env python

__author__ = 'hofmann'
__version__ = '0.1.2'

import io
from scripts.Archive.compress import Compress


class MetadataTable(Compress):
	"""Reading and writing a metadata table"""

	def __init__(self, separator="\t", logfile=None, verbose=True):
		"""
			Handle tab separated files

			@attention:

			@param separator: default character assumed to separate values in a file
			@type separator: str | unicode
			@param logfile: file handler or file path to a log file
			@type logfile: file | io.FileIO | StringIO.StringIO | str
			@param verbose: Not verbose means that only warnings and errors will be past to stream
			@type verbose: bool

			@return: None
			@rtype: None
		"""
		assert logfile is None or isinstance(logfile, str) or self.is_stream(logfile)
		assert isinstance(separator, str), "separator must be string"
		assert isinstance(verbose, bool), "verbose must be true or false"
		super(MetadataTable, self).__init__(label="MetadataReader", logfile=logfile, verbose=verbose)

		self._number_of_rows = 0
		self._meta_table = {}
		self._separator = separator
		self._list_of_column_names = []

	def clear(self):
		self._number_of_rows = 0
		self._meta_table = {}
		self._list_of_column_names = []

	def _has_unique_columns(self, list_of_column_names=None):
		if list_of_column_names is None:
			list_of_column_names = self._list_of_column_names
		return len(list_of_column_names) == len(set(list_of_column_names))

	def remove_empty_columns(self):
		for column_name in self.get_column_names():
			column = set(self.get_column(column_name))
			column = [value.strip() for value in column]
			if len(column) == 1 and '' in column:
				self._meta_table.pop(column_name)
				index = self._list_of_column_names.index(column_name)
				self._list_of_column_names.pop(index)

	def _parse_column_names(self, stream_input, separator):
			row = stream_input.readline().rstrip('\n').rstrip('\r')
			list_of_column_names = row.split(separator)
			assert self._has_unique_columns(list_of_column_names), "Column names must be unique!"
			return list_of_column_names

	def parse_file(self, file_path, separator=None, column_names=False, comment_line=None, as_list=True):
		"""
			Reading comma or tab separated values from a file

			@param file_path: path to file to be opened
			@type file_path: str | unicode
			@param separator: default character assumed to separate values in a file
			@type separator: str | unicode
			@param column_names: True if column names available
			@type column_names: bool
			@param comment_line: character or list of character indication comment lines
			@type comment_line: str | unicode | list[str|unicode]
			@param as_list: If true lists are returned, else dicts.
			@param as_list: bool

			@return: Generator of dictionary representing rows
			@rtype: generator[ dict[int|long|str|unicode, str|unicode] ]
			#
		"""
		with self.open(file_path) as file_handler:
			for row in self.parse_stream(file_handler, separator, column_names, comment_line, as_list):
				yield row

	def parse_stream(self, stream_input, separator=None, column_names=False, comment_line=None, as_list=True):
		"""
			Reading comma or tab separated values from a stream

			@param stream_input: stream
			@type stream_input: file | io.FileIO | StringIO.StringIO
			@param separator: default character assumed to separate values in a file
			@type separator: str | unicode
			@param column_names: True if column names available
			@type column_names: bool
			@param comment_line: character or list of character indication comment lines
			@type comment_line: str | unicode | list[str|unicode]
			@param as_list: If true lists are returned, else dicts.
			@param as_list: bool

			@return: Generator returning dictionary or list
			@rtype: generator[ dict[int|long|str|unicode, str|unicode] ]
			#
		"""
		if comment_line is None:
			comment_line = ['#']
		elif isinstance(comment_line, str):
			comment_line = [comment_line]

		if separator is None:
			separator = self._separator

		assert self.is_stream(stream_input)
		assert isinstance(separator, str)
		assert isinstance(comment_line, list)
		assert isinstance(column_names, bool)
		self.clear()

		# read column names
		number_of_columns = 0
		list_of_column_names = []
		if column_names:
			list_of_column_names = self._parse_column_names(stream_input, separator)
			number_of_columns = len(list_of_column_names)

		# read rows
		line_count = 0
		for line in stream_input:
			line_count += 1
			row = line.rstrip('\n').rstrip('\r')
			if line[0] in comment_line or len(row) == 0:
				continue

			row_cells = row.split(separator)
			if not column_names and number_of_columns == 0:
				number_of_columns = len(row_cells)

			if number_of_columns != len(row_cells):
				msg = "Format error. Bad number of values in line {}".format(line_count)
				self._logger.error(msg)
				raise ValueError(msg)

			if as_list:
				yield row_cells
				continue

			dict_row = dict()
			for index, value in enumerate(row_cells):
				if column_names:
					column_name = list_of_column_names[index]
				else:
					column_name = index
				dict_row[column_name] = row_cells[index].rstrip('\n').rstrip('\r')
			if number_of_columns == 0:
				list_of_column_names = sorted(dict_row.keys())
			yield dict_row

	def read(self, file_path, separator=None, column_names=False, comment_line=None):
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
			separator = self._separator

		assert isinstance(file_path, str)
		assert self.validate_file(file_path)
		assert isinstance(separator, str)
		assert isinstance(comment_line, list)
		assert isinstance(column_names, bool)

		self.clear()
		with self.open(file_path) as file_handler:
			self._logger.info("Reading file: '{}'".format(file_path))

			# read column names
			if column_names:
				self._list_of_column_names = self._parse_column_names(file_handler, separator)
				for column_name in self._list_of_column_names:
					self._meta_table[column_name] = []

			# read rows
			row_count = 0
			for line in file_handler:
				row_count += 1
				row = line.rstrip('\n').rstrip('\r')
				if line[0] in comment_line or len(row) == 0:
					continue
				self._number_of_rows += 1
				row_cells = row.split(separator)
				number_of_columns = len(self.get_column_names())
				if number_of_columns != 0 and number_of_columns != len(row_cells):
					msg = "Format error. Bad number of values in line {}".format(row_count)
					self._logger.error(msg)
					raise ValueError(msg)
				for index, value in enumerate(row_cells):
					if column_names:
						column_name = self._list_of_column_names[index]
					else:
						column_name = index
						if column_name not in self._meta_table:
							self._meta_table[column_name] = []
					self._meta_table[column_name].append(row_cells[index].rstrip('\n').rstrip('\r'))
				if number_of_columns == 0:
					self._list_of_column_names = sorted(self._meta_table.keys())

	def write(
		self, file_path, separator=None, column_names=False, compression_level=0,
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
			separator = self._separator

		assert isinstance(file_path, str)
		assert self.validate_dir(file_path, only_parent=True)
		assert isinstance(separator, str)
		assert isinstance(column_names, bool)
		assert isinstance(compression_level, int)
		assert 0 <= compression_level < 10
		assert exclude is None or isinstance(exclude, bool)
		assert value_list is None or isinstance(value_list, list)
		assert key_column_name is None or isinstance(key_column_name, str), "Invalid: {}".format(key_column_name)

		if compression_level > 0:
			file_handler = self.open(file_path, "w", compression_level)
		else:
			file_handler = open(file_path, "w")

		if column_names:
			if not isinstance(self._list_of_column_names[0], str):
				header = separator.join([str(index) for index in self._list_of_column_names])
			else:
				header = separator.join(self._list_of_column_names)
			file_handler.write(header + '\n')
		for row_number in range(0, self._number_of_rows):
			if exclude is not None:
				if not exclude and self._meta_table[key_column_name][row_number] not in value_list:
					continue
				if exclude and self._meta_table[key_column_name][row_number] in value_list:
					continue

			row = []
			for column_names in self._list_of_column_names:
				row.append(str(self._meta_table[column_names][row_number]))
			file_handler.write(separator.join(row) + '\n')
		file_handler.close()

	def get_column_names(self):
		"""
			Get list of column names

			@attention: returns list of indexes if no column names available

			@return: List of column names or indexes
			@rtype: list[str|int]
		"""
		return list(self._list_of_column_names)

	def get_number_of_rows(self):
		"""
			Get number of rows

			@attention:

			@return: Number of rows
			@rtype: int
		"""
		return self._number_of_rows

	def get_number_of_columns(self):
		"""
			Get number of columns

			@attention:

			@return: Number of rows
			@rtype: int
		"""
		return len(self.get_column_names())

	def get_row_index_of_value(self, value, column_name):
		"""
			Get index of value in a column

			@attention:

			@param value: value in column
			@type value: str | unicode
			@param column_name: column name
			@type column_name: int | long | str | unicode

			@return: index of value in a column, None if not there
			@rtype: None | int
		"""
		assert isinstance(column_name, (str, int))
		assert self.has_column(column_name), "Column '{}' not found!".format(column_name)

		if value in self._meta_table[column_name]:
			return self._meta_table[column_name].index(value)
		else:
			return None

	def has_column(self, column_name):
		"""
			Get index of value in a column

			@attention:

			@param column_name: column name
			@type column_name: int | long | str | unicode

			@return: True if column available
			@rtype: bool
		"""
		assert isinstance(column_name, (str, int))

		if column_name in self._meta_table:
			return True
		else:
			return False

	def get_column(self, column_name):
		"""
			Get a column

			@attention: use index if no name available

			@param column_name: column name
			@type column_name: int | long | str | unicode

			@return: Cell values of a column
			@rtype: list[str|unicode]
		"""
		assert isinstance(column_name, (str, int))
		assert self.has_column(column_name), "Column '{}' not found!".format(column_name)
		return list(self._meta_table[column_name])

	def get_empty_column(self, default_value=''):
		"""
			Get a empty column with the same number of rows as the current table

			@attention: empty list if number of rows is zero

			@param default_value: column name
			@type default_value: str | unicode

			@return: Column with cell values set to default value
			@rtype: list[str|unicode]
		"""
		assert isinstance(default_value, str)
		return [default_value] * self._number_of_rows

	def get_empty_row(self, default_value='', as_list=False):
		"""
			Get a empty column with the same number of rows as the current table

			@attention: empty list if number of rows is zero

			@param default_value: column name
			@type default_value: str | unicode
			@param as_list: return a list if true
			@type as_list: bool

			@return: Column with cell values set to default value
			@rtype: dict | list
		"""
		assert isinstance(default_value, str)
		assert isinstance(as_list, bool)
		if as_list:
			return [default_value] * len(self._list_of_column_names)
		row = {}
		for column_name in self._list_of_column_names:
			row[column_name] = default_value
		return row

	def insert_column(self, list_of_values=None, column_name=None):
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
			column_name = len(self._list_of_column_names)
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

	def insert_row(self, row):
		"""
			Insert a new row.

			@attention:

			@param row: Cell values of a row
			@type row: list[str|unicode] | dict

			@return: Nothing
			@rtype: None
		"""
		assert isinstance(row, (list, dict))
		assert len(row) == len(self._list_of_column_names)
		if isinstance(row, dict):
			diff = set(self._list_of_column_names).difference(set(row.keys()))
			if len(diff) != 0:
				msg = "Bad column names '{}', could not add row!".format(", ".join(diff))
				self._logger.error(msg)
				raise ValueError(msg)
			for column_name in self._list_of_column_names:
				self._meta_table[column_name].append(row[column_name])
		else:
			# assert len(row) == len(self._header)
			for index_column in range(len(row)):
				self._meta_table[self._list_of_column_names[index_column]].append(row[index_column])
		self._number_of_rows += 1

	def get_cell_value(self, key_column_name, key_value, value_column_name):
		"""
			Get the cell value at the index of a key in a key column

			@attention:

			@param key_column_name: column name
			@type key_column_name: str | unicode | int | long
			@param value_column_name: column name
			@type value_column_name: str | unicode | int | long
			@param key_value: key cell value
			@type key_value: str | unicode

			@return: None if key value is not there
			@rtype: str | unicode | None
		"""
		assert isinstance(key_column_name, (str, int))
		assert isinstance(value_column_name, (str, int))
		assert self.has_column(key_column_name), "Column '{}' not found!".format(key_column_name)
		assert self.has_column(value_column_name), "Column '{}' not found!".format(value_column_name)

		index = self.get_row_index_of_value(key_value, key_column_name)
		if index is not None:
			return self._meta_table[value_column_name][index]
		return None

	def validate_column_names(self, list_of_column_names):
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
			if not self.has_column(column_name):
				list_of_invalid_column_names.append(column_name)
				self._logger.info("Invalid columns: {}".format(", ".join(list_of_invalid_column_names)))
		if len(list_of_invalid_column_names) > 0:
			return False
		return True

	def concatenate(self, meta_table, strict=True):
		"""
			Concatenate two metadata tables

			@attention:

			@param meta_table: column name
			@type meta_table: MetadataTable
			@param strict: if true, both tables must have the same column names, else empty cells will be added where needed
			@type strict: bool

			@return: Nothing
			@rtype: None
		"""
		assert isinstance(meta_table, MetadataTable)
		assert isinstance(strict, bool)

		if len(self._list_of_column_names) == 0:
			strict = False
		if strict:
			valid_foreign_column_names = self.validate_column_names(meta_table.get_column_names())
			valid_own_column_names = meta_table.validate_column_names(self._list_of_column_names)
			if not valid_foreign_column_names or not valid_own_column_names:
				msg = "Column names are not identical!"
				self._logger.error(msg)
				raise ValueError(msg)
			for column_name in self._list_of_column_names:
				self._meta_table[column_name].extend(meta_table.get_column(column_name))
		else:
			for column_name in meta_table.get_column_names():
				if column_name not in self._list_of_column_names:
					self.insert_column(self.get_empty_column(), column_name)
				self._meta_table[column_name].extend(meta_table.get_column(column_name))

		self._number_of_rows += meta_table.get_number_of_rows()

		for column_name in self._list_of_column_names:
			if len(self._meta_table[column_name]) < self._number_of_rows:
				self._meta_table[column_name].extend([''] * (self._number_of_rows - len(self._meta_table[column_name])))

	def reduce_rows_to_subset(self, list_of_values, key_column_name):
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
		assert self.has_column(key_column_name), "Column '{}' not found!".format(key_column_name)

		new_meta_table = {}
		for column_name in self._list_of_column_names:
			new_meta_table[column_name] = []
		column = self.get_column(key_column_name)
		for index, value in enumerate(column):
			if value not in list_of_values:
				continue
			for column_name in self._list_of_column_names:
				new_meta_table[column_name].append(self._meta_table[column_name][index])
		self._meta_table = new_meta_table
		self._number_of_rows = len(self._meta_table[key_column_name])

	def get_map(self, key_column_name, value_column_name, unique_key=True):
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
		assert self.has_column(key_column_name), "Column '{}' not found!".format(key_column_name)
		assert self.has_column(value_column_name), "Column '{}' not found!".format(value_column_name)

		if key_column_name not in self._meta_table:
			self._logger.error("Column name '{}' not available!".format(key_column_name))
			return None
		if value_column_name not in self._meta_table:
			self._logger.error("Column name '{}' not available!".format(value_column_name))
			return None
		new_map = {}
		if len(self._meta_table) < 2:
			return new_map
		row_keys = self._meta_table[key_column_name]
		row_values = self._meta_table[value_column_name]
		for index, key in enumerate(row_keys):
			if unique_key and key in new_map:
				msg = "Key column is not unique! Key: '{}'".format(key)
				self._logger.error(msg)
				raise KeyError(msg)
			new_map[key] = row_values[index]
		return new_map

	def rename_column(self, old_column_name, new_column_name):
		"""
			Keep rows at key values of a column

			@attention:

			@param old_column_name: Column name
			@type old_column_name: str | unicode
			@param new_column_name: Column name
			@type new_column_name: str | unicode

			@return: Nothing
			@rtype: None
		"""
		assert isinstance(old_column_name, (str, int))
		assert isinstance(new_column_name, (str, int))
		assert self.has_column(old_column_name), "Column '{}' not found!".format(old_column_name)

		self._list_of_column_names[self._list_of_column_names.index(old_column_name)] = new_column_name
		self._meta_table[new_column_name] = self._meta_table.pop(old_column_name)
