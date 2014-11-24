#!/usr/bin/env python

__author__ = 'hofmann'

import os


class MetaTable:
	"""Reading and writing a meta table"""
	def __init__(self, separator="\t", logger=None):
		self._logger = logger
		self._number_of_rows = 0
		self._header = []
		self._meta_table = {}
		self._separator = separator

	def clear(self):
		self._header = []
		self._number_of_rows = 0
		self._meta_table = {}

	def remove_empty_columns(self):
		new_header = list(self._header)
		for title in self._header:
			column = set(self.get_column(title))
			column = [value.strip() for value in column]
			if len(column) == 1 and '' in column:
				new_header.remove(title)
		self._header = new_header

	def read(self, file_path, head=True):
		self.clear()
		if not os.path.isfile(file_path):
			if self._logger:
				self._logger.error("MetaTable: no file found at: '{}'".format(file_path))
			return
		with open(file_path) as file_handler:
			if head:
				self._header = file_handler.readline().strip().split(self._separator)
				for column_name in self._header:
					self._meta_table[column_name] = []
			for line in file_handler:
				if line.startswith('#') or len(line.strip().strip('\r')) == 0:
					continue
				self._number_of_rows += 1
				row = line.split(self._separator)
				for index in range(0, len(row)):
					if head:
						self._meta_table[self._header[index]].append(row[index].strip())
					else:
						if index not in self._meta_table:
							self._meta_table[index] = []
						self._meta_table[index].append(row[index].strip())

	def write(self, file_path, head=True):
		with open(file_path, "w") as file_handler:
			if head:
				header = self._separator.join(self._header)
				file_handler.write(header + '\n')
			for row_number in range(0, self._number_of_rows):
				row = []
				for head in self._header:
					if len(self._meta_table[head]) > row_number:
						row.append(str(self._meta_table[head][row_number]))
					else:
						row.append('')
				file_handler.write(self._separator.join(row) + '\n')

	def get_header(self):
		return self._header

	def get_number_of_rows(self):
		return self._number_of_rows

	def get_entry_index(self, column_name, entry_name):
		if entry_name in self._meta_table[column_name]:
			return self._meta_table[column_name].index(entry_name)
		else:
			return None

	def get_column(self, column_name):
		if column_name in self._meta_table:
			return self._meta_table[column_name]
		else:
			return None

	def get_empty_column(self, default_value=''):
		return [default_value] * self._number_of_rows

	def set_column(self, values=None, column_name=None):
		if column_name is None:
			column_name = len(self._header)
		if values is None:
			values = []
		#if self._number_of_rows == 0:
		#	self._number_of_rows = len(values)
		#elif self._number_of_rows < len(values):
		#	self._number_of_rows = len(values)
		#	if self._logger:
		#		self._logger.warning("MetaTable: Inconsistent length of columns! {} - {}".format(self._number_of_rows, len(values)))
		if column_name not in self._header:
			self._header.append(column_name)
		self._meta_table[column_name] = values

	def get_new_row(self):
		row = {}
		for head in self._header:
			row[head] = ''
		return row

	def add_row(self, row):
		diff = set(self._header).difference(set(row.keys()))
		if len(diff) != 0:
			if self._logger:
				self._logger.error("MetaTable: Bad header, could not add row!")
			return
		self._number_of_rows += 1
		for head in self._header:
			self._meta_table[head].append(row[head])

	def get_cell_value(self, key_header, key_value, value_header):
		if key_header not in self._header or value_header not in self._header:
			return None
		index = self.get_entry_index(key_header, key_value)
		if index is not None:
			return self._meta_table[value_header][index]
		return None

	def are_valid_header(self, list_of_header):
		for header in list_of_header:
			if header not in self._header:
				return False
		return True

	def concatenate(self, meta_table, strict=True):
		if len(self._header) == 0:
			strict = False
		if strict:
			if not self.are_valid_header(meta_table.get_header()) or not meta_table.are_valid_header(self._header):
				if self._logger:
					self._logger.error("MetaTable: header are not identical!")
				return
			for header in self._header:
				self._meta_table[header].extend(meta_table.get_column(header))
		else:
			for header in meta_table.get_header():
				if header in self._header:
					self._meta_table[header].extend(meta_table.get_column(header))
				else:
					new_column = self.get_empty_column()
					new_column.extend(meta_table.get_column(header))
					self.set_column(new_column, header)
		self._number_of_rows += meta_table.get_number_of_rows()
