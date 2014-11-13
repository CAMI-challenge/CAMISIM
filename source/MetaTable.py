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
		self._meta_table = {}
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
				if line.startswith('#'):
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
				header = ""
				for value in self._header:
					header += value + self._separator
				file_handler.write(header.rstrip('\t') + '\n')
			for row_number in range(0, self._number_of_rows):
				row = []
				for value in self._header:
					if len(self._meta_table[value]) > row_number:
						row.append(str(self._meta_table[value][row_number]))
					else:
						row.append('')
				file_handler.write(self._separator.join(row) + '\n')

	def get_number_of_rows(self):
		return self._number_of_rows

	def get_entry_index(self, column_name, entry_name):
		if entry_name in self._meta_table[column_name]:
			return self._meta_table[column_name].index(entry_name)
		else:
			return -1

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
		if self._number_of_rows == 0:
			self._number_of_rows = len(values)
		elif self._number_of_rows < len(values):
			self._number_of_rows = len(values)
			if self._logger:
				self._logger.warning("MetaTable: Inconsistent length of columns!")
		if column_name not in self._header:
			self._header.append(column_name)
		self._meta_table[column_name] = values
