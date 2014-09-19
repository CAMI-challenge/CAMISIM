#!/usr/bin/env python

__author__ = 'hofmann'

import sys
import os
import argparse
import re


class MetaTable:
	"""Reading and writing a meta table"""
	def __init__(self, file_path=None, separator="\t", head=True):
		self.number_of_rows = 0
		self.header = []
		self.clear()
		self.file_path = file_path
		self.separator = separator
		if file_path is not None and os.path.isfile(file_path):
			self.meta_table = self.read_meta_table_by_column(head)
		elif file_path is not None:
			print "Error (MetaTable): no file found at:", file_path

	def clear(self):
		self.file_path = None
		self.header = []
		self.number_of_rows = 0
		self.meta_table = {}

	def read_meta_table_by_column(self, head=True):
		meta_table = {}
		with open(self.file_path) as file_handler:
			if head:
				self.header = file_handler.readline().strip().split(self.separator)
				#print "headers", len(self.header), self.header
				for column_name in self.header:
					meta_table[column_name] = []
			for line in file_handler:
				self.number_of_rows += 1
				if line[0] != '#':
					row = line.split(self.separator)
					for index in range(0, len(row)):
						if head:
							meta_table[self.header[index]].append(row[index].strip())
						else:
							if index not in meta_table:
								meta_table[index] = []
							meta_table[index].append(row[index].strip())
		return meta_table

	def write_meta_table_by_column(self, file_path):
		with open(file_path, "w") as file_handler:
			#number_of_rows = len(self.meta_table[self.header[0]])
			header = ""
			for value in self.header:
				header += value + self.separator
			file_handler.write(header.strip() + "\n")
			for row_number in range(0, self.number_of_rows):
				row = []
				for value in self.header:
					if len(self.meta_table[value]) > row_number:
						row.append(str(self.meta_table[value][row_number]))
					else:
						row.append('')
				file_handler.write(self.separator.join(row) + "\n")

	def get_entry_index(self, column_name, entry_name):
		if entry_name in self.meta_table[column_name]:
			return self.meta_table[column_name].index(entry_name)
		else:
			return -1

	def get_column(self, column_name):
		if column_name in self.meta_table:
			return self.meta_table[column_name]
		else:
			return None

	def set_column(self, column_name, values=[]):
		if column_name not in self.header:
			self.header.append(column_name)
		if self.number_of_rows == 0:
			self.number_of_rows = len(values)
		self.meta_table[column_name] = values
