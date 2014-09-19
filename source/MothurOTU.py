__author__ = 'hofmann'

import sys
import os
import argparse
import re
import operator
from collections import Counter


class MothurOTU:
	"""Reading and writing a meta table"""
	def __init__(self, file_path, taxonomy=None, otu_separator="\t", strain_separator=","):
		self.file_path = file_path
		self.otu_separator = otu_separator
		self.strain_separator = strain_separator
		self.taxonomy = taxonomy
		#self.header=[]
		self.otu_cutoff = "unique"
		self.rank = "strain"
		self.otu_lists_by_cutoff_raw = {}
		self.otu_list = []
		self.read_mothur_clustering_file()
		#self.ranks = ['root', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
		self.ranks = ['strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom', 'root']
		self.tmp_tax_list_by_rank = {}

	def get_otus(self, cutoffs, otu_ids, unknown_genomes_id_column):
		result = {}
		number_otu_ids = len(otu_ids)
		for row_index in range(0, number_otu_ids):
			otu = self.get_otu(cutoffs[row_index], int(float(otu_ids[row_index])))
			if otu is None:
				result[unknown_genomes_id_column[row_index]] = None
				continue
			result[unknown_genomes_id_column[row_index]] = []
			for item in otu:
				if item not in otu_ids:
					result[unknown_genomes_id_column[row_index]].append(item)
		return result

	def get_otu(self, cutoff, otu_id):
		if cutoff not in self.otu_lists_by_cutoff_raw:
			print "MothurOTU: bad cutoff {}".format(cutoff)
			return None
		#if otu_id not in self.otu_lists_by_cutoff_raw[cutoff]:
		#	print "MothurOTU: bad otu_id".format(otu_id)
		#	return None
		if otu_id >= len(self.otu_lists_by_cutoff_raw[cutoff]["otu"]):
			print "MothurOTU: bad otu_id".format(otu_id)
			return None
		return self.otu_lists_by_cutoff_raw[cutoff]["otu"][otu_id]

	def read_mothur_clustering_file(self):
		self.rank = "strain"
		print "loading mothur file"
		with open(self.file_path) as file_handler:
			for line in file_handler:
				if line[0] != '#' and "label" not in line:
					line = line.strip()
					otu_list = []
					row = line.split(self.otu_separator)
					#if row[0] == "unique":
					cutoff_of_otu = row[0]
					number_of_otu = row[1]
					for otu_string in row[2:]:
						otu = otu_string.split(self.strain_separator)
						for index, element in enumerate(otu):
							otu[index] = element
						otu_list.append(otu)
					self.otu_lists_by_cutoff_raw[cutoff_of_otu] = {"count": number_of_otu, "otu": otu_list}
					print cutoff_of_otu  # , self.otu_lists_by_cutoff[cutoff_of_otu]
		self.otu_list = self.otu_lists_by_cutoff_raw[self.otu_cutoff]["otu"]
		print "loading finished"

	def otu_to_other_rank(self, otu, rank, excluded_genomes_id_list, unknown_id=None):
		if rank not in self.tmp_tax_list_by_rank:
			self.tmp_tax_list_by_rank[rank] = {}
		species_list = []
		species_list_tmp2 = []
		for item_id in otu:
			if item_id.split("_")[0] in excluded_genomes_id_list:
				continue
			species_id = item_id.split(".")[0]
			if species_id.isdigit():
				if species_id in self.tmp_tax_list_by_rank[rank]:
					species_ncbi = self.tmp_tax_list_by_rank[rank][species_id]
				else:
					species_ncbi = self.taxonomy.get_ncbi_of_rank(rank, int(float(species_id)))
					self.tmp_tax_list_by_rank[rank][species_id] = species_ncbi
				if species_ncbi is not None:
					species_list.append(str(species_ncbi))
					if rank == "genus" or rank == "species":
						species_list_tmp2.append(str(item_id) + " -> " + str(species_ncbi))
				#else:
				#	#species_list.append(str(species_id))
				#	print "MothurOTU: ncbi of rank", rank, " not found:", species_id, "  otu size:", len(otu)
			#else:
			#	species_list.append(str(species_id))
			#	print "MothurOTU: not a ncbi id:", species_id, "  otu size:", len(otu)
		#species_list = species_list_tmp
		if unknown_id is not None and len(species_list) > 0:
			print unknown_id, rank
			self.otu_list_to_stdout(Counter(species_list), True)
			if rank == "genus" or rank == "species":
				print unknown_id, rank
				self.otu_list_to_stdout(species_list_tmp2)
		return species_list

	def otu_to_other_rank_old(self, otu, rank, unknown_genomes_id_column, unknown_id=None):
		if rank not in self.tmp_tax_list_by_rank:
			self.tmp_tax_list_by_rank[rank] = {}
		species_list_tmp = []
		species_list_tmp2 = []
		for item_id in otu:
			if item_id in unknown_genomes_id_column:
				continue
			species_id = item_id.split(".")[0]
			if species_id.isdigit():
				if species_id in self.tmp_tax_list_by_rank[rank]:
					species_ncbi = self.tmp_tax_list_by_rank[rank][species_id]
				else:
					species_ncbi = self.taxonomy.get_ncbi_of_rank(rank, int(float(species_id)))
					self.tmp_tax_list_by_rank[rank][species_id] = species_ncbi
				if species_ncbi is not None:
					species_list_tmp.append(str(species_ncbi))
					if rank == "genus" or rank == "species":
						species_list_tmp2.append(str(item_id) + " -> " + str(species_ncbi))
				#else:
				#	#species_list.append(str(species_id))
				#	print "MothurOTU: ncbi of rank", rank, " not found:", species_id, "  otu size:", len(otu)
			#else:
			#	species_list.append(str(species_id))
			#	print "MothurOTU: not a ncbi id:", species_id, "  otu size:", len(otu)
		species_list = set(species_list_tmp)
		if unknown_id is not None and len(species_list) > 0:
			print unknown_id, rank
			self.otu_list_to_stdout(Counter(species_list_tmp), True)
			if rank == "genus" or rank == "species":
				print unknown_id, rank
				self.otu_list_to_stdout(species_list_tmp2)
		return species_list

	@staticmethod
	def otu_list_to_stdout(otu_list, is_dict=False, max_item_count=10):
		result = ""
		item_count = 0
		if not is_dict:
			for item in otu_list:
				item_count += 1
				if item_count % max_item_count == 0:
					print ""
				print item+",",
		else:
			for item in otu_list:
				item_count += 1
				if item_count % max_item_count == 0:
					print ""
				print "{}: {},".format(item, otu_list[item]),
		print ""
		return result

	def otus_to_other_rank_list(self, otu):
		if self.rank not in self.tmp_tax_list_by_rank:
			self.tmp_tax_list_by_rank[self.rank] = {}
		species_list = []
		for species_id in otu:
			species_id = species_id.split(".")[0]
			if species_id.isdigit():
				if species_id in self.tmp_tax_list_by_rank[self.rank]:
					species_ncbi = self.tmp_tax_list_by_rank[self.rank][species_id]
				else:
					species_ncbi = self.taxonomy.get_ncbi_of_rank(self.rank, species_id)
					self.tmp_tax_list_by_rank[self.rank][species_id] = species_ncbi
				if species_ncbi is not None:
					species_list.append(str(species_ncbi))
				else:
					species_list.append(str(species_id))
					print "MothurOTU: ncbi id not found:", species_id, "  otu size:", len(otu)
			else:
				species_list.append(str(species_id))
				print "MothurOTU: not a ncbi id:", species_id, "  otu size:", len(otu)
		return species_list

	def set_rank(self, rank, unknown_genomes_id_column=[]):
		self.rank = rank
		self.otu_list = []
		current_otu_number = 0
		#number_of_otus = len(self.otu_lists_by_cutoff_raw[self.otu_cutoff]["otu"])
		for otu in self.otu_lists_by_cutoff_raw[self.otu_cutoff]["otu"]:
			#if number_of_otus % 100 == 0:
			#print "{0} of {1} otu done".format(current_otu_number, number_of_otus)
			current_otu_number += 1
			new_otu = list(self.otu_to_other_rank(otu, rank, unknown_genomes_id_column))
			#new_otu = self.otus_to_other_rank_list(otu)
			if len(new_otu) > 0:
				self.otu_list.append(new_otu)
			else:
				self.otu_list.append([])
				#print "No rank for otu:", otu

	def set_cutoff_rank(self, otu_cutoff, rank, unknown_genomes_id_column=[]):
		self.otu_cutoff = otu_cutoff
		self.set_rank(rank, unknown_genomes_id_column)

	def set_cutoff_raw(self, otu_cutoff):
		self.otu_cutoff = otu_cutoff
		self.otu_list = self.otu_lists_by_cutoff_raw[self.otu_cutoff]["otu"]

	def get_otu_of(self, unique_id):
		#unique_id = unique_id.split(".")[0]
		#print unique_id
		#for otu in self.otu_list:
		for index, otu in enumerate(self.otu_lists_by_cutoff_raw[self.otu_cutoff]["otu"]):
			#print otu,
			if unique_id in otu:
				return index, self.otu_list[index]
		#print ""
		return None, None

	def get_otu_list_raw(self):
		return self.otu_lists_by_cutoff_raw[self.otu_cutoff]["otu"]

	def get_otu_list(self):
		return self.otu_list

	def get_number_of_otu(self):
		return self.otu_lists_by_cutoff_raw[self.otu_cutoff]["count"]

	def get_otu_majority_ncbi(self, otu):
		id_count = Counter(otu)
		id_set = set(otu)
		max_count = 0
		majority_id = None
		for unique_id in id_set:
			if unique_id.isdigit() and id_count[unique_id] > max_count:
				majority_id = unique_id
				max_count = id_count[unique_id]
		return majority_id

	def get_otu_ncbi_tax_prediction(self, otu, unknown_genomes_id_column, unknown_id=None):
		rank_index = 1
		otu = self.otu_to_other_rank(otu, self.ranks[rank_index], unknown_genomes_id_column, unknown_id)
		#otu = self.otu_to_other_rank(set(otu), self.ranks[rank_index], unknown_genomes_id_column, unknown_id)
		while len(otu) > 0 and float(max(Counter(otu).iteritems(), key=operator.itemgetter(1))[1])/len(otu) < .9 and (rank_index + 1) < len(self.ranks):
			#print self.ranks[rank_index], Counter(otu)
			rank_index += 1
			tmp = self.otu_to_other_rank(otu, self.ranks[rank_index], unknown_genomes_id_column, unknown_id)
			if len(tmp) == 0:
				break
			otu = tmp
		print ""
		if len(otu) == 0:
			return None, ""
		dominant_id = max(Counter(otu).iteritems(), key=operator.itemgetter(1))[0]
		return dominant_id, self.ranks[rank_index - 1]

	def get_difference_of_number_of_unique_id_to_ids_within_otus(self):
		total_number = 0
		total_set = set()
		for otu in self.otu_list:
			#print "otu:", otu
			total_number += len(otu)
			#print len(new_set), new_set, otu, len(otu), "\n",
			total_set = total_set.union(otu)
		return len(total_set), total_number

	def get_distinct_number_of_rank(self):
		major_otu_list = []
		for otu in self.otu_list:
			majority_id = self.get_otu_majority_ncbi(otu)
			if majority_id is not None:
				major_otu_list.append(majority_id)
		return len(set(major_otu_list))

	#new_set = Counter(otu_to_species_list(otu, taxonomy, "species"))
	def get_best_cutoff(self, rank="species"):
		best_distance = None
		best_cutoff = None
		for otu_cutoff in self.otu_lists_by_cutoff_raw:
			self.set_cutoff_rank(otu_cutoff, rank)
			number_of_unique_ids, difference_within_otus = self.get_difference_of_number_of_unique_id_to_ids_within_otus()
			distinct_number_of_rank = self.get_distinct_number_of_rank()
			#print " distinct_number_of_rank:", cutoff, distinct_number_of_rank
			new_distance = abs(number_of_unique_ids - distinct_number_of_rank) + 1.5 * abs(number_of_unique_ids - difference_within_otus)
			if best_distance is None or new_distance < best_distance:
				best_distance = new_distance
				best_cutoff = otu_cutoff
			print otu_cutoff, "number_of_unique_species:", number_of_unique_ids, " distinct_number_of_rank:", distinct_number_of_rank, " difference_within_otus:", difference_within_otus, " new_distance:", new_distance
		print "Best distance:", best_distance, "Cutoff:", best_cutoff
		return best_cutoff

	def to_string(self, cutoff="unique"):
		for otu_group in self.otu_lists_by_cutoff_raw[cutoff]["otu"]:
			for item in otu_group:
				print str(item),
			print ""
