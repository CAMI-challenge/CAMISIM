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
		self.cluster_separator = otu_separator
		self.strain_separator = strain_separator
		self.taxonomy = taxonomy
		#self.header=[]
		self.cluster_cutoff = "unique"
		self.rank = "strain"
		self.lists_of_cluster_by_cutoff_raw = {}
		self.genome_id_to_index_mapping = {}
		self.list_of_cluster = []
		self.read_mothur_clustering_file()
		#self.ranks = ['root', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
		self.ranks = ['strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom', 'root']
		self.tmp_tax_list_by_rank = {}

	def get_clusters(self, cutoffs, unpublished_genome_ids_column):
		result = {}
		number_of_rows = len(unpublished_genome_ids_column)
		for row_index in range(0, number_of_rows):
			unpublished_genome_id = unpublished_genome_ids_column[row_index]
			list_of_index, list_of_clusters = self.get_cluster_of(unpublished_genome_id, cutoffs[row_index])
			if unpublished_genome_id.strip() == "" or len(list_of_clusters) == 0:
				result[unpublished_genome_id] = None
				continue
			result[unpublished_genome_id] = []
			for cluster in list_of_clusters:
				for element in cluster:
					genome_id = ".".join(element.split("_")[0].split(".")[:2])
					#sys.stderr.write("{} vs {}\n".format(unpublished_genome_id, genome_id))
					if genome_id not in unpublished_genome_ids_column:
						result[unpublished_genome_id].append(genome_id)
		return result

	def read_mothur_clustering_file(self):
		self.rank = "strain"
		print "loading mothur file"
		with open(self.file_path) as file_handler:
			for line in file_handler:
				if line[0] != '#' and "label" not in line:
					line = line.strip()
					list_of_cluster = []
					row = line.split(self.cluster_separator)
					#if row[0] == "unique":
					cutoff_of_cluster = row[0]
					self.genome_id_to_index_mapping[cutoff_of_cluster] = {}
					number_of_cluster = row[1]
					cluster_index = 0
					for cluster_as_string in row[2:]:
						list_of_strains = cluster_as_string.split(self.strain_separator)
						for strain in list_of_strains:
							#genome_id = strain.split('_')[0]
							genome_id = ".".join(strain.split("_")[0].split(".")[:2])
							if genome_id not in self.genome_id_to_index_mapping[cutoff_of_cluster]:
								self.genome_id_to_index_mapping[cutoff_of_cluster][genome_id] = []
							self.genome_id_to_index_mapping[cutoff_of_cluster][genome_id].append(cluster_index)
						list_of_cluster.append(list_of_strains)
						cluster_index += 1
					self.lists_of_cluster_by_cutoff_raw[cutoff_of_cluster] = {"count": number_of_cluster, "otu": list_of_cluster}
					print cutoff_of_cluster  # , self.otu_lists_by_cutoff[cutoff_of_otu]
		self.list_of_cluster = self.lists_of_cluster_by_cutoff_raw[self.cluster_cutoff]["otu"]
		print "loading finished"

	def cluster_to_other_rank(self, cluster, rank, list_of_excluded_genomes_ids, unpublished_sequence_id=None, debug=False):
		if rank not in self.tmp_tax_list_by_rank:
			self.tmp_tax_list_by_rank[rank] = {}
		ncbi_id_list = []
		ncbi_id_list_details = []
		for element in cluster:
			#genome_id = element.split("_")[0]
			genome_id = ".".join(element.split("_")[0].split(".")[:2])
			if genome_id in list_of_excluded_genomes_ids:
				#sys.stderr.write("genome_id: {} {}\n".format(genome_id, element))
				continue
			ncbi_id = genome_id.split(".")[0]
			if not ncbi_id.isdigit():
				continue
			if ncbi_id in self.tmp_tax_list_by_rank[rank]:
				ncbi_higher_rank = self.tmp_tax_list_by_rank[rank][ncbi_id]
			else:
				ncbi_higher_rank = self.taxonomy.get_ncbi_of_rank(rank, int(float(ncbi_id)))
				# put into quick loopup list, for speedup
				self.tmp_tax_list_by_rank[rank][ncbi_id] = ncbi_higher_rank
			if ncbi_higher_rank is not None:
				ncbi_id_list.append(str(ncbi_higher_rank))
				if rank == "genus" or rank == "species":
					ncbi_id_list_details.append(str(element) + " -> " + str(ncbi_higher_rank))
		if debug and unpublished_sequence_id is not None and len(ncbi_id_list) > 0:
			print unpublished_sequence_id, rank
			self.cluster_list_to_stdout(Counter(ncbi_id_list), True)
			if rank == "genus" or rank == "species":
				print unpublished_sequence_id, rank
				self.cluster_list_to_stdout(ncbi_id_list_details)
		return ncbi_id_list

	@staticmethod
	def cluster_list_to_stdout(otu_list, is_dict=False, max_item_count=10):
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

	def cluster_to_other_rank_list_old(self, cluster):
		if self.rank not in self.tmp_tax_list_by_rank:
			self.tmp_tax_list_by_rank[self.rank] = {}
		tax_id_list = []
		for element in cluster:
			ncbi_id = element.split(".")[0]
			if ncbi_id.isdigit():
				if ncbi_id in self.tmp_tax_list_by_rank[self.rank]:
					species_ncbi = self.tmp_tax_list_by_rank[self.rank][ncbi_id]
				else:
					species_ncbi = self.taxonomy.get_ncbi_of_rank(self.rank, ncbi_id)
					self.tmp_tax_list_by_rank[self.rank][ncbi_id] = species_ncbi
				if species_ncbi is not None:
					tax_id_list.append(str(species_ncbi))
				else:
					tax_id_list.append(str(ncbi_id))
					print "MothurOTU: ncbi id not found:", ncbi_id, "  otu size:", len(cluster)
			else:
				tax_id_list.append(str(ncbi_id))
				print "MothurOTU: not a ncbi id:", ncbi_id, "  otu size:", len(cluster)
		return tax_id_list

	def set_rank(self, rank, unpublished_genomes_id_column=[]):
		self.rank = rank
		self.list_of_cluster = []
		number_of_clusters_done = 0
		#number_of_clusters = len(self.otu_lists_by_cutoff_raw[self.otu_cutoff]["otu"])
		for cluster in self.lists_of_cluster_by_cutoff_raw[self.cluster_cutoff]["otu"]:
			#if number_of_otus % 100 == 0:
			#print "{0} of {1} otu done".format(number_of_clusters_done, number_of_clusters)
			number_of_clusters_done += 1
			cluster_of_higher_rank = self.cluster_to_other_rank(cluster, rank, unpublished_genomes_id_column)
			#new_otu = self.otus_to_other_rank_list(otu)
			if len(cluster_of_higher_rank) > 0:
				self.list_of_cluster.append(cluster_of_higher_rank)
			else:
				self.list_of_cluster.append([])
				#print "No rank for otu:", otu

	def set_cutoff_rank(self, cluster_cutoff, rank, unpublished_genomes_id_column=None):
		if unpublished_genomes_id_column is None:
			unpublished_genomes_id_column = []
		self.cluster_cutoff = cluster_cutoff
		self.set_rank(rank, unpublished_genomes_id_column)

	def set_cutoff_raw(self, cluster_cutoff):
		self.cluster_cutoff = cluster_cutoff
		self.list_of_cluster = self.lists_of_cluster_by_cutoff_raw[self.cluster_cutoff]["otu"]

	def genome_id_exists(self, genome_id):
		genome_id = ".".join(genome_id.split("_")[0].split(".")[:2])
		if genome_id not in self.genome_id_to_index_mapping[self.cluster_cutoff]:
			print "Warning: genome_id not found in clusters", self.cluster_cutoff, genome_id
			return False
		return True

	def get_cluster(self, cutoff, cluster_index):
		if cutoff not in self.lists_of_cluster_by_cutoff_raw:
			print "MothurOTU: bad cutoff {}".format(cutoff)
			return None
		if cluster_index >= len(self.lists_of_cluster_by_cutoff_raw[cutoff]["otu"]):
			print "MothurOTU: bad cluster_id".format(cluster_index)
			return None
		return self.lists_of_cluster_by_cutoff_raw[cutoff]["otu"][cluster_index]

	def get_cluster_of(self, genome_id, cluster_cutoff=None):
		if cluster_cutoff is None:
			cluster_cutoff = self.cluster_cutoff
		if not self.genome_id_exists(genome_id) or genome_id.strip() == '' or cluster_cutoff.strip() == '':
			return [], []
		list_of_index = self.genome_id_to_index_mapping[cluster_cutoff][genome_id]
		if len(set(list_of_index)) > 1:
			print "Warning: multiple marker genes in different clusters", cluster_cutoff, genome_id, set(list_of_index)
		#index = list_of_index[0]
		return list_of_index, [self.list_of_cluster[index] for index in list_of_index]

	def get_cluster_list_raw(self):
		return self.lists_of_cluster_by_cutoff_raw[self.cluster_cutoff]["otu"]

	def get_cluster_list(self):
		return self.list_of_cluster

	def get_number_of_clusters(self):
		return self.lists_of_cluster_by_cutoff_raw[self.cluster_cutoff]["count"]

	def get_cluster_majority_ncbi(self, cluster):
		id_count = Counter(cluster)
		id_set = set(cluster)
		max_count = 0
		majority_id = None
		for unique_id in id_set:
			if unique_id.isdigit() and id_count[unique_id] > max_count:
				majority_id = unique_id
				max_count = id_count[unique_id]
		return majority_id

	def get_cluster_ncbi_tax_prediction(self, cluster, unpublished_genome_ids_column, unpublished_id=None):
		rank_index = 1
		cluster = self.cluster_to_other_rank(cluster, self.ranks[rank_index], unpublished_genome_ids_column, unpublished_id)
		while len(cluster) > 0 \
				and float(max(Counter(cluster).iteritems(), key=operator.itemgetter(1))[1])/len(cluster) < .9 \
				and (rank_index + 1) < len(self.ranks):
			#print self.ranks[rank_index], Counter(otu)
			rank_index += 1
			tmp = self.cluster_to_other_rank(cluster, self.ranks[rank_index], unpublished_genome_ids_column, unpublished_id)
			if len(tmp) == 0:
				break
			cluster = tmp
		print ""
		if len(cluster) == 0:
			return None, ""
		dominant_id = max(Counter(cluster).iteritems(), key=operator.itemgetter(1))[0]
		return dominant_id, self.ranks[rank_index - 1]

	def get_difference_of_number_of_unique_id_to_ids_within_otus(self):
		total_number = 0
		total_set = set()
		for otu in self.list_of_cluster:
			#print "otu:", otu
			total_number += len(otu)
			#print len(new_set), new_set, otu, len(otu), "\n",
			total_set = total_set.union(otu)
		return len(total_set), total_number

	def get_distinct_number_of_rank(self):
		major_otu_list = []
		for otu in self.list_of_cluster:
			majority_id = self.get_cluster_majority_ncbi(otu)
			if majority_id is not None:
				major_otu_list.append(majority_id)
		return len(set(major_otu_list))

	#new_set = Counter(otu_to_species_list(otu, taxonomy, "species"))
	def get_best_cutoff(self, rank="species"):
		best_distance = None
		best_cutoff = None
		for otu_cutoff in self.lists_of_cluster_by_cutoff_raw:
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
		for otu_group in self.lists_of_cluster_by_cutoff_raw[cutoff]["otu"]:
			for item in otu_group:
				print str(item),
			print ""
