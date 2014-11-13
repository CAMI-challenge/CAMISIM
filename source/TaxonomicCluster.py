__author__ = 'hofmann'

import sys
import operator
from collections import Counter


class TaxonomicCluster:
	"""Reading and writing a meta table"""
	def __init__(self, mothur_cluster, taxonomy, logger=None):
		self.mothur_cluster = mothur_cluster
		self.taxonomy = taxonomy
		self.logger = logger
		self.rank = "strain"
		self.ranks = ['strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom', 'root']
		self.taxids_by_element = {}

	def cluster_to_other_rank(self, cluster, rank, list_of_excluded_elements, unpublished_sequence_id=None, debug=False):
		ncbi_id_list = []
		index_of_rank = self.ranks.index(rank)
		#if index_of_rank is None:
		for element in cluster:
			if element in list_of_excluded_elements:
				continue
			ncbi_id = element.split('.')[0]
			if not ncbi_id.isdigit():
				continue
			if not element in self.taxids_by_element:
				#self.taxonomy.get_ncbi_of_rank(rank, int(float(ncbi_id)))
				self.taxids_by_element[element] = self.taxonomy.get_lineage_of_legal_ranks(ncbi_id, ranks=self.ranks, default_value=None)
			ncbi_higher_rank = self.taxids_by_element[element][index_of_rank]
			if ncbi_higher_rank is not None:
				ncbi_id_list.append(str(ncbi_higher_rank))
		if debug and unpublished_sequence_id is not None and len(ncbi_id_list) > 0 and self.logger:
			self.logger.debug("{id}\t{rank}".format(id=unpublished_sequence_id, rank=rank))
			self.mothur_cluster.cluster_list_to_handle(Counter(ncbi_id_list), sys.stderr)
		return ncbi_id_list

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
				del tmp
				break
			cluster = tmp
		if len(cluster) == 0:
			return None, ""
		dominant_id = max(Counter(cluster).iteritems(), key=operator.itemgetter(1))[0]
		del cluster
		return dominant_id, self.ranks[rank_index - 1]
