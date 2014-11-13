__author__ = 'hofmann'

import sys
import textwrap


class MothurCluster:
	"""Reading and writing a meta table"""
	def __init__(self, otu_separator="\t", element_separator=",", logger=None):
		self.logger = logger
		self.cluster_separator = otu_separator
		self.element_separator = element_separator
		self._cluster_by_cutoff = {}
		self.element_to_index_mapping = {}

	@staticmethod
	def element_to_genome_id(element):
		if '.' in element:
			return ".".join(element.split("_")[0].split(".")[:2])
		else:
			return element

	def get_clusters_of_elements(self, cutoffs, list_of_elements):
		result = {}
		number_of_elements = len(list_of_elements)
		for element_index in range(0, number_of_elements):
			element = list_of_elements[element_index]
			if len(element.strip()) == 0:
				result[element] = None
				continue
			list_of_index, list_of_clusters = self.get_cluster_of_cutoff_of_element(cutoffs[element_index], element)
			if len(list_of_clusters) == 0:
				result[element] = None
				continue
			result[element] = []
			for cluster in list_of_clusters:
				for element in cluster:
					if element in list_of_elements:
						continue
					#genome_id = self.element_to_genome_id(element)
					result[element].append(element)
		return result

	def read_mothur_clustering_file(self, file_handler):
		self.element_to_index_mapping = {}
		if self.logger:
			self.logger.info("loading mothur cluster file")

		for line in file_handler:
			if line.startswith('#') or line.startswith("label") or len(line) < 2:
				continue
			line = line.strip()
			list_of_cluster = []
			row = line.split(self.cluster_separator)
			cutoff = row[0]
			cluster_amount = row[1]
			if self.logger:
				self.logger.info("cutoff: {}".format(cutoff))
			self.element_to_index_mapping[cutoff] = {}
			cluster_index = 0
			for cluster_as_string in row[2:]:
				list_of_elements = cluster_as_string.split(self.element_separator)
				for element in list_of_elements:
					genome_id = self.element_to_genome_id(element)
					if genome_id not in self.element_to_index_mapping[cutoff]:
						self.element_to_index_mapping[cutoff][genome_id] = []
					self.element_to_index_mapping[cutoff][genome_id].append(cluster_index)
				list_of_cluster.append(list_of_elements)
				cluster_index += 1
			self._cluster_by_cutoff[cutoff] = {"count": cluster_amount, "cluster": list_of_cluster}
		if self.logger:
			self.logger.info("loading finished")

	def get_sorted_lists_of_cutoffs(self, precision=2, reverse=False):
		lists_of_cutoff = list(self._cluster_by_cutoff.keys())
		del lists_of_cutoff[lists_of_cutoff.index("unique")]
		lists_of_cutoff = sorted(set([str(round(float(cutoff), precision)) for cutoff in lists_of_cutoff]), reverse=reverse)
		if not reverse:
			tmp = lists_of_cutoff
			lists_of_cutoff = ["unique"]
			lists_of_cutoff.extend(tmp)
		else:
			lists_of_cutoff.append("unique")
		return lists_of_cutoff

	@staticmethod
	def cluster_list_to_handle(list_of_cluster, handle=sys.stdout, width=80):
		if not isinstance(list_of_cluster, dict):
			handle.write(textwrap.fill(", ".join(list_of_cluster)+'\n', width))
		else:
			line = ", ".join('{}: {}'.format(key, value) for key, value in list_of_cluster.items())
			handle.write(textwrap.fill(line, width)+'\n')

	def element_exists(self, cutoff, element):
		#element = ".".join(element.split("_")[0].split(".")[:2])
		if element not in self.element_to_index_mapping[cutoff]:
			#if self.logger:
			#	self.logger.warning("{} not found in {}".format(element, cutoff))
			return False
		return True

	def get_cluster_of_cutoff_of_index(self, cutoff, cluster_index):
		if cutoff not in self._cluster_by_cutoff:
			if self.logger:
				self.logger.error("MothurCluster: bad cutoff {}".format(cutoff))
			return None
		if cluster_index >= len(self._cluster_by_cutoff[cutoff]["cluster"]):
			if self.logger:
				self.logger.error("MothurCluster: bad cluster index".format(cluster_index))
			return None
		return self._cluster_by_cutoff[cutoff]["cluster"][cluster_index]

	def get_cluster_of_cutoff_of_element(self, cutoff, element):
		if cutoff not in self._cluster_by_cutoff:
			if self.logger:
				self.logger.error("Bad cutoff: {}".format(cutoff))
			return [], []
		if not self.element_exists(cutoff, element) or element.strip() == '' or cutoff.strip() == '':
			if self.logger:
				self.logger.warning("Bad element: {} in {}".format(element, cutoff))
			return [], []
		list_of_index = self.element_to_index_mapping[cutoff][element]
		if len(set(list_of_index)) > 1:
			if self.logger:
				self.logger.warning("{}: Multiple elements found. {}: {}".format(cutoff, element, ", ".join(set(list_of_index))))
				#print "Warning: multiple marker genes in different clusters", cutoff, element, set(list_of_index)
		return list_of_index, [self._cluster_by_cutoff[cutoff]["cluster"][index] for index in list_of_index]

	def get_cluster_of_cutoff(self, cutoff="unique"):
		if cutoff not in self._cluster_by_cutoff:
			if self.logger:
				self.logger.error("Bad cutoff: {}".format(cutoff))
			return None
		return self._cluster_by_cutoff[cutoff]["cluster"]

	def get_cluster_count_of_cutoff(self, cutoff="unique"):
		if cutoff not in self._cluster_by_cutoff:
			if self.logger:
				self.logger.error("Bad cutoff: {}".format(cutoff))
			return None
		return self._cluster_by_cutoff[cutoff]["count"]

	def to_string_cutoff(self, cutoff="unique"):
		result_string = "{}\n".format(cutoff)
		if cutoff not in self._cluster_by_cutoff:
			if self.logger:
				self.logger.error("Bad cutoff: {}".format(cutoff))
			return None
		for otu_group in self._cluster_by_cutoff[cutoff]["cluster"]:
			result_string += ", ".join(otu_group)+'\n'
		return result_string