__author__ = 'hofmann'

import sys
import os
import math
import textwrap


class MothurCluster:
	"""Reading and writing a meta table"""
	def __init__(self, precision, otu_separator="\t", element_separator=",", logger=None):
		assert isinstance(precision, int)

		self._precision = int(math.log10(precision))
		self._logger = logger
		self.cluster_separator = otu_separator
		self.element_separator = element_separator
		self._cluster_by_cutoff = {}
		self.element_to_index_mapping = {}
		self._unique_threshold = "unique"

	@staticmethod
	def element_to_genome_id(element):
		if '.' in element:
			prefix, suffix = element.split(".", 1)
			suffix = suffix.split("_", 1)[0]
			return "{pre}.{suf}".format(pre=prefix, suf=suffix)
		elif element.count('_') >= 5:
			return element.rsplit('_', 5)[0]
		else:
			return element

	def get_max_threshold(self):
		lists_of_thresholds = list(self._cluster_by_cutoff.keys())
		lists_of_thresholds.remove(self._unique_threshold)
		lists_of_thresholds = sorted(lists_of_thresholds)
		if len(lists_of_thresholds) == 0:
			return self._unique_threshold
		return lists_of_thresholds[-1]

	def get_clusters_of_elements(self, threshold, list_of_elements):
		if not threshold == "unique":
			assert isinstance(threshold, (int, float))
		assert isinstance(list_of_elements, list)

		result = {}
		number_of_elements = len(list_of_elements)
		for element_index in xrange(0, number_of_elements):
			element = list_of_elements[element_index]
			if len(element.strip()) == 0:
				result[element] = None
				continue
			list_of_index, list_of_clusters = self.get_cluster_of_cutoff_of_element(threshold, element)
			if len(list_of_clusters) == 0:
				result[element] = None
				continue
			result[element] = []
			for cluster in list_of_clusters:
				for c_element in cluster:
					if c_element in list_of_elements:
						continue
					#genome_id = self.element_to_genome_id(element)
					result[element].append(c_element)
		return result

	def read(self, file_path):
		if not os.path.isfile(file_path):
			if self._logger:
				self._logger.error("[MothurCluster] No file found at: '{}'".format(file_path))
			return
		if self._logger:
			self._logger.info("[MothurCluster] Reading cluster file '{}'".format(file_path))

		with open(file_path) as file_handler:
			self.element_to_index_mapping = {}
			for line in file_handler:
				if line.startswith('#') or line.startswith("label") or len(line) < 2:
					continue
				line = line.strip()
				list_of_cluster = []
				row = line.split(self.cluster_separator)
				cutoff = row[0]
				if cutoff.isdigit():
					cutoff = str(float(cutoff))
				cluster_amount = row[1]
				if self._logger:
					self._logger.info("[MothurCluster] Reading threshold: {}".format(cutoff))
				self.element_to_index_mapping[cutoff] = {}
				cluster_index = 0
				for cluster_as_string in row[2:]:
					list_of_elements = cluster_as_string.split(self.element_separator)
					new_list_of_elements = [self.element_to_genome_id(element) for element in list_of_elements]
					for element in new_list_of_elements:
						#genome_id = self.element_to_genome_id(element)
						if element not in self.element_to_index_mapping[cutoff]:
							self.element_to_index_mapping[cutoff][element] = []
						self.element_to_index_mapping[cutoff][element].append(cluster_index)
					list_of_cluster.append(new_list_of_elements)
					cluster_index += 1
				self._cluster_by_cutoff[cutoff] = {"count": cluster_amount, "cluster": list_of_cluster}
			if self._logger:
				self._logger.info("[MothurCluster] Reading finished")

	def get_prediction_thresholds(self, minimum=0):
		subset = set()
		list_of_cutoff = list(self._cluster_by_cutoff.keys())
		list_as_float = []
		for cutoff in list_of_cutoff:
			if not '.' in cutoff:
				continue
			list_as_float.append(float(cutoff))

		for cutoff in list_of_cutoff:
			if not '.' in cutoff:
				continue
			threshold = round(float(cutoff), self._precision)

			if threshold >= minimum and threshold in list_as_float:
				subset.add(threshold)

		return subset

	def get_sorted_lists_of_cutoffs(self, reverse=False):
		lists_of_cutoff = list(self._cluster_by_cutoff.keys())
		lists_of_cutoff.remove("unique")
		#lists_of_cutoff = sorted(set([str(round(float(cutoff), precision)) for cutoff in lists_of_cutoff]), reverse=reverse)
		lists_of_cutoff = sorted(lists_of_cutoff, reverse=reverse)
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

	def element_exists(self, threshold, element):
		if not threshold == "unique":
			assert isinstance(threshold, (int, float))
			threshold = "{th:.{pre}f}".format(th=threshold, pre=self._precision)

		if threshold not in self.element_to_index_mapping:
			if self._logger:
				self._logger.error("[MothurCluster] Cutoff key error: {}\nAvailable keys: '{}'".format(threshold, ','.join(self.element_to_index_mapping.keys())))
			return False

		#element = ".".join(element.split("_")[0].split(".")[:2])
		if element not in self.element_to_index_mapping[threshold]:
			#if self.logger:
			#	self.logger.warning("{} not found in {}".format(element, cutoff))
			return False
		return True

	def get_cluster_of_cutoff_of_index(self, cutoff, cluster_index):
		if cutoff not in self._cluster_by_cutoff:
			if self._logger:
				self._logger.error("[MothurCluster] Bad cutoff {}".format(cutoff))
			return None
		if cluster_index >= len(self._cluster_by_cutoff[cutoff]["cluster"]):
			if self._logger:
				self._logger.error("[MothurCluster] Bad cluster index".format(cluster_index))
			return None
		return self._cluster_by_cutoff[cutoff]["cluster"][cluster_index]

	def get_cluster_of_cutoff_of_element(self, threshold, element):
		if not threshold == "unique":
			assert isinstance(threshold, (int, float))
			threshold = "{th:.{pre}f}".format(th=threshold, pre=self._precision)

		if threshold not in self._cluster_by_cutoff:
			if self._logger:
				self._logger.error("[MothurCluster] Bad cutoff: {}".format(threshold))
			return [], []
		if not self.element_exists(float(threshold), element) or element.strip() == '' or threshold.strip() == '':
			if self._logger:
				self._logger.warning("[MothurCluster] Bad element: {} in {}".format(element, threshold))
			return [], []
		list_of_index = self.element_to_index_mapping[threshold][element]
		if len(set(list_of_index)) > 1:
			if self._logger:
				self._logger.debug("[MothurCluster] {}: Multiple elements found. {}: {}".format(threshold, element, ", ".join([str(item) for item in set(list_of_index)])))
				#print "Warning: multiple marker genes in different clusters", cutoff, element, set(list_of_index)
		return list_of_index, [self._cluster_by_cutoff[threshold]["cluster"][index] for index in list_of_index]

	def get_cluster_of_cutoff(self, threshold="unique"):
		if not threshold == "unique":
			assert isinstance(threshold, (int, float))
			threshold = "{th:.{pre}f}".format(th=threshold, pre=self._precision)

		if threshold not in self._cluster_by_cutoff:
			if self._logger:
				self._logger.error("[MothurCluster] Bad cutoff: {}".format(threshold))
			return None
		return self._cluster_by_cutoff[threshold]["cluster"]

	def get_cluster_count_of_cutoff(self, threshold="unique"):
		if not threshold == "unique":
			assert isinstance(threshold, (int, float))
			threshold = "{th:.{pre}f}".format(th=threshold, pre=self._precision)

		if threshold not in self._cluster_by_cutoff:
			if self._logger:
				self._logger.error("[MothurCluster] Bad cutoff: {}".format(threshold))
			return None
		return self._cluster_by_cutoff[threshold]["count"]

	def to_string_cutoff(self, threshold="unique"):
		if not threshold == "unique":
			assert isinstance(threshold, (int, float))
			threshold = "{th:.{pre}f}".format(th=threshold, pre=self._precision)

		result_string = "{}\n".format(threshold)
		if threshold not in self._cluster_by_cutoff:
			if self._logger:
				self._logger.error("[MothurCluster] Bad cutoff: {}".format(threshold))
			return None
		for otu_group in self._cluster_by_cutoff[threshold]["cluster"]:
			result_string += ", ".join(otu_group)+'\n'
		return result_string