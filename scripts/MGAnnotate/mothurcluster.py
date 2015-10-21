__author__ = 'hofmann'

import sys
import os
import math
import textwrap
from scripts.Validator.validator import Validator
from scripts.MetaDataTable.metadatatable import MetadataTable


class MothurCluster(Validator):
	"""Reading and writing a meta table"""

	_label = "MothurCluster"

	def __init__(
		self, precision, otu_separator="\t", element_separator=",", data_table_iid_mapping=None,
		logfile=None, verbose=False, debug=False):
		assert isinstance(precision, int)
		assert data_table_iid_mapping is None or isinstance(data_table_iid_mapping, MetadataTable)
		super(MothurCluster, self).__init__(logfile=logfile, verbose=verbose, debug=debug)
		self._precision = int(math.log10(precision))
		self.cluster_separator = otu_separator
		self.element_separator = element_separator
		self._cutoff_to_cluster = {}
		self._gid_to_cluster_index_list = {}
		self._unique_threshold = "unique"
		self._iid_gid = {}
		if data_table_iid_mapping is not None:
			self._iid_gid = data_table_iid_mapping.get_map(0, 1)

	def get_max_threshold(self):
		lists_of_thresholds = list(self._cutoff_to_cluster.keys())
		lists_of_thresholds.remove(self._unique_threshold)
		lists_of_thresholds = sorted(lists_of_thresholds)
		if len(lists_of_thresholds) == 0:
			return self._unique_threshold
		return lists_of_thresholds[-1]

	def get_clusters_of_elements(self, threshold, list_of_query_gid):
		if not threshold == "unique":
			assert isinstance(threshold, (int, float))
		assert isinstance(list_of_query_gid, list)

		result = {}
		for gid in list_of_query_gid:
			if len(gid.strip()) == 0:
				self._logger.warning("Empty gid in list of query gid!")
				result[gid] = None
				continue
			list_of_index, list_of_clusters = self.get_cluster_of_cutoff_of_gid(threshold, gid)
			if len(list_of_clusters) == 0:
				result[gid] = None
				continue
			result[gid] = []
			for cluster in list_of_clusters:
				for element in cluster:
					if self._iid_gid[element] in list_of_query_gid:
						continue
					result[gid].append(element)
		return result

	def read(self, file_path, list_of_query_gid=None):
		if not os.path.isfile(file_path):
			self._logger.error("No file found at: '{}'".format(file_path))
			return
		self._logger.info("Reading cluster file '{}'".format(file_path))

		with open(file_path) as file_handler:
			self._gid_to_cluster_index_list = {}
			for line in file_handler:
				line = line.strip()
				if line.startswith('#') or line.startswith("label") or len(line) == 0:
					continue
				list_of_cluster = []
				row = line.split(self.cluster_separator)
				cutoff = row[0]
				if cutoff.isdigit():
					cutoff = str(float(cutoff))
				cluster_amount = row[1]
				self._logger.info("Reading threshold: {}".format(cutoff))
				self._gid_to_cluster_index_list[cutoff] = {}
				cluster_index = 0
				for cluster_as_string in row[2:]:
					list_of_elements = cluster_as_string.split(self.element_separator)
					list_of_cluster.append(list_of_elements)
					if list_of_query_gid is None:
						continue
					list_of_gid = [self._iid_gid[gid] for gid in list_of_elements]
					for gid in list_of_gid:
						if gid not in list_of_query_gid:
							continue
						if gid not in self._gid_to_cluster_index_list[cutoff]:
							self._gid_to_cluster_index_list[cutoff][gid] = []
						self._gid_to_cluster_index_list[cutoff][gid].append(cluster_index)
					cluster_index += 1
				self._cutoff_to_cluster[cutoff] = {"count": cluster_amount, "cluster": list_of_cluster}
			self._logger.info("Reading finished")

	def get_prediction_thresholds(self, minimum=0):
		subset = set()
		list_of_cutoff = list(self._cutoff_to_cluster.keys())
		list_as_float = []
		for cutoff in list_of_cutoff:
			if '.' not in cutoff:
				continue
			list_as_float.append(float(cutoff))

		for cutoff in list_of_cutoff:
			if '.' not in cutoff:
				continue
			threshold = round(float(cutoff), self._precision)

			if threshold >= minimum and threshold in list_as_float:
				subset.add(threshold)

		return subset

	def get_sorted_lists_of_cutoffs(self, reverse=False):
		lists_of_cutoff = list(self._cutoff_to_cluster.keys())
		lists_of_cutoff.remove("unique")
		lists_of_cutoff = sorted(lists_of_cutoff, reverse=reverse)
		if not reverse:
			tmp = lists_of_cutoff
			lists_of_cutoff = ["unique"]
			lists_of_cutoff.extend(tmp)
		else:
			lists_of_cutoff.append("unique")
		return lists_of_cutoff

	@staticmethod
	def cluster_list_to_stream(list_of_cluster, stream=sys.stdout, width=80):
		if not isinstance(list_of_cluster, dict):
			stream.write(textwrap.fill(", ".join(list_of_cluster) + '\n', width))
		else:
			line = ", ".join('{}: {}'.format(key, value) for key, value in list_of_cluster.items())
			stream.write(textwrap.fill(line, width) + '\n')

	def element_exists(self, threshold, gid):
		if not threshold == "unique":
			assert isinstance(threshold, (int, float))
			threshold = "{th:.{pre}f}".format(th=threshold, pre=self._precision)

		if threshold not in self._gid_to_cluster_index_list:
			self._logger.error("Cutoff key error: {}\nAvailable keys: '{}'".format(threshold, ','.join(self._gid_to_cluster_index_list.keys())))
			return False

		if gid not in self._gid_to_cluster_index_list[threshold]:
			# 	self.logger.warning("{} not found in {}".format(element, cutoff))
			return False
		return True

	def get_cluster_of_cutoff_of_index(self, cutoff, cluster_index):
		if cutoff not in self._cutoff_to_cluster:
			self._logger.error("Bad cutoff {}".format(cutoff))
			return None
		if cluster_index >= len(self._cutoff_to_cluster[cutoff]["cluster"]):
			self._logger.error("Bad cluster index".format(cluster_index))
			return None
		return self._cutoff_to_cluster[cutoff]["cluster"][cluster_index]

	def get_cluster_of_cutoff_of_gid(self, threshold, gid):
		if not threshold == "unique":
			assert isinstance(threshold, (int, float))
			threshold = "{th:.{pre}f}".format(th=threshold, pre=self._precision)

		if threshold not in self._cutoff_to_cluster:
			self._logger.error("Bad cutoff: {}".format(threshold))
			return [], []
		if not self.element_exists(float(threshold), gid) or gid.strip() == '' or threshold.strip() == '':
			self._logger.warning("Bad element: {} in {}".format(gid, threshold))
			return [], []
		list_of_index = self._gid_to_cluster_index_list[threshold][gid]
		if len(set(list_of_index)) > 1:
			self._logger.debug("{}: Multiple elements found. {}: {}".format(threshold, gid, ", ".join([str(item) for item in set(list_of_index)])))
		return list_of_index, [self._cutoff_to_cluster[threshold]["cluster"][index] for index in list_of_index]

	def get_cluster_of_cutoff(self, threshold="unique"):
		if not threshold == "unique":
			assert isinstance(threshold, (int, float))
			threshold = "{th:.{pre}f}".format(th=threshold, pre=self._precision)

		if threshold not in self._cutoff_to_cluster:
			if self._logger:
				self._logger.error("Bad cutoff: {}".format(threshold))
			return None
		return self._cutoff_to_cluster[threshold]["cluster"]

	def get_cluster_count_of_cutoff(self, threshold="unique"):
		if not threshold == "unique":
			assert isinstance(threshold, (int, float))
			threshold = "{th:.{pre}f}".format(th=threshold, pre=self._precision)

		if threshold not in self._cutoff_to_cluster:
			if self._logger:
				self._logger.error("Bad cutoff: {}".format(threshold))
			return None
		return self._cutoff_to_cluster[threshold]["count"]

	def to_string_cutoff(self, threshold="unique"):
		if not threshold == "unique":
			assert isinstance(threshold, (int, float))
			threshold = "{th:.{pre}f}".format(th=threshold, pre=self._precision)

		result_string = "{}\n".format(threshold)
		if threshold not in self._cutoff_to_cluster:
			if self._logger:
				self._logger.error("Bad cutoff: {}".format(threshold))
			return None
		for otu_group in self._cutoff_to_cluster[threshold]["cluster"]:
			result_string += ", ".join(otu_group) + '\n'
		return result_string
