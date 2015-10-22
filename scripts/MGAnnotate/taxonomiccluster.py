__author__ = 'hofmann'

import operator
from collections import Counter
from scripts.Validator.validator import Validator
from scripts.MGAnnotate.mothurcluster import MothurCluster
from scripts.NcbiTaxonomy.ncbitaxonomy import NcbiTaxonomy


class TaxonomicCluster(Validator):
	"""Reading and writing a meta table"""

	_label = "TaxonomicCluster"

	def __init__(self, mothur_cluster, taxonomy, iid_tid_map, minimum_support=.9, logfile=None, verbose=True, debug=False):
		"""

		@param mothur_cluster:
		@type mothur_cluster: MothurCluster
		@param taxonomy:
		@type taxonomy: NcbiTaxonomy
		@param iid_tid_map:
		@type iid_tid_map: dict[str|unicode, str|unicode]
		@param minimum_support:
		@type minimum_support: float
		@param logfile: file handler or file path to a log file
		@type logfile: file | FileIO | StringIO | basestring
		@param verbose: Not verbose means that only warnings and errors will be past to stream
		@type verbose: bool
		@param debug: Display debug messages
		@type debug: bool
		"""
		assert isinstance(mothur_cluster, MothurCluster)
		assert isinstance(taxonomy, NcbiTaxonomy)
		assert isinstance(iid_tid_map, dict)
		assert isinstance(minimum_support, float)
		super(TaxonomicCluster, self).__init__(logfile=logfile, verbose=verbose, debug=debug)
		self._mothur_cluster = mothur_cluster
		self.taxonomy = taxonomy
		self._ranks = ['strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']  # , 'root'
		self._iid_tid_map = iid_tid_map
		self._iid_to_tid_lineage = {}
		self._minimum_support = minimum_support

	def predict_tax_id_of(self, cluster_raw, lowest_predicted_novelty=None):
		"""

		@param cluster_raw:
		@type cluster_raw: list[str|unicode]
		@param lowest_predicted_novelty:
		@type lowest_predicted_novelty: None|dict[str|unicode, str|unicode]

		@return:
		@rtype: tuple[None|str|unicode, str|unicode, int|long]
		"""
		assert isinstance(cluster_raw, list)
		assert isinstance(lowest_predicted_novelty, dict)
		list_of_valid_iid = self.load_lineages(cluster_raw)
		root = {"count": 0, "c": {}, 'p': None}

		total_count = [0] * len(self._ranks)
		for iid in list_of_valid_iid:
			node = root
			for rank_index in xrange(len(self._ranks) - 1, 0, -1):
				tax_id = self._iid_to_tid_lineage[iid][rank_index]
				if tax_id is None:
					break
				if tax_id not in node["c"]:
					node["c"][tax_id] = {"count": 0, "c": {}, 'p': node, 'r': rank_index, "id": tax_id}
				node["c"][tax_id]["count"] += 1
				node = node["c"][tax_id]
				total_count[rank_index] += 1
		# print 'root', root
		# print 'total_count', total_count
		node = root
		max_child_node = None
		while node is not None:
			max_node_count = 0
			for child_node in node["c"]:
				if max_child_node is not None and node["c"][child_node]["count"] == max_child_node["count"]:
					max_node_count += 1
				if max_child_node is None or node["c"][child_node]["count"] > max_child_node["count"]:
					max_child_node = node["c"][child_node]
					max_node_count = 1
			if max_child_node is None or len(max_child_node["c"]) == 0 or max_node_count > 1:
				# if max_node_count > 1:
				# impossible to get over 50% support
				# 	self.logger.warning("[TaxonomicCluster] max_node_count > 1: {id}, {count}".format(id=max_child_node["id"], count=max_node_count))
				node = None
			else:
				node = max_child_node
				max_child_node = None

		# print 'max_child_node', max_child_node
		# print 'total_count', total_count
		# sys.exit()
		if max_child_node is None:
			# print 'root', root
			return None, "", 0

		# print 'max_child_node', max_child_node
		list_of_candidate = [max_child_node]
		parent_node = max_child_node['p']
		while parent_node is not None:
			list_of_candidate.append(parent_node)
			parent_node = parent_node['p']

		for node in list_of_candidate:
			if float(node["count"]) / total_count[node['r']] >= self._minimum_support:
				novelty = self._ranks[node['r'] - 1]
				# if previous novelty lower
				if lowest_predicted_novelty["novelty"] in self._ranks and self._ranks.index(lowest_predicted_novelty["novelty"]) + 1 <= node['r']:
					novelty = lowest_predicted_novelty["novelty"]
				return node["id"], novelty, node["count"]
		return None, "", 0

	def cluster_to_other_rank(self, list_of_iid, index_of_rank, unpublished_sequence_id=None):
		"""

		@param list_of_iid:
		@type list_of_iid: list[str|unicode] | set[str|unicode]
		@param index_of_rank:
		@type index_of_rank: int | long
		@param unpublished_sequence_id:
		@type unpublished_sequence_id: None | str|unicode

		@return:
		@rtype: list[str|unicode]
		"""
		assert isinstance(list_of_iid, list)
		assert isinstance(index_of_rank, (int, long))
		assert unpublished_sequence_id is None or isinstance(unpublished_sequence_id, basestring)
		ncbi_id_list = []
		for iid in list_of_iid:
			ncbi_higher_rank = self._iid_to_tid_lineage[iid][index_of_rank]
			if ncbi_higher_rank is None:
				continue
			ncbi_id_list.append(str(ncbi_higher_rank))

		# if self._debug and unpublished_sequence_id is not None and len(ncbi_id_list) > 0:
		# 	self._logger.debug("{id}\t{rank}".format(id=unpublished_sequence_id, rank=self._ranks[index_of_rank]))
		# 	self._mothur_cluster.cluster_list_to_stream(Counter(ncbi_id_list), sys.stderr)
		return ncbi_id_list

	# no longer used
	def get_cluster_ncbi_tax_prediction(self, cluster, unpublished_id=None):
		"""

		@param cluster:
		@type cluster: list[str|unicode]
		@param unpublished_id:
		@type unpublished_id: None|str|unicode

		@return:
		@rtype: tuple[None|str|unicode, str|unicode]
		"""
		assert isinstance(cluster, list)
		assert unpublished_id is None or isinstance(unpublished_id, basestring)
		rank_index = 1
		list_of_valid_elements = self.load_lineages(cluster)
		ncbi_id_list = []
		for rank_index in range(1, len(self._ranks)):
			ncbi_id_list = self.cluster_to_other_rank(list_of_valid_elements, rank_index, unpublished_id)
			if len(ncbi_id_list) == 0:
				continue
			dominant_id, dominant_support = max(Counter(ncbi_id_list).iteritems(), key=operator.itemgetter(1))
			if dominant_id is None:
				continue
			if float(dominant_support) / len(ncbi_id_list) >= self._minimum_support:
				break

		if len(ncbi_id_list) == 0:
			return None, ""
		dominant_id = max(Counter(ncbi_id_list).iteritems(), key=operator.itemgetter(1))[0]
		del ncbi_id_list
		return dominant_id, self._ranks[rank_index - 1]

	def has_consistent_lineage(self, iid1, iid2):
		"""

		@param iid1:
		@type iid1: str|unicode
		@param iid2:
		@type iid2: str|unicode

		@return:
		@rtype: bool
		"""
		assert isinstance(iid1, basestring)
		assert isinstance(iid2, basestring)
		set1 = set(self._iid_to_tid_lineage[iid1])
		set2 = set(self._iid_to_tid_lineage[iid2])
		if None in set1:
			set1.remove(None)
		if None in set2:
			set2.remove(None)
		if len(set2) > len(set1):
			tmp = set2
			set2 = set1
			set1 = tmp
		if len(set1) == len(set1.union(set2)):
			return True
		return False

	def load_lineages(self, cluster):
		"""

		@param cluster:
		@type cluster: list[str|unicode]

		@return:
		@rtype: list[str|unicode]
		"""
		assert isinstance(cluster, list)
		list_of_valid_elements = set()
		for iid in cluster:
			ncbi_id = self._iid_tid_map[iid]
			if not ncbi_id.isdigit():
				self._logger.warning("Bad tax id '{}': {id}".format(iid, id=ncbi_id))
				continue
			if iid not in self._iid_to_tid_lineage:
				self._iid_to_tid_lineage[iid] = self.taxonomy.get_lineage_of_legal_ranks(ncbi_id, ranks=self._ranks, default_value=None)
			list_of_valid_elements.add(iid)
		return list_of_valid_elements
