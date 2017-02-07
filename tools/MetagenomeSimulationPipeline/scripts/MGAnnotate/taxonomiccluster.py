__author__ = 'hofmann'

from scripts.Validator.validator import Validator
from scripts.MGAnnotate.mothurcluster import MothurCluster
from scripts.NcbiTaxonomy.ncbitaxonomy import NcbiTaxonomy


class TaxonomicCluster(Validator):
	"""Reading and writing a meta table"""

	_label = "TaxonomicCluster"

	def __init__(self, mothur_cluster, taxonomy, iid_tid_map, set_reference_genome_ncbi, minimum_support=.9, logfile=None, verbose=True, debug=False):
		"""
		Constructor

		@param mothur_cluster: A handle to a MothurCluster object
		@type mothur_cluster: MothurCluster
		@param taxonomy: Handle to NcbiTaxonomy
		@type taxonomy: NcbiTaxonomy
		@param iid_tid_map: A map from internal id to the taxonomic id
		@type iid_tid_map: dict[str|unicode, str|unicode]
		@param set_reference_genome_ncbi: NCBI taxonomic ids for reference genomes
		@type set_reference_genome_ncbi: set[str|unicode]]
		@param minimum_support: Minimum percentage of elements that must support a specific taxid
		@type minimum_support: float
		@param logfile: File handler or file path to a log file
		@type logfile: file | FileIO | StringIO | basestring
		@param verbose: Not verbose means that only warnings and errors will be past to stream
		@type verbose: bool
		@param debug: Display debug messages
		@type debug: bool
		"""
		assert isinstance(mothur_cluster, MothurCluster)
		assert isinstance(taxonomy, NcbiTaxonomy)
		assert isinstance(iid_tid_map, dict)
		assert isinstance(set_reference_genome_ncbi, set)
		assert isinstance(minimum_support, float)
		super(TaxonomicCluster, self).__init__(logfile=logfile, verbose=verbose, debug=debug)
		self._mothur_cluster = mothur_cluster
		self.taxonomy = taxonomy
		self._ranks = ['strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']  # , 'root'
		self._iid_tid_map = iid_tid_map
		self._iid_to_tid_lineage = {}
		self._set_reference_genome_ncbi = set_reference_genome_ncbi
		self._minimum_support = minimum_support

	def predict_tax_id_of(self, cluster, lowest_predicted_novelty=None):
		"""
		Get the predicted taxid of a cluster

		@param cluster: List of internal ids within a otu
		@type cluster: list[str|unicode]
		@param lowest_predicted_novelty: Lowest predicted novelty so far
		@type lowest_predicted_novelty: None|dict[str|unicode, str|unicode]

		@return: taxonomic_id, novelty, amount_of_supporting_elements
		@rtype: tuple[None|str|unicode, str|unicode, int|long]
		"""
		assert isinstance(cluster, list)
		assert isinstance(lowest_predicted_novelty, dict)

		list_of_valid_iid, set_of_ncbi_taxid = self.load_lineages(cluster)
		root = {"count": 0, "c": {}, 'p': None}

		# Build tree based from lineage of known references
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

		# Get node with highest number of strains
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
			if 'r' not in node:
				if lowest_predicted_novelty["novelty"] == '':
					lowest_predicted_novelty["novelty"] = 'superkingdom'
				return "1", "", 0
				# continue  # root node
			if float(node["count"]) / total_count[node['r']] >= self._minimum_support:
				novelty = self._ranks[node['r'] - 1]  # novelty of all references
				# if previous novelty lower
				if lowest_predicted_novelty["novelty"] in self._ranks and self._ranks.index(lowest_predicted_novelty["novelty"]) + 1 <= node['r']:
					novelty = lowest_predicted_novelty["novelty"]
				if self.is_near_genome_reference(list_of_valid_iid):  # or float(threshold) == 0.1
					if lowest_predicted_novelty["novelty"] == '' or self._ranks.index(lowest_predicted_novelty["novelty"]) + 1 > node['r']:
						lowest_predicted_novelty["novelty"] = novelty  # novelty of genome references
				return node["id"], novelty, node["count"]
		return None, "", 0

	def cluster_to_ncbi_of_a_rank(self, list_of_iid, index_of_rank, query_gid=None):
		"""
		Get a list of each taxonomic classifications of the elements in a cluster

		@param list_of_iid: A list of internal ids (usualy those within a otu cluster)
		@type list_of_iid: list[str|unicode] | set[str|unicode]
		@param index_of_rank: Index referencing to a rank in self._ranks
		@type index_of_rank: int | long
		@param query_gid: genome id this request is made for
		@type query_gid: None | str|unicode

		@return: List of each taxonomic classifications of the elements in a cluster
		@rtype: list[str|unicode]
		"""
		assert isinstance(list_of_iid, list)
		assert isinstance(index_of_rank, (int, long))
		assert query_gid is None or isinstance(query_gid, basestring)
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

	def has_consistent_lineage(self, iid1, iid2):
		"""
		Compares the lineages of two ids for inconsistencies

		@param iid1: Internal id
		@type iid1: str|unicode
		@param iid2: Internal id
		@type iid2: str|unicode

		@return: True if consistent
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

	def is_near_genome_reference(self, set_iids):
		"""
		Compares the lineages of two ids for inconsistencies

		@param set_iids: Set of internal id
		@type set_iids: set[str|unicode]

		@return: True if consistent
		@rtype: bool
		"""
		assert isinstance(set_iids, set)
		if None in set_iids:
			set_iids.remove(None)
		if None in self._set_reference_genome_ncbi:
			self._set_reference_genome_ncbi.remove(None)
		set2 = self._set_reference_genome_ncbi
		for iid in set_iids:
			set1 = set(self._iid_to_tid_lineage[iid])
			if not set1.isdisjoint(set2):
				return True
		return False

	def load_lineages(self, cluster):
		"""
		Get list of valid internal ids and save lineages of those

		@param cluster: List of internal ids from an otu cluster
		@type cluster: list[str|unicode]

		@return: List of internal ids that correspond to a valid taxid
		@rtype: tuple[set[str|unicode], set[str|unicode]]
		"""
		assert isinstance(cluster, list)
		list_of_valid_elements = set()
		set_of_ncbi_taxid = set()
		for iid in cluster:
			ncbi_id = self._iid_tid_map[iid]
			if not ncbi_id.isdigit():
				self._logger.warning("Bad tax id '{}': {id}".format(iid, id=ncbi_id))
				continue
			if iid not in self._iid_to_tid_lineage:
				self._iid_to_tid_lineage[iid] = self.taxonomy.get_lineage_of_legal_ranks(ncbi_id, ranks=self._ranks, default_value=None)
			set_of_ncbi_taxid.add(ncbi_id)
			list_of_valid_elements.add(iid)
		return list_of_valid_elements, set_of_ncbi_taxid
