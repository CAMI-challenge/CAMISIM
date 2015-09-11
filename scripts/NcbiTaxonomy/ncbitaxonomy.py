# original from Dmitrij Turaev

__author__ = 'hofmann'
__version__ = '0.0.6'


import os
import time
from taxonomynode import TaxonomyNode
from scripts.Validator.validator import Validator


class NcbiTaxonomy(Validator):
	"""Loading NCBI from SQL dump into dictionary for fast processing"""

	# TODO: if list of ranks given, validate ranks

	_label = "NcbiTaxonomy"

	default_ordered_legal_ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
	name_to_taxids = {}
	taxid_to_parent_taxid = {}
	taxid_to_name = {}
	taxid_to_rank = {}
	taxid_old_to_taxid_new = {}
	_has_node_tree = False

	def __init__(self, taxonomy_directory="./", build_node_tree=False, verbose=True, logfile=None):
		"""
			Loading NCBI from SQL dump files into dictionary.

			@attention: building a node tree requires several gigabytes of RAM !!!

			@param taxonomy_directory: directory containing ncbi dump
			@type taxonomy_directory: str | unicode
			@param build_node_tree: Building a node tree, maybe useful if subtree is needed.
			@type build_node_tree: bool
			@param verbose: If False, messages are only written to the logfile, if given
			@type verbose: bool
			@param logfile: file stream or file path of logfile
			@type logfile: None | file | FileIO | StringIO | basestring

			@return: None
			@rtype: None
		"""
		super(NcbiTaxonomy, self).__init__(logfile=logfile, verbose=verbose)
		assert self.validate_dir(taxonomy_directory, file_names=["names.dmp", "merged.dmp", "nodes.dmp"])
		self._file_path_ncbi_names = os.path.join(taxonomy_directory, "names.dmp")
		self._file_path_ncbi_merged = os.path.join(taxonomy_directory, "merged.dmp")
		self._file_path_ncbi_nodes = os.path.join(taxonomy_directory, "nodes.dmp")
		# self._gi_taxid_file = os.path.join(taxonomy_directory, "gi_taxid_nucl.dmp")

		start = time.time()

		if len(NcbiTaxonomy.taxid_to_name) == 0:
			NcbiTaxonomy._has_node_tree = build_node_tree
			self._build_ncbi_taxonomy(build_node_tree)
			self._read_names_file()
			self._read_merged_file()
		elif not NcbiTaxonomy._has_node_tree and build_node_tree:
			self._build_ncbi_taxonomy(build_node_tree)
		else:
			self._logger.info("Using previously loaded Taxonomy")

		end = time.time()
		self._logger.info("Done ({}s)".format(round(end - start), 1))

	def get_updated_taxid(self, taxid):
		"""
			Return current taxid, in case it was merged

			@attention: taxid is not accepted as digit!!!

			@param taxid: ncbi taxonomic identifier
			@type taxid: basestring

			@return: ncbi taxonomic identifier
			@rtype: str | unicode
		"""
		assert isinstance(taxid, basestring)
		if taxid not in NcbiTaxonomy.taxid_to_rank:
			self._logger.error("Invalid taxid: '{}'".format(taxid))
			raise ValueError("Invalid taxid")
		if taxid not in NcbiTaxonomy.taxid_old_to_taxid_new:
			return taxid
		taxid_new = NcbiTaxonomy.taxid_old_to_taxid_new[taxid]
		self._logger.warning("Merged id: '{}' -> '{}'".format(taxid, taxid_new))
		return taxid_new

	def get_scientific_name(self, taxid):
		"""
			Return scientific name of ncbi taxonomic identifier

			@attention: taxid is not accepted as digit!!!

			@param taxid: ncbi taxonomic identifier
			@type taxid: basestring

			@return: ncbi scientific name
			@rtype: str | unicode
		"""
		assert isinstance(taxid, basestring)
		taxid = self.get_updated_taxid(taxid)
		if taxid in NcbiTaxonomy.taxid_to_name:
			return NcbiTaxonomy.taxid_to_name[taxid]
		self._logger.error("No name available for taxid: {}".format(taxid))
		raise ValueError("Invalid taxid")

	def get_taxids_by_scientific_name(self, scientific_name):
		"""
			Return all available taxid that fit the scientific name

			@attention: Several taxid might be a hit for one scientific name

			@param scientific_name: ncbi scientific name or synonym
			@type scientific_name: basestring

			@return: list of ncbi taxonomic identifiers
			@rtype: str | unicode
		"""
		assert isinstance(scientific_name, basestring)
		scientific_name = scientific_name.lower()
		if scientific_name in NcbiTaxonomy.name_to_taxids:
			return list(NcbiTaxonomy.name_to_taxids[scientific_name])
		self._logger.error("No taxid available for scientific_name: {}".format(scientific_name))
		raise ValueError("Invalid scientific name")

	def get_lineage_of_legal_ranks(self, taxid, ranks=None, default_value=None):
		"""
			Return lineage of a specific taxonomic identifier, filtered by a list of legal ranks

			@attention: The list of ranks determines the order of the returned taxonomic identifiers

			@param taxid: ncbi taxonomic identifier
			@type taxid: basestring
			@param ranks: List of ncbi ranks in lower case
			@type ranks: list[basestring]
			@param default_value: Value at rank indexes at which the taxid of that specific rank is undefined
			@type default_value: None | basestring

			@return: list of ncbi taxonomic identifiers
			@rtype: list[str|unicode|None]
		"""
		assert isinstance(taxid, basestring)
		taxid = self.get_updated_taxid(taxid)
		count = 0
		if ranks is None:
			ranks = NcbiTaxonomy.default_ordered_legal_ranks

		lineage = [default_value] * len(ranks)
		original_rank = self.get_rank_of_taxid(taxid)
		if original_rank is not None and original_rank in ranks:
			lineage[ranks.index(original_rank)] = taxid

		while taxid != "1" and count < 50:
			count += 1
			taxid = NcbiTaxonomy.taxid_to_parent_taxid[taxid]
			rank = NcbiTaxonomy.taxid_to_rank[taxid]
			if rank in ranks:
				lineage[ranks.index(rank)] = taxid
		if count == 50:
			self._logger.error("Bad lineage?: {}".format(lineage))
			raise Warning("Strange Error")
		return lineage

	def get_lineage(self, taxid):
		"""
			Return lineage of a specific taxonomic identifier, filtered by a list of legal ranks

			@param taxid: ncbi taxonomic identifier
			@type taxid: basestring

			@return: list of ncbi taxonomic identifiers
			@rtype: list[str|unicode]
		"""
		assert isinstance(taxid, basestring)
		taxid = self.get_updated_taxid(taxid)
		count = 0
		if NcbiTaxonomy._has_node_tree:
			return TaxonomyNode.by_name[taxid].get_lineage()

		lineage = [taxid]
		while taxid != "1" and count < 50:
			count += 1
			taxid = NcbiTaxonomy.taxid_to_parent_taxid[taxid]
			lineage.append(taxid)
		if count == 50:
			self._logger.error("Bad lineage?: {}".format(lineage))
			raise Warning("Strange Error")
		return lineage

	def get_parent_taxid_of_legal_ranks(self, taxid, ranks=None):
		"""
			Returns taxonomic identifier of the first parent of legal rank and its rank

			@param taxid: ncbi taxonomic identifier
			@type taxid: basestring
			@param ranks: List of ncbi ranks in lower case
			@type ranks: list[basestring]

			@return: tuple ncbi taxonomic identifiers and its rank
			@rtype: tuple
		"""
		assert isinstance(taxid, basestring)
		taxid = self.get_updated_taxid(taxid)
		if ranks is None:
			ranks = NcbiTaxonomy.default_ordered_legal_ranks
		if taxid not in NcbiTaxonomy.taxid_to_parent_taxid:
			self._logger.error("No parent taxid available for taxid: {}".format(taxid))
			raise ValueError("Invalid taxid")
		taxid = NcbiTaxonomy.taxid_to_parent_taxid[taxid]
		while taxid is not None and NcbiTaxonomy.taxid_to_rank[taxid] not in ranks:
			taxid = NcbiTaxonomy.taxid_to_parent_taxid[taxid]
		return taxid, NcbiTaxonomy.taxid_to_rank[taxid]

	def get_parent_taxid(self, taxid):
		"""
			Return taxonomic identifier of the parent node

			@param taxid: ncbi taxonomic identifier
			@type taxid: basestring

			@return: ncbi taxonomic identifiers
			@rtype: str | unicode
		"""
		assert isinstance(taxid, basestring)
		taxid = self.get_updated_taxid(taxid)
		if taxid in NcbiTaxonomy.taxid_to_parent_taxid:
			return NcbiTaxonomy.taxid_to_parent_taxid[taxid]
		self._logger.error("No parent taxid available for taxid: {}".format(taxid))
		raise ValueError("Invalid taxid")

	def get_rank_of_taxid(self, taxid):
		"""
			Return rank of ncbi taxonomic identifier

			@param taxid: ncbi taxonomic identifier
			@type taxid: basestring

			@return: ncbi rank of taxonomic identifiers
			@rtype: str | unicode
		"""
		assert isinstance(taxid, basestring)
		taxid = self.get_updated_taxid(taxid)
		if taxid in NcbiTaxonomy.taxid_to_rank:
			return NcbiTaxonomy.taxid_to_rank[taxid]
		self._logger.error("No rank available for taxid: {}".format(taxid))
		raise ValueError("Invalid taxid")

	def _add_nodes(self, taxid, parent_taxid='', rank='', name=''):
		"""insert nodes into taxonomy tree."""
		new_node = TaxonomyNode.by_name.get(taxid)

		if new_node is None:
			TaxonomyNode(taxid, parent_taxid, rank, name)
		# check rank
		if rank == 'no rank':
			return
		ind1 = TaxonomyNode.allranks.index(rank)
		try:
			if not TaxonomyNode.by_name[parent_taxid].rank == 'no rank':
				ind2 = TaxonomyNode.allranks.index(TaxonomyNode.by_name[parent_taxid].rank)
				assert ind1 >= ind2
				# e.g. Ovis aries platyura ('species'), Oves aries ('species')
		except KeyError:
			self._logger.debug("__add_nodes KeyError: {}".format(parent_taxid))
			pass
		# add new node to parent's all_child_nodes
		# while parent_taxid in Node.byname:
		#    Node.byname[parent_taxid].all_child_nodes.add(newnode)
		#    parent_taxid = Node.byname[parent_taxid].taxid

	@staticmethod
	def _insert_into_dict(taxid, name, my_dict):
		name = name.lower()
		assert int(taxid)
		if name not in my_dict:
			my_dict[name] = set()
		my_dict[name].add(taxid)

	def _build_ncbi_taxonomy(self, build_node_tree):
		""" parse NCBI taxonomy files."""
		self._logger.info("Building taxonomy tree...")
		if build_node_tree:
			TaxonomyNode.by_name.clear()

		# names.dmp (taxid, name, unique name, name class):
		# 521095	|	Atopobium parvulum ATCC 33793	|		|	synonym	|
		# 521095	|	Atopobium parvulum DSM 20469	|		|	scientific name	|
		# 521095	|	Atopobium parvulum str. DSM 20469	|		|	equivalent name	|
		# 521095	|	Atopobium parvulum strain DSM 20469	|		|	equivalent name	|
		# e.g. entries for "1382" in names.dmp:
		# 	1382	|	"Streptococcus parvulus" Weinberg et al. 1937	|		|	synonym	|
		# 	1382	|	Atopobium parvulum	|		|	scientific name	|
		# 	1382	|	Atopobium parvulum (Weinberg et al. 1937) Collins and Wallbanks 1993	|		|	synonym	|
		# 	1382	|	Peptostreptococcus parvulus	|		|	synonym	|
		# 	1382	|	Peptostreptococcus parvulus (Weinberg et al. 1937) Smith 1957 (Approved Lists 1980)	|	|synonym	|
		# 	1382	|	Streptococcus parvulus	|		|	synonym	|
		# 	1382	|	Streptococcus parvulus (Weinberg et al. 1937) Cato 1983	|		|	synonym	|
		# 	1382	|	not "Streptococcus parvulus" Levinthal 1928	|		|	synonym	|

		self._logger.info("Reading 'nodes' file:\t'{}'".format(self._file_path_ncbi_nodes))
		with open(self._file_path_ncbi_nodes) as file_handler:
			for line in file_handler:
				elements = [el.strip() for el in line.split('|')]
				taxid, parent_taxid, rank = elements[0:3]
				rank = rank.lower()  # should be lower-case in file, but can't be bad to doublecheck
				NcbiTaxonomy.taxid_to_parent_taxid[taxid] = parent_taxid
				NcbiTaxonomy.taxid_to_rank[taxid] = rank
				if not build_node_tree:
					continue
				assert taxid not in TaxonomyNode.by_name
				self._add_nodes(taxid, parent_taxid=parent_taxid, rank=rank)

		with open(self._file_path_ncbi_names) as file_handler:
			for line in file_handler:
				taxid, name, unique, name_class, sonst = [el.strip() for el in line.split('|')]
				self._insert_into_dict(taxid, name, NcbiTaxonomy.name_to_taxids)
				if not build_node_tree:
					continue
				try:
					my_node = TaxonomyNode.by_name[taxid]
					assert taxid == my_node.taxid
				except KeyError:
					self._logger.error("build_ncbi_taxonomy KeyError: {}".format(taxid))
					continue

				if name_class == 'scientific name':
					my_node.unique_name = unique
					my_node.scientific_name = name

				elif name_class == 'synonym':
					my_node.synonyms.append(name)
					# example: Bacteroides corrodens: Campylobacter ureolyticus (taxid 827), Eikenella corrodens (taxid 539)
					self._insert_into_dict(taxid, name, TaxonomyNode.by_synonym)

				elif name_class == 'equivalent name':
					my_node.equivalent_name.append(name)
					self._insert_into_dict(taxid, name, TaxonomyNode.by_equivalent)

				elif name_class == 'in-part' or name_class == 'includes' or \
					name_class == 'blast name' or name_class == 'genbank common name' or\
					name_class == 'misspelling' or name_class == 'authority':
					pass
		# update the taxonomy!
		TaxonomyNode.update()

	# read NCBI names file
	def _read_names_file(self):
		with open(self._file_path_ncbi_names) as fin:
			self._logger.info("Reading 'names' file:\t'{}'".format(self._file_path_ncbi_names))
			for line in fin:
				# 65      |       Herpetosiphon aurantiacus       |               |       scientific name |
				taxid, name, disambiguation, nametype, more = line.strip().split('|')
				if nametype.strip() == 'scientific name':
					NcbiTaxonomy.taxid_to_name[taxid.strip()] = name.strip()

	# read NCBI merged file
	def _read_merged_file(self):
		with open(self._file_path_ncbi_merged) as fin:
			self._logger.info("Reading 'merged' file:\t'{}'".format(self._file_path_ncbi_merged))
			for line in fin:
				# 5085       |       746128  |
				old_taxid, new_taxid, sonst = line.strip().split('|')
				NcbiTaxonomy.taxid_old_to_taxid_new[old_taxid.strip()] = new_taxid.strip()
