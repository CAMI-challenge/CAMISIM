__author__ = 'hofmann'
# original from Dmitrij Turaev

import os
from TaxonomyNode import TaxonomyNode


class NcbiTaxonomy(object):
	"""Loading NCBI from SQL dump into dictionary for fast processing, tons of RAM required"""
	default_ordered_legal_ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
	name_to_taxids = {}
	taxid_to_parent_taxid = {}
	taxid_to_name = {}
	taxid_to_rank = {}
	taxid_old_to_taxid_new = {}

	def __init__(self, taxonomy_directory='', build_node_tree=False, logger=None):
		self._logger = logger
		self._ncbi_names_file = os.path.join(taxonomy_directory, "names.dmp")
		self._ncbi_merged_file = os.path.join(taxonomy_directory, "merged.dmp")
		self._ncbi_nodes_file = os.path.join(taxonomy_directory, "nodes.dmp")
		#self._gi_taxid_file = os.path.join(taxonomy_directory, "gi_taxid_nucl.dmp")
		self.__build_ncbi_taxonomy(build_node_tree)
		self.__read_names_file()
		self.__read_merged_file()

	def get_scientific_name(self, taxid):
		if taxid in NcbiTaxonomy.taxid_old_to_taxid_new:
			taxid = NcbiTaxonomy.taxid_old_to_taxid_new[taxid]
		if taxid in NcbiTaxonomy.taxid_to_name:
			return NcbiTaxonomy.taxid_to_name[taxid]
		if self._logger:
			self._logger.error("get_scientific_name KeyError: {}".format(taxid))
		return None

	def get_taxids_by_scientific_name(self, scientific_name):
		if scientific_name in NcbiTaxonomy.name_to_taxids:
			return NcbiTaxonomy.name_to_taxids[scientific_name]
		if self._logger:
			self._logger.debug("get_taxids_by_name KeyError: {}".format(scientific_name))
		return None

	def get_lineage_of_legal_ranks(self, taxid, ranks=None, default_value=None):
		count = 0
		if taxid in NcbiTaxonomy.taxid_old_to_taxid_new:
			taxid = NcbiTaxonomy.taxid_old_to_taxid_new[taxid]
		if ranks is None:
			ranks = NcbiTaxonomy.default_ordered_legal_ranks
		lineage = [default_value] * len(ranks)
		#if self._has_node_tree:
		#	lineage = TaxonomyNode.by_name[taxid].get_lineage()
		#else:
		while taxid != "1" and count < 50:
			count += 1
			taxid = NcbiTaxonomy.taxid_to_parent_taxid[taxid]
			rank = NcbiTaxonomy.taxid_to_rank[taxid]
			if rank in ranks:
				lineage[ranks.index(rank)] = taxid
		if count == 50 and self._logger:
			self._logger.error("Bad lineage?: {}".format(lineage))
		return lineage

	def get_lineage(self, taxid):
		count = 0
		if taxid in NcbiTaxonomy.taxid_old_to_taxid_new:
			taxid = NcbiTaxonomy.taxid_old_to_taxid_new[taxid]
		if self._has_node_tree:
			lineage = TaxonomyNode.by_name[taxid].get_lineage()
		else:
			lineage = [taxid]
			while taxid != "1" and count < 50:
				count += 1
				taxid = NcbiTaxonomy.taxid_to_parent_taxid[taxid]
				lineage.append(taxid)
			if count == 50 and self._logger:
				self._logger.error("Bad lineage?: {}".format(lineage))
		return lineage

	@staticmethod
	def get_parent_taxid_of_legal_ranks(taxid, ranks=None):
		if taxid in NcbiTaxonomy.taxid_old_to_taxid_new:
			taxid = NcbiTaxonomy.taxid_old_to_taxid_new[taxid]
		if ranks is None:
			ranks = NcbiTaxonomy.default_ordered_legal_ranks
		taxid = NcbiTaxonomy.taxid_to_parent_taxid[taxid]
		while taxid is not None and NcbiTaxonomy.taxid_to_rank[taxid] not in ranks:
			taxid = NcbiTaxonomy.taxid_to_parent_taxid[taxid]
		return taxid, NcbiTaxonomy.taxid_to_rank[taxid]

	@staticmethod
	def get_parent_taxid(taxid):
		if taxid in NcbiTaxonomy.taxid_old_to_taxid_new:
			taxid = NcbiTaxonomy.taxid_old_to_taxid_new[taxid]
		try:
			return NcbiTaxonomy.taxid_to_parent_taxid[taxid]
		except KeyError:
			return None

	@staticmethod
	def get_rank_of_taxid(taxid):
		if taxid in NcbiTaxonomy.taxid_old_to_taxid_new:
			taxid = NcbiTaxonomy.taxid_old_to_taxid_new[taxid]
		if taxid in NcbiTaxonomy.taxid_to_rank:
			return NcbiTaxonomy.taxid_to_rank[taxid]
		return None

	def __add_nodes(self, taxid, parent_taxid='', rank='', name=''):
		"""insert nodes into taxonomy tree."""
		new_node = TaxonomyNode.by_name.get(taxid)

		if new_node is None:
			TaxonomyNode(taxid, parent_taxid, rank, name)
			#newnode = TaxonomyNode(taxid, parent_taxid, rank, name)
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
			if self._logger:
				self._logger.debug('__add_nodes KeyError: {}'.format(parent_taxid))
			pass
		# add new node to parent's all_child_nodes
		#while parent_taxid in Node.byname:
		#    Node.byname[parent_taxid].all_child_nodes.add(newnode)
		#    parent_taxid = Node.byname[parent_taxid].taxid

	@staticmethod
	def __insert_into_dict(taxid, name, my_dict):
		name = name.lower()
		assert int(taxid)
		if name not in my_dict:
			my_dict[name] = set()
		my_dict[name].add(taxid)

	def __build_ncbi_taxonomy(self, build_node_tree):
		""" parse NCBI taxonomy files."""
		self._has_node_tree = build_node_tree
		if self._logger:
			self._logger.info('Reading NCBI taxonomy files and building taxonomy tree...')
		if build_node_tree:
			TaxonomyNode.by_name.clear()

		# names.dmp (taxid, name, unique name, name class):
		# 521095	|	Atopobium parvulum ATCC 33793	|		|	synonym	|
		# 521095	|	Atopobium parvulum DSM 20469	|		|	scientific name	|
		# 521095	|	Atopobium parvulum str. DSM 20469	|		|	equivalent name	|
		# 521095	|	Atopobium parvulum strain DSM 20469	|		|	equivalent name	|
		# e.g. entries for "1382" in names.dmp:
		#	1382	|	"Streptococcus parvulus" Weinberg et al. 1937	|		|	synonym	|
		#	1382	|	Atopobium parvulum	|		|	scientific name	|
		#	1382	|	Atopobium parvulum (Weinberg et al. 1937) Collins and Wallbanks 1993	|		|	synonym	|
		#	1382	|	Peptostreptococcus parvulus	|		|	synonym	|
		#	1382	|	Peptostreptococcus parvulus (Weinberg et al. 1937) Smith 1957 (Approved Lists 1980)	|	|synonym	|
		#	1382	|	Streptococcus parvulus	|		|	synonym	|
		#	1382	|	Streptococcus parvulus (Weinberg et al. 1937) Cato 1983	|		|	synonym	|
		#	1382	|	not "Streptococcus parvulus" Levinthal 1928	|		|	synonym	|

		with open(self._ncbi_nodes_file) as file_handler:
			for line in file_handler:
				elements = [el.strip() for el in line.split('|')]
				taxid, parent_taxid, rank = elements[0:3]
				#import pdb; pdb.set_trace()
				rank = rank.lower()  # should be lower-case in file, but can't be bad to doublecheck
				NcbiTaxonomy.taxid_to_parent_taxid[taxid] = parent_taxid
				NcbiTaxonomy.taxid_to_rank[taxid] = rank
				if not build_node_tree:
					continue
				assert taxid not in TaxonomyNode.by_name
				self.__add_nodes(taxid, parent_taxid=parent_taxid, rank=rank)

		#file_handler = open(self._ncbi_names_file)
		with open(self._ncbi_names_file) as file_handler:
			for line in file_handler:
				taxid, name, unique, name_class, sonst = [el.strip() for el in line.split('|')]
				self.__insert_into_dict(taxid, name, NcbiTaxonomy.name_to_taxids)
				if not build_node_tree:
					continue
				try:
					my_node = TaxonomyNode.by_name[taxid]
					assert taxid == my_node.taxid
				except KeyError:
					if self._logger:
						self._logger.error('build_ncbi_taxonomy KeyError: {}'.format(taxid))
					continue

				if name_class == 'scientific name':
					my_node.unique_name = unique
					my_node.scientific_name = name

				elif name_class == 'synonym':
					my_node.synonyms.append(name)
					# example: Bacteroides corrodens: Campylobacter ureolyticus (taxid 827), Eikenella corrodens (taxid 539)
					self.__insert_into_dict(taxid, name, TaxonomyNode.by_synonym)

				elif name_class == 'equivalent name':
					my_node.equivalent_name.append(name)
					self.__insert_into_dict(taxid, name, TaxonomyNode.by_equivalent)

				elif name_class == 'in-part' or name_class == 'includes' or \
					name_class == 'blast name' or name_class == 'genbank common name' or\
					name_class == 'misspelling' or name_class == 'authority':
					pass
		# update the taxonomy!
		TaxonomyNode.update()

	# read NCBI names file
	def __read_names_file(self):
		with open(self._ncbi_names_file) as fin:
			if self._logger:
				self._logger.info("Reading NCBI namesfile: {}".format(self._ncbi_names_file))
			for line in fin:
				#65      |       Herpetosiphon aurantiacus       |               |       scientific name |
				taxid, name, disambiguation, nametype, more = line.strip().split('|')
				if nametype.strip() == 'scientific name':
					NcbiTaxonomy.taxid_to_name[taxid.strip()] = name.strip()

	# read NCBI merged file
	def __read_merged_file(self):
		with open(self._ncbi_merged_file) as fin:
			if self._logger:
				self._logger.info("Reading NCBI mergedfile (deprecated taxon IDs): {}".format(self._ncbi_merged_file))
			for line in fin:
				#5085       |       746128  |
				old_taxid, new_taxid, sonst = line.strip().split('|')
				NcbiTaxonomy.taxid_old_to_taxid_new[old_taxid.strip()] = new_taxid.strip()
