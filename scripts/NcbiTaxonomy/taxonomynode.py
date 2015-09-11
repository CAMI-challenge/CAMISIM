__author__ = 'hofmann'
# original from Dmitrij Turaev

import sys


class TaxonomyNode(object):
	"""class to describe NCBI taxonomy tree.
	
	This class has to know: parent, children, taxid, rank, gi(?), 
	scientific name, synonyms(?), equivalent name(?)
	"""
	allranks = [
		'root', 'superkingdom', 'kingdom', 'subkingdom', 'superphylum', 'phylum',
		'subphylum', 'superclass', 'class', 'subclass', 'infraclass', 'superorder', 'order', 'suborder',
		'infraorder', 'parvorder', 'superfamily', 'family', 'subfamily', 'tribe', 'subtribe', 'genus',
		'subgenus', 'species group', 'species subgroup', 'species', 'subspecies', 'varietas', 'forma']
	by_name = {}
	by_rank = {}
	by_synonym = {}
	by_equivalent = {}
	by_scientific_name = {}
	ambiguous = {}
	inactive_top_nodes = []
	user_provided = []
	# counter = 0
	
	def __set_scientific_name(self, name):
		# ~ print name
		self._scientific_name = name
		key = name.lower()
		if not key:
			return
		if key not in TaxonomyNode.by_scientific_name:
			TaxonomyNode.by_scientific_name[key] = self.taxid
			return
		# problem: two identical scientific names; examples:
		# Bironella / Bironella <subgenus>
		# Bacillus <stick insect> / Bacillus <bacterium>
		# #newkey = Node.bysciname[key].unique_name.lower()
		# #if newkey:
		# #	Node.bysciname[newkey] = Node.bysciname[key]
		# #	del Node.bysciname[key]
		# #newnewkey = self.unique_name.lower()
		# assert that there was a solution for the problem
		# #assert newkey or newnewkey
		# #if newnewkey:
		# #	assert newkey != newnewkey
		# #	Node.bysciname[newnewkey] = self
		# #else:
		# #	Node.bysciname[key] = self
		# ## check if the unique name works correctly, e.g. for carma !!!!!!
		# dictionary for ambiguous cases
		# #Node.ambiguous[key] = [k for k in [key, newkey, newnewkey] if Node.bysciname.get(k)]
		if isinstance(TaxonomyNode.by_scientific_name[key], list):
			TaxonomyNode.by_scientific_name[key].append(self.taxid)
		else:
			assert int(TaxonomyNode.by_scientific_name[key])
			TaxonomyNode.by_scientific_name[key] = [TaxonomyNode.by_scientific_name[key], self.taxid]

	def __get_scientific_name(self):
		return self._scientific_name
	scientific_name = property(__get_scientific_name, __set_scientific_name)

	def __init__(self, taxid, parent_taxid, rank='', name='', unique_name=''):
		""" constructor. """
		self.taxid = taxid
		self.rank = rank
		# ~ self.sciname = name
		self.scientific_name = name
		self.synonyms = []
		self.equivalent_name = []
		self.gi = ''
		self.parent = None
		self.parent_taxid = parent_taxid
		self.children = set()
		self.leafs = set()
		self.lineage = []
		# ~ self.all_child_nodes = []
		self.all_child_nodes = set()  # higher performance?
		self.node_active = True
		self.unique_name = unique_name

		assert taxid not in TaxonomyNode.by_name
		TaxonomyNode.by_name[taxid] = self
		TaxonomyNode.by_rank.setdefault(rank, []).append(self)
		if taxid.startswith('[u]'):
			TaxonomyNode.user_provided.append(self)

	def update_node(self):
		""" update nodes with some information which depends on
		their neighborhood and should be done after all nodes are created. """
		# replace parent taxids by the parent objects themselves
		self.parent = TaxonomyNode.by_name[self.parent_taxid]
		# tell parents that they have children
		# ~ import pdb; pdb.set_trace()
		# ~ if self not in self.parent.children:
		# ~ Node.counter += 1
		# ~ if Node.counter % 1000 == 0: print self.counter
		self.parent.children.add(self)
		# update the "bysciname" dictionary using "parent:child" as key
		# ~ assert self.parent.sciname
		# ~ key = ':'.join([self.parent.sciname.lower(), self.sciname.lower()])
		# ~ key = self.sciname.lower()
		# ~ self.bysciname[key] = self
		# ~ print key; import pdb; pdb.set_trace()

	def get_leafs(self, leafs=None):
		""" get the terminal leafs for a particular node."""
		if leafs is None:
			if self.leafs:
				return	  # leafs are already known, nothing to do
			leafs = self.leafs  # starting point
		if self.children:
			for child in self.children:
				child.get_leafs(leafs)
		else:
			leafs.add(self)

	def get_child_nodes(self, child_nodes=None):
		""" get all child nodes (not only the direct children) for a particular node. """
		if child_nodes is None:  # starting point
			# ~ del self.all_child_nodes[:]
			# ~ self.all_child_nodes.clear()
			if self.all_child_nodes:
				return	 # nothing to do
			child_nodes = self.all_child_nodes
		else:
			# #~ if self not in child_nodes:
			# ~ child_nodes.append(self)
			child_nodes.add(self.taxid)
		if self.children:
			for child in self.children:
				if child.taxid == '1':
					continue  # don't get caught in infinite loop
				child.get_child_nodes(child_nodes)

	def get_all_descendant_taxids(self):
		""" return taxids of all descendants """
		self.get_child_nodes()
		# return set(node.taxid for node in self.all_child_nodes)
		return self.all_child_nodes

	def get_lineage(self, lineage=None):
		""" get taxonomy lineage up to the root for this node """
		if lineage is None:
			lineage = []
		lineage.append(self.taxid)
		if self.taxid == '1':
			return list(reversed(lineage))
		else:
			return self.parent.get_lineage(lineage)

	@staticmethod
	def update():
		for key, my_node in TaxonomyNode.by_name.items():
			my_node.update_node()

	@staticmethod
	def active_parent_nodes_consistency(oCurrNode):
		""" active_parent_nodes_consistency(oCurrNode)
			oCurrNode ... class Node
			------------------------------------------------------------------------
			Checks, if the parent node of an ncbi taxonomy node has at least one active
			child node. Every parent node of the ncbi taxonomy node that has no active
			child nodes will be inactivated.
			Returns the last node that has a parent with at least one active child.
		"""
		try:
			bCheckActiveChildNodes = False
			oChildrenNodes = oCurrNode.parent.children
			for oChildNode in oChildrenNodes:
				bCheckActiveChildNodes = bCheckActiveChildNodes or oChildNode.node_active

			if not bCheckActiveChildNodes:
				oCurrNode.parent.node_active = False
				oCurrNode = TaxonomyNode.active_parent_nodes_consistency(oCurrNode.parent)

		except Exception as e:
			print >> sys.stderr, str(e)
			# logging.error(str(e))

		return oCurrNode

	@staticmethod
	def inactivate_branch(taxid):
		""" inactivate_branch(taxid)
			taxid ... ncbi taxonomy id
			------------------------------------------------------------------------
			Inactivates an ncbi taxonomy node with id = taxid including all 
			children by changing node attribute "node_active" to False. 
			Furthermore inactivates all successive parent nodes, which 
			have no active child nodes. Adds the top inactivated node within 
			the inactivated branch to class dictionary Node.inactive_top_nodes 
			to permit a reactivation of this branch later on.
		"""
		try:
			oCurrNode = TaxonomyNode.by_name.get(str(taxid))
			if not oCurrNode.all_child_nodes:
				oCurrNode.get_child_nodes()
			oCurrNode.node_active = False
			for oChildNode in oCurrNode.all_child_nodes:
				oChildNode.node_active = False

			oCurrNode = TaxonomyNode.active_parent_nodes_consistency(oCurrNode)

			TaxonomyNode.inactive_top_nodes.append(oCurrNode)

		except Exception as e:
			print >> sys.stderr, str(e)
			# logging.error(str(e))

	@staticmethod
	def activate_branch(taxid):
		""" "activate" a branch of the taxonomy tree. """
		oCurrNode = TaxonomyNode.by_name.get(str(taxid))
		TaxonomyNode.inactive_top_nodes.remove(oCurrNode)
		oCurrNode.get_child_nodes()
		oCurrNode.node_active = True
		for oChildNode in oCurrNode.all_child_nodes:
			oChildNode.node_active = True

	@staticmethod
	def find_parent_by_rank(oCurrNode, vFindParentRank, origNode=None, takeLowerIfMissing=False):
		""" find parent node for given taxonomic rank. examples:
		root (no rank); viruses (superkingdom); unclassified phages (no rank); Streptococcus phage YMC-2011 (species)
		cellular organisms; Bacteria; Cyanobacteria (phylum); Chroococcales (order); Cyanothece; Cyanothece sp. (rank "class" missing)
		"""
		if oCurrNode.sciname == 'root':  # last resort
			return False
		vFindParentRank = vFindParentRank.lower().strip()
		if not origNode:
			origNode = oCurrNode   # remember starting point
		target_ind = TaxonomyNode.allranks.index(vFindParentRank)
		# if origNode has rank 'no rank', go on with its parent
		try:
			my_ind = TaxonomyNode.allranks.index(oCurrNode.rank)
		except ValueError:
			assert oCurrNode.rank == 'no rank'
			return TaxonomyNode.find_parent_by_rank(oCurrNode.parent, vFindParentRank, origNode)
		# if parent node has rank 'no rank', go on with its parent
		parent_ind = None
		while not parent_ind:
			try:
				parent_ind = TaxonomyNode.allranks.index(oCurrNode.parent.rank)
			except ValueError:
				assert oCurrNode.parent.rank == 'no rank'
				oCurrNode = oCurrNode.parent
		# check the rank positions
		if my_ind == target_ind:		# it's me
			oReturnParentNode = oCurrNode
		elif parent_ind == target_ind:  # it's my parent
			oReturnParentNode = oCurrNode.parent
		elif parent_ind < target_ind:   # missed it
			if takeLowerIfMissing:
				oReturnParentNode = oCurrNode
			else:
				oReturnParentNode = False
		else:
			oReturnParentNode = TaxonomyNode.find_parent_by_rank(oCurrNode.parent, vFindParentRank, origNode)
		return oReturnParentNode