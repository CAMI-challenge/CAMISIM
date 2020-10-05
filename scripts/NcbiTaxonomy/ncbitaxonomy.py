# original from Dmitrij Turaev

__author__ = 'Peter Hofmann'
__version__ = '0.1.5'


import os
import time
import fnmatch
import tempfile
from .taxonomynode import TaxonomyNode
from scripts.Validator.validator import Validator
from scripts.Archive.archive import Archive


class NcbiTaxonomy(Validator):
    """
    Loading NCBI from SQL dump into dictionary for fast processing

    @type name_to_taxids: dict[str, set[str]]
    @type taxid_to_parent_taxid: dict[str, str]
    @type taxid_to_name: dict[str, str]
    @type taxid_to_rank: dict[str, str]
    @type taxid_old_to_taxid_new: dict[str, str]
    @type _has_node_tree: bool
    """

    # TODO: if list of ranks given, validate ranks

    default_ordered_legal_ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
    name_to_taxids = {}
    taxid_to_parent_taxid = {}
    taxid_to_name = {}
    taxid_to_rank = {}
    taxid_old_to_taxid_new = {}
    _has_node_tree = False

    def __init__(self, taxonomy_path="./", temporary_directory=None, build_node_tree=False, verbose=True, logfile=None):
        """
        Loading NCBI from SQL dump files into dictionary.

        @attention: building a node tree requires several gigabytes of RAM !!!

        @param taxonomy_path: directory containing ncbi dump
        @type taxonomy_path: str | unicode
        @param build_node_tree: Building a node tree, maybe useful if subtree is needed.
        @type build_node_tree: bool
        @param verbose: If False, messages are only written to the logfile, if given
        @type verbose: bool
        @param logfile: file stream or file path of logfile
        @type logfile: None | file | FileIO | StringIO | str

        @return: None
        @rtype: None
        """
        super(NcbiTaxonomy, self).__init__(label="NcbiTaxonomy", logfile=logfile, verbose=verbose)
        assert isinstance(taxonomy_path, str), "Invalid taxonomy directory."
        assert temporary_directory is None or self.validate_dir(temporary_directory)
        assert isinstance(build_node_tree, bool)

        assert os.path.exists(taxonomy_path), "Invalid taxonomy directory."
        self._tmp_dir = None

        if not self.validate_dir(taxonomy_path, silent=True):
            archive = Archive()
            assert archive.is_archive(taxonomy_path), "Can not read taxonomy. Unknown archive."
            if temporary_directory is None:
                self._tmp_dir = tempfile.mkdtemp()
            else:
                self._tmp_dir = tempfile.mkdtemp(dir=temporary_directory)
            archive.extract_all(taxonomy_path, self._tmp_dir)
            folder_name = os.listdir(self._tmp_dir)[0]
            taxonomy_path = os.path.join(self._tmp_dir, folder_name)

        assert self.validate_dir(taxonomy_path, file_names=["names.dmp", "merged.dmp", "nodes.dmp"])

        taxonomy_path = self.get_full_path(taxonomy_path)
        self._file_path_ncbi_names = os.path.join(taxonomy_path, "names.dmp")
        self._file_path_ncbi_merged = os.path.join(taxonomy_path, "merged.dmp")
        self._file_path_ncbi_nodes = os.path.join(taxonomy_path, "nodes.dmp")
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

    def __exit__(self, type, value, traceback):
        super(NcbiTaxonomy, self).__exit__(type, value, traceback)
        if self.validate_dir(self._tmp_dir, silent=True):
            import shutil
            shutil.rmtree(self._tmp_dir)
        self.tmp_dir = None

    def __del__(self):
        super(NcbiTaxonomy, self).__del__()
        if self.validate_dir(self._tmp_dir, silent=True):
            import shutil
            shutil.rmtree(self._tmp_dir)
        self.tmp_dir = None

    def has_taxid(self, taxid):
        """
        Return current taxid, in case it was merged

        @attention: taxid is not accepted as digit!!!

        @param taxid: ncbi taxonomic identifier
        @type taxid: str

        @return: True if taxid exists in taxdump
        @rtype: bool
        """
        assert isinstance(taxid, str)
        if taxid in NcbiTaxonomy.taxid_to_rank:
            return True
        return False

    def get_updated_taxid(self, taxid):
        """
        Return current taxid, in case it was merged

        @attention: taxid is not accepted as digit!!!

        @param taxid: ncbi taxonomic identifier
        @type taxid: str

        @return: ncbi taxonomic identifier
        @rtype: str | unicode
        """
        assert isinstance(taxid, str)
        if taxid in NcbiTaxonomy.taxid_to_rank:
            return taxid
        if taxid not in NcbiTaxonomy.taxid_old_to_taxid_new:
            self._logger.error("Invalid taxid: '{}'".format(taxid))
            raise ValueError("Invalid taxid")

        taxid_new = NcbiTaxonomy.taxid_old_to_taxid_new[taxid]
        self._logger.warning("Merged id: '{}' -> '{}'".format(taxid, taxid_new))
        return taxid_new

    def get_scientific_name(self, taxid):
        """
        Return scientific name of ncbi taxonomic identifier

        @attention: taxid is not accepted as digit!!!

        @param taxid: ncbi taxonomic identifier
        @type taxid: str

        @return: ncbi scientific name
        @rtype: str | unicode
        """
        assert isinstance(taxid, str)
        taxid = self.get_updated_taxid(taxid)
        if taxid in NcbiTaxonomy.taxid_to_name:
            return NcbiTaxonomy.taxid_to_name[taxid]
        self._logger.error("No name available for taxid: {}".format(taxid))
        raise ValueError("Invalid taxid")

    def get_taxids_by_scientific_name(self, scientific_name, silent=False):
        """
        Return all available taxid that fit the scientific name

        @attention: Several taxid might be a hit for one scientific name

        @param scientific_name: ncbi scientific name or synonym
        @type scientific_name: str

        @return: list of ncbi taxonomic identifiers
        @rtype: set[str | unicode] | None
        """
        assert isinstance(scientific_name, str)
        scientific_name = scientific_name.lower()
        if scientific_name in NcbiTaxonomy.name_to_taxids:
            return set(NcbiTaxonomy.name_to_taxids[scientific_name])
        if not silent:
            self._logger.error("No taxid available for scientific_name: {}".format(scientific_name))
            raise ValueError("Invalid scientific name")
        return None

    def get_taxids_by_scientific_name_wildcard(self, scientific_name):
        """
        Return all available taxid that fit the scientific name

        @attention: Several taxid might be a hit for one scientific name

        @param scientific_name: ncbi scientific name or synonym
        @type scientific_name: str

        @return: set of ncbi taxonomic identifiers
        @rtype: set[str | unicode] | None
        """
        assert isinstance(scientific_name, str)
        scientific_name = scientific_name.lower()
        matches = fnmatch.filter(self.name_to_taxids.keys(), scientific_name)
        set_of_tax_id = set()
        for match in matches:
            set_of_tax_id.update(set(self.name_to_taxids[match]))
        if len(set_of_tax_id) > 1:
            self._logger.warning(
                "Several matches '{}' found for scientific_name: '{}'".format(", ".join(matches), scientific_name))
            return set_of_tax_id
        elif len(set_of_tax_id) == 0:
            return None
        return set_of_tax_id

    def get_lineage_of_legal_ranks(self, taxid, ranks=None, default_value=None, as_name=False, inherit_rank=False):
        """
        Return lineage of a specific taxonomic identifier, filtered by a list of legal ranks

        @attention: The list of ranks determines the order of the returned taxonomic identifiers

        @param taxid: ncbi taxonomic identifier
        @type taxid: str
        @param ranks: List of ncbi ranks in lower case
        @type ranks: list[str]
        @param default_value: Value at rank indexes at which the taxid of that specific rank is undefined
        @type default_value: None | str
        @param as_name: return scientific name if true, not taxonomic id
        @type as_name: bool
        @param inherit_rank: name unnamed rank names by known ones, species -> root
        @type inherit_rank: bool

        @return: list of ncbi taxonomic identifiers
        @rtype: list[str|unicode|None]
        """
        assert isinstance(taxid, str)
        taxid = self.get_updated_taxid(taxid)
        if ranks is None:
            ranks = NcbiTaxonomy.default_ordered_legal_ranks

        lineage = [default_value] * len(ranks)
        original_rank = self.get_rank_of_taxid(taxid)
        if original_rank is not None and original_rank in ranks:
            if as_name:
                lineage[ranks.index(original_rank)] = NcbiTaxonomy.taxid_to_name[taxid]
            else:
                lineage[ranks.index(original_rank)] = taxid
        try:
            rank_counter = ranks.index(NcbiTaxonomy.taxid_to_rank[taxid]) # starting at rank of original tax id
        except ValueError: # rank is not in ranks
            rank_counter = ranks.index(ranks[-1]) # choose lowest rank then
        while taxid != "1":
            taxid = NcbiTaxonomy.taxid_to_parent_taxid[taxid]
            rank = NcbiTaxonomy.taxid_to_rank[taxid]
            if rank in ranks:
                current_rank_counter = ranks.index(NcbiTaxonomy.taxid_to_rank[taxid])
                rank_difference = rank_counter - current_rank_counter
                if rank_difference > 1:
                    for i in range(current_rank_counter, rank_counter - 1):
                        lineage[i] = "" # add empty name to list if name is missing in the taxonomy
                rank_counter = current_rank_counter
                if as_name:
                    lineage[ranks.index(rank)] = NcbiTaxonomy.taxid_to_name[taxid]
                else:
                    lineage[ranks.index(rank)] = taxid

        # todo: sort ranks
        if inherit_rank:
            rank_previous = default_value
            tmp_list = enumerate(lineage)
            if self.default_ordered_legal_ranks.index(ranks[0]) < self.default_ordered_legal_ranks.index(ranks[-1]):
                tmp_list = reversed(list(enumerate(lineage)))
            for index, value in tmp_list:
                if value == default_value:
                    lineage[index] = rank_previous
                else:
                    rank_previous = value
        return lineage

    def get_lineage(self, taxid):
        """
        Return lineage of a specific taxonomic identifier, filtered by a list of legal ranks

        @param taxid: ncbi taxonomic identifier
        @type taxid: str

        @return: list of ncbi taxonomic identifiers
        @rtype: list[str|unicode]
        """
        assert isinstance(taxid, str)
        taxid = self.get_updated_taxid(taxid)
        if NcbiTaxonomy._has_node_tree:
            return TaxonomyNode.by_name[taxid].get_lineage()

        lineage = [taxid]
        while taxid != "1":
            taxid = NcbiTaxonomy.taxid_to_parent_taxid[taxid]
            lineage.append(taxid)
        return lineage

    def get_parent_taxid_of_legal_ranks(self, taxid, ranks=None):
        """
        Returns taxonomic identifier of the first parent of legal rank and its rank

        @param taxid: ncbi taxonomic identifier
        @type taxid: str
        @param ranks: List of ncbi ranks in lower case
        @type ranks: list[str]

        @return: tuple ncbi taxonomic identifiers and its rank
        @rtype: tuple
        """
        assert isinstance(taxid, str)
        taxid = self.get_updated_taxid(taxid)
        if ranks is None:
            ranks = NcbiTaxonomy.default_ordered_legal_ranks
        if taxid not in NcbiTaxonomy.taxid_to_parent_taxid:
            self._logger.error("No parent taxid available for taxid: {}".format(taxid))
            raise ValueError("Invalid taxid")
        taxid = NcbiTaxonomy.taxid_to_parent_taxid[taxid]
        while taxid is not None and taxid != "1" and NcbiTaxonomy.taxid_to_rank[taxid] not in ranks:
            taxid = NcbiTaxonomy.taxid_to_parent_taxid[taxid]
        if NcbiTaxonomy.taxid_to_rank[taxid] not in ranks:
            return None, None
        return taxid, NcbiTaxonomy.taxid_to_rank[taxid]

    def get_parent_taxid(self, taxid):
        """
        Return taxonomic identifier of the parent node

        @param taxid: ncbi taxonomic identifier
        @type taxid: str

        @return: ncbi taxonomic identifiers
        @rtype: str | unicode
        """
        assert isinstance(taxid, str)
        taxid = self.get_updated_taxid(taxid)
        if taxid in NcbiTaxonomy.taxid_to_parent_taxid:
            return NcbiTaxonomy.taxid_to_parent_taxid[taxid]
        self._logger.error("No parent taxid available for taxid: {}".format(taxid))
        raise ValueError("Invalid taxid")

    def get_rank_of_taxid(self, taxid):
        """
        Return rank of ncbi taxonomic identifier

        @param taxid: ncbi taxonomic identifier
        @type taxid: str

        @return: ncbi rank of taxonomic identifiers
        @rtype: str | unicode
        """
        assert isinstance(taxid, str)
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
        # 521095    |    Atopobium parvulum ATCC 33793    |        |    synonym    |
        # 521095    |    Atopobium parvulum DSM 20469    |        |    scientific name    |
        # 521095    |    Atopobium parvulum str. DSM 20469    |        |    equivalent name    |
        # 521095    |    Atopobium parvulum strain DSM 20469    |        |    equivalent name    |
        # e.g. entries for "1382" in names.dmp:
        #     1382    |    "Streptococcus parvulus" Weinberg et al. 1937    |        |    synonym    |
        #     1382    |    Atopobium parvulum    |        |    scientific name    |
        #     1382    |    Atopobium parvulum (Weinberg et al. 1937) Collins and Wallbanks 1993    |        |    synonym    |
        #     1382    |    Peptostreptococcus parvulus    |        |    synonym    |
        #     1382    |    Peptostreptococcus parvulus (Weinberg et al. 1937) Smith 1957 (Approved Lists 1980)    |    |synonym    |
        #     1382    |    Streptococcus parvulus    |        |    synonym    |
        #     1382    |    Streptococcus parvulus (Weinberg et al. 1937) Cato 1983    |        |    synonym    |
        #     1382    |    not "Streptococcus parvulus" Levinthal 1928    |        |    synonym    |

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

    # ###############
    # newick
    # ###############

    # ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
    def _add_lineage_to_tree(self, root, lineage):
        """
        Adding a lineage to a dictionary based tree

        @param root: Root node
        @type root: dict[str,dict]
        @param lineage: A lineage
        @type lineage: list[str]

        @rtype: None
        """
        node = root
        for taxid in lineage:
            if taxid is None:
                continue
            if taxid not in node:
                node[taxid] = {}
            node = node[taxid]

    # (A,B,(C,D)E)F;
    def _node_to_newick(self, node, node_name):
        """
        Create a newick sting based on a tree

        @param node:
        @type node: dict[str,dict]
        @param node_name:
        @type node_name: str

        @return: newick string
        @rtype: str
        """
        if len(node) == 0:
            return node_name
        child_nodes = []
        for name in sorted(node.keys()):
            child_nodes.append(self._node_to_newick(node[name], name))
        return "({}){}".format(",".join(child_nodes), node_name)

    def to_newick(self, stream, ranks=None):
        """
        Export taxonomy as newick formated string.

        @attention: Always rooted with id '1'

        @param stream: Output stream
        @type stream: file | FileIO | StringIO
        @param ranks: List of legal ranks
        @type ranks: list[str]

        @rtype: None
        """
        # build tree
        if ranks is None:
            ranks = self.default_ordered_legal_ranks
        root = {}
        for taxid in sorted(self.taxid_to_rank.keys()):
            lineage = self.get_lineage_of_legal_ranks(taxid, ranks=ranks)
            self._add_lineage_to_tree(root, lineage)

        # build newick string
        stream.write("{};\n".format(self._node_to_newick(root, '1')))

    def to_map(self, stream):
        """
        Exporting a map of all taxonomic ids to its respective taxonomic name.

        @param stream: Output stream
        @type stream: file | FileIO | StringIO

        @rtype: None
        """
        # for taxid in set_of_strains:
        for taxid, name in self.taxid_to_name.items():
            stream.write("{}\t{}\n".format(taxid, name))

    def lca(self, tax_id1, tax_id2):
        """

        @param tax_id1: ncbi taxonomic identifier
        @type tax_id1: str
        @param tax_id2: ncbi taxonomic identifier
        @type tax_id2: str

        @return: ncbi taxonomic identifier
        @rtype: str
        """
        ranks = self.default_ordered_legal_ranks
        ranks.reverse()
        consistent_lineage = True
        lineage1 = self.get_lineage_of_legal_ranks(tax_id1, ranks=ranks)
        lineage2 = self.get_lineage_of_legal_ranks(tax_id2, ranks=ranks)
        for index, value in enumerate(lineage1):
            if value is None:
                continue
            if lineage2[index] is None:
                continue
            if value != lineage2[index]:
                consistent_lineage = False
                continue
            if not consistent_lineage:
                self._logger.info("Inconsitent lineage: {} vs {}".format(tax_id1, tax_id2))
            return value
        if not consistent_lineage:
            self._logger.info("Inconsitent lineage: {} vs {}".format(tax_id1, tax_id2))
        return "1"
