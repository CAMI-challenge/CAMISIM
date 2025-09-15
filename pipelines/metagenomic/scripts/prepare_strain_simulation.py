import random
import os
import argparse
import io
import sys
from collections import Counter
import numpy.random as np_random
#from Bio import Phylo

class PrepareStrainSimulation():

    directory_template_filenames = ["simujobparams.pm", "template.tree"]
    column_name_gid = "genome_ID"
    column_name_ncbi="NCBI_ID"
    column_name_source="source"
    column_name_novelty_category="novelty_category"
    column_name_otu="OTU"
    list_of_column_names = [column_name_gid, column_name_ncbi, column_name_novelty_category, column_name_otu]
    separator='\t'
    filename_prefix="simulated_"
    keep_original=True
    max_processors=1
    cats = 0
    draw = 0
    per_cat = 0
    rest = 0
    per_otu = 0
    #tmp_dir=self._tmp_dir,
    #logfile=self._logfile, verbose=self._verbose, debug=self._debug,

    def get_filenames_strains(self, file_path_template_newick_tree):
        """
        Get list of file names of simulated genomes by reading template newick tree

        @attention: 'ancestor' is assumed to be part of tree as original sequence and will not be included

        @param file_path_template_newick_tree: File path to newick file
        @type file_path_template_newick_tree: str | unicode

        @return: list of file names of simulated genomes
        @rtype: list[str|unicode]
        """
        #assert self.validate_file(file_path_template_newick_tree)
        list_of_filenames_strains = []
        tree = Phylo.read(file_path_template_newick_tree, 'newick')
        for leaf in tree.get_terminals():
            prefix = leaf.name
            if prefix.lower() == "ancestor":
                continue
            list_of_filenames_strains.append("{prefix}.fasta".format(prefix=prefix))
        return list_of_filenames_strains

    #def init(directory_template=None, seed=None):
    def init(self, seed=None):
        """
            Initialize with seed

            @attention:

            @param executable_sim: filepath to 'simujobrun.pl', default is 'StrainSimulationWrapper/sgEvolver/simujobrun.pl'
            @type executable_sim: str | unicode
            @param directory_template: directory with 'simujobparams.pm', 'template.tree'
            @type directory_template: str | unicode
            @param column_name_gid: Name of genomic ID column
            @type column_name_gid: str | unicode
            @param column_name_ncbi: Name of NCBI taxid column
            @type column_name_ncbi: str | unicode
            @param column_name_source: Name of genomic ID column
            @type column_name_source: str | unicode
            @param separator: separator used in metadata file
            @type separator: str | unicode
            @param filename_prefix: filename prefix of simulated genomes
            @type filename_prefix: str | unicode
            @param keep_original: If true, original genomes will be kept, else only simulated genomes are used
            @type keep_original: bool
            @param max_processors: maximum number of processors available to be used
            @type max_processors: int
            @param tmp_dir: working directory or place temporary files can be stored
            @type tmp_dir: str | unicode
            @param logfile: file handler or file path to a log file
            @type logfile: str | file | io.FileIO | StringIO.StringIO
            @param verbose: Not verbose means that only warnings and errors will be past to stream
            @type verbose: bool
            @param debug: If True logger will output DEBUG messages
            @type debug: bool
            @param seed: The seed used for initiation of the 'random' module
            @type seed: int | float | str | unicode

            @return: None
            @rtype: None
        """
        assert isinstance(self.keep_original, bool)
        assert isinstance(self.separator, str)
        assert isinstance(self.column_name_gid, str)
        assert isinstance(self.column_name_ncbi, str)
        assert isinstance(self.column_name_source, str)
        assert isinstance(self.filename_prefix, str)

        #self._debug = debug
        #if debug:
        #	self._logger.set_level(self._logger.DEBUG)

        if seed is not None:
            random.seed(seed)
            np_random.seed(abs(hash(seed)) % 4294967295)  # numpy accepts only 32 bit integers

        assert isinstance(self.max_processors, int)
        #self._max_processors = max_processors#

        #directory_sgevolver = self.get_full_path(os.path.join(os.path.dirname(__file__), "sgEvolver"))
        #self._executable_sim = executable_sim
        #if self._executable_sim is None:
        #	self._executable_sim = os.path.join(directory_sgevolver, "simujobrun.pl")
        #assert self.validate_file(self._executable_sim, executable=True)

        #if self._directory_template is None:
        #	self._directory_template = self.get_full_path(os.path.join(os.path.dirname(__file__), "sgEvolver", "simulation_dir"))
        #assert self.validate_dir(self._directory_template, file_names=[self._filename_tree, self._filename_parameter])

        #self._tmp_dir = tmp_dir
        #assert self.validate_dir(self._tmp_dir)

        #self._directory_strain = self.get_full_path(os.path.join(self._tmp_dir, "{gid}.strains"))
        #file_path_template_newick_tree = os.path.join(directory_template, directory_template_filenames[1])
        #filenames_strains = get_filenames_strains(file_path_template_newick_tree)
        #assert len(filenames_strains) > 0

    def get_genome_amounts_geometric(self, probability, max_genome_amount, geometric_probability=0.3):
        """
        Get amounts of genomes by original genome

        @param probability: Proportion of simulated original genomes
        @type probability: int  | float
        @param max_genome_amount: Total number of genomes
        @type max_genome_amount: int 

        @return: List of integers representing amount of strains
        @rtype: list[int]
        """
        assert isinstance(probability, (int, float))
        assert 0 <= probability <= 1
        assert isinstance(max_genome_amount, int)

        final_amounts = []
        while sum(final_amounts) < max_genome_amount:
            if random.uniform(0, 1) < probability:
                final_amounts.append(1)
            else:
                amount = 1 + np_random.geometric(geometric_probability)
                final_amounts.append(amount)

        final_amounts[-1] -= sum(final_amounts) - max_genome_amount
        return final_amounts

    def get_genome_amounts(self, probability, max_genome_amount):
        """
        Get amounts of genomes by original genome

        @param probability: Proportion of simulated original genomes
        @type probability: int  | float
        @param max_genome_amount: Total number of genomes
        @type max_genome_amount: int 

        @return: List of integers representing amount of strains
        @rtype: list[int]
        """
        assert isinstance(probability, (int, float))
        assert 0 <= probability <= 1
        assert isinstance(max_genome_amount, int)

        genome_amounts = self.get_genome_amounts_geometric(probability, max_genome_amount)
        diverence = Counter(genome_amounts)[1] / float(len(genome_amounts))
        if max_genome_amount >= 10:
            while abs(diverence - probability) > 0.05:
                # print "need: {}, gotten: {}".format(probability, diverence)
                genome_amounts = self.get_genome_amounts_geometric(probability, max_genome_amount)
                diverence = Counter(genome_amounts)[1] / float(len(genome_amounts))
        return genome_amounts

    def get_genome_amounts_geometric_fix(self, num_real_genomes, max_genome_amount, geometric_probability=0.3):
        """
        Get amounts of genomes by original genome

        @param num_real_genomes: exact number of real genomes
        @type num_real_genomes: int 
        @param max_genome_amount: Total number of genomes
        @type max_genome_amount: int 

        @return: List of integers representing amount of strains
        @rtype: list[int]
        """
        assert isinstance(num_real_genomes, int)
        assert isinstance(max_genome_amount, int)

        final_amounts = [1] * num_real_genomes
        index = 0
        while index < len(final_amounts):
            if sum(final_amounts) >= max_genome_amount:
                break
            final_amounts[index] += 1 + np_random.geometric(geometric_probability)
            index += 1

        final_amounts[index-1] -= sum(final_amounts) - max_genome_amount
        return final_amounts

    def get_genome_amounts(self, probability, max_genome_amount, num_real_genomes=None, silent=True):
        """
        Get amounts of genomes by original genome

        @param probability: Proportion of simulated original genomes
        @type probability: int  | float
        @param max_genome_amount: Total number of genomes
        @type max_genome_amount: int 
        @param num_real_genomes: exact number of real genomes
        @type num_real_genomes: int 

        @return:
        @rtype: list[int]
        """
        assert probability is None or isinstance(probability, (int, float))
        if probability:
            assert 0 <= probability <= 1
        assert isinstance(max_genome_amount, int)
        assert isinstance(num_real_genomes, int)
        assert isinstance(silent, bool)

        if num_real_genomes is not None:
            genome_amounts = self.get_genome_amounts_geometric_fix(num_real_genomes, max_genome_amount)
        else:
            genome_amounts = self.get_genome_amounts(probability, max_genome_amount)

        #if not silent:
        #	print_distribution(genome_amounts)
        #	message = "Do you accept this distribution? [y/n]"
        #	while not get_confirmation(message):
        #		if num_real_genomes is not None:
        #			genome_amounts = get_genome_amounts_geometric_fix(num_real_genomes, max_genome_amount)
        #		else:
        #			genome_amounts = get_genome_amounts(probability, max_genome_amount)
        #		print_distribution(genome_amounts)
        return genome_amounts
    
    def has_unique_columns(self, list_of_column_names=None):
        if list_of_column_names is None:
            list_of_column_names = self._list_of_column_names
        return len(list_of_column_names) == len(set(list_of_column_names))
    
    def parse_column_names(self, stream_input, separator):
        row = stream_input.readline().rstrip('\n').rstrip('\r')
        list_of_column_names = row.split(separator)
        assert self.has_unique_columns(list_of_column_names), "Column names must be unique!"
        return list_of_column_names

    def read(self, file_path, separator=None, column_names=False, comment_line=None):
        """
            Reading comma or tab separated values in a file as table

            @param file_path: path to file to be opened
            @type file_path: str | unicode
            @param separator: default character assumed to separate values in a file
            @type separator: str | unicode
            @param column_names: True if column names available
            @type column_names: bool
            @param comment_line: character or list of character indication comment lines
            @type comment_line: str | unicode | list[str|unicode]

            @return: None
            @rtype: None
        """
        if comment_line is None:
            comment_line = ['#']
        elif isinstance(comment_line, str):
            comment_line = [comment_line]

        if separator is None:
            separator="\t"

        #assert isinstance(file_path, str)
        #assert self.validate_file(file_path)
        assert isinstance(separator, str)
        assert isinstance(comment_line, list)
        assert isinstance(column_names, bool)

        if isinstance(file_path, (str, bytes, os.PathLike)):
            file_handler = open(file_path, 'r')
        elif isinstance(file_path, io.TextIOWrapper):
            file_handler = file_path
        else:
            raise TypeError("Invalid type for file_path")

        #with open(file_path, 'r') as file_handler:
            #self._logger.info("Reading file: '{}'".format(file_path))

        meta_table = {}
        number_of_rows=0

        # read column names
        if column_names:
            list_of_column_names = self.parse_column_names(file_handler, separator)
            for column_name in list_of_column_names:
                meta_table[column_name] = []
            
        #global list_of_column_names

        # read rows
        row_count = 0
        for line in file_handler:

            row_count += 1
            row = line.rstrip('\n').rstrip('\r')
            if line[0] in comment_line or len(row) == 0:
                continue
            number_of_rows += 1
            row_cells = row.split(separator)
            number_of_columns = len(list(self.list_of_column_names))
            if number_of_columns != 0 and number_of_columns != len(row_cells):
                msg = "Format error. Bad number of values in line {}".format(row_count)
                #self._logger.error(msg)
                raise ValueError(msg)
            for index, value in enumerate(row_cells):
                if column_names:
                    column_name = self.list_of_column_names[index]
                else:
                    column_name = index
                    if column_name not in meta_table:
                        meta_table[column_name] = []
                meta_table[column_name].append(row_cells[index].rstrip('\n').rstrip('\r'))
            if number_of_columns == 0:
                self.list_of_column_names = sorted(meta_table.keys())

        return meta_table, number_of_rows
        
    def has_column(self, meta_table, column_name):
        """
            Get index of value in a column

            @attention:

            @param column_name: column name
            @type column_name: int | long | str | unicode

            @return: True if column available
            @rtype: bool
        """

        assert isinstance(column_name, (str, int))

        if column_name in meta_table:
            return True
        else:
            return False
        
    def validate_column_names(self, meta_table, list_of_column_names):
        """
            Validate that a list of column names exists in the loaded table

            @attention:

            @param list_of_column_names: column name
            @type list_of_column_names: list[str|unicode]

            @return: True if all column names exist
            @rtype: bool
        """
        assert isinstance(list_of_column_names, list)

        list_of_invalid_column_names = []
        for column_name in list_of_column_names:
            if not self.has_column(meta_table, column_name):
                list_of_invalid_column_names.append(column_name)
                #self._logger.info("Invalid columns: {}".format(", ".join(list_of_invalid_column_names)))
        if len(list_of_invalid_column_names) > 0:
            return False
        return True

    def get_column(self, metadata_table, column_name):
        """
            Get a column

            @attention: use index if no name available

            @param column_name: column name
            @type column_name: int | long | str | unicode

            @return: Cell values of a column
            @rtype: list[str|unicode]
        """
        assert isinstance(column_name, (str, int))
        assert self.has_column(metadata_table, column_name), "Column '{}' not found!".format(column_name)
        return list(metadata_table[column_name])

    def get_number_of_rows(self, metadata_table):
        """
        Get number of rows

        @attention:

        @return: Number of rows
        @rtype: int
        """
        first_value = next(iter(metadata_table.values()))
        return len(first_value)

    def add_strain(self, novelty_category_dict, novelty_category_name, otu_id, strain_id):
        """
        Add a strain

        @param otu_id: Identifier of otu.
        @type otu_id: str | unicode
        @param strain_id: Identifier of strain.
        @type strain_id: str | unicode

        @rtype: None
        """
        assert isinstance(otu_id, str)
        assert isinstance(strain_id, str)
        if otu_id not in novelty_category_dict[novelty_category_name][1]:
            novelty_category_dict[novelty_category_name][1][otu_id] = []
        novelty_category_dict[novelty_category_name][1][otu_id].append(strain_id)
        novelty_category_dict[novelty_category_name][0] += 1

        return novelty_category_dict

    def per_category(self, per_cat, rest):
        """
        Get maximum number of strains valid to be drawn from a category.

        @rtype: int
        """
        return per_cat + (rest > 0)

    def get_all_strains(self, name, novelty_category):
        """
        Return list of all strains from this category

        @return: List of drawn strain id
        @rtype: list[str|unicode]
        """
        out = []
        for otu_id in novelty_category[name][1]:
            for strain_id in novelty_category[name][1][otu_id]:
                out.append(strain_id)
        return out

    def recalc(self, draw, cats, per_cat, rest, number_of_drawn):
        """
        Adjust values based on the number of strains drawn from the previous category

        @param number_of_drawn: Number of drawn strains from the previous category
        @rtype: None
        """
        assert isinstance(number_of_drawn, int)
        draw -= number_of_drawn
        cats -= 1
        per_cat = int(draw/cats)
        rest = draw % cats
     
        return draw, cats, per_cat, rest

    def draw_strains_subset(self, categories, name, total, limit_per_otu):
        """
        Draw subset of strains from this category

        @param total: Amount to be drawn.
        @type total: int
        @param limit_per_otu: Maximum number of strains to be drawn from an otu
        @type limit_per_otu: int

        @return: List of drawn strain id
        @rtype: list[str|unicode]
        """
        assert isinstance(total, int)
        assert isinstance(limit_per_otu, int)
        assert total <= categories[name][0]
        assert limit_per_otu > 0
        drawn_strain = []
        overhead = []
        drawn_strain_count_overall = 0
        
        # Bug?
        # for otu_id in random.sample(categories[name][1].keys(), len(categories[name][1])):
        
        for otu_id in random.sample(list(categories[name][1].keys()), len(categories[name][1])):
            drawn_strain_count_otu = 0
            for strain_id in random.sample(categories[name][1][otu_id], len(categories[name][1][otu_id])):
                if drawn_strain_count_otu < limit_per_otu and drawn_strain_count_overall < total:
                    drawn_strain.append(strain_id)
                    drawn_strain_count_otu += 1
                    drawn_strain_count_overall += 1
            overhead += list(set(categories[name][1][otu_id])-set(drawn_strain))

        if drawn_strain_count_overall < total:
            # out += overhead[0:total-drawn_strain_cound_overall]
            #print "\n", total-drawn_strain_count_overall, len(overhead), len(drawn_strain), "\n"
            drawn_strain += random.sample(overhead, total-drawn_strain_count_overall)
        return drawn_strain

    def substract(self, draw, per_cat, rest):
        """
        Adjust values by subtracting the maximum number of strains that can be drawn from a category.

        @rtype: None
        """
        draw -= per_cat
        if rest > 0:
            rest -= 1
        return draw, rest

    def draw_strains(self, categories, number_of_strains_to_draw, max_amount_of_strains_per_otu):
        """
        Get list of drawn strains

        @param categories:
        @type: dict[str, NoveltyCategory]
        @param number_of_strains_to_draw: Amount of strains to be drawn.
        @type: int
        @param max_amount_of_strains_per_otu: Maximum number of strains to be drawn from one category.
        @type: int

        @return: @raise ValueError:
        @rtype: list[str|unicode]
        """
        assert isinstance(categories, dict)
        assert isinstance(number_of_strains_to_draw, int)
        assert isinstance(max_amount_of_strains_per_otu, int)
        number_of_novelty_categories = len(categories)

        cats = number_of_novelty_categories
        draw = number_of_strains_to_draw
        per_cat = int(number_of_strains_to_draw/number_of_novelty_categories)
        rest = number_of_strains_to_draw % number_of_novelty_categories
        per_otu = max_amount_of_strains_per_otu

        drawn_genome_id = []
        for x in sorted(categories, key=lambda name: categories[name][0]):
            if categories[x][0] < self.per_category(per_cat, rest):
                # draw all strains
                drawn_genome_id += categories[x].get_all_strains(x, categories)
                draw, cats, per_cat, rest = self.recalc(draw, cats, per_cat, rest, categories[x][0])
            else:
                # draw subset of strains
                drawn_genome_id += self.draw_strains_subset(categories, x, self.per_category(per_cat, rest), per_otu)
                draw, rest = self.substract(draw, per_cat, rest)

        if len(drawn_genome_id) < number_of_strains_to_draw:
            msg = "Could only draw {} samples instead of {}".format(len(drawn_genome_id), number_of_strains_to_draw)
            #self._logger.error(msg)
            raise ValueError(msg)
        return drawn_genome_id

    def get_drawn_genome_id(self, metadata_table, number_of_strains, number_of_strains_per_otu):
        """
        Get a list of drawn genome ids.

        @param metadata_table: Metadata containing genome id, otu and a novelty category
        @type metadata_table: MetadataTable
        @param number_of_strains: Amount of strains to be drawn
        @type number_of_strains: int
        @param number_of_strains_per_otu:
        @type number_of_strains_per_otu: int

        @return: list of drawn genome ids @raise ValueError:
        @rtype: list[str|unicode]
        """
        #assert isinstance(metadata_table, MetadataTable)
        assert isinstance(number_of_strains, int)
        assert isinstance(number_of_strains_per_otu, int)
        genomes_read_from_metadata = 0
        lost_lines = 0

        novelty_categories = dict()

        required_headers = [self.column_name_gid, self.column_name_otu, self.column_name_novelty_category]
        if not self.validate_column_names(metadata_table, required_headers):
            msg = 'A header is missing, needed: {}'.format(', '.join(required_headers))
            #self._logger.error(msg)
            raise ValueError(msg)

        column_genome_id = self.get_column(metadata_table, self.column_name_gid)
        column_otu = self.get_column(metadata_table, self.column_name_otu)
        column_novelty = self.get_column(metadata_table, self.column_name_novelty_category)
        row_count = self.get_number_of_rows(metadata_table)

        # read metadata from provided file
        for row_index in range(0, row_count):
            if column_genome_id[row_index] == '' or column_otu[row_index] == '' or column_novelty[row_index] == '':
                lost_lines += 1
                continue
            genomes_read_from_metadata += 1
            novelty = column_novelty[row_index]
            if novelty not in novelty_categories:
                novelty_categories[novelty] = [0, {}] # first value is the number of strains and the second the OTUs
            novelty_categories = self.add_strain(novelty_categories, novelty, column_otu[row_index], column_genome_id[row_index])

        if lost_lines > 0:
            msg = "Invalid lines/ unusable genomes: {}".format(lost_lines)
            #self._logger.debug(msg)

        if genomes_read_from_metadata < number_of_strains:
            print(genomes_read_from_metadata)
            print(number_of_strains)
            msg = 'Not enough data to draw.'
            #self._logger.error(msg)
            raise ValueError(msg)

        # sample OTUs from the pool of available OTUs
        drawn_genome_id = self.draw_strains(novelty_categories, number_of_strains, number_of_strains_per_otu)
        return drawn_genome_id

    def write(self, metadata_table_community, file_path, number_of_rows, column_names=False, exclude=None, value_list=None, key_column_name=None):
        """
        Write tab separated files

        @attention: No comments will be written

        @param file_path: path to file to be opened
        @type file_path: str | unicode
        @param separator: file handler or file path to a log file
        @type separator: str | unicode
        @param column_names: True if column names should be written
        @type column_names: bool
        @param compression_level: any value above 0 will compress files
        @type compression_level: int | long
        @param exclude: If True, rows with a value in the value_list at the key_column_names are removed, False: all others are removed
        @type exclude: None | bool
        @param value_list:
        @type value_list: list[str|unicode]
        @param key_column_name: column name of excluded or included rows
        @type key_column_name: str | unicode

        @return: None
        @rtype: None
        """
        assert isinstance(file_path, str)
        assert isinstance(self.separator, str)
        assert isinstance(column_names, bool)
        #assert isinstance(compression_level, int)
        #assert 0 <= compression_level < 10
        assert exclude is None or isinstance(exclude, bool)
        assert value_list is None or isinstance(value_list, list)
        assert key_column_name is None or isinstance(key_column_name, str), "Invalid: {}".format(key_column_name)

        #if compression_level > 0:
            #file_handler = self.open(file_path, "w", compression_level)
        #else:
        file_handler = open(file_path, "w")

        if column_names:
            if not isinstance(self.list_of_column_names[0], str):
                header = self.separator.join([str(index) for index in self.list_of_column_names])
            else:
                header = self.separator.join(self.list_of_column_names)
            file_handler.write(header + '\n')
        for row_number in range(0, number_of_rows):
            if exclude is not None:
                if not exclude and metadata_table_community[key_column_name][row_number] not in value_list:
                    continue
                if exclude and metadata_table_community[key_column_name][row_number] in value_list:
                    continue

            row = []
            for column_names in self.list_of_column_names:
                row.append(str(metadata_table_community[column_names][row_number]))
            file_handler.write(self.separator.join(row) + '\n')
        file_handler.close()

    def get_genome_id_to_path_map(self, dict_mapping_genome_id_to_paths, list_of_drawn_genome_id):
        """
        Get a dictionary mapping genome id to the path of their genome

        @param file_path_of_file_mapping_genome_id_to_paths: File path to file with format 'id \t path'
        @type file_path_of_file_mapping_genome_id_to_paths: str | unicode
        @param list_of_drawn_genome_id: List of genome identifiers
        @type list_of_drawn_genome_id: list[str|unicode]

        @return: genome ids mapped to their gnome file path
        @rtype: dict[str|unicode, str|unicode]
        """
        msg = "Dict mapping genome id to paths is missing one or more genome id."
        assert set(dict_mapping_genome_id_to_paths.keys()).issuperset(list_of_drawn_genome_id), msg
        return {genome_id: dict_mapping_genome_id_to_paths[genome_id] for genome_id in list_of_drawn_genome_id}

    def reduce_rows_to_subset(self, metadata_table_community, list_of_values, key_column_name):
        """
        Keep rows at key values of a column

        @attention:

        @param list_of_values: Cell values of table column
        @type list_of_values: list[str|unicode]
        @param key_column_name: Column name
        @type key_column_name: str | unicode

        @return: Nothing
        @rtype: None
        """

        assert isinstance(key_column_name, (str, int))
        assert isinstance(list_of_values, list)
        assert self.has_column(metadata_table_community, key_column_name), "Column '{}' not found!".format(key_column_name)

        new_meta_table = {}
        for column_name in self.list_of_column_names:
            new_meta_table[column_name] = []
        column = self.get_column(metadata_table_community, key_column_name)
        for index, value in enumerate(column):
            if value not in list_of_values:
                continue
            for column_name in self.list_of_column_names:
                new_meta_table[column_name].append(metadata_table_community[column_name][index])
        metadata_table_community = new_meta_table
        #self._number_of_rows = len(self._meta_table[key_column_name])

        return metadata_table_community

    def get_genome_id_to_amounts(self, list_of_drawn_genome_id, genome_amounts):
        """
        Assign genome IDs to genome amounts

        @param list_of_drawn_genome_id:
        @type list_of_drawn_genome_id: list[str | unicode]
        @param genome_amounts: List of integers representing amount of strains
        @type genome_amounts: list[int]

        @return: Mapping from genome id to the amount of strains
        @rtype : dict[str | unicode, int]
        """
        assert isinstance(list_of_drawn_genome_id, list)
        assert isinstance(genome_amounts, list)
        genome_id_to_amounts = {}
        for index, genome_id in enumerate(list_of_drawn_genome_id):
            genome_id_to_amounts[genome_id] = genome_amounts[index]
        return genome_id_to_amounts

    def simulate_strains(self, meta_table, genome_id_to_amounts, genome_id_to_file_path_genome, genome_id_to_file_path_gff=None):
        """
        Uses sgEvolver to generate strain-level diversity around an isolate assembly
        and add randomly picked strains to genome_id_to_file_path_genome and metadata table.

        @attention genome_id_to_file_path_genome: Will be extended with IDs and file paths to the strains

        @param meta_table: Metadata table containing genome information
        @type meta_table: MetadataTable
        @param genome_id_to_amounts: Mapping from genome id to the amount of strains
        @type genome_id_to_amounts: dict[str, int]
        @param genome_id_to_file_path_genome: Mapping from genome id to the file path of the genome
        @type genome_id_to_file_path_genome: dict[str, str]
        @param genome_id_to_file_path_gff: Mapping from genome id to the file path of the gene annotations of a genome
        @type genome_id_to_file_path_gff: dict[str, str]

        @return: Nothing
        @rtype: None
        """
        assert isinstance(meta_table, dict)
        assert isinstance(genome_id_to_amounts, dict)
        assert isinstance(genome_id_to_file_path_genome, dict)
        assert genome_id_to_file_path_gff is None or isinstance(genome_id_to_file_path_gff, dict)
        if genome_id_to_file_path_gff is None:
            msg = "No gff file (gene annotation) was given. Simulating strains without such a file can break genes."
            #self._logger.warning(msg)
        #for file_path in genome_id_to_file_path_genome.values():
            #self.validate_file(file_path)
        #if genome_id_to_file_path_gff is not None:
            #for file_path in genome_id_to_file_path_gff.values():
                #self.validate_file(file_path)

        self.write_to_tsv(genome_id_to_amounts, genome_id_to_file_path_genome, genome_id_to_file_path_gff)
        #self._pick_random_strains(meta_table, genome_id_to_amounts, genome_id_to_file_path_genome)

        # read file and generate strain diversity for each assembly
        # then subsample the strains

    def touch(self, file_path):
        file_handle = open(file_path, 'w')
        file_handle.close()

    def write_to_tsv(self, genome_id_to_amounts, genome_id_to_file_path_genome, genome_id_to_file_path_gff=None):
        """
        Use sgEvolver to generate strain-level diversity around an isolate assembly.

        @attention genome_id_to_file_path_genome: Will be extended with IDs and file paths to the strains

        @param genome_id_to_amounts: Mapping from genome id to the amount of strains
        @type genome_id_to_amounts: dict[str, int]
        @param genome_id_to_file_path_genome: Mapping from genome id to the file path of the genome
        @type genome_id_to_file_path_genome: dict[str, str]
        @param genome_id_to_file_path_gff: Mapping from genome id to the file path of the gene annotations of a genome
        @type genome_id_to_file_path_gff: dict[str, str]

        @return: Nothing
        @rtype: None
        """
        with open("genome_id_to_file_amount_gff.tsv", 'w') as genome_id_to_file_path_genome_file:

            for genome_id in genome_id_to_file_path_genome.keys():

                seed =  random.randint(0, sys.maxsize)

                if genome_id_to_file_path_gff is None:
                    genome_id_to_file_path_genome_file.write("{id}\t{path}\t{amount}\t{seed}\n".format(id=genome_id, path=genome_id_to_file_path_genome[genome_id], 
                                                                                                  amount=genome_id_to_amounts[genome_id], seed=seed))
                else:
                    file_path_gff = genome_id_to_file_path_gff[genome_id]
                    genome_id_to_file_path_genome_file.write("{id}\t{path}\t{amount}\t{seed}\t{gff}\n".format(id=genome_id, path=genome_id_to_file_path_genome[genome_id], 
                                                                                                  amount=genome_id_to_amounts[genome_id], seed=seed, gff=file_path_gff))
    
    def prepare(self, genomes_total, genomes_real, seed, metadata_file, max_strains_per_otu, id_to_genome_dict, id_to_gff_dict):
        
        verbose = False
    
        # pick how much a strain will be simulated
        genome_amounts = []

        self.init(#executable_sim=None,
            #directory_template=directory_simulation_template,
            seed=seed)

        probability = None  # 1-options.communities[community_id]["evolve"]
        genome_amounts = self.get_genome_amounts(
            probability=probability,
            max_genome_amount=genomes_total,
            num_real_genomes=genomes_real,
            silent=not verbose
            )
        number_of_strains = len(genome_amounts)

        # draw strains
        #self._logger.info("Drawing strains.")
        metadata_table_community, number_of_rows = self.read(metadata_file, separator=None, column_names=True)

        list_of_drawn_genome_id = self.get_drawn_genome_id(
            metadata_table=metadata_table_community,
            number_of_strains=number_of_strains,
            number_of_strains_per_otu=max_strains_per_otu)
    
        # write unused data to separate file
        if isinstance(metadata_file, (str, bytes, os.PathLike)):
            old_base_name = os.path.basename(metadata_file)
        elif isinstance(metadata_file, io.TextIOWrapper):
            old_base_name = os.path.basename(metadata_file.name)
        else:
            raise TypeError("Invalid type for file_path")
        #old_base_name = os.path.basename(metadata_file)

        file_prefix, extention = os.path.splitext(old_base_name)
        new_file_name = "unused_c{index}_{prefix}{ext}".format(
            prefix=file_prefix,
            #index=community.id,
            index=0,
            ext=extention)
        metadata_new_file_path = os.path.join("./", new_file_name)
        self.write(metadata_table_community=metadata_table_community,
            file_path=metadata_new_file_path,
            number_of_rows=number_of_rows,
            exclude=True,
            value_list=list_of_drawn_genome_id,
            key_column_name=self.column_name_gid,
            column_names=True)

        # get path for every genome
        genome_id_to_file_path_gff = None
        if id_to_gff_dict:
            genome_id_to_file_path_gff = self.get_genome_id_to_path_map(
                id_to_gff_dict, list_of_drawn_genome_id)
        genome_id_to_path_map = self.get_genome_id_to_path_map(
            id_to_genome_dict, list_of_drawn_genome_id)
    
        # concatenate
        metadata_table_community = self.reduce_rows_to_subset(metadata_table_community, list_of_drawn_genome_id, self.column_name_gid)

        # Is this necessary? In metagenomesimulation.py metadata_table seems to bee an empty table?
        #metadata_table.concatenate(metadata_table_community, strict=False)

        # validate correct format of files
        #self._logger.info("Validating raw sequence files!")
        #assert self.validate_format(
        #    list_of_file_paths=genome_id_to_path_map.values(),
        #    file_format="fasta",
        #    sequence_type="dna",
        #    ambiguous=True
        #    ), "Validation of file format failed!"

        # simulate diversity around strains
        #if community.simulate_strains:
        genome_id_to_amounts = self.get_genome_id_to_amounts(list_of_drawn_genome_id, genome_amounts)
        self.simulate_strains(
            meta_table={},  # is this dict empty? see metagenomesimulation.py line 232
            genome_id_to_amounts=genome_id_to_amounts,
            genome_id_to_file_path_genome=genome_id_to_path_map,
            genome_id_to_file_path_gff=genome_id_to_file_path_gff)
        # adopt new list that includes simulated strains
        #self._logger.info("Validating simulated sequence files!")
        # also need to be checked after sgevolver???
        #for genome_id, file_path in genome_id_to_path_map.items():
        #    if genome_id in list_of_drawn_genome_id:
        #        continue
        #    assert self.validate_sequence_file(
        #        file_path,
        #        file_format="fasta",
        #        sequence_type="dna",
        #        ambiguous=True)
        #list_of_drawn_genome_id = genome_id_to_path_map.keys()

def parse_tsv_to_dict(file):
    result_dict = {}
    with open(file, 'r') as f:
        for line in f:
            key, value = line.strip().split('\t')
            result_dict[key] = value
    return result_dict

# TODO mehrere communities
# main method and entry point of this script
if __name__ == "__main__":

    prepareStrainSimulation = PrepareStrainSimulation()

    parser = argparse.ArgumentParser()
    #parser.add_argument(
    #    "-template",
    #    help="The directory containing the strain simulation template",
    #    action='store',
    #    default="")
    parser.add_argument(
        "-genomes_total",
        type=int,
        help="Total number of simulated genomes")
    parser.add_argument(
        "-genomes_real",
        type=int,
        help="Number of genomes used from the input genomes")
    parser.add_argument(
        '-seed', 
        type=int, 
        help='The seed to use')
    parser.add_argument(
        "-metadata",
        help="The metadata table",
        action='store',
        type=argparse.FileType('r'),
        default="")
    parser.add_argument(
        '-max_strains_per_otu', 
        type=int, 
        help='Maximum number of strains drawn from genomes belonging to a single OTU')
    #parser.add_argument(
    #    "-genome_locations",
    #    help="The metadata table",
    #    action='store',
    #    type=argparse.FileType('r'),
    #    default="")
    parser.add_argument(
        "-id_to_genome_file",
        help="File containing tuples mapping genome IDs to reference genomes with key<tab>file_path pairs",
        type=str,
        required=True)
    parser.add_argument(
        "-id_to_gff_file",
        help="File containing tuples mapping genome IDs to annotation files with key<tab>file_path pairs",
        type=str,
        required=False)
    
    options = parser.parse_args()
    
    #directory_simulation_template = options.template
    genomes_total = options.genomes_total
    genomes_real = options.genomes_real
    seed = options.seed
    metadata_file = options.metadata
    max_strains_per_otu = options.max_strains_per_otu

    # Parse the tab-separated files into dictionaries
    id_to_genome_dict = parse_tsv_to_dict(options.id_to_genome_file)
    id_to_gff_dict = parse_tsv_to_dict(options.id_to_gff_file) if options.id_to_gff_file else None

    prepareStrainSimulation.prepare(genomes_total, genomes_real, seed, metadata_file, max_strains_per_otu, id_to_genome_dict, id_to_gff_dict)