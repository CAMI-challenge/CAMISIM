# !/usr/bin/env python3

__author__ = 'Peter Hofmann'
__original_author__ = "Aaron Daring"
__version__ = '0.0.5'


from scripts.parallel import TaskCmd, runCmdParallel, reportFailedCmd
from scripts.Validator.validator import Validator
from scripts.MetaDataTable.metadatatable import MetadataTable
import sys
import os
import random
import tempfile
import shutil
import numpy.random as np_random
from collections import Counter
from Bio import Phylo


class GenomeOrganizer(Validator):
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
			genome_amounts = self._get_genome_amounts_geometric_fix(num_real_genomes, max_genome_amount)
		else:
			genome_amounts = self._get_genome_amounts(probability, max_genome_amount)

		if not silent:
			self.print_distribution(genome_amounts)
			message = "Do you accept this distribution? [y/n]"
			while not self.get_confirmation(message):
				if num_real_genomes is not None:
					genome_amounts = self._get_genome_amounts_geometric_fix(num_real_genomes, max_genome_amount)
				else:
					genome_amounts = self._get_genome_amounts(probability, max_genome_amount)
				self.print_distribution(genome_amounts)
		return genome_amounts

	@staticmethod
	def print_distribution(genome_amounts):
		"""
		Print genome amounts to console

		@param genome_amounts: number of genomes for each original genome
		@type genome_amounts: list[int]

		@return: Nothing
		@rtype: None
		"""
		assert isinstance(genome_amounts, list)
		counter = Counter(genome_amounts)
		text = "{sep}".join(["{}: {}".format(counter[k], k) for k in counter]).format(sep=os.linesep)
		print("{sep}Using {genoms} original genomes.{sep}<#genomes>: <#strains>{sep}{counter}".format(
			genoms=len(genome_amounts), counter=text, sep=os.linesep))

	def _get_genome_amounts(self, probability, max_genome_amount):
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

		genome_amounts = self._get_genome_amounts_geometric(probability, max_genome_amount)
		diverence = Counter(genome_amounts)[1] / float(len(genome_amounts))
		if max_genome_amount >= 10:
			while abs(diverence - probability) > 0.05:
				# print "need: {}, gotten: {}".format(probability, diverence)
				genome_amounts = self._get_genome_amounts_geometric(probability, max_genome_amount)
				diverence = Counter(genome_amounts)[1] / float(len(genome_amounts))
		return genome_amounts

	@staticmethod
	def _get_genome_amounts_exponential(probability, max_genome_amount):
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
			amount = np_random.geometric(probability)
			final_amounts.append(amount)

		final_amounts[-1] -= sum(final_amounts) - max_genome_amount
		return final_amounts

	@staticmethod
	def _get_genome_amounts_geometric_fix(num_real_genomes, max_genome_amount, geometric_probability=0.3):
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

	@staticmethod
	def _get_genome_amounts_geometric(probability, max_genome_amount, geometric_probability=0.3):
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

	@staticmethod
	def _get_genome_amounts_uniform(probability, max_genome_amount):
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
				amount = 1 + random.randint(1, 3)
				final_amounts.append(amount)

		final_amounts[-1] -= sum(final_amounts) - max_genome_amount
		return final_amounts

	def get_confirmation(self, message):
		"""
			Confirm something with user

			@attention:

			@param message: Question to be confirmed
			@type message: str

			@return: Nothing
			@rtype: bool
		"""
		assert isinstance(message, str)

		if not message:
			raise AssertionError("asd")
		user_input = raw_input("{}\n>".format(message)).lower()
		while True:
			if self.is_boolean_state(user_input):
				return self.get_boolean_state(user_input)
			user_input = raw_input("Please type 'n' for no, or 'y' for yes:\n>").lower()


class StrainSimulationWrapper(GenomeOrganizer):
	"""
	StrainSimulationWrapper
	Generates additional substrain-level diversity around a draft assembly
	"""
	_label = "StrainSimulationWrapper"
	_filename_parameter = "simujobparams.pm"
	_filename_tree = "template.tree"

	_directory_template_filenames = ["simujobparams.pm", "template.tree"]

	def __init__(
		self, executable_sim=None, directory_template=None,
		column_name_gid="genome_ID", column_name_ncbi="NCBI_ID", column_name_source="source", separator='\t',
		filename_prefix="simulated_", keep_original=True,
		max_processors=1, tmp_dir=None, logfile=None, verbose=True, debug=False, seed=None):
		"""
			Initialize instance with seed

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
		super(StrainSimulationWrapper, self).__init__(logfile, verbose)
		assert isinstance(keep_original, bool)
		assert isinstance(separator, str)
		assert isinstance(column_name_gid, str)
		assert isinstance(column_name_ncbi, str)
		assert isinstance(column_name_source, str)
		assert isinstance(filename_prefix, str)
		assert isinstance(debug, bool)

		if tmp_dir is None:
			tmp_dir = tempfile.gettempdir()

		self._debug = debug
		if debug:
			self._logger.set_level(self._logger.DEBUG)

		if seed is not None:
			random.seed(seed)
			np_random.seed(abs(hash(seed)) % 4294967295)  # numpy accepts only 32 bit integers

		assert isinstance(max_processors, int)
		self._max_processors = max_processors

		self._separator = separator
		self._column_name_gid = column_name_gid
		self._column_name_ncbi = column_name_ncbi
		self._column_name_source = column_name_source
		self._filename_prefix = filename_prefix
		self._keep_original = keep_original
		self._directory_template = directory_template

		directory_sgevolver = self.get_full_path(os.path.join(os.path.dirname(__file__), "sgEvolver"))
		self._executable_sim = executable_sim
		if self._executable_sim is None:
			self._executable_sim = os.path.join(directory_sgevolver, "simujobrun.pl")
		assert self.validate_file(self._executable_sim, executable=True)

		if self._directory_template is None:
			self._directory_template = self.get_full_path(os.path.join(os.path.dirname(__file__), "sgEvolver", "simulation_dir"))
		assert self.validate_dir(self._directory_template, file_names=[self._filename_tree, self._filename_parameter])

		self._tmp_dir = tmp_dir
		assert self.validate_dir(self._tmp_dir)

		self._directory_strain = self.get_full_path(os.path.join(self._tmp_dir, "{gid}.strains"))
		file_path_template_newick_tree = os.path.join(self._directory_template, self._directory_template_filenames[1])
		self._filenames_strains = self.get_filenames_strains(file_path_template_newick_tree)
		assert len(self._filenames_strains) > 0

	@staticmethod
	def _get_seed():
		return random.randint(0, sys.maxsize)

	def _get_simulate_cmd(self, directory_strains, filepath_genome, filepath_gff):
		"""
		Get system command to start simulation. Change directory to the strain directory and start simulating strains.

		@param directory_strains: Directory for the simulated strains
		@type directory_strains: str | unicode
		@param filepath_genome: Genome to get simulated strains of
		@type filepath_genome: str | unicode
		@param filepath_gff: gff file with gene annotations
		@type filepath_gff: str | unicode

		@return: System command line
		@rtype: str
		"""
		cmd_run_simujobrun = "cd {dir}; {executable} {filepath_genome} {filepath_gff} {seed}" + " >> {log}"
		cmd = cmd_run_simujobrun.format(
			dir=directory_strains,
			executable=self._executable_sim,
			filepath_genome=filepath_genome,
			filepath_gff=filepath_gff,
			seed=self._get_seed(),
			log=os.path.join(directory_strains, os.path.basename(filepath_genome) + ".sim.log")
		)
		return cmd

	def _prepare_simulation_subfolder(self, directory_strains):
		"""
		Create strain directory and copy templates and parameter file into it.

		@param directory_strains: Directory for the simulated strains
		@type directory_strains: str | unicode

		@return: Nothing
		@rtype: None
		"""
		if not os.path.exists(directory_strains):
			os.mkdir(directory_strains)
		for filename in self._directory_template_filenames:
			src = os.path.join(self._directory_template, filename)
			dst = os.path.join(directory_strains, filename)
			shutil.copy(src, dst)

	@staticmethod
	def get_genome_id_to_amounts(list_of_drawn_genome_id, genome_amounts):
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

	def simulate_strains(
		self, meta_table, genome_id_to_amounts, genome_id_to_file_path_genome, genome_id_to_file_path_gff=None):
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
		assert isinstance(meta_table, MetadataTable)
		assert isinstance(genome_id_to_amounts, dict)
		assert isinstance(genome_id_to_file_path_genome, dict)
		assert genome_id_to_file_path_gff is None or isinstance(genome_id_to_file_path_gff, dict)
		if genome_id_to_file_path_gff is None:
			msg = "No gff file (gene annotation) was given. Simulating strains without such a file can break genes."
			self._logger.warning(msg)
		for file_path in genome_id_to_file_path_genome.values():
			self.validate_file(file_path)
		if genome_id_to_file_path_gff is not None:
			for file_path in genome_id_to_file_path_gff.values():
				self.validate_file(file_path)
		self._simulate_strains(genome_id_to_amounts, genome_id_to_file_path_genome, genome_id_to_file_path_gff)
		self._pick_random_strains(meta_table, genome_id_to_amounts, genome_id_to_file_path_genome)

		# read file and generate strain diversity for each assembly
		# then subsample the strains
	def _simulate_strains(self, genome_id_to_amounts, genome_id_to_file_path_genome, genome_id_to_file_path_gff=None):
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
		tasks = []
		file_path_empty_file = None
		if genome_id_to_file_path_gff is None:
			file_path_empty_file = self.get_full_path(tempfile.mktemp(dir=self._tmp_dir))
			touch(file_path_empty_file)

		genome_id_to_file_path_genome_copy = genome_id_to_file_path_genome.copy()
		for genome_id in genome_id_to_file_path_genome_copy.keys():
			if self._keep_original and genome_id_to_amounts[genome_id] == 1:
				continue
			directory_strain = self._directory_strain.format(gid=genome_id)
			self._prepare_simulation_subfolder(directory_strain)
			file_path_genome = genome_id_to_file_path_genome[genome_id]
			if genome_id_to_file_path_gff is None:
				file_path_gff = file_path_empty_file
			else:
				file_path_gff = genome_id_to_file_path_gff[genome_id]
			self._logger.info("Simulating strain evolution of '{}'".format(genome_id))
			tasks.append(
				TaskCmd(self._get_simulate_cmd(
					directory_strains=directory_strain,
					filepath_genome=file_path_genome,
					filepath_gff=file_path_gff)))
		list_of_fails = runCmdParallel(tasks, maxProc=self._max_processors)

		if file_path_empty_file is not None:
			if os.path.exists(file_path_empty_file):
				os.remove(file_path_empty_file)

		if list_of_fails is not None:
			for message in reportFailedCmd(list_of_fails):
				self._logger.error(message)
			msg = "Simulation of strains failed."
			self._logger.error(msg)
			raise OSError(msg)

	def _pick_random_strains(self, meta_table, genome_id_to_amounts, genome_id_to_file_path_genome):
		"""
		Add randomly picked strains to genome_id_to_file_path_genome and metadata table.

		@param meta_table: Metadata table containing genome information
		@type meta_table: MetadataTable
		@param genome_id_to_file_path_genome:
		@type genome_id_to_file_path_genome: dict[str, str]
		@param genome_id_to_amounts:
		@type genome_id_to_amounts: dict[str, int]

		@return: Nothing
		@rtype: None
		"""
		assert isinstance(meta_table, MetadataTable)

		genome_id_to_file_path_genome_copy2 = genome_id_to_file_path_genome.copy()
		for genome_id in genome_id_to_file_path_genome_copy2.keys():
			if self._keep_original and genome_id_to_amounts[genome_id] == 1:
				continue
			directory_strain = self._directory_strain.format(gid=genome_id)

			amount = genome_id_to_amounts[genome_id]
			if self._keep_original and amount == 1:
				return

			if self._keep_original:
				amount -= 1
			else:
				genome_id_to_file_path_genome.pop(genome_id)

			genome_taxid = meta_table.get_cell_value(self._column_name_gid, genome_id, self._column_name_ncbi)
			sample = random.sample(range(0, len(self._filenames_strains)), amount)
			for index in sample:
				filename = self._filenames_strains[index]
				name, ext = os.path.splitext(filename)
				# index = name.split("Taxon")[1]
				new_id = "{prefix}{id}.{index}".format(prefix=self._filename_prefix, id=genome_id, index=name)
				# new_id = "{prefix}{id}.{index}".format(prefix=self._filename_prefix, id=genome_id, index=index)
				source = os.path.join(directory_strain, filename)
				destination = os.path.join(directory_strain, new_id + ".fna")
				os.rename(source, destination)
				genome_id_to_file_path_genome[new_id] = destination

				row = meta_table.get_empty_row()
				row[self._column_name_gid] = new_id
				if self._column_name_source in row:
					row[self._column_name_source] = "simulated"
				if genome_taxid is None:
					sys.stderr.write("Bad genome_ID: {}\n".format(genome_id))
					genome_taxid = 1
				row[self._column_name_ncbi] = genome_taxid
				meta_table.insert_row(row)

	def get_filenames_strains(self, file_path_template_newick_tree):
		"""
		Get list of file names of simulated genomes by reading template newick tree

		@attention: 'ancestor' is assumed to be part of tree as original sequence and will not be included

		@param file_path_template_newick_tree: File path to newick file
		@type file_path_template_newick_tree: str | unicode

		@return: list of file names of simulated genomes
		@rtype: list[str|unicode]
		"""
		assert self.validate_file(file_path_template_newick_tree)
		list_of_filenames_strains = []
		tree = Phylo.read(file_path_template_newick_tree, 'newick')
		for leaf in tree.get_terminals():
			prefix = leaf.name
			if prefix.lower() == "ancestor":
				continue
			list_of_filenames_strains.append("{prefix}.fasta".format(prefix=prefix))
		return list_of_filenames_strains


def touch(file_path):
	file_handle = open(file_path, 'w')
	file_handle.close()
