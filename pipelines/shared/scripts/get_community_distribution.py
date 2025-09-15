import sys
import random
import os
import zipfile
from numbers import Number
import gzip
import bz2
import tempfile

filename_distribution_comunity_joint = "distribution_{sample_index}.txt"
boolean_states = {
		'yes': True, 'true': True, 'on': True,
		'no': False, 'false': False, 'off': False,
		'y': True, 't': True, 'n': False, 'f': False}
number_of_rows = 0
file_extensions_compression = {
		".zip": "zip",
		# ".7z": "7z",
		".gz": "gz",
		".bz2": "bz2",
		}
modes = ['r', 'w']
_open = {
		"gz": gzip.open,
		"bz2": bz2.BZ2File,
		"zip": zipfile.ZipFile,
		# "7z": tarfile.open,
		None: open,
		}
list_of_column_names = []


def get_initial_list(size_of_population, number_of_samples):
		"""
			Get initial list with zero initialized

			@attention: Each list in the list contains the distribution value for each sample

			@param size_of_population: Amount of genomes or individuals
			@type size_of_population: int | long
			@param number_of_samples: Number of samples
			@type number_of_samples: int | long

			@return: A list of lists.
			@rtype: list[list[float]]
		"""
		assert isinstance(size_of_population, int)
		assert isinstance(number_of_samples, int)
		return [[0.0] * number_of_samples for _ in range(size_of_population)]

def add_initial_log_distribution(list_population, mu, sigma):
		"""
			Adding a initial distribution

			@attention: Values for first sample

			@param list_population: Main list for all distributions
			@type : list[list[float]]
			@param mu: Mean
			@type mu: float
			@param sigma: standard deviation
			@type sigma: float

			@return: Nothing
			@rtype: None
		"""
		assert isinstance(list_population, list)
		assert isinstance(mu, (float, int))
		assert isinstance(sigma, (float, int))
		for index in range(len(list_population)):
			list_population[index][0] = random.lognormvariate(mu, sigma)

def lt_zero(value):
		"""
			Prevent values <= 0

			@attention:

			@param value:
			@type value: float | int | long

			@return: value if > 0, else 0.001
			@rtype: float | int | long
		"""
		if value <= 0:
			# > 0 to prevent extinction
			return 0.001
		else:
			return value

def add_replicates(list_population, mu, sigma):
		"""
			Adding gaussian noise to the first drawn abundances

			@attention:

			@param list_population: Main list for all distributions
			@type : list[list[float]]
			@param mu: Mean
			@type mu: float
			@param sigma: standard deviation
			@type sigma: float

			@return: Nothing
			@rtype: None
		"""
		assert isinstance(list_population, list)
		assert isinstance(mu, (float, int))
		assert isinstance(sigma, (float, int))
		for index_p in range(len(list_population)):
			initial_log_distribution = list_population[index_p][0]
			for index_i in range(len(list_population[index_p])-1):
				list_population[index_p][index_i+1] = lt_zero(initial_log_distribution + random.gauss(mu, sigma))

def add_timeseries_gauss(list_population, mu, sigma):
		"""
			Adding gaussian noise sequentially to the previous sample

			@attention:

			@param list_population: Main list for all distributions
			@type : list[list[float]]
			@param mu: Mean
			@type mu: float
			@param sigma: standard deviation
			@type sigma: float

			@return: Nothing
			@rtype: None
		"""
		assert isinstance(list_population, list)
		assert isinstance(mu, (float, int))
		assert isinstance(sigma, (float, int))
		for index_p in range(len(list_population)):
			for index_i in range(len(list_population[index_p])-1):
				if list_population[index_p][index_i] > 0:
					list_population[index_p][index_i+1] = lt_zero(list_population[index_p][index_i] + random.gauss(mu, sigma))
				else:
					# extinction
					list_population[index_p][index_i+1] = 0.0

def add_timeseries_lognorm(list_population, mu, sigma):
		"""
			each abundance profile is produced by
			- draw new value from lognorm distribution
			- add old and new value and divide by 2

			@attention:

			@param list_population: Main list for all distributions
			@type : list[list[float]]
			@param mu: Mean
			@type mu: float
			@param sigma: standard deviation
			@type sigma: float

			@return: Nothing
			@rtype: None
		"""
		assert isinstance(list_population, list)
		assert isinstance(mu, (float, int))
		assert isinstance(sigma, (float, int))
		for index_p in range(len(list_population)):
			for index_i in range(len(list_population[index_p])-1):
				list_population[index_p][index_i+1] = (list_population[index_p][index_i] + random.lognormvariate(mu, sigma))/2

def add_differential(list_population, mu, sigma):
		"""
			Abundance is drawn independently from previous lognorm distributions

			@attention:

			@param list_population: Main list for all distributions
			@type : list[list[float]]
			@param mu: Mean
			@type mu: float
			@param sigma: standard deviation
			@type sigma: float

			@return: Nothing
			@rtype: None
		"""
		assert isinstance(list_population, list)
		assert isinstance(mu, (float, int))
		assert isinstance(sigma, (float, int))
		for index_p in range(len(list_population)):
			for index_i in range(len(list_population[index_p])-1):
				list_population[index_p][index_i+1] = random.lognormvariate(mu, sigma)

def is_boolean_state(word):
		"""
			Test for boolean state

			@param word: A word
			@type word: str | unicode

			@return: True if word is identified as an word equivalent to true or false
			@rtype: bool
		"""
		return str(word) in boolean_states

def get_boolean_state(word):
		"""
			Get boolean from word

			@param word: A word
			@type word: str | unicode

			@return: True if word is identified as an word equivalent to true
			@rtype: bool
		"""
		assert str(word) in boolean_states
		return boolean_states[str(word)]

def random_distribution_to_relative_abundance(list_population, precision=10):
		"""
			Replace random distributions with relative abundances

			@attention: limited to first 20

			@param list_population: Main list for all distributions
			@type list_population: list[list[float]]
			@param precision: Precision, numbers after decimal point
			@type precision: int
		"""
		number_of_samples = len(list_population[0])
		for index_i in range(number_of_samples):
			total_abundance = 0.0
			for index_p in range(len(list_population)):
				total_abundance += list_population[index_p][index_i]
			for index_p in range(len(list_population)):
				list_population[index_p][index_i] = round(list_population[index_p][index_i] / float(total_abundance), precision)


def get_lists_of_distributions(
		size_of_population, number_of_samples, modus, log_mu, log_sigma, gauss_mu=None, gauss_sigma=None,
		view_distribution=False):
		"""
			Get list of distributions of all samples

			@attention:

			@param size_of_population: Amount of genomes or individuals
			@type size_of_population: int | long
			@param number_of_samples: Number of samples
			@type number_of_samples: int | long
			@param modus: 'replicates', 'timeseries_normal','timeseries_lognormal', 'differential'
			@type modus: str
			@param log_mu: Mean for log
			@type log_mu: float
			@param log_sigma: standard deviation for log
			@type log_sigma: float
			@param gauss_mu: Mean for gauss
			@type gauss_mu: float
			@param gauss_sigma: standard deviation for gauss
			@type gauss_sigma: float

			@return: Main list for all distributions
			@rtype: list[list[float]]
		"""
		assert isinstance(size_of_population, int)
		assert isinstance(number_of_samples, int)
		assert isinstance(modus, str)
		assert isinstance(log_mu, (float, int))
		assert isinstance(log_sigma, (float, int))
		assert isinstance(gauss_mu, (float, int))
		assert isinstance(gauss_sigma, (float, int))
		if gauss_mu is None:
			gauss_mu = 0
		if gauss_sigma is None:
			# TODO: gauss sigma needs proper dependence of log sigma
			gauss_sigma = 3 * log_sigma

		list_population = get_initial_list(size_of_population, number_of_samples)
		while True:
			add_initial_log_distribution(list_population, log_mu, log_sigma)

			if modus == 'replicates':
				add_replicates(list_population, gauss_mu, gauss_sigma)
			elif modus == 'timeseries_normal':
				add_timeseries_gauss(list_population, gauss_mu, gauss_sigma)
			elif modus == 'timeseries_lognormal':
				add_timeseries_lognorm(list_population, log_mu, log_sigma)
			elif modus == 'differential':
				add_differential(list_population, log_mu, log_sigma)

			if not view_distribution:
				break

		random_distribution_to_relative_abundance(list_population)

		return list_population

def has_unique_columns(list_of_column_names):
		return len(list_of_column_names) == len(set(list_of_column_names))

def parse_column_names(stream_input, separator):
			row = stream_input.readline().rstrip('\n').rstrip('\r')
			list_of_column_names_from_file = row.split(separator)
			assert has_unique_columns(list_of_column_names_from_file), "Column names must be unique!"

def get_compression_type(file_path):
		"""
		Return compression type assumed by filename

		@param file_path: Path to file
		@type file_path: str | unicode

		@return: compression type, None if no compression
		@rtype: str | None
		"""

		assert isinstance(file_path, str)
		filename, extension = os.path.splitext(file_path)

		if extension == ".zip" and not zipfile.is_zipfile(file_path):
			return None

		if extension in file_extensions_compression:
			return file_extensions_compression[extension]
		else:
			return None

def validate_number(digit, minimum=None, maximum=None, zero=True, key=None, silent=False):
		"""
			Validate that a variable is a number within a specific range if given.

			@attention: valid minimum <= digit <= maximum

			@param digit: Any number such as int, float, long
			@type digit: Number
			@param minimum: valid minimum <= digit
			@type minimum: Number
			@param maximum: valid digit <= maximum
			@type maximum: Number
			@param zero: If 0 is to be excluded
			@type zero: bool
			@param silent: If True, no error message will be made
			@type silent: bool

			@return: bool
			@rtype: bool
		"""
		# TODO: digit >= -1 can not be properly tested yet
		assert isinstance(digit, Number), type(digit)

		prefix = ""
		if key:
			prefix = "'{}' ".format(key)

		if minimum and digit < minimum:
			#if not silent:
				#self._logger.error("{}Invalid digit, must be bigger than {}, but was {}".format(prefix, minimum, digit))
			return False

		if maximum and digit > maximum:
			#if not silent:
				#self._logger.error("{}Invalid digit, must be smaller than {}, but was {}".format(prefix, maximum, digit))
			return False

		if not zero and digit == 0:
			#if not silent:
				#self._logger.error("{}Invalid digit, must not be {}".format(prefix, digit))
			return False
		return True

def openy(file_path, mode='r', compresslevel=5, compression_type=None):
		"""
		Open a file for reading or writing

		@attention: When reading file and compression_type None, type will be guessed.

		@param file_path: Path to file
		@type file_path: str | unicode
		@param mode: mode a file is opened with. 'r' or 'w'
		@type mode: str | unicode
		@param compresslevel: Higher level is slower but likely smaller. 0-9, except zip 0-8.
		@type compresslevel: int
		@param compression_type: "zip", "gz", "bz2",
		@type compression_type: str | unicode

		@return: Return a file object
		@rtype: file
		"""
		assert mode in modes, "Unsupported mode '{}'.".format(mode)
		if compression_type is None:
			compression_type = get_compression_type(file_path)
		if mode == 'r':
			return _open[compression_type](file_path, mode=mode)
		elif compression_type == "gz":
			assert validate_number(compresslevel, minimum=0, maximum=9)
			return _open[compression_type](file_path, mode='w', compresslevel=compresslevel)
		elif compression_type == "bz2":
			assert validate_number(compresslevel, minimum=0, maximum=9)
			return _open[compression_type](file_path, mode='w', compresslevel=compresslevel)
		elif compression_type == "zip":
			assert validate_number(compresslevel, minimum=0, maximum=8)
			return _open[compression_type](file_path, mode='w', compression=compresslevel)

def read(file_path, separator=None, column_names=False, comment_line=None):
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
		global number_of_rows

		if comment_line is None:
			comment_line = ['#']
		elif isinstance(comment_line, str):
			comment_line = [comment_line]

		if separator is None:
			separator="\t"

		assert isinstance(file_path, str)
		#assert self.validate_file(file_path)
		assert isinstance(separator, str)
		assert isinstance(comment_line, list)
		assert isinstance(column_names, bool)

		with open(file_path, 'r') as file_handler:
			#self._logger.info("Reading file: '{}'".format(file_path))

			meta_table = {}

			# read column names
			#if column_names:
			#	list_of_column_names = parse_column_names(file_handler, separator)
			#	for column_name in list_of_column_names:
			#		meta_table[column_name] = []
			
			global list_of_column_names

			# read rows
			row_count = 0
			for line in file_handler:

				row_count += 1
				row = line.rstrip('\n').rstrip('\r')
				if line[0] in comment_line or len(row) == 0:
					continue
				number_of_rows += 1
				row_cells = row.split(separator)
				number_of_columns = len(list(list_of_column_names))
				if number_of_columns != 0 and number_of_columns != len(row_cells):
					msg = "Format error. Bad number of values in line {}".format(row_count)
					#self._logger.error(msg)
					raise ValueError(msg)
				for index, value in enumerate(row_cells):
					if column_names:
						column_name = list_of_column_names[index]
					else:
						column_name = index
						if column_name not in meta_table:
							meta_table[column_name] = []
					meta_table[column_name].append(row_cells[index].rstrip('\n').rstrip('\r'))
				if number_of_columns == 0:
					list_of_column_names = sorted(meta_table.keys())

			return meta_table

def has_column(meta_table, column_name):
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

def get_map(metadata_table, key_column_name, value_column_name, unique_key=True):
		"""
			Keep rows at key values of a column

			@attention:

			@param key_column_name: Column name
			@type key_column_name: str | unicode | int | long
			@param value_column_name: Column name
			@type value_column_name: str | unicode | int | long

			@return: map
			@rtype: dict[str|unicode, str|unicode]

			@raises: KeyError
		"""

		assert isinstance(key_column_name, (str, int))
		assert isinstance(value_column_name, (str, int))
		assert has_column(metadata_table, key_column_name), "Column '{}' not found!".format(key_column_name)
		assert has_column(metadata_table, value_column_name), "Column '{}' not found!".format(value_column_name)

		if key_column_name not in metadata_table:
			# ??? self._logger.error("Column name '{}' not available!".format(key_column_name))
			return None
		if value_column_name not in metadata_table:
			# ??? self._logger.error("Column name '{}' not available!".format(value_column_name))
			return None
		new_map = {}
		if len(metadata_table) < 2:
			return new_map
		row_keys = metadata_table[key_column_name]
		row_values = metadata_table[value_column_name]
		for index, key in enumerate(row_keys):
			if unique_key and key in new_map:
				msg = "Key column is not unique! Key: '{}'".format(key)
				# ??? self._logger.error(msg)
				raise KeyError(msg)
			new_map[key] = row_values[index]
		return new_map

def get_genome_id_to_path_map(file_path_of_file_mapping_genome_id_to_paths):
		"""
		Get a dictionary mapping genome id to the path of their genome

		@param file_path_of_file_mapping_genome_id_to_paths: File path to file with format 'id \t path'
		@type file_path_of_file_mapping_genome_id_to_paths: str | unicode
		@param list_of_drawn_genome_id: List of genome identifiers
		@type list_of_drawn_genome_id: list[str|unicode]

		@return: genome ids mapped to their gnome file path
		@rtype: dict[str|unicode, str|unicode]
		"""
		genome_id_to_path_map = {}
		metadat_table = read(file_path_of_file_mapping_genome_id_to_paths)
		
		#if mdt.get_number_of_rows() > 0:
		if metadat_table:
			genome_id_to_path_map = get_map(metadat_table, 0, 1, unique_key=True)
		#msg = "'{}' is missing one or more genome id".format(os.path.basename(file_path_of_file_mapping_genome_id_to_paths))
		#assert set(genome_id_to_path_map.keys()).issuperset(list_of_drawn_genome_id), msg
		#return {genome_id: genome_id_to_path_map[genome_id] for genome_id in list_of_drawn_genome_id}
		return genome_id_to_path_map

def write_distribution_file(stream_out, genome_id_to_abundance, sample_index):
		"""
		Write abundance file for each sample

		@param stream_out: Output stream
		@type stream_out: file | FileIO | StringIO
		@param genome_id_to_abundance: Drawn distribution for each genome id
		@type genome_id_to_abundance: dict[str|unicode, list[float]]
		"""
		
		for genome_id in genome_id_to_abundance:
			distribution = [str(abundance) for abundance in genome_id_to_abundance[genome_id]][sample_index-1]
			stream_out.write("{id}\t{distr}\n".format(id=genome_id, distr=distribution))

def str_to_bool(s):
	if s == 'True':
		return True
	elif s == 'False':
		return False
	else:
		raise ValueError("Flag 'verbose' must either be set to True or False.")

def get_distribution_file_paths(directory, number_of_samples):
        """
        Generate directory paths for each sample

        @param directory: Output stream
        @type directory: str | unicode
        @param number_of_samples: Number of samples
        @type number_of_samples: int | long

        @return: list of directories
        @rtype: list[str | unicode]
        """
        file_path = os.path.join(directory, "distribution_{sample_index}.txt")
        return [file_path.format(sample_index=sample_index) for sample_index in range(number_of_samples)]

# TODO mehrere communities
# main method and entry point of this script
if __name__ == "__main__":

	number_of_samples = int(sys.argv[1])
	file_path_of_drawn_genome_location = sys.argv[2]
	mode = sys.argv[3]
	log_mu = int(sys.argv[4])
	log_sigma = int(sys.argv[5])
	gauss_mu = int(sys.argv[6])
	gauss_sigma = int(sys.argv[7])
	verbose = str_to_bool(sys.argv[8])
	seed = int(sys.argv[9])

	if seed is not None:
			random.seed(seed)

	genome_id_to_path_map = get_genome_id_to_path_map(file_path_of_drawn_genome_location)

	list_of_drawn_genome_id = genome_id_to_path_map.keys()

	directory_out_distributions = "./"
	list_of_file_paths_distribution = get_distribution_file_paths(
            directory_out_distributions, number_of_samples)

	list_of_distributions = get_lists_of_distributions(
			size_of_population=len(list_of_drawn_genome_id),
			number_of_samples=number_of_samples,
			modus=mode,
			log_mu=log_mu, log_sigma=log_sigma,
			gauss_mu=gauss_mu, gauss_sigma=gauss_sigma,
			view_distribution=verbose
		)

	assert len(list_of_drawn_genome_id) == len(list_of_distributions)
	genome_id_to_distributions = dict(zip(list_of_drawn_genome_id, list_of_distributions))

	file_path_distributions = tempfile.mktemp(tempfile.gettempdir())
	#file_path_distributions = "./file.txt"
	for sample_index, file_path_distribution in enumerate(list_of_file_paths_distribution):
		with open(file_path_distribution, 'w') as stream_out:
			write_distribution_file(stream_out=stream_out, genome_id_to_abundance=genome_id_to_distributions, sample_index=sample_index)