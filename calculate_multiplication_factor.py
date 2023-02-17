#!/usr/bin/env python
 
import os
import tempfile
import sys
import decimal
from Bio import SeqIO
from Bio.Seq import Seq
from numbers import Number
import string

formats = ["fasta", "fastq"]

alphabets = {
		"rna": ["GAUC", "GAUCRYWSMKHBVDN"],
		"dna": ["GATC", "GATCRYWSMKHBVDN", "GATCRYWSMKHBVDN"],
		"protein": ["ACDEFGHIKLMNPQRSTVWYBXZJUO", "ACDEFGHIKLMNPQRSTVWYBXZJUO"]
		}

legal_text_characters = string.printable

qformats = {
		"sanger": [0, 40, 33],
		"solexa": [-5, 40, 64],
		"illumina": [0, 41, 33]  # 1.8+
	}

sequence_indicators = {
		"fasta": ">",
		"fastq": "@"
		}

def validate_characters(text, legal_alphabet=string.printable, key=None, silent=False):
    """
    Validate that only legal characters are contained in a text

    @attention:

    @param text: Some string
    @type text: str | unicode
    @param legal_alphabet: String of legal characters
    @type legal_alphabet: str | unicode
    @param key: If True, no error message will be made
    @type key: str | None
    @param silent: If True, no error message will be made
    @type silent: bool

    @return: bool
    @rtype: bool
    """
    prefix = ""
    if key:
        prefix = "'{}' ".format(key)

    set_legal_alphabet = set(legal_alphabet)
    set_text = set(text)
    if not set_legal_alphabet.issuperset(set_text):
        if not silent:
            difference = set_text.difference(set_legal_alphabet)
            difference.discard(set_legal_alphabet)
            #self._logger.error("{}Invalid characters: '{}'".format(prefix, ", ".join(difference)))
        return False
    return True


def validate_sequence(sequence, key=None, silent=False):
    """
    Validate that the sequence has only valid characters

    @attention:

	@param sequence: sequence
	@type sequence: Seq
	@param key: If True, no error message will be made
	@type key: str | None
	@param silent: If True, no error message will be made
	@type silent: bool

	@return: True if valid
	@rtype: bool
	"""
    assert isinstance(sequence, Seq)
    assert isinstance(silent, bool)

    prefix = ""
    if key:
        prefix = "'{}' ".format(key)

    if not len(sequence) > 0:
        #if not silent:
            #self._logger.error("{}Empty sequence".format(prefix))
        return False

    if not validate_characters(
        sequence.upper(), key=key, silent=silent):
        return False
    return True

def validate_sequence_id(identifier, used_ids=None, key=None, silent=False):
		"""
			Validate that the sequence identifier has only valid characters

			@attention:

			@param identifier: sequence
			@type identifier: str | unicode
			@param used_ids: Set of used up ids, that should not be repeated
			@type used_ids: set
			@param key: If True, no error message will be made
			@type key: str | unicode | None
			@param silent: If True, no error message will be made
			@type silent: bool

			@return: True if valid
			@rtype: bool
		"""
		set_of_seq_id = used_ids
		if set_of_seq_id is None:
			set_of_seq_id = set()
		assert isinstance(set_of_seq_id, set)
		assert isinstance(identifier, str)
		assert isinstance(silent, bool)

		prefix = ""
		if key:
			prefix = "'{}' ".format(key)

		if not len(identifier) > 0:
			#if not silent:
				#self._logger.error("{}Missing sequence id".format(prefix))
			return False

		if not validate_characters(identifier, legal_alphabet=legal_text_characters, key=key, silent=silent):
			return False

		if identifier in set_of_seq_id:
			#if not silent:
				#self._logger.error("{}Repeated sequence id '{}'".format(
					#prefix, identifier))
			return False
		return True

def validate_sequence_description(description, key=None, silent=False):
	"""
	Validate that the sequence description has only valid characters

	@attention:

	@param description: sequence
	@type description: str | unicode
	@param key: If True, no error message will be made
	@type key: str | unicode | None
	@param silent: If True, no error message will be made
	@type silent: bool

	@return: True if valid
	@rtype: bool
"""
	assert isinstance(description, str)
	assert key is None or isinstance(key, str)
	assert isinstance(silent, bool)

	if not validate_characters(description, legal_alphabet=legal_text_characters, key=key, silent=silent):
		return False
	return True

def validate_sequence_quality(quality, qformat="Illumina", key=None, silent=False):
		"""
			Validate that the sequence description has only valid characters

			@attention:

			@param quality: quality of each letter
			@type quality: list[int]
			@param qformat: 'Illumina', 'Sanger', 'Solexa'
			@type qformat: str | unicode
			@param key: If True, no error message will be made
			@type key: str | unicode | None
			@param silent: If True, no error message will be made
			@type silent: bool

			@return: True if valid, else False
			@rtype: bool
		"""
		assert isinstance(quality, list)
		assert isinstance(silent, bool)
		assert isinstance(qformat, str)
		qformat = qformat.lower()
		assert qformat in qformats, "{} not in {}".format(qformat, qformats.keys())

		prefix = ""
		if key:
			prefix = "'{}' ".format(key)

		# offset = self._qformats[qformat][2]
		minimum = qformats[qformat][0]  # +offset
		maximum = qformats[qformat][1]  # +offset

		invalid_indexes = [
			"{}: '{}'".format(index, value) for index, value in enumerate(quality) if not minimum <= value <= maximum]
		if len(invalid_indexes) > 0:
			#if not silent:
				#self._logger.error("{}Invalid quality at: {}.".format(prefix, ", ".join(invalid_indexes)))
			return False
		return True


def validate_sequence_record(seq_record, set_of_seq_id, file_format, key=None, silent=False):
			result = True
			if not validate_sequence(seq_record.seq, key=key, silent=silent):
				result = False
			if not validate_sequence_id(seq_record.id, used_ids=set_of_seq_id, key=key, silent=silent):
				result = False
			set_of_seq_id.add(seq_record.id)
			if not validate_sequence_description(seq_record.description, key=key, silent=silent):
				result = False
			if file_format == "fastq":
				if not validate_sequence_quality(
					seq_record.letter_annotations["phred_quality"], key=key, silent=silent):
					result = False
			return result

def validate_file_start(stream_handle, file_format):
		"""
			Validate that a stream with sequences starts with the correct character

			@param stream_handle: Any kind of stream type
			@type stream_handle: file | io.FileIO | StringIO.StringIO
			@param file_format: Format of the file at the file_path provided. Valid: 'fasta', 'fastq'
			@type file_format: str | unicode

			@return: True if the first character is correct
			@rtype: bool
		"""
		#assert is_stream(stream_handle)
		assert isinstance(file_format, str)
		file_format = file_format.lower()
		assert file_format in formats

		sequence_indicator = sequence_indicators[file_format]

		head = stream_handle.read(1)
		stream_handle.seek(0)
		if not head:
			return False
		if not head.startswith(sequence_indicator):
			return False
		return True

def get_sequence_lengths(
	file_path, file_format, sequence_type, ambiguous, key=None, silent=False):
	"""
	Validate a file to be correctly formatted and have sequences of a minimum length

	@attentioget_multiplication_factorn: Currently only phred quality for fastq files

	@param file_path: Path to file containing sequences
	@type file_path: str | unicode
	@param file_format: Format of the file at the file_path provided. Valid: 'fasta', 'fastq'
	@type file_format: str | unicode
	@param sequence_type: Are the sequences DNA or RNA? Valid: 'rna', 'dna', 'protein'
	@type sequence_type: str | unicode
	@param ambiguous: True or False, DNA example for strict 'GATC',  ambiguous example 'GATCRYWSMKHBVDN'
	@type ambiguous: bool
	@param key: If True, no error message will be made
	@type key: str | None
	@param silent: If True, no error message will be made
	@type silent: bool

	@return: True if the file is correctly formatted
	@rtype: tuple[int|long, int|long]
	"""
	#assert self.validate_file(file_path)
	assert isinstance(file_format, str)
	file_format = file_format.lower()
	assert file_format in formats
	assert isinstance(sequence_type, str)
	sequence_type = sequence_type.lower()
	assert sequence_type in alphabets

	prefix = ""
	if key:
		prefix = "'{}' ".format(key)

	if ambiguous:
		alphabet = alphabets[sequence_type][1]
	else:
		alphabet = alphabets[sequence_type][0]

	set_of_seq_id = set()
	total_length = 0
	sequence_count = 0
	min_sequence_length = None
	with open(file_path) as file_handle:
		if not validate_file_start(file_handle, file_format):
			msg = "{}Invalid beginning of file '{}'.".format(prefix, os.path.basename(file_path))
			#self._logger.error(msg)
			raise IOError(msg)
		try:
			for seq_record in SeqIO.parse(file_handle, file_format):
				sequence_count += 1
				if not validate_sequence_record(seq_record, set_of_seq_id, file_format, key=None, silent=False):
					msg = "{}{}. sequence '{}' is invalid.".format(prefix, sequence_count, seq_record.id)
					#self._logger.error(msg)
					raise IOError(msg)
				if not min_sequence_length:
					min_sequence_length = len(seq_record.seq)
				total_length += len(seq_record.seq)
				if len(seq_record.seq) < min_sequence_length:
					min_sequence_length = len(seq_record.seq)
		except Exception as e:
			if not silent:
				msg = "{}Corrupt sequence in file '{}'.\nException: {}".format(prefix, os.path.basename(file_path), e.message)
				#self._logger.error(msg)
				raise IOError(msg)
			return False
	if sequence_count == 0:
		return 0, 0
	return min_sequence_length, total_length


def stream_sequences_of_min_length(
	stream_input, stream_output, sequence_min_length, file_format="fasta"):
	"""
	Stream sequences of a minimum length

	@attention file_format: Anything but 'fasta' is not supported, yet

	@param stream_input: input stream of sequence file
	@type stream_input: file | FileIO | StringIO
	@param stream_output: Output stream
	@type stream_output: file | FileIO | StringIO
	@param stream_output: Output stream mapping
	@type stream_output: file | FileIO | StringIO
	@param sequence_min_length: Minimum length of sequences
	@type sequence_min_length: int | long
	@param file_format: 'fasta' format by default.
	@type file_format: str | unicode

	@return: Total length of all sequences (base pairs)
	@rtype: int | long
	"""
	total_base_pairs = 0
	for seq_record in SeqIO.parse(stream_input, file_format):
		# remove description, else art illumina messes up sam format
		seq_record.description = ''
		if len(seq_record.seq) < sequence_min_length:
			#self._logger.debug("'{}', Removing short sequence '{}', length: {}".format(
				#os.path.basename(stream_input.name), seq_record.id, len(seq_record.seq)))
			continue
		stream_output.write(seq_record.format(file_format))
		total_base_pairs += len(seq_record.seq)
	return total_base_pairs


def remove_short_sequences(file_path, min_sequence_length, file_format="fasta"):
    """"
    Copies a genome with sequences shorter than a minimum removed.

    @param file_path: File genome id associated with the abundance of a genome
    @type file_path: str | unicode
    @param min_sequence_length: Minimum length of a sequence in base pairs
    @type min_sequence_length: int 
    @param file_format:
    @type file_format: str | unicode

    @return: File path of the genome with removed short sequences
    @rtype: str | unicode
    """
    #assert self.validate_file(file_path)
    assert isinstance(min_sequence_length, int), "Expected natural digit"
    assert isinstance(file_format, str), "Expected file format 'fasta'"

    file_path_output = tempfile.mktemp(tempfile.gettempdir())
    with open(file_path) as stream_input, open(file_path_output, 'w') as stream_output:
        total_base_pairs = stream_sequences_of_min_length(
            stream_input, stream_output,
            sequence_min_length=min_sequence_length,
            file_format=file_format
            )
        if total_base_pairs == 0:
            msg = "No valid sequences > {} found!".format(min_sequence_length)
            #self._logger.error(msg)
            raise IOError(msg)
    return file_path_output

def has_unique_columns(list_of_column_names):
		return len(list_of_column_names) == len(set(list_of_column_names))

def parse_column_names(stream_input, separator):
			row = stream_input.readline().rstrip('\n').rstrip('\r')
			list_of_column_names_from_file = row.split(separator)
			assert has_unique_columns(list_of_column_names_from_file), "Column names must be unique!"

def parse_stream(stream_input, separator, column_names=False, comment_line=None, as_list=True):
	"""
	Reading comma or tab separated values from a stream

	@param stream_input: stream
	@type stream_input: file | io.FileIO | StringIO.StringIO
	@param separator: default character assumed to separate values in a file
	@type separator: str | unicode
	@param column_names: True if column names available
	@type column_names: bool
	@param comment_line: character or list of character indication comment lines
	@type comment_line: str | unicode | list[str|unicode]
	@param as_list: If true lists are returned, else dicts.
	@param as_list: bool

	@return: Generator returning dictionary or list
	@rtype: generator[ dict[int|long|str|unicode, str|unicode] ]
	#
	"""
	if comment_line is None:
		comment_line = ['#']
	elif isinstance(comment_line, str):
		comment_line = [comment_line]

	#assert self.is_stream(stream_input)
	assert isinstance(separator, str)
	assert isinstance(comment_line, list)
	assert isinstance(column_names, bool)

	# read column names
	number_of_columns = 0
	list_of_column_names = []
	if column_names:
		list_of_column_names = parse_column_names(stream_input, separator)
		number_of_columns = len(list_of_column_names)

	# read rows
	line_count = 0
	for line in stream_input:
		line_count += 1
		row = line.rstrip('\n').rstrip('\r')
		if line[0] in comment_line or len(row) == 0:
			continue

		row_cells = row.split(separator)
		if not column_names and number_of_columns == 0:
			number_of_columns = len(row_cells)

		if number_of_columns != len(row_cells):
			msg = "Format error. Bad number of values in line {}".format(line_count)
			#self._logger.error(msg)
			raise ValueError(msg)

		if as_list:
			yield row_cells
			continue

		dict_row = dict()
		for index, value in enumerate(row_cells):
			if column_names:
				column_name = list_of_column_names[index]
			else:
				column_name = index
			dict_row[column_name] = row_cells[index].rstrip('\n').rstrip('\r')
		if number_of_columns == 0:
			list_of_column_names = sorted(dict_row.keys())
		yield dict_row


def parse_file(file_path, separator=None, column_names=False, comment_line=None, as_list=True):
    """
    Reading comma or tab separated values from a file

    @param file_path: path to file to be opened
    @type file_path: str | unicode
    @param separator: default character assumed to separate values in a file
    @type separator: str | unicode
    @param column_names: True if column names available
    @type column_names: bool
    @param comment_line: character or list of character indication comment lines
    @type comment_line: str | unicode | list[str|unicode]
    @param as_list: If true lists are returned, else dicts.
    @param as_list: bool

    @return: Generator of dictionary representing rows
    @rtype: generator[ dict[int|long|str|unicode, str|unicode] ]
    #
    """
    with open(file_path) as file_handler:
        for row in parse_stream(file_handler, '\t', column_names, comment_line, as_list):
            yield row


# read genome location file
def read_genome_location_file(file_path):
    """
    Read file with the file paths of gnomes

    @param file_path: File genome id associated with the file path of a genome
    @type file_path: str | unicode

    @return: Dictionary of genome id to file path
    @rtype: dict[str|unicode, str|unicode]
    """
    #self._logger.info('Reading genome location file')
    #assert self.validate_file(file_path)
    dict_id_file_path = {}

    iterator_distributions = parse_file(file_path, as_list=True)
    for genome_id, file_path_genome in iterator_distributions:
        assert genome_id != '', "Invalid genomid: '{}'".format(genome_id)
        assert file_path_genome != '', "Invalid file path: '{}'".format(genome_id)
        #assert self.validate_file(file_path_genome), "Invalid file path: '{}'".format(genome_id)

        # check uniqueness
        assert genome_id not in dict_id_file_path, "Genome '{}' not unique in the distribution file!".format(genome_id)
        dict_id_file_path[genome_id] = file_path_genome
    return dict_id_file_path


def get_multiplication_factor(
    dict_id_file_path, dict_id_abundance, total_size, min_sequence_length,
    file_format="fasta", sequence_type="dna", ambiguous=True):
    """
    A factor is calculated based on total size of a sample to calculate the required covered
    # coverage = abundance * factor
    Files will be validated while the length of sequences are determined.

    @attention min_sequence_length: Sequences that are shorter than the expected fragment size are removed

    @param dict_id_file_path: Dictionary of genome id to file path
    @type dict_id_file_path: dict[str|unicode, str|unicode]
    @param dict_id_abundance: Dictionary of genome id to abundance
    @type dict_id_abundance: dict[str|unicode, float]
    @param total_size: Size of sample in base pairs
    @type total_size: int 
    @param min_sequence_length: Minimum length of a sequence in base pairs
    @type min_sequence_length: int 
    @param file_format: fasta or fastq
    @type file_format: str | unicode
    @param sequence_type: dna or rna or protein
    @type sequence_type: str | unicode
    @param ambiguous: DNA example for strict 'GATC',  ambiguous example 'GATCRYWSMKHBVDN'
    @type ambiguous: bool

    @return: Factor abundances will be multiplied by
    @rtype: float
    """
    assert isinstance(dict_id_file_path, dict), "Expected dictionary, genome id as key, file path as value"
    assert isinstance(dict_id_abundance, dict), "Expected dictionary, genome id as key, abundance as value"
    assert isinstance(total_size, (float, int)), "Expected natural digit"
    assert isinstance(min_sequence_length, int), "Expected natural digit"
    assert isinstance(file_format, str), "Expected file format 'fasta'"
    assert isinstance(sequence_type, str), "Expected sequence type 'rna' or 'dna' or 'protein'"
    assert isinstance(ambiguous, bool)

    relative_size_total = 0
    for genome_id, abundance in dict_id_abundance.items():
        try:
            min_seq_length, genome_length = get_sequence_lengths(
                file_path=dict_id_file_path[genome_id],
                file_format=file_format,
                sequence_type=sequence_type,
                ambiguous=ambiguous,
                key=None,
                silent=False)

            if min_seq_length < min_sequence_length:
                #self._logger.info("Genome '{}' has sequences below minimum, creating filtered copy.".format(genome_id))
                new_file_path = remove_short_sequences(
                    dict_id_file_path[genome_id], min_sequence_length, file_format="fasta")
                dict_id_file_path[genome_id] = new_file_path
                #self._temporary_files.add(new_file_path)

        except IOError as e:
            #self._remove_temporary_files()
            raise e

        relative_size = abundance * genome_length
        relative_size_total += relative_size
    return total_size / float(relative_size_total)

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

def read_distribution_file(file_path):
    """
    Read file with the distribution of a sample

    @param file_path: File genome id associated with the abundance of a genome
    @type file_path: str | unicode

    @return: Dictionary of genome id to file path
    @rtype: dict[str|unicode, float]
    """
    #self._logger.info('Reading distribution file')
    #assert self.validate_file(file_path)
    dict_id_abundance = {}
    # dict_id_file_path = {}

    iterator_distributions = parse_file(file_path, as_list=True)
    # for genome_id, abundance, genome_length, file_path_genome in iterator_distributions:
    abundance_sum = 0.
    for genome_id, abundance in iterator_distributions:
        assert genome_id != '', "Invalid genom id: '{}'".format(genome_id)
        assert abundance != '', "Invalid abundance: '{}'".format(genome_id)
        abundance = float(abundance)
        assert validate_number(abundance, zero=True), "Invalid abundance: '{}'".format(genome_id)

        assert genome_id not in dict_id_abundance, "Genome '{}' not unique in the distribution file!".format(genome_id)
        dict_id_abundance[genome_id] = abundance
        abundance_sum += abundance
    dict_id_abundance = {x : dict_id_abundance[x]/abundance_sum for x in dict_id_abundance} # normalise to 1
    return dict_id_abundance

# main method and entry point of this script
if __name__ == "__main__":

    fragment_size_mean = int(sys.argv[1])
    fragment_size_standard_deviation = int(sys.argv[2])
    total_size = float(sys.argv[3])
    file_path_genome_locations = sys.argv[4]
    file_path_distribution = sys.argv[5]

    min_sequence_length = fragment_size_mean - fragment_size_standard_deviation

    dict_id_file_path = read_genome_location_file(file_path_genome_locations)

    dict_id_abundance = read_distribution_file(file_path_distribution)

    factor = get_multiplication_factor(
            dict_id_file_path, dict_id_abundance, total_size, min_sequence_length,
            file_format="fasta", sequence_type="dna", ambiguous=True)

    print(factor)