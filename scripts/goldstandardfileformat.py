__original_author__ = 'majda'
__author__ = 'hofmann'
__version__ = '0.0.2'


import re
import argparse
import sys
import os
import io
from Bio import SeqIO


class GoldStandardFileFormat():

    _label = "GoldStandardFileFormat"
    _separator="\t"
    _column_name_gid = ""
    _column_name_ncbi = ""
    fixed_name = {}

    def __init__(
        self, column_name_gid="genome_ID", column_name_ncbi="NCBI_ID", separator='\t', logfile=None, verbose=True):
        """
            Constructor

            @attention:

            @param column_name_gid: Column name of genome identifier
            @type column_name_gid: str | unicode
            @param column_name_ncbi: Column name of taxonomy identifier
            @type column_name_ncbi: str | unicode
            @param separator: character separating values in file
            @type separator: str | unicode
            @param logfile: file handler or file path to a log file
            @type logfile: file | io.FileIO | StringIO.StringIO | str | unicode
            @param verbose: Not verbose means that only warnings and errors will be past to stream
            @type verbose: bool

            @return: None
            @rtype: None
        """
        self._column_name_gid = column_name_gid
        self._column_name_ncbi = column_name_ncbi
        self._separator = separator

    # ###############
    # genome location
    # ###############

    def get_dict_sequence_to_genome_id(self, file_path_genome_locations, project_dir, set_of_genome_id=None, nanosim_real_fastq=False, wgsim=False):
        """
            Get a map, sequence id to genome id from an abundance file.

            @attention: Assuming genome id in first column, sequence id in second column

            @param file_path_genome_locations: List of file paths of abundance files
            @type file_path_genome_locations: str | unicode
            @param set_of_genome_id: Limit sequences to this set of genome id
            @type set_of_genome_id: set[str | unicode]

            @return: Mapping of sequence id to genome id
            @rtype: dict[str | unicode, str | unicode]
        """
        #assert isinstance(file_path_genome_locations, str)
        #assert self.validate_file(file_path_genome_locations)
        unique_id_to_genome_file_path = self.get_dict_unique_id_to_genome_file_path(file_path_genome_locations)
        if set_of_genome_id is None:
            set_of_genome_id = set(unique_id_to_genome_file_path.keys())
        else:
            assert set(unique_id_to_genome_file_path.keys()).issuperset(set_of_genome_id)

        sequence_id_to_genome_id = {}
        for genome_id in set_of_genome_id:
            file_path_genome = unique_id_to_genome_file_path[genome_id]
            #assert self.validate_file(file_path_genome)

            if not os.path.isabs(file_path_genome):
                file_path_genome = os.path.join(project_dir,file_path_genome)

            for seq_record in SeqIO.parse(file_path_genome, "fasta"):

                # If fastq files are generated directly with nanosim, the sequence id does not contain the version of the sequence record anymore.
                # To still be able to print the version of the sequence record to the read mapping file, we create a dict here, that holds
                # the sequence id without the version as key and the whole id as value to be accessed later on.
                if(nanosim_real_fastq):
                    self.fixed_name[seq_record.id.split('.',1)[0].replace("_","-")] = seq_record.id

                sequence_id_to_genome_id[seq_record.id] = genome_id

        return sequence_id_to_genome_id

    def get_dict_unique_id_to_genome_file_path(self, file_path_mapping):
        """
            Get a map, original sequence name to anonymous sequence name from a mapping file.

            @attention: anonymous name in second column

            @param file_path_mapping: File path to mapping file
            @type file_path_mapping: str | unicode

            @return: Mapping of anonymous sequence name to original sequence name
            @rtype: dict[str | unicode, str | unicode]
        """
        #assert isinstance(file_path_mapping, str)

        table = self.read(file_path_mapping, separator=self._separator)
        dict_mapping = self.get_map(table, 0, 1)
        return dict_mapping

    # ###############
    # metadata file
    # ###############

    def get_dict_genome_id_to_tax_id(
        self, file_path_metadata):
        """
            Get a map, genome id to taxonomic id from a metadata file.

            @attention: "genome_ID" and "NCBI_ID" assumed default column names.

            @param file_path_metadata: File path to metadata file
            @type file_path_metadata: str | unicode

            @return: Mapping of  genome id to taxonomic id
            @rtype: dict[str | unicode, str | unicode]
        """
        #assert isinstance(file_path_metadata, str)
        assert isinstance(self._column_name_gid, str)
        assert isinstance(self._column_name_ncbi, str)
        assert isinstance(self._separator, str)

        table = self.read(file_path_metadata, separator=self._separator, column_names=True)
        dict_genome_id_to_tax_id = self.get_map(table, self._column_name_gid, self._column_name_ncbi)
        return dict_genome_id_to_tax_id

    # ###############
    # anonymisation mapping file
    # ###############

    def get_dict_anonymous_to_original_id(self, file_path_mapping):
        """
            Get a map, anonymous sequence name to original sequence name from a mapping file.

            @attention: anonymous name in second column

            @param file_path_mapping: File path to mapping file
            @type file_path_mapping: str | unicode

            @return: Mapping of anonymous sequence name to original sequence name
            @rtype: dict[str | unicode, str | unicode]
        """
        #assert isinstance(file_path_mapping, str)

        table = self.read(file_path_mapping, separator=self._separator)
        dict_mapping = self.get_map(table, 1, 0)
        return dict_mapping

    def get_dict_sequence_name_to_anonymous(self, file_path_mapping):
        """
            Get a map, original sequence name to anonymous sequence name from a mapping file.

            @attention: anonymous name in second column

            @param file_path_mapping: File path to mapping file
            @type file_path_mapping: str | unicode

            @return: Mapping of anonymous sequence name to original sequence name
            @rtype: dict[str | unicode, str | unicode]
        """
        #assert isinstance(file_path_mapping, str)

        table = self.read(file_path_mapping, separator=self._separator)
        dict_mapping = self.get_map(table, 0, 1)
        return dict_mapping

    # ###############
    # position file
    # ###############

    def get_dict_sequence_name_to_positions(self, list_of_sam_position_files):
        """
            Get a map, sequence name to list of starting position from a mapping file.

            @attention: First column are sequence names, second column the start position

            @param list_of_sam_position_files: List of sam  position files
            @type list_of_sam_position_files: list[str|unicode]

            @return: Mapping of sequence name to list of starting position
            @rtype: dict[str | unicode, list[long]]
        """
        assert isinstance(list_of_sam_position_files, list)

        dict_original_seq_pos = {}

        for sam_position_file in list_of_sam_position_files:
            table = self.read(sam_position_file, separator=self._separator)
            column_key = table.get_column(0)
            column_values = table.get_column(1)

            for index_row in range(len(column_key)):
                key = column_key[index_row]
                value = column_values[index_row]
                seq_without_index = key.split("-")[0]
                if seq_without_index not in dict_original_seq_pos:
                    dict_original_seq_pos[seq_without_index] = []
                dict_original_seq_pos[seq_without_index].append(value)
        return dict_original_seq_pos

    # ###############
    # write gold standard assembly mapping
    # ###############

    def write_gs_read_mapping(
        self, stream_output, dict_anonymous_to_read_id, dict_sequence_to_genome_id, dict_genome_id_to_tax_id, nanosim_real_fastq, wgsim):
        """
            Write the gold standard for every read

            @attention:

            @param stream_output: File path for the output to be written to
            @type stream_output: file | FileIO | StringIO
            @param dict_anonymous_to_read_id: Mapping of anonymous sequence name to original sequence name
            @type dict_anonymous_to_read_id: dict[str | unicode, str | unicode]
            @param dict_sequence_to_genome_id: Mapping of sequence id to genome id
            @type dict_sequence_to_genome_id: dict[str | unicode, str | unicode]
            @param dict_genome_id_to_tax_id: Mapping of  genome id to taxonomic id
            @type dict_genome_id_to_tax_id: dict[str | unicode, str | unicode]

            @return: Nothing
            @rtype: None
        """
        #assert self.is_stream(stream_output)
        assert isinstance(dict_anonymous_to_read_id, dict)
        assert isinstance(dict_sequence_to_genome_id, dict)
        assert isinstance(dict_genome_id_to_tax_id, dict)

        row_format = "{aid}\t{gid}\t{tid}\t{sid}\n"
        line = '#' + row_format.format(
            aid="anonymous_read_id",
            gid="genome_id",
            tid="tax_id",
            sid="read_id")
        stream_output.write(line)
        for anonymous_id in sorted(dict_anonymous_to_read_id):
            read_id = dict_anonymous_to_read_id[anonymous_id]

            if nanosim_real_fastq:
                # If fastq files are generated directly with nanosim, the sequence id does not any "-" but "_".
                tmp = read_id.split('_', 1)[0]
                # If fastq files are generated directly with nanosim, the sequence id does not contain the version of the sequence record anymore.
                # To still be able to print the version of the sequence record to the read mapping file, we retrieve the whole sequence id from the dict.
                sequence_id = self.fixed_name[tmp]
            elif wgsim:
                # Change in CAMISIM 2:
                # When using wgsim 1.0 installed via conda instead of wgsim 0.3.0 delivered with CAMISIM 1, this error occured:
                # ValueError: missing '-' reads2anonymous: CP001958.1_986000_986310_0:0:0_0:0:0_1167/1
                # This is because the read ID in the simulated wgsim reads do not contain any '-' anymore.
                if '_' in read_id:
                    sequence_id = read_id.split('_', 1)[0]
                else:
                    msg = "missing '_' reads2anonymous: {}\n".format(read_id)
                    #self._logger.error(msg)
                    raise ValueError(msg)
            else:    
                if '-' not in read_id:
                    msg = "missing '-' reads2anonymous: {}\n".format(read_id)
                    #self._logger.error(msg)
                    raise ValueError(msg)  
                sequence_id = read_id.rsplit('-', 1)[0]
            if sequence_id not in dict_sequence_to_genome_id:
                msg = "sequence_id '{}' not found in mapping\n".format(sequence_id)
                #self._logger.error(msg)
                raise KeyError(msg)
            genome_id = dict_sequence_to_genome_id[sequence_id]
            tax_id = dict_genome_id_to_tax_id[genome_id]
            # final_dict[anonymous_id]= (genome_id,meta_tax_id,seq_id)

            # For nanosim only print the sequence id with version number and index of the read.
            # For fastq files converted from fasta files generated with nanosim, the read id is already formatted in this way.
            if nanosim_real_fastq:
                read_id = sequence_id + "-" + read_id.split("_")[2]

            line = row_format.format(
                aid=anonymous_id,
                gid=genome_id,
                tid=tax_id,
                sid=read_id,
            )
            stream_output.write(line)

    def write_gsa_contig_mapping(
        self, stream_output,
        dict_sequence_name_to_anonymous, dict_original_seq_pos, dict_sequence_to_genome_id, dict_genome_id_to_tax_id):
        """
            Write the gold standard for every read

            @attention:

            @param stream_output: stream output to be written to
            @type stream_output: file | FileIO | StringIO
            @param dict_sequence_name_to_anonymous:
            @type dict_sequence_name_to_anonymous: dict[str | unicode, str | unicode]
            @param dict_original_seq_pos: Mapping of sequence name to list of starting position
            @type dict_original_seq_pos: dict[str | unicode, list[long]]
            @param dict_sequence_to_genome_id: Mapping of sequence id to genome id
            @type dict_sequence_to_genome_id: dict[str | unicode, str | unicode]
            @param dict_genome_id_to_tax_id: Mapping of  genome id to taxonomic id
            @type dict_genome_id_to_tax_id: dict[str | unicode, str | unicode]

            @return: Nothing
            @rtype: None
        """
        assert self.is_stream(stream_output)
        assert isinstance(dict_sequence_name_to_anonymous, dict)
        assert isinstance(dict_original_seq_pos, dict)
        assert isinstance(dict_sequence_to_genome_id, dict)
        assert isinstance(dict_genome_id_to_tax_id, dict)
        row_format = "{name}\t{genome_id}\t{tax_id}\t{seq_id}\t{count}\t{position_0}\t{position_1}\n"

        stream_output.write("#anonymous_contig_id\tgenome_id\ttax_id\tcontig_id\tnumber_reads\tstart_position\tend_position\n")

        for original_contig_id, anonymous_contig_id in dict_sequence_name_to_anonymous.items():
            seq_info = original_contig_id.strip().rsplit("_from_", 1)
            # print(seq_info)
            sequence_id = seq_info[0]
            # pos_start, pos_end = re.findall(r'\d+', seq_info[1])[:2]
            pos_start = int(seq_info[1].split("_", 1)[0])
            pos_end = int(seq_info[1].split("_to_", 1)[1].split("_", 1)[0])

            # check if read is in contig
            count = 0
            for number in dict_original_seq_pos[sequence_id]:
                if pos_start <= int(number) <= pos_end:
                    count += 1

            # write output
            if sequence_id not in dict_sequence_to_genome_id:
                msg = "Bad sequence_id '{}'\n".format(sequence_id)
                self._logger.error(msg)
                raise KeyError(msg)
            genome_id = dict_sequence_to_genome_id[sequence_id]
            tax_id = dict_genome_id_to_tax_id[genome_id]
            stream_output.write(row_format.format(
                name=anonymous_contig_id,
                genome_id=genome_id,
                tax_id=tax_id,
                seq_id=sequence_id,
                count=count,
                position_0=pos_start,
                position_1=pos_end)
                )

    def gs_contig_mapping(
        self, file_path_genome_locations, file_path_metadata, file_path_id_map, list_file_paths_read_positions,
        stream_output):
        """
            Write the gold standard for every read

            @attention:

            @param stream_output: Stream output to be written to
            @type stream_output: file | FileIO | StringIO
            @param file_path_genome_locations:
            @type file_path_genome_locations: str | unicode
            @param file_path_metadata: Metadata file path, "genome_ID" and "NCBI_ID" assumed default column names.
            @type file_path_metadata: str | unicode
            @param file_path_id_map:
            @type file_path_id_map: str | unicode
            @param list_file_paths_read_positions:
            @type list_file_paths_read_positions: list[str | unicode]

            @return: Nothing
            @rtype: None
        """
        dict_sequence_to_genome_id = self.get_dict_sequence_to_genome_id(file_path_genome_locations)
        dict_genome_id_to_tax_id = self.get_dict_genome_id_to_tax_id(file_path_metadata)
        dict_original_seq_pos = self.get_dict_sequence_name_to_positions(list_file_paths_read_positions)
        dict_sequence_name_to_anonymous = self.get_dict_sequence_name_to_anonymous(file_path_id_map)
        self.write_gsa_contig_mapping(
            stream_output, dict_sequence_name_to_anonymous, dict_original_seq_pos,
            dict_sequence_to_genome_id, dict_genome_id_to_tax_id)

    def gs_read_mapping(self, file_path_genome_locations, file_path_metadata, file_path_id_map, stream_output, project_dir, nanosim_real_fastq, wgsim):
        """
            Write the gold standard for every read

            @attention:

            @param stream_output: Stream output to be written to
            @type stream_output: file | FileIO | StringIO
            @param file_path_genome_locations:
            @type file_path_genome_locations: str | unicode
            @param file_path_metadata: Metadata file path, "genome_ID" and "NCBI_ID" assumed default column names.
            @type file_path_metadata: str | unicode
            @param file_path_id_map:get_dict_sequence_to_genome_id
            @type file_path_id_map: str | unicode

            @return: Nothing
            @rtype: None
        """
        dict_sequence_to_genome_id = self.get_dict_sequence_to_genome_id(file_path_genome_locations, project_dir, nanosim_real_fastq=nanosim_real_fastq, wgsim=wgsim)
        dict_genome_id_to_tax_id = self.get_dict_genome_id_to_tax_id(file_path_metadata)
        dict_anonymous_to_read_id = self.get_dict_anonymous_to_original_id(file_path_id_map)
        self.write_gs_read_mapping(
            stream_output, dict_anonymous_to_read_id, dict_sequence_to_genome_id, dict_genome_id_to_tax_id, nanosim_real_fastq=nanosim_real_fastq, wgsim=wgsim)


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
            separator = self._separator

        #assert isinstance(file_path, str)
        #assert self.validate_file(file_path)
        assert isinstance(separator, str)
        assert isinstance(comment_line, list)
        assert isinstance(column_names, bool)

        meta_table = {}

        if isinstance(file_path, (str, bytes, os.PathLike)):
            file_handler = open(file_path)

            with file_handler as file_handler:
                meta_table = self.process_lines(file_handler, meta_table, comment_line, separator, column_names)

        elif isinstance(file_path, io.TextIOWrapper):
            file_handler = file_path
            meta_table = self.process_lines(file_handler, meta_table, comment_line, separator, column_names)
        else:
            raise TypeError("Invalid type for file_path")

        return meta_table    

    def process_lines(self, file_handler, meta_table, comment_line, separator, column_names):

        #self._logger.info("Reading file: '{}'".format(file_path))

        # read column names
        if column_names:
            list_of_column_names = self._parse_column_names(file_handler, separator)
            for column_name in list_of_column_names:
                meta_table[column_name] = []

        # read rows
        row_count = 0
        for line in file_handler:
            row_count += 1
            row = line.rstrip('\n').rstrip('\r')
            if line[0] in comment_line or len(row) == 0:
                continue
            row_cells = row.split(separator)

            if column_names:
                number_of_columns = len(list_of_column_names)
            else:
                number_of_columns = 0

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

    def _parse_column_names(self, stream_input, separator):
            row = stream_input.readline().rstrip('\n').rstrip('\r')
            list_of_column_names = row.split(separator)
            assert (len(list_of_column_names) == len(set(list_of_column_names))), "Column names must be unique!"
            return list_of_column_names

    def get_map(self, meta_table, key_column_name, value_column_name, unique_key=True):
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
        assert self.has_column(meta_table, key_column_name), "Column '{}' not found!".format(key_column_name)
        assert self.has_column(meta_table, value_column_name), "Column '{}' not found!".format(value_column_name)

        if key_column_name not in meta_table:
            #self._logger.error("Column name '{}' not available!".format(key_column_name))
            return None
        if value_column_name not in meta_table:
            #self._logger.error("Column name '{}' not available!".format(value_column_name))
            return None
        new_map = {}
        if len(meta_table) < 2:
            return new_map
        row_keys = meta_table[key_column_name]
        row_values = meta_table[value_column_name]
        for index, key in enumerate(row_keys):
            if unique_key and key in new_map:
                msg = "Key column is not unique! Key: '{}'".format(key)
                #self._logger.error(msg)
                raise KeyError(msg)
            new_map[key] = row_values[index]
        return new_map     

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
		"-contig",
		help="map the contigs for the anonymous gsa",
		action="store_true",
		default=False)
    parser.add_argument(
		"-input",
		help="input file (e.g. file path to the temporary anonymous mapping), reads from std.in by default",
		action='store',
		type=argparse.FileType('r'),
		default=sys.stdin)
    parser.add_argument(
		"-genomes",
		help="file with genome locations",
		action='store',
		type=argparse.FileType('r'),
		default="")
    parser.add_argument(
		"-metadata",
		help="metadata file",
		action='store',
		type=argparse.FileType('r'),
		default="")
    parser.add_argument(
		"-out",
		help="output file",
		action='store',
		type=argparse.FileType('w'),
		default=sys.stdout)
    parser.add_argument(
		"-projectDir",
		help="the path to the project directory",
		action='store',
		default="")
    parser.add_argument(
		"-nanosim_real_fastq",
		help="read files in fastq format generated directly with nanosim 3",
		action="store_true",
		default=False)
    parser.add_argument(
		"-wgsim",
		help="read files in fastq format generated directly with wgsim",
		action="store_true",
		default=False)     
    options = parser.parse_args()

    input_file_stream = options.input
    file_path_genome_locations = options.genomes
    file_path_metadata = options.metadata
    stream_output = options.out
    project_dir = options.projectDir
    nanosim_real_fastq = options.nanosim_real_fastq
    wgsim = options.wgsim

    goldStandardFileFormat = GoldStandardFileFormat()

    goldStandardFileFormat.gs_read_mapping(file_path_genome_locations, file_path_metadata, input_file_stream, stream_output, project_dir, nanosim_real_fastq, wgsim)
