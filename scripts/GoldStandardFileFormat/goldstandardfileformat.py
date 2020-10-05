__original_author__ = 'majda'
__author__ = 'hofmann'
__version__ = '0.0.2'


import re
from Bio import SeqIO
from scripts.MetaDataTable.metadatatable import MetadataTable
from scripts.Validator.validator import Validator


class GoldStandardFileFormat(Validator):

    _label = "GoldStandardFileFormat"

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
        super(GoldStandardFileFormat, self).__init__(logfile=logfile, verbose=verbose)
        self._column_name_gid = column_name_gid
        self._column_name_ncbi = column_name_ncbi
        self._separator = separator

    # ###############
    # genome location
    # ###############

    def get_dict_sequence_to_genome_id(self, file_path_genome_locations, set_of_genome_id=None):
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
        assert isinstance(file_path_genome_locations, str)
        assert self.validate_file(file_path_genome_locations)
        unique_id_to_genome_file_path = self.get_dict_unique_id_to_genome_file_path(file_path_genome_locations)
        if set_of_genome_id is None:
            set_of_genome_id = set(unique_id_to_genome_file_path.keys())
        else:
            assert set(unique_id_to_genome_file_path.keys()).issuperset(set_of_genome_id)

        sequence_id_to_genome_id = {}
        for genome_id in set_of_genome_id:
            file_path_genome = unique_id_to_genome_file_path[genome_id]
            assert self.validate_file(file_path_genome)
            for seq_record in SeqIO.parse(file_path_genome, "fasta"):
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
        assert isinstance(file_path_mapping, str)

        table = MetadataTable(logfile=self._logfile, verbose=self._verbose)
        table.read(file_path_mapping, separator=self._separator)
        dict_mapping = table.get_map(0, 1)
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
        assert isinstance(file_path_metadata, str)
        assert isinstance(self._column_name_gid, str)
        assert isinstance(self._column_name_ncbi, str)
        assert isinstance(self._separator, str)

        table = MetadataTable(logfile=self._logfile, verbose=self._verbose)
        table.read(file_path_metadata, separator=self._separator, column_names=True)
        dict_genome_id_to_tax_id = table.get_map(self._column_name_gid, self._column_name_ncbi)
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
        assert isinstance(file_path_mapping, str)

        table = MetadataTable(logfile=self._logfile, verbose=self._verbose)
        table.read(file_path_mapping, separator=self._separator)
        dict_mapping = table.get_map(1, 0)
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
        assert isinstance(file_path_mapping, str)

        table = MetadataTable(logfile=self._logfile, verbose=self._verbose)
        table.read(file_path_mapping, separator=self._separator)
        dict_mapping = table.get_map(0, 1)
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

        table = MetadataTable(logfile=self._logfile, verbose=self._verbose)
        for sam_position_file in list_of_sam_position_files:
            table.read(sam_position_file, separator=self._separator)
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
        self, stream_output, dict_anonymous_to_read_id, dict_sequence_to_genome_id, dict_genome_id_to_tax_id):
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
        assert self.is_stream(stream_output)
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
            if '-' not in read_id:
                msg = "missing '-' reads2anonymous: {}\n".format(read_id)
                self._logger.error(msg)
                raise ValueError(msg)
            sequence_id = read_id.rsplit('-', 1)[0]
            if sequence_id not in dict_sequence_to_genome_id:
                msg = "sequence_id '{}' not found in mapping\n".format(sequence_id)
                self._logger.error(msg)
                raise KeyError(msg)
            genome_id = dict_sequence_to_genome_id[sequence_id]
            tax_id = dict_genome_id_to_tax_id[genome_id]
            # final_dict[anonymous_id]= (genome_id,meta_tax_id,seq_id)
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

    def gs_read_mapping(self, file_path_genome_locations, file_path_metadata, file_path_id_map, stream_output):
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

            @return: Nothing
            @rtype: None
        """
        dict_sequence_to_genome_id = self.get_dict_sequence_to_genome_id(file_path_genome_locations)
        dict_genome_id_to_tax_id = self.get_dict_genome_id_to_tax_id(file_path_metadata)
        dict_anonymous_to_read_id = self.get_dict_anonymous_to_original_id(file_path_id_map)
        self.write_gs_read_mapping(
            stream_output, dict_anonymous_to_read_id, dict_sequence_to_genome_id, dict_genome_id_to_tax_id)
