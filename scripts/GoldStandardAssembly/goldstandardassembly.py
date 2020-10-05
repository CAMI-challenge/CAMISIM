__author__ = 'hofmann'
__version__ = '0.0.5'


import os
# import time
import tempfile
import shutil
import scripts
from .samtoolswrapper import SamtoolsWrapper


class GoldStandardAssembly(SamtoolsWrapper):

    _label = "GoldStandardAssembly"

    _list_of_reference_file_extension = [".fna", ".fasta"]

    def __init__(
        self, file_path_samtools="samtools", max_processes=1, tmp_dir=None, logfile=None, verbose=True, debug=False):
        """
            Collection of Methods related to gold standard assemblies

            @attention:

            @param file_path_samtools: path to the samtools executable
            @type file_path_samtools: str | unicode
            @param max_processes: Maximum number of processes used in parallel
            @type max_processes: int | long
            @param tmp_dir: Temp directory for temporary data if needed
            @type tmp_dir: str | unicode
            @param logfile: file handler or file path to a log file
            @type logfile: file | io.FileIO | StringIO.StringIO | str
            @param verbose: Not verbose means that only warnings and errors will be past to stream
            @type verbose: bool
            @param debug: Display debug messages
            @type debug: bool

            @return: None
            @rtype: None
        """
        super(GoldStandardAssembly, self).__init__(
            file_path_samtools=file_path_samtools,
            max_processes=max_processes,
            tmp_dir=tmp_dir,
            logfile=logfile, verbose=verbose, debug=debug
        )
        self._temp_merges_bam_directory = tempfile.mkdtemp(dir=self._tmp_dir)
        self._bamToGold = os.path.join(os.path.dirname(scripts.__file__), "bamToGold.pl")
        assert self.validate_file(self._bamToGold)

    def _close(self):
        if self.validate_dir(self._temp_merges_bam_directory, silent=True):
            shutil.rmtree(self._temp_merges_bam_directory)
        self._logger = None

    def bam_reads_to_contigs(
        self, file_path_bam, file_path_fasta_ref, file_path_output,
        min_length=1, min_coverage=1, file_path_log=None):
        """
            Convert a bam file to a gold standard assembly.

            @attention:

            @param file_path_fasta_ref: path to reference fasta file
            @type file_path_fasta_ref: str | unicode
            @param file_path_bam: path to bam files containing reads based on the reference fasta
            @type file_path_bam: str | unicode
            @param file_path_output: output fasta file path
            @type file_path_output: str | unicode
            @param min_length: Minimum length of sequence.
            @type min_length: int
            @param min_coverage: Minimum coverage at a site of a sequence.
            @type min_coverage: int
            @param file_path_log: output path for a log file
            @type file_path_log: str | unicode

            @return: None
            @rtype: None
        """
        assert self.validate_file(file_path_fasta_ref)
        assert self.validate_file(file_path_bam)
        assert self.validate_dir(file_path_output, only_parent=True)
        if file_path_log is None:
            file_path_log = os.path.splitext(file_path_output)[0] + ".log"
        assert self.validate_dir(file_path_log, only_parent=True)
        cmd = "{exe} -st '{samtools}' -r '{ref}' -b '{bam}' -l '{min_length}' -c '{min_coverage}' >> '{output}' 2>> '{logfile}'".format(
            exe=self._bamToGold,
            samtools=self._file_path_samtools,
            ref=file_path_fasta_ref,
            bam=file_path_bam,
            min_length=min_length,
            min_coverage=min_coverage,
            output=file_path_output,
            logfile=file_path_log,
            )
        if self._debug:
            self._logger.debug(cmd)
        exit_status = os.system(cmd)
        if not exit_status == 0:
            msg = "Error occurred converting '{}'\n{}".format(os.path.basename(file_path_bam), cmd)
            self._logger.error(msg)
            raise OSError(msg)

    def get_dict_id_to_file_path_bam_from_dir(self, directory):
        """
            Get a dictionary with unique id to bam file path, parsed from files in a directory

            @attention: assumes bam files to be named after unique genome id

            @param directory: directory path containing bam files
            @type directory: str | unicode

            @return: dictionary of list of bam files by a key based on filename
            @rtype: dict[str|unicode, str|unicode]
        """
        dict_id_to_filepath_bam = {}
        list_of_folder_file_paths = self.get_files_in_directory(directory, self._bam_file_extension)
        for file_path_bam in list_of_folder_file_paths:
            key = os.path.splitext(os.path.basename(file_path_bam))[0]
            if key in dict_id_to_filepath_bam:
                msg = "Can not establish unique key '{}'".format(key)
                self._logger.error(msg)
                raise KeyError(msg)
            dict_id_to_filepath_bam[key] = file_path_bam
        return dict_id_to_filepath_bam

    def get_dict_id_to_file_path_reference_from_dir(self, directory):
        """
            Get a dictionary with unique id to '.fna', '.fasta' file path, parsed from files in a directory

            @attention:

            @param directory: directory path containing reference files
            @type directory: str | unicode

            @return: dictionary of list of reference files by a key based on filename
            @rtype: dict[str|unicode, str|unicode]
        """
        dict_id_to_filepath_bam = {}
        list_of_folder_file_paths = []
        for file_extension in self._list_of_reference_file_extension:
            list_of_folder_file_paths.extend(self.get_files_in_directory(directory, file_extension))

        for file_path_bam in list_of_folder_file_paths:
            key = os.path.splitext(os.path.basename(file_path_bam))[0]
            if key in dict_id_to_filepath_bam:
                msg = "Can not establish unique key '{}'".format(key)
                self._logger.error(msg)
                raise KeyError(msg)
            dict_id_to_filepath_bam[key] = file_path_bam

        return dict_id_to_filepath_bam

    def pooled_gold_standard_by_dir(self, list_of_directory_bam, dict_id_to_file_path_fasta, file_path_output=None):
        """
            Make a gold standard assembly merging bam files of several samples

            @attention: bam files must have same name to be merged

            @param list_of_directory_bam: list of directories containing bam files
            @type list_of_directory_bam: list[str|unicode]
            @param dict_id_to_file_path_fasta: path to reference files by key
            @type dict_id_to_file_path_fasta: dict[str|unicode, str|unicode]

            @return: output file path
            @rtype: str | unicode
        """
        if file_path_output is None:
            file_path_output = tempfile.mktemp(dir=self._tmp_dir)
        self._logger.info("Creating pooled gold standard")
        self.merge_bam_files_by_list_of_dir(list_of_directory_bam, output_dir=self._temp_merges_bam_directory)
        dict_id_to_file_path_bam = self.get_dict_id_to_file_path_bam_from_dir(self._temp_merges_bam_directory)
        return self.gold_standard_assembly(dict_id_to_file_path_bam, dict_id_to_file_path_fasta, file_path_output)

    def gold_standard_assembly(self, dict_id_to_file_path_bam, dict_id_to_file_path_fasta, file_path_output=None):
        """
            Make a gold standard assembly using bam files of a samples

            @attention:

            @param dict_id_to_file_path_bam: path to reference fasta files by key
            @type dict_id_to_file_path_bam: dict[str|unicode, str|unicode]
            @param dict_id_to_file_path_fasta: path to reference files by key
            @type dict_id_to_file_path_fasta: dict[str|unicode, str|unicode]
            @param file_path_output: output fasta file path
            @type file_path_output: str | unicode

            @return: output file path
            @rtype: str | unicode
        """
        if file_path_output is None:
            file_path_output = tempfile.mktemp(dir=self._tmp_dir)
        assert isinstance(dict_id_to_file_path_bam, dict)
        assert isinstance(dict_id_to_file_path_fasta, dict)
        assert len(dict_id_to_file_path_bam) != 0, "Empty bam file list"
        assert len(dict_id_to_file_path_fasta) != 0, "Empty reference fasta file list"
        set_fasta = set(dict_id_to_file_path_fasta.keys())
        set_bam = set(dict_id_to_file_path_bam.keys())
        is_superset = set_fasta.issuperset(set_bam)
        set_bam.discard(set_fasta)
        assert is_superset, "Key error: '{}'".format(", ".join(set_bam))
        file_path_output = self.get_full_path(file_path_output)
        file_path_output = self.get_available_file_path(file_path_output)

        # create all contigs
        for key, file_path_bam in dict_id_to_file_path_bam.items():
            file_path_fasta_ref = dict_id_to_file_path_fasta[key]
            self.bam_reads_to_contigs(
                file_path_fasta_ref=file_path_fasta_ref,
                file_path_bam=file_path_bam,
                file_path_output=file_path_output)
        return file_path_output
