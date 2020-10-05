__author__ = 'hofmann'
__version__ = '0.0.3.1'

import os
import subprocess
import shutil
import tempfile
from scripts.parallel import TaskCmd, runCmdParallel, reportFailedCmd
from scripts.Validator.validator import Validator


class SamtoolsWrapper(Validator):

	_label = "SamtoolsWrapper"

	_sam_file_extension = ".sam"
	_bam_file_extension = ".bam"

	def __init__(self, file_path_samtools="samtools", max_processes=1, max_memory=1, compression_level=5,tmp_dir=None, logfile=None, verbose=True, debug=False):
		"""
			Collection of Methods to accomplish samtools tasks

			@attention: samtools version >= 1.0 expected

			@param file_path_samtools: path to the samtools executable
			@type file_path_samtools: str | unicode
			@param max_processes: Maximum number of processes used in parallel
			@type max_processes: int | long
			@param max_memory: Maximum available Memory (RAM) in gigabyte
			@type max_memory: int | long
			@param compression_level: Compression level used. 0-9 (0 means no compression)
			@type compression_level: int | long
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
		super(SamtoolsWrapper, self).__init__(logfile, verbose, debug)

		if tmp_dir is None:
			tmp_dir = tempfile.gettempdir()

		assert isinstance(tmp_dir, str)
		assert isinstance(verbose, bool), "Verbose must be true or false"
		assert isinstance(max_processes, int), "'max_processes' must be a digit"
		assert isinstance(max_memory, int), "'max_memory' must be a digit"
		assert isinstance(compression_level, int), "'compression_level' must be a digit"
		self.validate_number(max_processes, zero=False, minimum=1)
		self.validate_number(max_memory, zero=False, minimum=1)
		self.validate_number(compression_level, zero=False, minimum=0, maximum=9)
		assert self.validate_dir(tmp_dir)
		assert self.validate_file(file_path=file_path_samtools, executable=True)

		self._tmp_dir = tmp_dir
		self._file_path_samtools = file_path_samtools
		self._max_processes = max_processes
		self._max_memory = max_memory
		self._compression_level = compression_level

	# #######################################################
	#
	# 				Convert sam to bam
	#
	# #######################################################

	def _get_sam_to_bam_cmd(self, file_path_sam, output_dir, max_memory=-1):
		"""
			Return system command as string.
			Command will create a sorted by position and indexed bam file from a sam file.

			@attention:

			@param file_path_sam: file path
			@type file_path_sam: str | unicode
			@param output_dir: output directory
			@type output_dir: str | unicode
			@param max_memory: maximum available memory in gigabyte
			@type max_memory: int | long

			@return: system command
			@rtype: str
		"""
		if max_memory == -1:
			max_memory = self._max_memory
		file_name = os.path.splitext(os.path.basename(file_path_sam))[0]
		file_path_bam = os.path.join(output_dir, file_name)
		# cmd = "{samtools} view -bS {input} | {samtools} sort - {output}; {samtools} index {output}.bam"
		prefix_temp_files = tempfile.mktemp(dir=self._tmp_dir, prefix="temp_sam_to_sorted_bam")

		cmd_stream_sam_file = "{samtools} view -bS {input}"
		cmd_sort_bam_file = "{samtools} sort -l {compression} -m {memory}G -o {output}.bam -O bam -T {prefix}"
		cmd_index_bam_file = "{samtools} index {output}.bam"

		cmd = cmd_stream_sam_file + " | " + cmd_sort_bam_file + "; " + cmd_index_bam_file
		return cmd.format(
			samtools=self._file_path_samtools,
			input=file_path_sam,
			compression=self._compression_level,
			memory=max_memory,
			output=file_path_bam,
			prefix=prefix_temp_files
			)

	def convert_sam_to_bam(self, directory_sam, output_dir="./"):
		"""
			Converts all SAM-files in current directory to BAM-Format

			@attention:

			@param directory_sam: directory or file path
			@type directory_sam: str | unicode
			@param output_dir: output directory
			@type output_dir: str | unicode

			@return: None
			@rtype: None

			@raises: AssertionError
		"""
		directory_sam = self.get_full_path(directory_sam)
		output_dir = self.get_full_path(output_dir)
		sam_is_file = self.validate_file(directory_sam, silent=True)
		sam_is_folder = self.validate_dir(directory_sam, silent=True)
		bam_is_folder = self.validate_dir(output_dir, silent=True)
		assert sam_is_file or sam_is_folder, "Invalid file or directory: '{}'".format(directory_sam)
		assert bam_is_folder, "Invalid file or directory: '{}'".format(output_dir)

		sam_argument = [directory_sam]
		if sam_is_folder:
			sam_argument = self.get_files_in_directory(directory_sam, self._sam_file_extension)

		self.convert_sam_to_bam_by_list(sam_argument, output_dir)

	def convert_sam_to_bam_by_list(self, list_of_sam_files, output_dir="./"):
		"""
			Converts all SAM-files in current directory to BAM-Format

			@attention:

			@param list_of_sam_files: list of sam file paths
			@type list_of_sam_files: list[str|unicode]
			@param output_dir: output directory
			@type output_dir: str | unicode

			@return: None
			@rtype: None

			@raises: OSError | AssertionError
		"""
		bam_is_folder = self.validate_dir(output_dir, silent=True)
		assert isinstance(list_of_sam_files, list), "Expected list of file paths"
		assert bam_is_folder, "Invalid file or directory: '{}'".format(output_dir)
		# add commands to a list of tasks to run them in parallel
		tasks = []
		# cmd = "{exe} {sam} {name}"
		for sam_file_path in list_of_sam_files:
			cmd = self._get_sam_to_bam_cmd(sam_file_path, output_dir)
			tasks.append(TaskCmd(cmd))

		fail_list = runCmdParallel(tasks, maxProc=self._max_processes)
		if fail_list is not None:
			for message in reportFailedCmd(fail_list):
				self._logger.error(message)
			msg = "Converting sam files to bam files failed."
			self._logger.error(msg)
			raise OSError(msg)

	# #######################################################
	#
	# 				Merge bam files
	#
	# #######################################################

	def _get_merge_bam_cmd(self, list_of_file_paths, file_name_output, max_memory=-1):
		"""
			Return system command as string.
			Command will create a sorted and indexed bam file from a sam file.

			@attention:

			@param file_name_output: file path with no file extension '.bam'
			@type file_name_output: str | unicode
			@param list_of_file_paths: list of bam file paths
			@type list_of_file_paths: list[str|unicode]

			@return: system command
			@rtype: str
		"""
		if max_memory == -1:
			max_memory = self._max_memory
		prefix_temp_files = tempfile.mktemp(dir=self._tmp_dir, prefix="temp_merged_to_sorted_bam")

		cmd_merge_bam_files = "{samtools} merge -u - '{input_list}'"
		cmd_sort_bam_file = "{samtools} sort -l {compression} -m {memory}G -o {output}.bam -O bam -T {prefix}"
		cmd_index_bam_file = "{samtools} index {output}.bam"
		cmd = cmd_merge_bam_files + " | " + cmd_sort_bam_file + "; " + cmd_index_bam_file

		# cmd = "{samtools} merge - '{input_list}' | {samtools} sort - {output}; {samtools} index {output}.bam"
		return cmd.format(
			samtools=self._file_path_samtools,
			input_list="' '".join(list_of_file_paths),
			compression=self._compression_level,
			memory=max_memory,
			output=file_name_output,
			prefix=prefix_temp_files
			)

	def merge_bam_files_by_dict(self, dict_of_bam_files, output_dir):
		"""
			Merge lists of bam files into one.

			@attention: dictionary keys used as file names

			@param dict_of_bam_files: dictionary list of bam file paths as value
			@type dict_of_bam_files: dict[str|unicode, list[str|unicode]]
			@param output_dir: output directory
			@type output_dir: str | unicode

			@return: None
			@rtype: None

			@raises: OSError | AssertionError
		"""
		output_dir = self.get_full_path(output_dir)
		bam_is_folder = self.validate_dir(output_dir, silent=True)
		assert isinstance(dict_of_bam_files, dict), "Expected dictionary of file paths"
		assert bam_is_folder, "Invalid file or directory: '{}'".format(output_dir)
		for key, list_of_bam_paths in dict_of_bam_files.items():
			for file_path in list_of_bam_paths:
				assert self.validate_file(file_path), "Invalid file: '{}'".format(file_path)

		# add commands to a list of tasks to run them in parallel
		tasks = []
		for filename, list_of_bam_paths in dict_of_bam_files.items():
			if len(list_of_bam_paths) == 1:
				# move bam instead of merge, if only one
				file_path = list_of_bam_paths[0]
				self._logger.warning("List contains only one file: '{}'".format(file_path))
				out_file_path = os.path.join(output_dir, filename+self._bam_file_extension)
				shutil.copy2(file_path, out_file_path)
				continue
			cmd = self._get_merge_bam_cmd(list_of_bam_paths, os.path.join(output_dir, filename))
			tasks.append(TaskCmd(cmd))
		fail_list = runCmdParallel(tasks, maxProc=self._max_processes)
		if fail_list is not None:
			for message in reportFailedCmd(fail_list):
				self._logger.error(message)
			msg = "Converting sam files to bam files failed."
			self._logger.error(msg)
			raise OSError(msg)

	def merge_bam_files_by_list_of_dir(self, list_of_dir, output_dir):
		"""
			Merge lists of bam files with same prefix found in several directories.

			@attention:

			@param list_of_dir: list of directories containing bam files
			@type list_of_dir: list[str|unicode]
			@param output_dir: output directory
			@type output_dir: str | unicode

			@return: None
			@rtype: None

			@raises: AssertionError
		"""
		assert isinstance(list_of_dir, list)
		assert isinstance(output_dir, str)
		assert self.validate_dir(output_dir)
		for directory in list_of_dir:
			assert self.validate_dir(directory)

		self._logger.info("Merging bam files.")
		dict_of_bam_files = dict()
		for directory in list_of_dir:
			list_of_file_paths = self.get_files_in_directory(directory, ".bam")
			for file_path in list_of_file_paths:
				filename, extension = os.path.splitext(os.path.basename(file_path))
				if filename not in dict_of_bam_files:
					dict_of_bam_files[filename] = []
				dict_of_bam_files[filename].append(file_path)
		self.merge_bam_files_by_dict(dict_of_bam_files, output_dir)

	# #######################################################
	#
	# 		Parse 'read' start positions from sam/bam
	#
	# #######################################################

	def read_start_positions_from_list_of_sam(self, list_of_file_paths, output_file=None):
		"""
			Parse 'read' start positions from sam files.

			@attention:

			@param list_of_file_paths: list of sam file paths
			@type list_of_file_paths: list[str|unicode]
			@param output_file: output file path
			@type output_file: str | unicode

			@return: output file path
			@rtype: str | unicode

			@raises: AssertionError
		"""
		if output_file is None:
			output_file = tempfile.mktemp(dir=self._tmp_dir, prefix="start_positions")
		assert isinstance(list_of_file_paths, list)
		assert isinstance(output_file, str)
		assert self.validate_dir(output_file, only_parent=True)
		for file_path in list_of_file_paths:
			assert self.validate_file(file_path)

		with open(output_file, 'w') as file_handle_write:
			for file_path in list_of_file_paths:
				with open(file_path, 'r') as file_handle_read:
					for line in file_handle_read:
						line = line.rstrip()
						if line.startswith('@'):
							# header is ignored
							continue
						line_explode = line.split('\t')
						# seq_ID = line_explode[0]
						# positon = line_explode[3]
						file_handle_write.write(line_explode[0]+'\t'+line_explode[3]+'\n')
		return output_file

	def read_start_positions_from_list_of_bam(self, list_of_file_paths, output_file=None):
		"""
			Parse 'read' start positions from bam files.

			@attention:

			@param list_of_file_paths: list of sam file paths
			@type list_of_file_paths: list[str|unicode]
			@param output_file: output file path
			@type output_file: str | unicode

			@return: output file path
			@rtype: str | unicode

			@raises: AssertionError | OSError
		"""
		if output_file is None:
			output_file = tempfile.mktemp(dir=self._tmp_dir, prefix="read_start_positions")
		assert isinstance(list_of_file_paths, list)
		assert isinstance(output_file, str)
		assert self.validate_dir(output_file, only_parent=True)
		for file_path in list_of_file_paths:
			assert self.validate_file(file_path)

		cmd = "set -o pipefail; {samtools} view '{bamfile}' | awk '{{print $1 \"\\t\" $4}}' >> '{output}'"
		for file_path in list_of_file_paths:
			# exit_status = os.system(
			exit_status = subprocess.call(
				cmd.format(samtools=self._file_path_samtools, bamfile=file_path, output=output_file),
				shell=True,
				executable="bash")
			if exit_status != 0:
				msg = "Error occurred parsing '{}'".format(file_path)
				self._logger.error(msg)
				raise OSError(msg)
		return output_file

	def read_start_positions_from_dir_of_sam(self, directory, output_file=None):
		"""
			Parse 'read' start positions from sam files in a directory.

			@attention: caller of method needs to take care of output_file

			@param directory: directory containing sam files
			@type directory: str | unicode
			@param output_file: output file path
			@type output_file: str | unicode

			@return: output file path
			@rtype: str | unicode

			@raises: AssertionError
		"""
		assert output_file is None or isinstance(output_file, str)
		assert isinstance(directory, str)
		assert self.validate_dir(directory)

		list_of_file_paths = self.get_files_in_directory(directory, self._sam_file_extension)
		return self.read_start_positions_from_list_of_sam(list_of_file_paths, output_file)

	def read_start_positions_from_dir_of_bam(self, directory, output_file=None):
		"""
			Parse 'read' start positions from bam files in a directory.

			@attention: caller of method needs to take care of output_file

			@param directory: directory containing bam files
			@type directory: str | unicode
			@param output_file: output file path
			@type output_file: str | unicode

			@return: output file path
			@rtype: str | unicode

			@raises: AssertionError
		"""
		assert output_file is None or isinstance(output_file, str)
		assert isinstance(directory, str)
		assert self.validate_dir(directory)

		list_of_file_paths = self.get_files_in_directory(directory, self._bam_file_extension)
		return self.read_start_positions_from_list_of_bam(list_of_file_paths, output_file)
