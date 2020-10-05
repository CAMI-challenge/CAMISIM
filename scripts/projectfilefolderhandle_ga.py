__author__ = 'Peter Hofmann'

import os
import tempfile
import shutil
import datetime
import time
from scripts.Validator.validator import Validator


class ProjectFileFolderHandle(Validator):
	"""
	Dealing with file and folder locations related to the data produced
	"""

	_label = "ProjectFileFolderHandle"

	_TMP = True
	_HardDrive = False

	# ###################
	# sub folder names
	# ###################

	_folder_name_source = "source"
	# name of folder containing all log files
	# folder_name_logfiles = "logfiles"

	# ###################
	#   file names
	# ###################

	_filename_internal_id_map = "id_mapping.tsv"
	_suffix_16S = "16S_rRNA"
	# file_mg_05s = "5S_rRNA.fna"
	_filename_mg = "{}.fna"
	# file_mg_23s = "23S_rRNA.fna"
	_filename_cluster_mg = "mothur_cluster_{}.list"
	_filename_metadata = "meta_data.tsv"

	def __init__(self, tmp_dir, output_dir, time_stamp=None, logfile=None, verbose=True, debug=False):
		"""
		Constructor

		@param tmp_dir: Directory for temporary data
		@type tmp_dir: str | unicode
		@param output_dir: Directory where final data will be placed
		@type output_dir: str | unicode
		@param time_stamp: timestamp as string
		@type time_stamp: str | unicode
		@param logfile: file | FileIO | StringIO | str
		@param verbose: Not verbose means that only warnings and errors will be past to stream
		@type verbose: bool
		@param debug: Display debug messages
		@type debug: bool
		"""
		assert isinstance(tmp_dir, str)
		assert isinstance(output_dir, str)
		assert time_stamp is None or isinstance(time_stamp, str)
		assert self.validate_dir(tmp_dir)

		self._tmp_dir = tempfile.mkdtemp(dir=tmp_dir)
		self._directory_output = output_dir
		self._time_stamp = time_stamp
		if time_stamp is None:
			self._time_stamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y.%m.%d_%H.%M.%S')
		super(ProjectFileFolderHandle, self).__init__(logfile, verbose, debug)
		self._make_dir(output_dir)

	def get_time_stamp(self):
		return self._time_stamp

	def get_output_directory(self):
		"""
		Get directory where final data will be placed

		@return: Directory where final data will be placed
		@rtype: str | unicode
		"""
		return self._directory_output

	def remove_directory_temp(self):
		"""
		Delete temporary data

		@return: Nothing
		@rtype: None
		"""
		if os.path.exists(self._tmp_dir):
			assert os.path.isdir(self._tmp_dir)
			shutil.rmtree(self._tmp_dir)

	def _make_dir(self, directory):
		"""
		Create folder at given location, it it does not exists already.

		@param directory: Number of samples.
		@type directory: str | unicode

		@return: Nothing
		@rtype: None
		"""
		assert self.validate_dir(directory, only_parent=True)
		if os.path.exists(directory):
			assert os.path.isdir(directory)
		else:
			os.mkdir(directory)

	def _get_root_directory(self, is_tmp):
		"""
		Get root directory based on whether it is at a temporary location or output location.

		@type is_tmp: bool

		@return: temporary working directory
		@rtype: str | unicode
		"""
		if is_tmp:
			return self._tmp_dir
		else:
			return self._directory_output

	# ###################
	#   directories
	# ###################

	def get_tmp_wd(self):
		"""
		Get location of temporary working directory.

		@return: temporary working directory
		@rtype: str | unicode
		"""
		return self._tmp_dir

	# ###################
	#   file paths
	# ###################

	def get_file_path_internal_id_map(self):
		"""
		Get file location of the mapping for internal ids to original ids and taxonomic assignments.

		@return: File location of pooled gold standard assembly
		@rtype: str | unicode
		"""
		root_dir = self._get_root_directory(self._HardDrive)
		return os.path.join(root_dir, self._filename_internal_id_map)

	def get_file_path_meta_data_table(self):
		"""
		Get file location of the metadata for genomes

		@return: File location of pooled gold standard assembly
		@rtype: str | unicode
		"""
		root_dir = self._get_root_directory(self._HardDrive)
		return os.path.join(root_dir, self._filename_metadata)

	def get_file_path_mg_16s(self):
		"""
		Get file location of 16s sequences

		@return: File location of pooled gold standard assembly
		@rtype: str | unicode
		"""
		root_dir = self._get_root_directory(self._HardDrive)
		return os.path.join(root_dir, self._filename_mg.format(self._suffix_16S))

	def get_file_path_cluster_mg_16s(self):
		"""
		Get file location of 16s sequences

		@return: File location of pooled gold standard assembly
		@rtype: str | unicode
		"""
		root_dir = self._get_root_directory(self._HardDrive)
		return os.path.join(root_dir, self._filename_cluster_mg.format(self._suffix_16S))
