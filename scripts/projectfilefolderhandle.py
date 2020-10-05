__author__ = 'hofmann'

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

	# Validate raw genomes
	# 0 Community design	TRUE
	# 1 Move Genomes		TRUE
	# 2 Simulate reads.
	# 		Copy distribution file / genome location file to project folder if directly given
	# 		Input project folder, unless previous steps done
	# 		Output tmp, to be archived
	# 3 Gold standard assembly
	# 	!!	input tmp, unless previous steps done
	# 		output tmp, to be archived
	# 4 Anonymization
	# 	!!	input tmp, unless previous steps done
	# 		output tmp, to be archived
	# 5 Archive
	# 		- Compression strength option?
	# 		- requires list of all that is to be archived
	# 	!!	input tmp, unless previous steps done
	# 		output project folder

	# 0 disabled
	# 1 discard
	# 2 keep
	# 3 compress

	_TMP = True
	_HardDrive = False
	# (a) reads / (a) gsa / (a) pgsa
	_location_reads = [_HardDrive, _HardDrive]
	_location_gsa = [_HardDrive, _HardDrive]
	_location_pgsa = [_HardDrive, _HardDrive]

	# ###################
	#   sub folder_names
	# ###################

	_folder_name_internal = "internal"
	# _folder_name_comunity_design = "comunity_design"
	_folder_name_distribution = "distributions"
	_folder_name_genomes = "source_genomes"
	# _folder_name_meta_data = "meta_data"
	# folder_name_simulated = "simulated_genomes"
	_folder_name_bam = "bam"
	# _folder_name_sam = "sam"
	_folder_name_reads = "reads"
	_folder_name_contigs = "contigs"
	# _folder_name_logfiles = "logfiles"
	_folder_name_sample = "sample_{id}"

	_sub_folders_sample = [_folder_name_bam, _folder_name_reads, _folder_name_contigs]
	_sub_folders_output = [_folder_name_internal, _folder_name_distribution, _folder_name_genomes]

	# ###################
	#   file names
	# ###################

	_filename_genome_locations = "genome_locations.tsv"
	_filename_distribution = "distribution.txt"

	_filename_anonymous_reads = "anonymous_reads.fq"
	# filename_reads_anonymous_mapping = "reads_anonymous_mapping.tsv"

	_filename_gsa = "gsa.fasta"
	_filename_anonymous_gsa = "anonymous_gsa.fasta"
	# filename_gsa_anonymous_mapping = "gsa_anonymous_mapping.tsv"

	_filename_gsa_pooled = "gsa_pooled.fasta"
	_filename_anonymous_gsa_pooled = "anonymous_gsa_pooled.fasta"
	# filename_pooled_gsa_mapping = "pooled_" + filename_gsa_anonymous_mapping

	# filename_gsa = "gsa.fasta"
	# filename_pooled_gsa = "pooled_" + filename_gsa
	_filename_reads_mapping = "reads_mapping.tsv"
	_filename_gsa_mapping = "gsa_mapping.tsv"
	_filename_pooled_gsa_mapping = "gsa_pooled_mapping.tsv"

	_filename_log = "pipeline.log"
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
		self._tmp_dir = tempfile.mkdtemp(dir=tmp_dir)
		self._directory_output = output_dir
		self._time_stamp = time_stamp
		if time_stamp is None:
			self._time_stamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y.%m.%d_%H.%M.%S')
		super(ProjectFileFolderHandle, self).__init__(logfile, verbose, debug)

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

	def make_directory_structure(self, number_of_samples):
		"""
		Create folder structure at output and temporary location

		@param number_of_samples: Number of samples.
		@type number_of_samples: int | long

		@return: Nothing
		@rtype: None
		"""
		assert isinstance(number_of_samples, int)
		self.make_directory_temp_structure(number_of_samples)
		self.make_directory_output_structure(number_of_samples)

	def make_directory_temp_structure(self, number_of_samples):
		"""
		Create folder structure at temporary location

		@param number_of_samples: Number of samples.
		@type number_of_samples: int | long

		@return: Nothing
		@rtype: None
		"""
		assert isinstance(number_of_samples, int)
		self._make_directory_structure(self._TMP, number_of_samples)

	def make_directory_output_structure(self, number_of_samples):
		"""
		Create folder structure at output location

		@param number_of_samples: Number of samples.
		@type number_of_samples: int | long

		@return: Nothing
		@rtype: None
		"""
		assert isinstance(number_of_samples, int)
		self._make_directory_structure(self._HardDrive, number_of_samples)

	def _make_directory_structure(self, is_tmp, number_of_samples):
		"""
		Create folder structure at temporary location

		@param is_tmp: Location where the directory structure is to be created.
		@type is_tmp: bool
		@param number_of_samples: Number of samples.
		@type number_of_samples: int | long

		@return: Nothing
		@rtype: None
		"""
		assert isinstance(is_tmp, bool)
		assert isinstance(number_of_samples, int)
		dir_main = self._get_root_directory(is_tmp)
		self._make_dir(dir_main)
		for sub_folder in self._sub_folders_output:
			directory = os.path.join(dir_main, sub_folder)
			self._make_dir(directory)
		for sample_index in range(number_of_samples):
			dir_sample = self.get_sample_dir(is_tmp, str(sample_index))
			self._make_dir(dir_sample)
			for sub_folder in self._sub_folders_sample:
				sub_directory = os.path.join(dir_sample, sub_folder)
				self._make_dir(sub_directory)

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

	def get_tmp_wd(self):
		"""
		Get location of temporary working directory.

		@return: temporary working directory
		@rtype: str | unicode
		"""
		return self._tmp_dir

	def _get_root_directory(self, is_tmp):
		"""
		Get root directory baseed on whether it is at a temporary location or output location.

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

	def get_bam_dirs(self):
		"""
		Get list of bam directories of all samples

		@attention: The list includes previous runs!

		@return: List of bam directories
		@rtype: list[str|unicode]
		"""
		out_dir = self.get_output_directory()
		list_of_dirs = [
			os.path.join(out_dir, folder_name) for folder_name in os.listdir(out_dir)
			if os.path.isdir(os.path.join(out_dir, folder_name))]
		sample_dirs = sorted([
			directory for directory in list_of_dirs
			if self.validate_dir(directory, sub_directories=self._sub_folders_sample, silent=True)])
		return [os.path.join(sample_dir, self._folder_name_bam) for sample_dir in sample_dirs]

	def get_distribution_dir(self):
		"""
		Get directory where distribution files are located.

		@return: distribution directory
		@rtype: str | unicode
		"""
		root_dir = self._directory_output
		return os.path.join(root_dir, self._folder_name_distribution)

	def get_genome_dir(self):
		"""
		Get directory where genome files are located.

		@return: distribution directory
		@rtype: str | unicode
		"""
		root_dir = self._directory_output
		return os.path.join(root_dir, self._folder_name_genomes)

	def get_meta_data_dir(self):
		"""
		Get directory where metadata files are located.

		@return: metadata directory
		@rtype: str | unicode
		"""
		root_dir = self._directory_output
		return os.path.join(root_dir, self._folder_name_internal)

	def get_bam_dir(self, sample_id):
		"""
		Get directory where bam files are located.

		@type sample_id: str | unicode

		@return: bam directory
		@rtype: str | unicode
		"""
		assert isinstance(sample_id, str)
		sample_dir = self.get_sample_dir(self._HardDrive, sample_id)
		return os.path.join(sample_dir, self._folder_name_bam)

	def get_reads_dir(self, is_input, sample_id):
		"""
		Get directory where fastq files are located.

		@type is_input: bool
		@type sample_id: str | unicode

		@return: fastq directory
		@rtype: str | unicode
		"""
		assert isinstance(is_input, bool)
		assert isinstance(sample_id, str)

		if is_input:
			sample_dir = self.get_sample_dir(self._location_reads[0], sample_id)
		else:
			sample_dir = self.get_sample_dir(self._HardDrive, sample_id)
		return os.path.join(sample_dir, self._folder_name_reads)

	def get_contigs_dir(self, is_input, sample_id):
		"""
		Get directory where fastq files are located.

		@type is_input: bool
		@type sample_id: str | unicode

		@return: fastq directory
		@rtype: str | unicode
		"""
		assert isinstance(is_input, bool)
		assert isinstance(sample_id, str)

		if is_input:
			sample_dir = self.get_sample_dir(self._location_reads[0], sample_id)
		else:
			sample_dir = self.get_sample_dir(self._HardDrive, sample_id)
		return os.path.join(sample_dir, self._folder_name_contigs)

	def get_logfile_dir(self):
		"""
		Get directory where log files are located.

		@return: logfile directory
		@rtype: str | unicode
		"""
		root_dir = self._directory_output
		return os.path.join(root_dir, self._folder_name_internal)

	def get_sample_dir(self, is_tmp, sample_id):
		"""
		Get directory where sample files are located.

		@type is_tmp: bool
		@type sample_id: str | unicode

		@return: sample directory
		@rtype: str | unicode
		"""
		assert isinstance(is_tmp, bool)
		assert isinstance(sample_id, str)
		root_dir = self._get_root_directory(is_tmp)
		folder_name = "{}_{}".format(self._time_stamp, self._folder_name_sample.format(id=sample_id))
		return os.path.join(root_dir, folder_name)

	# ###################
	#   file paths
	# ###################

	def get_anonymous_gsa_pooled_file_path(self):
		"""
		Get file location of the gold standard assembly based on pooled sample reads.

		@return: file location of pooled gold standard assembly
		@rtype: str | unicode
		"""
		root_dir = self._get_root_directory(self._HardDrive)
		return os.path.join(
			root_dir, self._filename_anonymous_gsa_pooled)

	def get_gsa_pooled_file_path(self):
		"""
		Get file location of the gold standard assembly based on pooled sample reads.

		@return: file location of pooled gold standard assembly
		@rtype: str | unicode
		"""
		root_dir = self._get_root_directory(self._HardDrive)
		return os.path.join(
			root_dir, self._filename_gsa_pooled)

	def get_anonymous_gsa_pooled_map_file_path(self):
		"""
		Get file location of the anonymous gold standard assembly based on pooled sample reads.

		@return: file location of anonymous pooled gold standard assembly
		@rtype: str | unicode
		"""
		root_dir = self._get_root_directory(self._HardDrive)
		return os.path.join(
			root_dir, self._filename_pooled_gsa_mapping)

	def get_gsa_file_path(self, sample_id):
		"""
		Get file location of the anonymous gold standard assembly.

		@param sample_id: sample id of a sample
		@type sample_id: str | unicode

		@return: file location of anonymous gold standard assembly
		@rtype: str | unicode
		"""
		assert isinstance(sample_id, str)
		output_dir = self.get_contigs_dir(self._HardDrive, sample_id)
		return os.path.join(
			output_dir, self._filename_gsa)

	def get_anonymous_gsa_file_path(self, sample_id):
		"""
		Get file location of the anonymous gold standard assembly.

		@param sample_id: sample id of a sample
		@type sample_id: str | unicode

		@return: file location of anonymous gold standard assembly
		@rtype: str | unicode
		"""
		assert isinstance(sample_id, str)
		output_dir = self.get_contigs_dir(self._HardDrive, sample_id)
		return os.path.join(output_dir, self._filename_anonymous_gsa)

	def get_anonymous_gsa_map_file_path(self, sample_id):
		"""
		Get file location of the anonymous gold standard assembly mapping.

		@param sample_id: sample id of a sample
		@type sample_id: str | unicode

		@return: file location of anonymous gold standard assembly mapping
		@rtype: str | unicode
		"""
		assert isinstance(sample_id, str)
		output_dir = self.get_contigs_dir(self._HardDrive, sample_id)
		return os.path.join(output_dir, self._filename_gsa_mapping)

	def get_anonymous_reads_file_path(self, sample_id):
		"""
		Get file location of the anonymous gold standard assembly mapping.

		@param sample_id: sample id of a sample
		@type sample_id: str | unicode

		@return: file location of anonymous gold standard assembly mapping
		@rtype: str | unicode
		"""
		assert isinstance(sample_id, str)
		fastq_dir = self.get_reads_dir(self._HardDrive, sample_id)
		return os.path.join(fastq_dir, self._filename_anonymous_reads)

	def get_anonymous_reads_map_file_path(self, sample_id):
		"""
		Get file location of the anonymous reads mapping.

		@param sample_id: sample id of a sample
		@type sample_id: str | unicode

		@return: file location of anonymous reads mapping
		@rtype: str | unicode
		"""
		assert isinstance(sample_id, str)
		fastq_dir = self.get_reads_dir(self._HardDrive, sample_id)
		return os.path.join(fastq_dir, self._filename_reads_mapping)

	def get_distribution_file_path(self, sample_id):
		"""
		Get file location of a distribution file of a specific sample.

		@param sample_id: sample id of a sample
		@type sample_id: str | unicode

		@return: file location of distribution file
		@rtype: str | unicode
		"""
		assert isinstance(sample_id, str)
		return os.path.join(
			self.get_sample_dir(self._HardDrive, sample_id), self._filename_distribution)

	def get_distribution_file_path_list(self, number_of_samples):
		"""
		Get file locations of all distribution files.

		@param number_of_samples: Number of samples.
		@type number_of_samples: int | long

		@return: file location of distribution file
		@rtype: str | unicode
		"""
		assert isinstance(number_of_samples, int)
		return [self.get_distribution_file_path(str(sample_index)) for sample_index in range(number_of_samples)]

	def get_genome_location_file_path(self):
		"""
		Get file location of file containing genome locations by genome ids.

		@return: file location of file containing genome locations by genome ids.
		@rtype: str | unicode
		"""
		root_dir = self._directory_output
		return os.path.join(
			root_dir, self._folder_name_internal, self._filename_genome_locations)

	def get_log_file_path(self):
		"""
		Get logfile location.

		@return: logfile location.
		@rtype: str | unicode
		"""
		root_dir = self._directory_output
		return os.path.join(
			root_dir, self._folder_name_internal, self._filename_log)

	def get_genome_metadata_file_path(self):
		"""
		Get metadata file location.

		@return: metadata file location.
		@rtype: str | unicode
		"""
		root_dir = self._directory_output
		return os.path.join(
			root_dir, self._folder_name_internal, self._filename_metadata)
