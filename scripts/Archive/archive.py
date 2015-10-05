__author__ = 'hofmann'
__version__ = '0.0.2'

import io
import StringIO
from compress import Compress
import tarfile


class Archive(Compress):
	"""Reading and writing archived and/or compressed files"""

	_label = "Archive"

	_modes = {
		'r': {
			"gz": "r:gz",
			"bz2": "r:bz2",
			"tar": "r:",
			},
		'w': {
			"gz": "w:gz",
			"bz2": "w:bz2",
			"tar": "w:",
			}
		}

	_archive_mode_read_stream = {
		"gz": "r|gz",
		"bz2": "r|bz2",
		"tar": "r|",
		}

	_archive_mode_write_stream = {
		"gz": "w|gz",
		"bz2": "w|bz2",
		"tar": "w|",
		}

	_file_extensions_tar = ".tar"

	def __init__(self, default_compression="gz", logfile=None, verbose=True):
		"""
			Constructor

			@attention:

			@param default_compression: default compression used for files
			@type default_compression: str | unicode
			@param logfile: file handler or file path to a log file
			@type logfile: file | io.FileIO | StringIO.StringIO | basestring
			@param verbose: Not verbose means that only warnings and errors will be past to stream
			@type verbose: bool

			@return: None
			@rtype: None
		"""
		assert logfile is None or isinstance(logfile, basestring) or self.is_stream(logfile)
		assert isinstance(default_compression, basestring), "separator must be string"
		assert isinstance(verbose, bool), "verbose must be true or false"
		assert default_compression.lower() in self._open, "Unknown compression: '{}'".format(default_compression)

		super(Archive, self).__init__(default_compression, logfile, verbose)

		self._open['tar'] = tarfile.open
		self._default_compression = default_compression

	@staticmethod
	def is_archive(file_path):
		"""
			Test if archive can be assumed by filename

			@param file_path: Path to file
			@type file_path: str | unicode

			@return: True if file is archive
			@rtype: str | None
		"""
		return tarfile.is_tarfile(file_path)

	def open_archive(self, file_path, compression_type=None, mode='r'):
		"""
			Test if archive can be assumed by filename

			@param file_path: Path to file
			@type file_path: str | unicode

			@return: True if stream
			@rtype: tarfile.TarFile
		"""
		assert mode in self._modes, "Unsupported mode".format(mode)
		if compression_type is None:
			compression_type = self.get_compression_type(file_path)
		assert compression_type in self._modes[mode], "Unsupported compression '{}' for archive files.".format(
			compression_type)
		assert self.is_archive(file_path)

		if compression_type is None:
			compression_type = 'tar'

		mode = self._modes[mode][compression_type]
		return self._open[compression_type](file_path, mode=mode)
