__author__ = 'hofmann'
__version__ = '0.0.6'

import os
import io
import time
import datetime
from scripts.Validator.validator import Validator
import gzip
import bz2
import zipfile
from scripts.parallel import TaskThread, runThreadParallel


class Compress(Validator):
    """Reading and writing compressed files"""

    _label = "Compress"

    _open = {
        "gz": gzip.open,
        "bz2": bz2.BZ2File,
        "zip": zipfile.ZipFile,
        # "7z": tarfile.open,
        None: open,
        }

    _file_extensions_compression = {
        ".zip": "zip",
        # ".7z": "7z",
        ".gz": "gz",
        ".bz2": "bz2",
        }

    _modes = ['r', 'w']

    def __init__(self, default_compression="gz", label="Compress", logfile=None, verbose=True, debug=False):
        """
        Constructor

        @attention:

        @param default_compression: default compression used for files
        @type default_compression: str | unicode
        @param logfile: file handler or file path to a log file
        @type logfile: file | io.FileIO | StringIO.StringIO | str
        @param verbose: Not verbose means that only warnings and errors will be past to stream
        @type verbose: bool
        @param debug: Display debug messages
        @type debug: bool

        @return: None
        @rtype: None
        """
        assert logfile is None or isinstance(logfile, str) or self.is_stream(logfile)
        assert isinstance(default_compression, str), "separator must be string"
        assert isinstance(verbose, bool), "verbose must be true or false"
        assert default_compression.lower() in self._open, "Unknown compression: '{}'".format(default_compression)

        super(Compress, self).__init__(label=label, logfile=logfile, verbose=verbose, debug=debug)

        # self._logger = LoggingWrapper(self._label, verbose=verbose)
        # if logfile is not None:
        #     self._logger.set_log_file(logfile)

        self._default_compression = default_compression

    def get_compression_type(self, file_path):
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

        if extension in self._file_extensions_compression:
            return self._file_extensions_compression[extension]
        else:
            return None

    def open(self, file_path, mode='r', compresslevel=5, compression_type=None):
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
        assert mode in self._modes, "Unsupported mode '{}'.".format(mode)
        if compression_type is None:
            compression_type = self.get_compression_type(file_path)
        if mode == 'r':
            return self._open[compression_type](file_path, mode=mode)
        elif compression_type == "gz":
            assert self.validate_number(compresslevel, minimum=0, maximum=9)
            return self._open[compression_type](file_path, mode='w', compresslevel=compresslevel)
        elif compression_type == "bz2":
            assert self.validate_number(compresslevel, minimum=0, maximum=9)
            return self._open[compression_type](file_path, mode='w', compresslevel=compresslevel)
        elif compression_type == "zip":
            assert self.validate_number(compresslevel, minimum=0, maximum=8)
            return self._open[compression_type](file_path, mode='w', compression=compresslevel)

    def compress_file(self, src, dst='./', compresslevel=5, compression_type=None, overwrite=False):
        """
        Compress a file

        @attention: When reading file and compression_type None, type will be guessed.

        @param src: Path to file
        @type src: str | unicode
        @param dst: Destination path, a directory or file path
        @type dst: str | unicode
        @param compresslevel: Higher level is slower but likely smaller. 0-9, except zip 0-8.
        @type compresslevel: int
        @param compression_type: "zip", "gz", "bz2",
        @type compression_type: str | unicode
        @param overwrite: If false, a path will renamed if not available
        @type overwrite: bool

        @return: True if stream
        @rtype: None
        """
        if compression_type is None:
            compression_type = self.get_compression_type(dst)
        if compression_type is None:
            compression_type = self._default_compression
        compression_type = compression_type.lower()
        assert compression_type in self._open, "Unknown compression type: '{}'".format(compression_type)

        self._logger.info("Compressing '{}'".format(os.path.basename(src)))
        time_start = time.time()

        src = self.get_full_path(src)
        dst = self.get_full_path(dst)
        if not self.validate_file(src) or not self.validate_dir(dst, only_parent=True):
            msg = "Failed compressing '{}'!".format(src)
            self._logger.error(msg)
            raise IOError(msg)

        if self.validate_dir(dst, silent=True):
            extension = ".{}".format(compression_type)
            dst = os.path.join(dst, os.path.basename(src) + extension)

        if not overwrite:
            dst = self.get_available_file_path(dst)

        with open(src, 'rb') as read_handler, self.open(dst, 'w', compresslevel, compression_type) as write_handler:
            write_handler.writelines(read_handler)

        time_end = time.time()
        time_elapsed = str(datetime.timedelta(seconds=round(time_end - time_start)))
        self._logger.info("Done compressing '{file}' in {time}s.".format(time=time_elapsed, file=os.path.basename(dst)))

    def compress_list_of_files(
        self, list_of_file_paths, dst, compresslevel=5,
        compression_type=None, overwrite=False, max_processors=1):
        """
        Compress list of files

        @attention: When reading file and compression_type None, type will be guessed.

        @param list_of_file_paths: Path to file
        @type list_of_file_paths: list[str|unicode]
        @param dst: Destination path, a directory or file path
        @type dst: str | unicode
        @param compresslevel: Higher level is slower but likely better. 0-9, except zip 0-8.
        @type compresslevel: int
        @param compression_type: "zip", "gz", "bz2",
        @type compression_type: str | unicode
        @param overwrite: If false, a path will renamed if not available
        @type overwrite: bool
        @param max_processors: Maximum number processors used for compressing files simultaneously
        @type max_processors: int

        @return: True if stream
        @rtype: None
        """
        if not self.validate_dir(dst, silent=True):
            assert self.validate_dir(dst, only_parent=True), "Bad destination: '{}'".format(dst)
        task_list = []
        for file_path in list_of_file_paths:
            if self.validate_file(file_path):
                args = (file_path, dst, compresslevel, compression_type, overwrite)
                task_list.append(TaskThread(_compress_file, args))
            else:
                msg = "File not found '{}'".format(file_path)
                self._logger.error(msg)
                raise IOError(msg)
        list_of_return_values = runThreadParallel(task_list, maxThreads=max_processors)
        for index, return_value in enumerate(list_of_return_values):
            assert return_value is None, "Compressing of '{}' failed. '{}'".format(list_of_file_paths[index], return_value)

    def compress_list_tuples(
        self, list_of_tuples, compresslevel=5,
        compression_type=None, overwrite=False, max_processors=1):
        """
        Compress list of files

        @attention: When reading file and compression_type None, type will be guessed.

        @param list_of_tuples: Path to file and destination folder
        @type list_of_tuples: list[tuple[str|unicode, str|unicode]]
        @param compresslevel: Higher level is slower but likely better. 0-9, except zip 0-8.
        @type compresslevel: int
        @param compression_type: "zip", "gz", "bz2",
        @type compression_type: str | unicode
        @param overwrite: If false, a path will renamed if not available
        @type overwrite: bool
        @param max_processors: Maximum number processors used for compressing files simultaneously
        @type max_processors: int

        @return: True if stream
        @rtype: None
        """
        task_list = []
        for file_path, dst in list_of_tuples:
            self._logger.debug("Compressing '{file}' to '{dst}'".format(file=file_path, dst=dst))
            if not self.validate_dir(dst, silent=True):
                assert self.validate_dir(dst, only_parent=True), "Bad destination: '{}'.".format(dst)
            if self.validate_file(file_path):
                args = (file_path, dst, compresslevel, compression_type, overwrite)
                task_list.append(TaskThread(_compress_file, args))
            else:
                msg = "File not found '{}'".format(file_path)
                self._logger.error(msg)
                raise IOError(msg)
        list_of_return_values = runThreadParallel(task_list, maxThreads=max_processors)
        for index, return_value in enumerate(list_of_return_values):
            assert return_value is None, "Compressing of '{}' failed. '{}'".format(list_of_tuples[index][0], return_value)


def _compress_file(src, dst='./', compresslevel=5, compression_type=None, overwrite=False):
    # TODO: make this unnecessary
    # workaround since pickling a method is a pain
    """
    Compress a file

    @attention: When reading file and compression_type None, type will be guessed.

    @param src: Path to file
    @type src: str | unicode
    @param dst: Destination path, a directory or file path
    @type dst: str | unicode
    @param compresslevel: Higher level is slower but likely smaller. 0-9, except zip 0-8.
    @type compresslevel: int
    @param compression_type: "zip", "gz", "bz2",
    @type compression_type: str | unicode
    @param overwrite: If false, a path will renamed if not available
    @type overwrite: bool

    @return: True if stream
    @rtype: None
    """
    try:
        compressor = Compress(compression_type)
        compressor.compress_file(src, dst, compresslevel, compression_type, overwrite)
    except AssertionError as e:
        return e.message
    except IOError as e:
        return e.message
    return None
