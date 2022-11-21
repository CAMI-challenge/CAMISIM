__author__ = 'hofmann'
__version__ = '0.0.5'

import os
import io
import zipfile
import tarfile
from .compress import Compress


class Archive(Compress):
    """Reading and writing archived and/or compressed files"""

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
        @type logfile: file | io.FileIO | StringIO.StringIO | str
        @param verbose: Not verbose means that only warnings and errors will be past to stream
        @type verbose: bool

        @return: None
        @rtype: None
        """
        assert logfile is None or isinstance(logfile, str) or self.is_stream(logfile)
        assert isinstance(default_compression, str), "separator must be string"
        assert isinstance(verbose, bool), "verbose must be true or false"
        assert default_compression.lower() in self._open, "Unknown compression: '{}'".format(default_compression)

        super(Archive, self).__init__(label="Archive", default_compression=default_compression, logfile=logfile, verbose=verbose)

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
        return tarfile.is_tarfile(file_path) or zipfile.is_zipfile(file_path)

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

    @staticmethod
    def zip_directory(src_dir, dst):
        assert os.path.isdir(src_dir), "Invalid, not a directory: '{}'".format(src_dir)
        with Archive._open["zip"](dst, 'w', zipfile.ZIP_DEFLATED) as write_handler:
            Archive.zip_stream(src_dir, write_handler)

    @staticmethod
    def zip_stream(src_dir, output_stream):
        """

        @param src_dir:
        @type src_dir: str
        @param output_stream:
        @type output_stream: zipfile.ZipFile
        @return:
        """
        root_path = os.path.dirname(src_dir)
        assert os.path.isdir(src_dir), "Invalid, not a directory: '{}'".format(src_dir)
        for root, directories, files in os.walk(src_dir):
            for file_name in files:
                file_path = os.path.join(root, file_name)
                relative_path = os.path.relpath(file_path, root_path)
                output_stream.write(file_path, arcname=relative_path)

    def zip_decompress_all(self, file_path, output_directory):
        assert self.validate_dir(output_directory, only_parent=True)
        with Archive._open["zip"](file_path, "r") as read_handler:
            read_handler.extractall(output_directory)

    def tar_decompress_all(self, file_path, output_directory):
        assert self.validate_dir(output_directory, only_parent=True)
        with tarfile.open(file_path) as read_handler:
            def is_within_directory(directory, target):
                
                abs_directory = os.path.abspath(directory)
                abs_target = os.path.abspath(target)
            
                prefix = os.path.commonprefix([abs_directory, abs_target])
                
                return prefix == abs_directory
            
            def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
            
                for member in tar.getmembers():
                    member_path = os.path.join(path, member.name)
                    if not is_within_directory(path, member_path):
                        raise Exception("Attempted Path Traversal in Tar File")
            
                tar.extractall(path, members, numeric_owner=numeric_owner) 
                
            
            safe_extract(read_handler, output_directory)

    def extract_all(self, file_path, output_directory, compression_type=None, mode='r'):
        if compression_type is None:
            compression_type = self.get_compression_type(file_path)
        if compression_type == "zip":
            self.zip_decompress_all(file_path, output_directory)
            return
        assert compression_type in Archive._archive_mode_read_stream, "Unknown compression type: '{}'".format(compression_type)
        self.tar_decompress_all(file_path, output_directory)
