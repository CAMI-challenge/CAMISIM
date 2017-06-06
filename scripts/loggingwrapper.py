__author__ = 'Peter Hofmann'
__version__ = '0.1.4'

import sys
import io
import logging
import random
from io import StringIO

logger = logging.getLogger(__name__)


class LoggingWrapper(object):
    CRITICAL = logging.CRITICAL
    FATAL = logging.CRITICAL
    ERROR = logging.ERROR
    WARNING = logging.WARNING
    WARN = logging.WARN
    INFO = logging.INFO
    DEBUG = logging.DEBUG
    NOTSET = logging.NOTSET

    if sys.version_info < (3,):
        _levelNames = logging._levelNames
    else:
        _levelNames = logging._levelToName  # python 3

    _map_logfile_handler = dict()

    def __init__(self, label="", verbose=True, message_format=None, date_format=None, stream=sys.stderr):
        """
        Wrapper for the logging module for easy use.

        @attention: 'labels' are unique, LoggingWrapper with the same label will have the same streams!

        @param label: unique label for a LoggingWrapper
        @type label: str
        @param verbose: Not verbose means that only warnings and errors will be past to stream
        @type verbose: bool
        @param message_format: "%(asctime)s %(levelname)s: [%(name)s] %(message)s"
        @type message_format: str
        @param date_format: "%Y-%m-%d %H:%M:%S"
        @type date_format: str
        @param stream: To have no output at all, use "stream=None", stderr by default
        @type stream: file | FileIO | StringIO | None

        @return: None
        @rtype: None
        """
        assert isinstance(label, str)
        assert isinstance(verbose, bool)
        assert message_format is None or isinstance(message_format, str)
        assert message_format is None or isinstance(date_format, str)
        assert stream is None or self.is_stream(stream), stream

        if message_format is None:
            message_format = "%(asctime)s %(levelname)s: [%(name)s] %(message)s"
        if date_format is None:
            date_format = "%Y-%m-%d %H:%M:%S"
        self.message_formatter = logging.Formatter(message_format, date_format)
        old_label = label
        index = 0
        while label in logging.Logger.manager.loggerDict:
            index = random.randint(0, 99999999999)
            label = old_label
            label = label + " {}".format(index)

        self._label = label
        self._logger = logging.getLogger(label)

        if label in LoggingWrapper._map_logfile_handler:
            return

        LoggingWrapper._map_logfile_handler[label] = None
        self._logger.setLevel(logging.DEBUG)
        if stream is not None:
            if verbose:
                self.add_log_stream(stream=stream, level=logging.INFO)
            else:
                self.add_log_stream(stream=stream, level=logging.WARNING)

    def __exit__(self, type, value, traceback):
        self._close()

    def __enter__(self):
        return self

    def __del__(self):
        self._close()

    @staticmethod
    def is_stream(stream):
        """
        Test for streams

        @param stream: Any kind of stream type
        @type stream: file | io.FileIO | StringIO.StringIO

        @return: True if stream
        @rtype: bool
        """
        return hasattr(stream, 'read') and hasattr(stream, 'write')

    def get_label(self):
        return self._label

    def _close(self):
        """
        Close all logfile handler, unless given as stream.
        Remove all stream handler from handler list, stopping the log service

        @attention: only files opened by LoggingWrapper will be closed!
        If given as stream, logfiles will be kept open!

        @return: None
        @rtype: None
        """
        if hasattr(self, '_logger'):
            list_of_handlers = list(self._logger.handlers)
            for item in list_of_handlers:
                self._logger.removeHandler(item)
            if self._label not in LoggingWrapper._map_logfile_handler:
                return

            logfile_handler = LoggingWrapper._map_logfile_handler.pop(self._label)
            if logfile_handler is not None:
                logfile_handler.close()
        else:
            logger.warning('no attribute named "_logger"')

    def info(self, message):
        """
        Log general informative messages, that might be useful for the user.

        @param message: Message to be logged
        @type message: str

        @return: None
        @rtype: None
        """
        self._logger.info(message)

    def error(self, message):
        """
        Log an significant error that occured.

        @param message: Message to be logged
        @type message: str

        @return: None
        @rtype: None
        """
        self._logger.error(message)

    def debug(self, message):
        """
        Log a message for debugging puposes only.

        @param message: Message to be logged
        @type message: str

        @return: None
        @rtype: None
        """
        self._logger.debug(message)

    def critical(self, message):
        """
        Log a catastrophic error!

        @param message: Message to be logged
        @type message: str

        @return: None
        @rtype: None
        """
        self._logger.critical(message)

    def exception(self, message):
        """
        Log a exception with messages, that might be useful for the user.

        @attention: Call this only after an exception occurred, like in a "try..except.."!

        @param message: Message to be logged
        @type message: str

        @return: None
        @rtype: None
        """
        self._logger.exception(message)

    def warning(self, message):
        """
        Log warning messages, that the user should pay attention to.

        @param message: Message to be logged
        @type message: str

        @return: None
        @rtype: None
        """
        self._logger.warning(message)

    def set_level(self, level):
        """
        Set the minimum level of messages to be logged.

        Level of Log Messages
        CRITICAL    50
        ERROR    40
        WARNING    30
        INFO    20
        DEBUG    10
        NOTSET    0

        @param level: minimum level of messages to be logged
        @type level: int or long

        @return: None
        @rtype: None
        """
        assert level in self._levelNames

        list_of_handlers = self._logger.handlers
        for handler in list_of_handlers:
            handler.setLevel(level)

    def add_log_stream(self, stream=sys.stderr, level=logging.INFO):
        """
        Add a stream where messages are outputted to.

        @param stream: stderr/stdout or a file stream
        @type stream: file | FileIO | StringIO
        @param level: minimum level of messages to be logged
        @type level: int | long

        @return: None
        @rtype: None
        """
        assert self.is_stream(stream)
        # assert isinstance(stream, (file, io.FileIO))
        assert level in self._levelNames

        err_handler = logging.StreamHandler(stream)
        err_handler.setFormatter(self.message_formatter)
        err_handler.setLevel(level)
        self._logger.addHandler(err_handler)

    def set_log_file(self, log_file, mode='a', level=logging.INFO):
        """
        Add a stream where messages are outputted to.

        @attention: file stream will only be closed if a file path is given!

        @param log_file: file stream or file path of logfile
        @type log_file: file | FileIO | StringIO | str
        @param mode: opening mode for logfile, if a file path is given
        @type mode: str
        @param level: minimum level of messages to be logged
        @type level: int or long

        @return: None
        @rtype: None
        """
        assert isinstance(log_file, str) or self.is_stream(log_file)
        assert level in self._levelNames

        if LoggingWrapper._map_logfile_handler[self._label] is not None:
            self._logger.removeHandler(LoggingWrapper._map_logfile_handler[self._label])
            LoggingWrapper._map_logfile_handler[self._label].close()
            LoggingWrapper._map_logfile_handler[self._label] = None

        if self.is_stream(log_file):
            self.add_log_stream(stream=log_file, level=level)
            return

        try:
            err_handler_file = logging.FileHandler(log_file, mode)
            err_handler_file.setFormatter(self.message_formatter)
            err_handler_file.setLevel(level)
            self._logger.addHandler(err_handler_file)
            LoggingWrapper._map_logfile_handler[self._label] = err_handler_file
        except Exception:
            sys.stderr.write("[LoggingWrapper] Could not open '{}' for logging\n".format(log_file))
            return


class DefaultLogging(object):

    def __init__(self, label="DefaultLogging", logfile=None, verbose=False, debug=False):
        """
        Prototype class for any class needing a logger

        @attention:

        @param logfile: file handler or file path to a log file
        @type logfile: file | FileIO | StringIO | str
        @param verbose: Not verbose means that only warnings and errors will be past to stream
        @type verbose: bool
        @param debug: Display debug messages
        @type debug: bool

        @return: None
        @rtype: None
        """
        assert isinstance(debug, bool)

        self._logger = LoggingWrapper(label, verbose=verbose)
        if logfile:
            self._logger.set_log_file(logfile, mode='a')

        self._debug = debug
        if debug:
            self._logger.set_level(self._logger.DEBUG)

        self._logfile = None
        if isinstance(logfile, str):
            self._logfile = logfile
        else:
            if hasattr(logfile, 'name'):
                self._logfile = logfile.name
        self._verbose = verbose

    def __exit__(self, type, value, traceback):
        self._close()

    def __enter__(self):
        return self

    def __del__(self):
        self._close()

    def _close(self):
        self._logger = None

    def set_log_level(self, verbose, debug):
        """
        Simplified way to set log level.

        @attention verbose: Ignored if 'debug' true

        @param verbose: Display info messages and higher
        @type verbose: bool
        @param debug: Display debug messages and higher
        @type debug: bool

        @return: Nothing
        @rtype: None
        """
        if debug:
            self._logger.set_level(self._logger.DEBUG)
        elif verbose:
            self._logger.set_level(self._logger.INFO)
        else:
            self._logger.set_level(self._logger.WARNING)

    @staticmethod
    def is_stream(stream):
        """
        Test for streams

        @param stream: Any kind of stream type
        @type stream: file | io.FileIO | StringIO.StringIO

        @return: True if stream
        @rtype: bool
        """
        return hasattr(stream, 'read') and hasattr(stream, 'write')
