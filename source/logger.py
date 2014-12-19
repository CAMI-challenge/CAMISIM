__author__ = 'hofmann'

import sys
import logging


class Logger(object):
    message_formatter = logging.Formatter('%(levelname)s: %(message)s')
    _logger = None

    def __init__(self, label="", verbose=True, logfile=None):
        if Logger._logger is not None:
            return
        Logger._logger = logging.getLogger(label)
        Logger._logger.setLevel(logging.DEBUG)
        self._add_log_stderr(verbose)
        if logfile is not None:
            self._add_log_file(logfile)

    def add_log_file(self, logfile):
        self._add_log_file(logfile)

    @staticmethod
    def info(message):
        Logger._logger.info(message)

    @staticmethod
    def error(message):
        Logger._logger.error(message)

    @staticmethod
    def debug(message):
        Logger._logger.debug(message)

    @staticmethod
    def critical(message):
        Logger._logger.critical(message)

    @staticmethod
    def exception(message):
        Logger._logger.exception(message)

    @staticmethod
    def warning(message):
        Logger._logger.warning(message)

    @staticmethod
    def set_debug():
        Logger._logger.setLevel(logging.DEBUG)

    @staticmethod
    def _add_log_stderr(verbose=True):
        err_handler = logging.StreamHandler(sys.stderr)
        err_handler.setFormatter(Logger.message_formatter)
        if verbose:
            err_handler.setLevel(logging.INFO)
        else:
            err_handler.setLevel(logging.WARNING)
        Logger._logger.addHandler(err_handler)

    @staticmethod
    def _add_log_file(filename):
        if filename is None:
            Logger._logger.error("Could not add logfile!".format(filename))
            return
        try:
            err_handler_file = logging.FileHandler(filename, 'w')
            err_handler_file.setFormatter(Logger.message_formatter)
            err_handler_file.setLevel(logging.INFO)
            Logger._logger.addHandler(err_handler_file)
        except Exception:
            Logger._logger.error("Could not open %s for logging".format(filename))
            return