__author__ = 'hofmann'

import sys
import logging


class Logger(object):
    _message_formatter = logging.Formatter('%(levelname)s: %(message)s')

    def __init__(self, label="", verbose=True, logfile=None):
        self.logger = logging.getLogger(label)
        self.logger.setLevel(logging.DEBUG)
        self._add_log_stderr(verbose)
        if logfile is not None:
            self._add_log_file(logfile)

    def add_log_file(self, logfile):
        self._add_log_file(logfile)

    def info(self, message):
        self.logger.info(message)

    def error(self, message):
        self.logger.error(message)

    def debug(self, message):
        self.logger.debug(message)

    def critical(self, message):
        self.logger.critical(message)

    def exception(self, message):
        self.logger.exception(message)

    def warning(self, message):
        self.logger.warning(message)

    def set_debug(self):
        self.logger.setLevel(logging.DEBUG)

    def _add_log_stderr(self, verbose=True):
        err_handler = logging.StreamHandler(sys.stderr)
        err_handler.setFormatter(Logger._message_formatter)
        if verbose:
            err_handler.setLevel(logging.INFO)
        else:
            err_handler.setLevel(logging.WARNING)
        self.logger.addHandler(err_handler)

    def _add_log_file(self, filename):
        if filename is None:
            self.logger.error("Could not add logfile!".format(filename))
            return
        try:
            err_handler_file = logging.FileHandler(filename, 'w')
            err_handler_file.setFormatter(Logger._message_formatter)
            err_handler_file.setLevel(logging.INFO)
            self.logger.addHandler(err_handler_file)
        except Exception:
            self.logger.error("Could not open %s for logging".format(filename))
            return