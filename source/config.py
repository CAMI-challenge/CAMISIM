__author__ = 'hofmann'

from ConfigParser import SafeConfigParser
from logger import Logger


class Config(object):
	_boolean_states = {'1': True, 'yes': True, 'true': True, 'on': True,
					   '0': False, 'no': False, 'false': False, 'off': False,
					   'y': True, 't': True, 'n': False, 'f': False}

	def __init__(self, config_file, logger=None):
		self._logger = logger
		if not self._logger:
			self._logger = Logger("Config")
		self._config = SafeConfigParser()
		if isinstance(config_file, str):
			self._config.read(config_file)
			self._config_file_path = config_file
		elif isinstance(config_file, file):
			self._config.readfp(config_file)
			self._config_file_path = config_file.name
		else:
			self._logger.error("[Config] Bad config file type '{}'".format(config_file))
			return

	def has_missing_section(self, list_sections):
		for section in list_sections:
			if not self._config.has_section(section):
				return section
		return False

	def get_value(self, section, option, is_digit=False, is_boolean=False, verbose=True):
		if not self._config.has_section(section):
			if verbose:
				self._logger.error("[Config] Bad section '{}'".format(section))
			return None
		if not self._config.has_option(section, option):
			if verbose:
				self._logger.error("[Config] Bad option in '{}': {}".format(section, option))
			return None

		value = self._config.get(section, option)
		if value == '':
			return None

		if is_digit:
			return self._string_to_digit(value)

		if is_boolean:
			return self._is_true(value)
		return value

	def _string_to_digit(self, value):
		try:
			if '.' in value:
				return float(value)
			return int(value)
		except ValueError:
			self._logger.error("[Config] Bad value '{}'".format(value))
			return None

	def _is_true(self, value=''):
		if value is None or not isinstance(value, str):
			return None

		if value.lower() not in Config._boolean_states:
			self._logger.error("[Config] Bad value '{}'".format(value))
			return None
		return Config._boolean_states[value.lower()]
