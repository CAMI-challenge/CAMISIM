#!/usr/bin/env python

"""
    Copyright (C) 2014  Ivan Gregor

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Note that we could have written some parts of this code in a nicer way,
    but didn't have time. Be careful when reusing the source code.
"""

import os
import ConfigParser


class Config():
    """
        Initializes a config object to read elements from a section.
    """
    def __init__(self, openedConfigFile, section, defaultDict=None):
        """

        @param openedConfigFile:
        @param section:
        @param defaultDict: if mapping not present in the config file, this mapping is used
        @return:
        """
        self._config = ConfigParser.ConfigParser()
        self._config.readfp(openedConfigFile)
        self._section = section
        self._configFilePath = openedConfigFile.name
        if defaultDict is not None:
            self._defaultDict = defaultDict
        else:
            self._defaultDict = dict()

    def getConfigFile(self):
        return self._configFilePath

    def get(self, option):
        """
            Gets an element 'option' from a given section.
        """
        elem = self._config.get(self._section, option)
        if elem != '':
            return elem
        else:
            return self._defaultDict.get(option, None)


#class Config2():
#    def __init__(self, config, section):
#        """
#            Initializes a config object from an existing Config object to read from a different section.
#        """
#        self._config = config._config
#        self._section = section
#
#    def get(self, option):
#        """
#            Gets an element 'option' from a given section
#        """
#        elem = self._config.get(self._section, option)
#        if elem == '':
#            return None
#        else:
#            return elem


def _test():
    config = Config(open(os.path.normpath('/Users/ivan/Documents/work/python/PyCharm/PPSplus/config.cfg')), 'pPPS')
    print "|", config.get('inputFastaScaffoldsFile'), '|'
    print config.getConfigFile()

    #pattern = config.get('scaffoldPattern')
    #print pattern
    #sn = re.findall(pattern, 'Scaffold_1_4082367.lucy.pga.C2');
    #print sn

    #print config.get('fastaLineMaxChar')
    #print config.get('databaseFile')
    #print config.get('taxonomicRanks')

    #print config.get('rankIdAll')
    #print config.get('rankIdCut')
    #print config.get('rankIdCutMinBp')


if __name__ == "__main__":
    _test()