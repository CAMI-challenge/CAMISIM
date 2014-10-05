"""
    Generates test scenarios for the simulated tests.
"""
import os
from algbioi.com import csv


class FileHandler():
    def __init__(self, referenceConfig):
        self.referenceConfig = referenceConfig
        self.workingDir = os.path.dirname(referenceConfig)
        self.runTestSh = os.path.join(self.workingDir, 'runTest.sh')


class SubstDict():
    def __init__(self):
        self._propertyMap = {}
        self._featureMap = {}

    def setProperty(self, key, val):
        assert key not in self._propertyMap
        self._propertyMap[key] = val

    def getProperty(self, key):
        return self._propertyMap.get(key, None)

    def setFeature(self, key, val):
        assert key not in self._featureMap
        self._featureMap[key] = val

    def getFeature(self, key):
        return self._featureMap.get(key, None)


def substituteProperties(inContigFile, outContigFile, substDict):
    """
        Generates a new configuration file such that properties in the input configuration file are substituted
        according to the entries in the substitution dictionary

        @param inContigFile: input reference config
        @param outContigFile: output config with altered entries according to the change dict
        @param substDict: contain substitutions
        @type substDict: SubstDict
    """
    lines = csv.getColumnAsList(inContigFile, colNum=0, sep='\n', comment='_')
    out = csv.OutFileBuffer(outContigFile)
    for line in lines:
        key = line.split('=')[0]
        val = substDict.getProperty(key)
        if val is None:
            out.writeText(line + '\n')
        else:
            out.writeText(str(key) + '=' + str(val) + '\n')
    out.close()


def setConfigNameLogNamePipelineDir(inConfigFile, substDict):
    """

        @param inConfigFile: reference config input path
        @param substDict: pipeline dir entry will be added
        @type substDict: SubstDict
        @return: configName
    """
    rs = substDict.getProperty('excludeRefSeqRank')
    mg = substDict.getProperty('excludeRefMgRank')
    if rs is None:
        rs = 'no'
    if mg is None:
        mg = 'no'
    prefix = os.path.basename(inConfigFile).split('_', 1)[0]
    dirName = prefix + '_rs_' + rs + '_mg_' + mg
    configName = dirName + '.cfg'
    logName = dirName + '.log'
    workingDir = os.path.dirname(inConfigFile)
    substDict.setProperty('pipelineDir', os.path.join(workingDir, dirName))
    substDict.setFeature('configName', os.path.join(workingDir, configName))
    substDict.setFeature('logName', os.path.join(workingDir, logName))


def generateTestScenarios(fh, ppsMasterScript, ppspArgs, scenario='ppsp'):
    """
        Generates test scenarios, i.e. configuration files, pipeline directories and the master shell script.
        All configuration files and directories will be generated in the same directory as the reference configuration
        file is.

        @param fh: contains path to the reference configuration file and to the shell test script
        @type fh: FileHandler
        @param ppsMasterScript: path to the PPSP master script
        @param ppspArgs: arguments of the PPSP master script except for the configuration
        @param scenario: ppsp or pps or test
    """
    outSh = csv.OutFileBuffer(fh.runTestSh)
    outSh.writeText('#!/bin/bash\n\n')

    substDictList = getSubstDictList(scenario)

    count = 0
    for substDict in substDictList:
        count += 1
        # set configuration file name
        # set log file name
        # set pipeline dir
        setConfigNameLogNamePipelineDir(fh.referenceConfig, substDict)

        # generate configuration files
        substituteProperties(fh.referenceConfig, substDict.getFeature('configName'), substDict)

        # create pipeline directory
        os.mkdir(os.path.normpath(substDict.getProperty('pipelineDir')))

        # add an entry to the shell test script
        outSh.writeText('time python ' + ppsMasterScript + ' -c ' + substDict.getFeature('configName') + ' ' + ppspArgs +
                        ' > ' + substDict.getFeature('logName') + '\n\necho "' + str(count) + '"\n\n')

    outSh.close()


def getSubstDictList(scenario='ppsp'):
    """
        Gets the list of substitution entries

        @param scenario:
        @return:
    """
    substDictList = []
    #
    if scenario == 'ppsp':
        entry = SubstDict()
        entry.setProperty('excludeRefSeqRank', 'strain')
        substDictList.append(entry)
        #
        entry = SubstDict()
        entry.setProperty('excludeRefSeqRank', 'species')
        substDictList.append(entry)
        #
        entry = SubstDict()
        entry.setProperty('excludeRefSeqRank', 'genus')
        substDictList.append(entry)
        #
        entry = SubstDict()
        entry.setProperty('excludeRefMgRank', 'strain')
        entry.setProperty('excludeRefSeqRank', 'strain')
        substDictList.append(entry)
        #
        entry = SubstDict()
        entry.setProperty('excludeRefMgRank', 'strain')
        entry.setProperty('excludeRefSeqRank', 'species')
        substDictList.append(entry)
        #
        entry = SubstDict()
        entry.setProperty('excludeRefMgRank', 'strain')
        entry.setProperty('excludeRefSeqRank', 'genus')
        substDictList.append(entry)
        #
        entry = SubstDict()
        entry.setProperty('excludeRefMgRank', 'species')
        entry.setProperty('excludeRefSeqRank', 'species')
        substDictList.append(entry)
        #
        entry = SubstDict()
        entry.setProperty('excludeRefMgRank', 'genus')
        entry.setProperty('excludeRefSeqRank', 'genus')
        substDictList.append(entry)
    elif scenario == 'pps':
        entry = SubstDict()
        entry.setProperty('excludeRefSeqRank', 'strain')
        substDictList.append(entry)
        #
        entry = SubstDict()
        entry.setProperty('excludeRefSeqRank', 'species')
        substDictList.append(entry)
        #
        entry = SubstDict()
        entry.setProperty('excludeRefSeqRank', 'genus')
        substDictList.append(entry)
    elif scenario == 'test':
        entry = SubstDict()
        entry.setProperty('excludeRefMgRank', 'strain')
        entry.setProperty('excludeRefSeqRank', 'strain')
        substDictList.append(entry)
        #
        entry = SubstDict()
        entry.setProperty('excludeRefMgRank', 'strain')
        entry.setProperty('excludeRefSeqRank', 'species')
        substDictList.append(entry)
        #
        entry = SubstDict()
        entry.setProperty('excludeRefMgRank', 'strain')
        entry.setProperty('excludeRefSeqRank', 'genus')
        substDictList.append(entry)
    else:
        print('getSubstDictList: Unknown scenario:' + scenario)

    return substDictList


def _main():
    fh = FileHandler('/net/metagenomics/projects/PPSmg/tests/mercier042013/uniform/config_ex_rs_no_mg_no.cfg')
    scenario = 'test'
    ppspArgs = '-n -g -o s16 mg -t -a -p c -r -s'
    ppsMasterScript = '/net/metagenomics/projects/PPSmg/scripts/scriptsR29/algbioi/core/run.py'
    generateTestScenarios(fh, ppsMasterScript, ppspArgs, scenario)


if __name__ == "__main__":
    _main()