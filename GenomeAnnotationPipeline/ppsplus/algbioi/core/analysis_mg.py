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
import sys
import re
import subprocess
from Bio import SeqIO

from algbioi.com.csv import forEachLine
from algbioi.com.csv import isComment
from algbioi.com.csv import OutFileBuffer
from algbioi.com import common
from algbioi.core.dna_to_prot import dnaToProt


class _MgFiles():
    """
        Helper class to store all paths to the marker gene files.
    """
    def __init__(self, mgDir):
        self.mgDict = {}
        self.mgDir = mgDir # contains all marker genes

    def parse(self, line):
        val = line.split('\t')
        if len(val) == 3:
            self.addFilePath(val[0], val[1], val[2])

    def addFilePath(self, geneName, fileType, filePath):
        if geneName not in self.mgDict:
            self.mgDict[geneName] = {}
        self.mgDict[geneName][fileType] = os.path.join(self.mgDir, os.path.normpath(filePath))

    def getFilePath(self, geneName, fileType):
        if (geneName in self.mgDict) and (fileType in self.mgDict[geneName]):
            return self.mgDict[geneName][fileType]
        else:
            return None

    def getGeneNameList(self):
        resultList = []
        for name in self.mgDict:
            resultList.append(name)
        return resultList


class _MgRegions():
    """
        Helper class to read the hmmer dom? file to get the regions that correspond to the Amphora marker genes.
    """
    def __init__(self):
        self.entryDict = {}

    def parse(self, line):
        if not isComment(line, '#'):
            entryArray = line.split()
            if len(entryArray) == 23:
                aliFrom = int(entryArray[17])
                aliTo = int(entryArray[18])
                targetName = entryArray[0]
                #if targetName in self.entryDict:
                #    print str('MarkerGeneAnalysis: regions for the targetName were already set old: ' + str(self.entryDict[targetName])
                #    + ' new: ' + str([aliFrom, aliTo]))
                if targetName not in self.entryDict:
                    self.entryDict[targetName] = []
                self.entryDict[targetName].append([aliFrom, aliTo])
                #entry['target_name'] = entryArray[0]
                #entry['query_name'] = entryArray[3]
                #entry['dom_i_evalue'] = entryArray[12]
                #entry['ali_from'] = int(entryArray[17])
                #entry['ali_to'] = int(entryArray[18])
            else:
                print 'Line not considered: ', line

    def getEntryDict(self):
        return self.entryDict


class _MothurOutFileParser():
    """
        Helper class to parse the mothur output file and creates a standardized output file for each marker gene.
    """
    def __init__(self, outBuffer, source):
        self.outBuffer = outBuffer
        self.source = source

    def parse(self, line):
        lineArray = line.split()
        if len(lineArray) != 2:
            print '_MothurOutFileParser: wrong line', line
            return
        name = re.sub(r'^([0-9]+_[0-9]+)_[0-9]+_[0-9]+_[pr]+[0-2]$',r'\1', lineArray[0])
        tag = re.sub(r'^[0-9]+_[0-9]+_([0-9]+_[0-9]+_[pr]+[0-2])$',r'\1', lineArray[0])
        placementList = lineArray[1].replace('unclassified;', '').rsplit(';')
        if len(placementList) < 2:
            #print '_MothurOutFileParser: skip line', line
            return

        placement = placementList[-2]
        try:
            clade = int(re.sub('([0-9]+)\(.*', r'\1' , placement))
        except ValueError:
            return
        weight = float(re.sub('[0-9]+\(([0-9\.]+)\)', r'\1' , placement))

        entry = str(str(name) + '\t' + str(clade) + '\t' + str(weight) + '\t' + str(self.source) + '\t' + str(tag))
        if self.outBuffer.isEmpty():
            self.outBuffer.writeText(entry)
        else:
            self.outBuffer.writeText(str('\n' + entry))

    def finalize(self):
        self.outBuffer.close()


class MarkerGeneAnalysis():
    """
        Main class to perform the marker gene analysis based on the Amphora marker genes.
    """
    def __init__(self, config, mgDatabase, workingDir, mgWorkingDir):
        self.markerGeneListFileDir = os.path.normpath(mgDatabase)
        self.markerGeneListFile = os.path.join(self.markerGeneListFileDir, 'content.csv')
        #self.markerGeneListFile = os.path.normpath(configMG.get('markerGeneListFile'))
        self.markerGeneWorkingDir = mgWorkingDir #os.path.normpath(configMG.get('markerGeneWorkingDir'))
        self.hmmInstallDir = os.path.normpath(config.get('rnaHmmInstallDir'))
        self.hmmerBinDir = os.path.normpath(config.get('hmmerBinDir'))
        self.mothurParam = config.get('mothurClassifyParamOther')
        self.workingDir = workingDir
        self.hmmerBinDir = os.path.normpath(config.get('hmmerBinDir'))
        self.mothur = os.path.join(os.path.normpath(config.get('mothurInstallDir')), 'mothur')


    def runMarkerGeneAnalysis(self, fastaFileDNA, outLog=None):
        """
            Run hmmer HMM and mothur classify (bayesian), same param as for the 16S analysis.
        """
        #read list of marker genes
        mgFiles = forEachLine(self.markerGeneListFile, _MgFiles(self.markerGeneListFileDir))

        #translate DNA to protein sequences
        fastaFileProt = os.path.join(self.markerGeneWorkingDir, str(os.path.basename(fastaFileDNA) + '.PROT'))
        dnaToProt(fastaFileDNA, fastaFileProt)

        #read DNA fasta file
        try:
            handle = open(fastaFileDNA, "rU")
            dnaSeqDict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
            handle.close()
        except Exception:
            sys.stderr.write(str('Cannot read file: ' + str(fastaFileDNA)))
            raise

        #to output all predictions in one file
        outPredAllFileName = os.path.join(self.markerGeneWorkingDir,
                                           str(os.path.basename(fastaFileDNA) + '_all.mP'))
        outAllBuffer = OutFileBuffer(outPredAllFileName)

        #run HMM search
        mgList = mgFiles.getGeneNameList()

        if outLog is not None:
            stdoutLog = open(outLog,'w')
        else:
            stdoutLog = subprocess.STDOUT

        #for each gene perform the analysis separately
        for geneName in mgList:

            domFileArray = [os.path.join(self.markerGeneWorkingDir, str(geneName + '_1.dom')),
                            os.path.join(self.markerGeneWorkingDir, str(geneName + '_2.dom'))]
            outFileArray = [os.path.join(self.markerGeneWorkingDir, str(geneName + '_1.out')),
                            os.path.join(self.markerGeneWorkingDir, str(geneName + '_2.out'))]
            hmmFileArray = [mgFiles.getFilePath(geneName, 'hmmPROTPrim'),
                            mgFiles.getFilePath(geneName, 'hmmPROTSec')]
            cmdArray = list([])

            #define cmd
            for i in range(2):
                if hmmFileArray[i] is not None:
                    cmdArray.append(str(os.path.join(self.hmmerBinDir, 'hmmsearch') + ' --domtblout ' + domFileArray[i] + ' -E 0.01'
                               + ' -o ' + outFileArray[i] + ' ' + hmmFileArray[i] + ' ' + fastaFileProt))
                else:
                    cmdArray.append(None)

            #run cmd
            for cmd in cmdArray:
                if cmd is not None and os.name == 'posix':
                    hmmProc = subprocess.Popen(cmd, shell=True, bufsize=-1, cwd=self.hmmInstallDir, stdout=stdoutLog)
                    print 'run cmd:', cmd
                    hmmProc.wait()
                    print 'HMM  return code:', hmmProc.returncode
                    if hmmProc.returncode != 0:
                        raise Exception("Command returned with non-zero %s status: %s" % (hmmProc.returncode, cmd))
                else:
                    print 'Marker genes analysis, doesn`t run (no posix): ', cmd


            #get regions that match to the HMM profile ()
            entryDictList = []
            for i in range(2):
                if cmdArray[i] is not None:
                    entryDictList.append(forEachLine(domFileArray[i], _MgRegions()).getEntryDict())
                else:
                    entryDictList.append(None)

            entryDict1 = entryDictList[0]
            entryDict2 = entryDictList[1]

            #extract regions found in the protein sequences that were found by the HMM and generate corresponding DNA sequences
            regionDnaFasta = os.path.join(self.markerGeneWorkingDir, str(geneName + '_dna.gff'))
            outFileBuffer = OutFileBuffer(regionDnaFasta)

            for seqName in entryDict1:
                i = -1
                for e in entryDict1[seqName]:
                    i += 1
                    from1 = entryDict1[seqName][i][0]
                    to1 = entryDict1[seqName][i][1]
                    assert ((from1 != None) and (to1 != None))
                    #compare the results found by the primary and secondary HMM profiles
                    if (entryDict2 != None) and (seqName in entryDict2):
                        if len(entryDict2[seqName]) >= (i+1):
                            from2 = entryDict2[seqName][i][0]
                            to2 = entryDict2[seqName][i][1]
                            #if from1 != from2 or to1 != to2:
                            #    print str('Different positions in' + seqName + ' from1:' + str(from1) + ' from2:' + str(from2)
                            #                + ' to1:' + str(to1) + ' to2:' + str(to2))

                    #extract regions from the DNA sequences (consider 3 ORF and reverse complements)

                    #name of the whole sequence
                    dnaSeqName = re.sub(r'([0-9]+_[0-9]+)_[pr]+[012]', r'\1', seqName)
                    #whole DNA sequence
                    dnaSeq = dnaSeqDict[dnaSeqName].seq

                    #reverse complement (contains "pr")
                    tagRev = 'p'
                    if re.match(r'[0-9]+_[0-9]+_pr[012]', seqName):
                        dnaSeq = dnaSeq.reverse_complement()
                        tagRev = 'pr'

                    #shift "0"
                    if re.match(r'[0-9]+_[0-9]+_[pr]+0', seqName):
                        tagFrom = ((from1 - 1)*3)
                        tagTo = (to1*3)
                        tagRev += '0'
                        dnaSeq = dnaSeq[tagFrom:tagTo]

                    #shift "1"
                    elif re.match(r'[0-9]+_[0-9]+_[pr]+1', seqName):
                        tagFrom = (((from1 - 1)*3) + 1)
                        tagTo = ((to1*3) + 1)
                        tagRev += '1'
                        dnaSeq = dnaSeq[tagFrom:tagTo]

                    #shift "2"
                    elif re.match(r'[0-9]+_[0-9]+_[pr]+2', seqName):
                        tagFrom = (((from1 - 1)*3) + 2)
                        tagTo = ((to1*3) + 2)
                        tagRev += '2'
                        dnaSeq = dnaSeq[tagFrom:tagTo]

                    #error
                    else:
                        sys.stderr.write('Wrong seq name: ' + seqName + ' \n')
                        dnaSeq = None

                    tag = str(str(tagFrom) + '_' + str(tagTo) + '_' + tagRev)
                    outFileBuffer.writeText(str('>' + dnaSeqName + '_' + tag + '\n' + dnaSeq + '\n'))

            outFileBuffer.close()

            #if no marker gene found
            if outFileBuffer.isEmpty():
                continue

            #run mothur classify (bayesian? the same as for the 16S analysis)
            templateFile = mgFiles.getFilePath(geneName, 'templateDNA')
            taxonomyFile = mgFiles.getFilePath(geneName, 'taxonomyDNA')
            assert ((templateFile is not None) and (taxonomyFile is not None))
            cmd = str('time ' + self.mothur + ' "#classify.seqs(fasta=' + regionDnaFasta + ', template=' + templateFile
                + ', taxonomy=' +  taxonomyFile + ', ' + self.mothurParam + ')"')
            if os.name == 'posix':
                mothurProc = subprocess.Popen(cmd, shell=True, bufsize=-1, cwd=self.markerGeneWorkingDir, stdout=stdoutLog)
                print 'run cmd:', cmd
                mothurProc.wait()
                print 'mothur return code:', mothurProc.returncode
                if mothurProc.returncode != 0:
                    raise Exception("Command returned with non-zero %s status: %s" % (mothurProc.returncode, cmd))
            else:
                print 'Cannot run mothur since your system is not "posix" but', str('"' + os.name + '"'), '\n', cmd

            #transform the mothur output to a simple output (name, ncbid, weight)

            #mothurPredFileName = os.path.join(self.markerGeneWorkingDir,
            #                                  str(geneName + '_dna.' + os.path.basename(taxonomyFile) + 'onomy'))  # taxonomy
            #!!!!!!!!!!!!!
            mothurPredFileName = common.getMothurOutputFilePath(regionDnaFasta, taxonomyFile)
            if not os.path.isfile(mothurPredFileName):
                mothurPredFileName = common.getMothurOutputFilePath(regionDnaFasta, taxonomyFile, suffix='.bayesian.taxonomy')
                if not os.path.isfile(mothurPredFileName):
                    print("Can't open file: %s" % mothurPredFileName)

            outPredFileName = os.path.join(self.markerGeneWorkingDir,
                                           str(os.path.basename(fastaFileDNA) + '_' + geneName + '.mP'))
            outBuffer = OutFileBuffer(outPredFileName, bufferText=True)
            forEachLine(mothurPredFileName, _MothurOutFileParser(outBuffer, geneName))

            if not outAllBuffer.isEmpty():
                outAllBuffer.writeText('\n')
            outAllBuffer.writeText(outBuffer.getTextBuffer())

        if outLog is not None:
            stdoutLog.close()
        outAllBuffer.close()


    def setCandidatePlacement(self, sequences, taxonomy, fastaFileDNA):
        """
            Set candidate placement according to the marker gene analysis !!!
        """
        outPredAllFileName = os.path.join(self.markerGeneWorkingDir,
                                          str(os.path.basename(fastaFileDNA) + '_all.mP'))
        return forEachLine(outPredAllFileName, _SetCandidatePlacement(sequences, taxonomy)).getAssignedSeqCount()


class _SetCandidatePlacement():
    """
        Helper class to set the candidate placements.
    """
    def __init__(self, sequences, taxonomy):
        self.sequences = sequences
        self.taxonomy = taxonomy
        self.assignedIdList = []

    def parse(self, line):
        if line.strip() == '':
            return

        if re.match(r'^[0-9]+_[0-9]+\t[0-9]+\t[0-9\.]+\t[^\t]+\t[^\t]+$', line):
            scaffoldId = int(re.sub(r'^([0-9]+)_[0-9]+\t[0-9]+\t[0-9\.]+\t[^\t]+\t[^\t]+$',r'\1' ,line))
            contigId = int(re.sub(r'^[0-9]+_([0-9]+)\t[0-9]+\t[0-9\.]+\t[^\t]+\t[^\t]+$',r'\1' ,line))
            ncbid = int(re.sub(r'^[0-9]+_[0-9]+\t([0-9]+)\t[0-9\.]+\t[^\t]+\t[^\t]+$',r'\1' ,line))
            weight = float(re.sub(r'^[0-9]+_[0-9]+\t[0-9]+\t([0-9\.]+)\t[^\t]+\t[^\t]+$',r'\1' ,line))
            source = str(re.sub(r'^[0-9]+_[0-9]+\t[0-9]+\t[0-9\.]+\t([^\t]+)\t[^\t]+$',r'\1' ,line))
            tag = str(re.sub(r'^[0-9]+_[0-9]+\t[0-9]+\t[0-9\.]+\t[^\t]+\t([^\t]+)$',r'\1' ,line))

        if ncbid != 1:
            taxPathDict = self.taxonomy.getPathToRoot(ncbid)
            if taxPathDict is not None and taxPathDict.keys() >= 1:
                self.sequences.setCandidateTaxonomyPath(contigId, scaffoldId, taxPathDict, weight, source, tag)
                self.assignedIdList.append(contigId)
            else:
                sys.stderr.write(str('No taxonomic path found for ncbid: ' + str(ncbid)))

    def getAssignedSeqCount(self):
        return len(set(self.assignedIdList))



def test():
    print 'Marker gene analysis'

if __name__ == "__main__":
    test()