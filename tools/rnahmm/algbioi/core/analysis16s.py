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

from algbioi.com import common, csv
from scripts.configparserwrapper import ConfigParserWrapper


class RRNA16S():
    """
        A class that handels the rRNA 16S analysis.
    """
    def __init__(self, config, s16Database, workingDir):
        assert isinstance(config, ConfigParserWrapper)
        self._config = config
        self._workingDir = workingDir
        # self._refDir = os.path.normpath(s16Database)
        # self._refDict = csv.getMappingTuple(os.path.join(self._refDir, 'content.csv'), (0,1), (2,), '\t')

    def runHMM(self, inputFastaFile, outLog=None, hmmer=3, moltypes="ssu", kingdoms="arc,bac"):
        """
            Run the hidden markov model to get regions in the input sequences where the 16S and 23S genes are located.
        """
        if hmmer == 2:
            self.runHMM_2(inputFastaFile, outLog, moltypes=moltypes, kingdoms=kingdoms)
        elif hmmer == 3:
            self.runHMM_3(inputFastaFile, outLog, moltypes=moltypes, kingdoms=kingdoms)

    def runHMM_2(self, inputFastaFile, outLog=None, moltypes="ssu", kingdoms="arc,bac"):
        hmmInstallDir = self._config.get_value("MarkerGeneExtraction", 'rnaHmmInstallDir', is_path=True)
        rnammer_executable = self._config.get_value("MarkerGeneExtraction", 'rnammer', is_path=True)
        assert isinstance(hmmInstallDir, basestring)
        assert isinstance(rnammer_executable, basestring)

        out_file_name_prefix = os.path.join(self._workingDir, os.path.basename(inputFastaFile))

        hmmer_args = [
            "-i '{}'".format(inputFastaFile),
            "-o '{}'".format(out_file_name_prefix),
            "-r '{}'".format(rnammer_executable),
            "-m '{}'".format(moltypes),
            "-k '{}'".format(kingdoms),
            "-T '{}'".format(self._workingDir)
        ]
        cmd = "{script} {args}".format(
            script=os.path.join(hmmInstallDir, 'rna_hmm2.py'),
            args=" ".join(hmmer_args),
        )

        if os.name != 'posix':
            print 'Cannot run HMM since your system is not "posix" but', str('"' + os.name + '"'), '\n', cmd
            return

        if outLog is not None:
            stdoutLog = open(outLog, 'a')
        else:
            stdoutLog = subprocess.PIPE  # stdout=subprocess.STDOUT
        hmmProc = subprocess.Popen(cmd, shell=True, bufsize=-1, cwd=hmmInstallDir, stdout=stdoutLog)
        print 'run cmd:', cmd
        hmmProc.wait()
        if outLog is not None:
            stdoutLog.close()
        print 'HMM return code:', hmmProc.returncode
        if hmmProc.returncode != 0:
            raise Exception("Command returned with non-zero %s status: %s" % (hmmProc.returncode, cmd))

    def runHMM_3(self, inputFastaFile, outLog=None, moltypes="ssu", kingdoms="arc,bac"):
        hmmInstallDir = self._config.get_value("MarkerGeneExtraction", 'rnaHmmInstallDir', is_path=True)
        hmmerBinDir = self._config.get_value("MarkerGeneExtraction", 'hmmerBinDir', is_path=True)
        assert isinstance(hmmInstallDir, basestring)
        assert isinstance(hmmerBinDir, basestring)

        regionsFile = os.path.join(self._workingDir, str(os.path.basename(inputFastaFile) + '.gff'))
        hmmer_args = [
            "-i '{}'".format(inputFastaFile),
            "-o '{}'".format(regionsFile),
            "-r '{}'".format(os.path.join(hmmerBinDir, "hmmsearch")),
            "-m '{}'".format(moltypes),
            "-k '{}'".format(kingdoms)
        ]
        cmd = "{wrapper} {args}".format(
            hmmerBinDir=hmmerBinDir,
            wrapper=os.path.join(hmmInstallDir, 'rna_hmm3.py'),
            args=" ".join(hmmer_args)
        )
        if os.name != 'posix':
            print 'Cannot run HMM since your system is not "posix" but', str('"' + os.name + '"'), '\n', cmd
            return

        if outLog is not None:
            stdoutLog = open(outLog, 'a')
        else:
            stdoutLog = subprocess.PIPE  # stdout=subprocess.STDOUT
        hmmProc = subprocess.Popen(cmd, shell=True, bufsize=-1, cwd=hmmInstallDir, stdout=stdoutLog)
        print 'run cmd:', cmd
        hmmProc.wait()
        if outLog is not None:
            stdoutLog.close()
        print 'HMM return code:', hmmProc.returncode
        if hmmProc.returncode != 0:
            raise Exception("Command returned with non-zero %s status: %s" % (hmmProc.returncode, cmd))

        handle = open(inputFastaFile, "rU")
        record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        handle.close()
        # trunkoutputfilename = inputFastaFile.split( "/" )[-1]
        trunkoutputfilename = os.path.join(self._workingDir, os.path.basename(inputFastaFile))
        # parse results file line by line
        for line in open(regionsFile, "rU"):
            if line.startswith("#"):
                continue
            line = line.split()
            ident = line[0]
            start = int(line[3])
            stop = int(line[4])
            strand = line[6]
            gene = line[8]
            seq = record_dict[ident].seq
            if strand == "+":
                subseq = seq[start-1:stop]
            elif strand == "-":
                subseq = seq[start-1:stop].reverse_complement()
            else:
                sys.stderr.write(" analysis16s: invalid strand symbol")
                exit(1)

            outfile = open(trunkoutputfilename + "." + gene + ".fna", "a")
            print >> outfile, ">%s_%i_%i_%s" % (ident, start, stop, strand)
            print >> outfile, subseq
            outfile.close()

    def classify16S(self, inputFastaFile, outLog=None):
        """
            Run mothur to classify the sequences.
        """
        try:
            self._classify(16, inputFastaFile, outLog)
        except Exception:
            print('Mothur was not able to classify 16S sequences.')

    def classify23S(self, inputFastaFile, outLog=None):
        try:
            self._classify(23, inputFastaFile, outLog)
        except Exception:
            print('Mothur was not able to classify 23S sequences.')
    def classify5S(self, inputFastaFile, outLog=None):
        try:
            self._classify(5, inputFastaFile, outLog)
        except Exception:
            print('Mothur was not able to classify 5S sequences.')

    def _classify(self, mode, inputFastaFile, outLog=None):

        mothur = os.path.join(os.path.normpath(self._config.get('mothurInstallDir')), 'mothur')

        if mode == 16:
            extractedRegionsFasta = os.path.join(self._workingDir, str(os.path.basename(inputFastaFile) + '.16S_rRNA.fna'))
            taxonomyFile = os.path.join(self._refDir, os.path.normpath(self._refDict[('16S_rRNA','taxonomyDNA')][0]))
            templateFile = os.path.join(self._refDir, os.path.normpath(self._refDict[('16S_rRNA','templateDNA')][0]))
            #mothurPredFileName = str(extractedRegionsFasta[0:extractedRegionsFasta.rindex('.')] + '.taxonomy')
            mothurPredFileName = common.getMothurOutputFilePath(extractedRegionsFasta, taxonomyFile)
            predFileName = os.path.join(self._workingDir, str(os.path.basename(inputFastaFile) + '.16P'))

            #extractedRegionsFasta = str(inputFastaFile + '.16S_rRNA.fna')
            #templateFile = os.path.normpath(self._configRRNA16S.get('mothurClassifyParam16STemplate'))
            #taxonomyFile = os.path.normpath(self._configRRNA16S.get('mothurClassifyParam16STaxonomy'))
            #mothurPredFileName = str(inputFastaFile + '.16S_rRNA.bacteria+archaea.taxonomy')
            #mothurPredFileName = os.path.join(self._workingDir, str(os.path.basename(inputFastaFile) + '.16S_rRNA.bacteria+archaea.taxonomy'))
            #mothurPredFileName = os.path.join(self._workingDir, str(os.path.basename(inputFastaFile) + '.16S_rRNA.fasta.taxonomy'))
            #predFileName = str(inputFastaFile + '.16P')
        elif mode == 23:
            extractedRegionsFasta = os.path.join(self._workingDir, str(os.path.basename(inputFastaFile) + '.23S_rRNA.fna'))
            taxonomyFile = os.path.join(self._refDir, os.path.normpath(self._refDict[('23S_rRNA','taxonomyDNA')][0]))
            templateFile = os.path.join(self._refDir, os.path.normpath(self._refDict[('23S_rRNA','templateDNA')][0]))
            #mothurPredFileName = str(extractedRegionsFasta[0:extractedRegionsFasta.rindex('.')] + '.taxonomy')
            mothurPredFileName = common.getMothurOutputFilePath(extractedRegionsFasta, taxonomyFile)
            predFileName = os.path.join(self._workingDir, str(os.path.basename(inputFastaFile) + '.23P'))

            #extractedRegionsFasta = str(inputFastaFile + '.23S_rRNA.fna')
            #templateFile = os.path.normpath(self._configRRNA16S.get('mothurClassifyParam23STemplate'))
            #taxonomyFile = os.path.normpath(self._configRRNA16S.get('mothurClassifyParam23STaxonomy'))
            #mothurPredFileName = str(inputFastaFile + '.23S_rRNA.bacteria+archaea.taxonomy')
            #mothurPredFileName = os.path.join(self._workingDir, str(os.path.basename(inputFastaFile) + '.23S_rRNA.bacteria+archaea.taxonomy'))
            #mothurPredFileName = os.path.join(self._workingDir, str(os.path.basename(inputFastaFile) + '.23S_rRNA.fasta.taxonomy'))
            #predFileName = str(inputFastaFile + '.23P')
        elif mode == 5:
            #extractedRegionsFasta = str(inputFastaFile + '.5S_rRNA.fna')
            extractedRegionsFasta = os.path.join(self._workingDir, str(os.path.basename(inputFastaFile) + '.5S_rRNA.fna'))
            taxonomyFile = os.path.join(self._refDir, os.path.normpath(self._refDict[('5S_rRNA','taxonomyDNA')][0]))
            templateFile = os.path.join(self._refDir, os.path.normpath(self._refDict[('5S_rRNA','templateDNA')][0]))
            mothurPredFileName = common.getMothurOutputFilePath(extractedRegionsFasta, taxonomyFile)
            predFileName = os.path.join(self._workingDir, str(os.path.basename(inputFastaFile) + '.5P'))

            #templateFile = os.path.normpath(self._configRRNA16S.get('mothurClassifyParam5STemplate'))
            #taxonomyFile = os.path.normpath(self._configRRNA16S.get('mothurClassifyParam5STaxonomy'))
            #mothurPredFileName = os.path.join(self._workingDir,
            #                                  str(os.path.basename(inputFastaFile) + '.5S_rRNA.' + os.path.basename(taxonomyFile) + 'onomy'))#.taxonomy
            #predFileName = str(inputFastaFile + '.5P')

        else:
            raise Exception('Wrong branch')

        if not os.path.isfile(mothurPredFileName):
            mothurPredFileName = common.getMothurOutputFilePath(extractedRegionsFasta, taxonomyFile, suffix='.bayesian.taxonomy')

        param = self._config.get('mothurClassifyParamOther')

        cmd = str('time ' + mothur + ' "#classify.seqs(fasta=' + extractedRegionsFasta + ', template=' + templateFile
                + ', taxonomy=' + taxonomyFile + ', ' + param + ')"')

        if os.name == 'posix':
            if outLog is not None:
                stdoutLog = open(outLog, 'w')
            else:
                stdoutLog = subprocess.STDOUT
            mothurProc = subprocess.Popen(cmd, shell=True, bufsize=-1, cwd=self._workingDir, stdout=stdoutLog)
            print 'run cmd:', cmd
            mothurProc.wait()
            if outLog is not None:
                stdoutLog.close()
            print 'mothur return code:', mothurProc.returncode
            if mothurProc.returncode != 0:
                raise Exception("Command returned with non-zero %s status: %s" % (mothurProc.returncode, cmd))
        else:
            print 'Cannot run mothur since your system is not "posix" but', str('"' + os.name + '"'), '\n', cmd



        #transform mothur prediction files to the tab separated files
        self.mothurPredToTabSepPred(mothurPredFileName, predFileName)

    def mothurPredToTabSepPred(self, mothurPredFileName, outPredFileName):
        """
            Transforms the mothur output prediction file (*.taxonomy) to the tab separated prediction file seqName tab ncbid tab weight.
        """
        try:
            fr = open(os.path.normpath(mothurPredFileName),'r')
        except Exception:
            sys.stderr.write("Cannot open file:" + mothurPredFileName + '\n')
            raise
        else:
            try:
                fw = open(os.path.normpath(outPredFileName), 'w')
                lineCount = 0
                for line in fr:
                    line = common.noNewLine(line)
                    try:
                        if re.match(r'^[0-9]+_[0-9]+_[0-9]+_[0-9]+.*', line):
                            name = re.sub(r'([0-9]+_[0-9]+)_[0-9]+_[0-9]+_[\+\-\t ]+.*', r'\1' , line)
                            tag = re.sub(r'[0-9]+_[0-9]+_([0-9]+_[0-9]+_[\+\-]+)[\t ]+.*', r'\1' , line)
                            placementList = re.sub(r'[0-9]+_[0-9]+_[0-9]+_[0-9]+_[\+\-\t ]+(.*)', r'\1' , line.replace('unclassified;', '')).rsplit(';')
                            if len(placementList) < 2:
                                continue
                            placement = placementList[-2]
                            try:
                                clade = int(re.sub('([0-9]+)\(.*', r'\1' , placement))
                            except ValueError:
                                continue
                            weight = float(re.sub('[0-9]+\(([0-9\.]+)\)', r'\1' , placement))
                            lineCount += 1
                            if lineCount == 1:
                                fw.write(name + '\t' + str(clade) + '\t' + str(weight) + '\t' + str(tag))
                            else:
                                fw.write('\n' + name + '\t' + str(clade) + '\t' + str(weight) + '\t' + str(tag))
                    except Exception:
                        sys.stderr.write('Cannot parse line: ' + str(lineCount) +  'in file: ' + mothurPredFileName + '\n')
                        raise
            except Exception:
                sys.stderr.write("Cannot write to file:" + outPredFileName + '\n')
                raise
            finally:
                fw.close()
            fr.close()

    #
    def setCandidatePlacementFrom16S23S5S(self, sequences, taxonomy, inputFastaFile):
        #set assigned by the 16S countList[0]
        #set assigned by the 23S countList[1]
        #intersection =  both countList[3]
        #predFileName16S = str(inputFastaFile + '.16P')
        predFileName16S = os.path.join(self._workingDir, str(os.path.basename(inputFastaFile) + '.16P'))
        #predFileName23S = str(inputFastaFile + '.23P')
        predFileName23S = os.path.join(self._workingDir, str(os.path.basename(inputFastaFile) + '.23P'))
        #predFileName5S = str(inputFastaFile + '.5P')
        predFileName5S = os.path.join(self._workingDir, str(os.path.basename(inputFastaFile) + '.5P'))
        try:
            seqIdSet16S = self._setCandidatePlacement(sequences, taxonomy, predFileName16S, '16S_rRNA')
        except Exception:
            seqIdSet16S = set()
            print("Can't set candidate placement from 16S classification.")
        try:
            seqIdSet23S = self._setCandidatePlacement(sequences, taxonomy, predFileName23S, '23S_rRNA')
        except Exception:
            seqIdSet23S = set()
            print("Can't set candidate placement from 23S sequences.")
        try:
            seqIdSet5S = self._setCandidatePlacement(sequences, taxonomy, predFileName5S, '5S_rRNA')
        except Exception:
            seqIdSet5S = set()
            print("Can't set candidate placement from 5S sequences.")
        intersectSet = seqIdSet16S | seqIdSet23S | seqIdSet5S
        return [len(seqIdSet16S),len(seqIdSet23S),len(seqIdSet5S),len(intersectSet)]

    #
    def _setCandidatePlacement(self, sequences, taxonomy, predFileName, source):
        assignedIdList = []
        try:
            f = open(os.path.normpath(predFileName),'r')
        except Exception:
            print "Cannot open file:", predFileName
            raise
        else:
            for line in f:
                line = common.noNewLine(line)
                if re.match(r'^[0-9]+_[0-9]+\t[0-9]+\t[0-9\.]+\t[^\t]+$', line):
                    scaffoldId = int(re.sub(r'^([0-9]+)_[0-9]+\t[0-9]+\t[0-9\.]+\t[^\t]+$',r'\1' ,line))
                    contigId = int(re.sub(r'^[0-9]+_([0-9]+)\t[0-9]+\t[0-9\.]+\t[^\t]+$',r'\1' ,line))
                    ncbid = int(re.sub(r'^[0-9]+_[0-9]+\t([0-9]+)\t[0-9\.]+\t[^\t]+$',r'\1' ,line))
                    weight = float(re.sub(r'^[0-9]+_[0-9]+\t[0-9]+\t([0-9\.]+)\t[^\t]+$',r'\1' ,line))
                    tag = str(re.sub(r'^[0-9]+_[0-9]+\t[0-9]+\t[0-9\.]+\t([^\t]+)$',r'\1' ,line))
                    if ncbid != 1:
                        taxPathDict = taxonomy.getPathToRoot(ncbid)
                        if taxPathDict is not None and taxPathDict.keys() >= 1:
                            sequences.setCandidateTaxonomyPath(contigId, scaffoldId, taxPathDict, weight, source, tag)
                            assignedIdList.append(contigId)
                        else:
                            sys.stderr.write(str('No taxonomic path found for ncbid: ' + str(ncbid)))
        finally:
            f.close()

        return set(assignedIdList)



if __name__ == "__main__":
    #print '16S analysis'
    #line = '125_528_-    2(100);unclassified;unclassified;unclassified;unclassified;unclassified;unclassified;unclassified;'
    line = '125_528_-    2157(100);28890(95.3333);183925(95.3333);2158(95.3333);2159(95.3333);2172(95.3333);unclassified;unclassified;'

    try:
        name = re.sub('([0-9]+_[0-9]+)_[\+\-\t ]+.*', r'\1' , line)
        placement = re.sub('[0-9]+_[0-9]+_[\+\-\t ]+(.*)', r'\1' , line.replace('unclassified;', '')).rsplit(';')[-2]
        clade = int(re.sub('([0-9]+)\(.*', r'\1' , placement))
        weight = float(re.sub('[0-9]+\(([0-9\.]+)\)', r'\1' , placement))
    except Exception:
        print sys.stderr.write('Cannot parse line in file x' + '\n')

    print name
    print placement
    print clade
    print weight


