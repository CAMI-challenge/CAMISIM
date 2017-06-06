#!/usr/bin/env python

import os
import sys
import subprocess
import argparse
import string
import gzip

from algbioi.com.config import Config
from algbioi.com.csv import getMapping
from algbioi.com.csv import OutFileBuffer
from algbioi.com.fasta import fastaFileToDict


def main():
    """
        Wraps pIRS read simulator to simulate Illumina paired end reads.

        Sample config: /Users/ivan/Documents/work/binning/data/V35/simMetagenome/configMetagenome01.cfg
    """
    if os.name != 'posix':
        print 'runs only on posix systems'
        return

    #parse arguments
    parser = argparse.ArgumentParser(description='''A simple Metagenome Illumina read simulator that wraps pIRS''',
                                 epilog='''''')

    parser.add_argument('-c', '--config', nargs=1, type=file, required=True,
                        help='configuration file of the simulator', metavar='configMetagenome.cfg',
                        dest='config')

    parser.add_argument('-p', '--pIRS-param', action='store', nargs='+',
                        help='parameters of the pIRS simulator, e.g. "-Q 64 -E 1"',
                        dest='p')

    args = parser.parse_args()
    config = Config(args.config[0], 'Sim')

    pirsParam = ''
    if args.p:
        pirsParam = args.p[0]

    #reads configuration
    workingDir = config.get('workingDir')
    referenceSeq = config.get('referenceSeq')
    frequenciesInfo = config.get('frequenciesInfo')
    coverageFrequencyMultiplier = float(config.get('coverageFrequencyMultiplier'))
    pirsInstallDir = config.get('pirsInstallDir')
    insertSizeMean = int(config.get('insertSizeMean'))
    insertSizeSd = int(config.get('insertSizeSd'))
    readLength = int(config.get('readLength'))

    #check whether the pIRS optional parameters doesn`t contain those predefined elsewhere (e.g. in the config)
    if (string.count(pirsParam,'-m') != 0 or string.count(pirsParam,'-v') != 0 or string.count(pirsParam,'-l') != 0
        or string.count(pirsParam,'-x') != 0 or string.count(pirsParam,'-i') != 0 or string.count(pirsParam,'-o') != 0):
        print 'pIRS parameters -m -v -l (-x) must be set in the configuration file, parameters -i -o cannot be set '
        return

    #check working directory, create temporary directory
    tmpDir = os.path.join(workingDir,'tmp')
    if not os.path.isdir(workingDir):
        print str('The working directory does not exists, create it! (' + str(workingDir) + ')')
        return
    if not os.path.isdir(tmpDir):
        os.mkdir(tmpDir)

    seqNameToSeq = fastaFileToDict(referenceSeq)
    seqNameToFreq = getMapping(frequenciesInfo, 0, 1, sep='\t', comment = '#')

    outReads1Merged = OutFileBuffer(os.path.join(workingDir,'reads_1.fq'))
    outReads2Merged = OutFileBuffer(os.path.join(workingDir,'reads_2.fq'))

    for seqName in seqNameToFreq:
        seq = seqNameToSeq[seqName]
        coverage = float(seqNameToFreq[seqName][0])*coverageFrequencyMultiplier

        fastaFile = os.path.join(tmpDir,str(seqName + '.fna'))
        outBuffer = OutFileBuffer(fastaFile)
        outBuffer.writeText(str('>' + seqName + '\n' + seq + '\n'))
        outBuffer.close()

        cmd = str(os.path.join(pirsInstallDir,'pirs') + ' simulate -i ' + fastaFile + ' -x ' + str(coverage) +
                  ' -m ' + str(insertSizeMean) + ' -v ' + str(insertSizeSd) + ' -l ' + str(readLength)
                  + ' -o ' + seqName + ' ' + pirsParam)
        #print cmd
        proc = subprocess.Popen(cmd, shell=True, bufsize=-1, cwd=tmpDir)# stdout=subprocess.STDOUT, stderr=subprocess.STDOUT)
        proc.wait()
        if proc.returncode != 0:
            sys.stderr.write(str('command failed: ' + cmd))

        #append generated reads to the merged files
        reads1 = gzip.open(os.path.join(tmpDir, str(seqName + '_' + str(readLength) + '_' + str(insertSizeMean) + '_1.fq.gz')), 'rb')
        file1Content = reads1.read()
        outReads1Merged.writeText(str(file1Content.replace('@read_',str('@read_' + seqName + '_')) + '\n'))
        reads1.close()

        reads2 = gzip.open(os.path.join(tmpDir, str(seqName + '_' + str(readLength) + '_' + str(insertSizeMean) + '_2.fq.gz')), 'rb')
        file2Content = reads2.read()
        outReads2Merged.writeText(str(file2Content.replace('@read_',str('@read_' + seqName + '_')) + '\n'))
        reads2.close()

    outReads1Merged.close()
    outReads2Merged.close()


if __name__ == "__main__":
    main()