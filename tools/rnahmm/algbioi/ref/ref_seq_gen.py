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


    After the resources has been downloaded, i.e. paths to the following data are available:

    ncbi-refseq-microbial_64:
    ncbi-genomes-bacteria_20140428:
    ncbi-draftgenomes-bacteria_20140508:
    ncbi-hmp_20131125:

    Run the script in these steps:

    1) In method _main, set the following variables to the lists containing the downloaded reference data
    mapFilePathList ~ list of mapping files
    fastaFilePathList ~ list of the corresponding fasta files

    2) create three directories and accordingly modify variables in the _main:
    mergedDir
    sortedDir
    centroidsDir

    3) Run the script using:
    python algbioi/ref/ref_seq_gen.py

    4) The output non-redundant reference is contained in directory defined by variable centroidsDir in the _main,
    where the name of a fasta file denotes the taxon id of contained sequences


    Note that we know that this script could be written in a nicer way!
"""

import os
import sys
import signal
import string
import subprocess
import shutil

from Bio import SeqIO
from Bio.Seq import Seq

from algbioi.com import csv
from algbioi.com import fasta
from algbioi.com import common


def printStatDbk():
    """
        Print statistics of a DBK file.
    """
    seqIdSet = set()
    taxonSet = set()
    cumulativeLen = 0
    recordCount = 0
    zeros = 0
    #
    for record in SeqIO.parse(sys.stdin, "genbank"):
        recordCount += 1
        seqId = record.id

        if seqId in seqIdSet:
            print seqId, 'already in set', seqId
        else:
            seqIdSet.add(seqId)

        seq = str(record.seq)
        cumulativeLen += len(seq)

        if len(string.replace(common.noNewLine(seq), 'N', '')) == 0:
            zeros += 1

        taxonId = None

        for feature in record.features:
            if feature.type == "source":
                for xrefentry in feature.qualifiers["db_xref"]:
                    (key, val) = xrefentry.split(":")
                    if key == "taxon":
                        taxonId = int(val)
                        break
            if taxonId is not None:
                break

        if taxonId is None:
            print 'could not find taxonId for', seqId
        else:
            taxonSet.add(taxonId)

    print 'record count', recordCount
    print 'seq count', len(seqIdSet)
    print 'taxon id count', len(taxonSet)
    if len(seqIdSet) > 0:
        print 'avg. seq. len', cumulativeLen / len(seqIdSet)
    print 'zeros', zeros
    # Find out which sequences are just full of zeros!?


def findPlasmids(outPlasmidFilePath):
    """
        Read sequence descriptions from a DBK files (stdin), output sequence ids (record.id) if the corresponding
        description contain "plasmid".
        Plasmids can be also within the sequences!
    """
    # append to a file if it already exists
    if os.path.isfile(outPlasmidFilePath):
        outFileMode = 'a'
    else:
        outFileMode = 'w'
    outBuffer = csv.OutFileBuffer(outPlasmidFilePath, bufferText=False, fileOpenMode=outFileMode)

    recordCount = 0
    plasmidCount = 0

    for record in SeqIO.parse(sys.stdin, "genbank"):
        recordCount += 1

        if string.find(record.description, 'plasmid') != -1:
            outBuffer.writeText(str(str(record.id) + '\n'))
            plasmidCount += 1

    outBuffer.close()
    print 'file, records, plasmids:', outPlasmidFilePath, recordCount, plasmidCount
    #for feature in record.features:
    #    if feature.type == "source":


def mergeSequences(mapFilePathList, fastaFilePathList, outputDir):
    """
        Reads all sequences. For each taxonId creates a file that contain all sequences
        mapped to this taxonId. If a seqId appears more than one it is ignored since
        acession numbers are unique.

        @param mapFilePathList: list of files where each contain mapping: seqId -> taxonId
        @param fastaFilePathList: list of fasta files that contain mapping: seqId -> seq
    """
    taxonIdToOutBuffer = {}
    seqIdSet = set()

    totalSeqCount = 0
    totalStoredSeqCount = 0
    totalIdenticalSeqCount = 0

    for mapFilePath, fastaFilePath in zip(mapFilePathList, fastaFilePathList):
        print 'processing', mapFilePath, fastaFilePath
        seqCount = 0
        storedSeqCount = 0

        seqIdToSeq = fasta.fastaFileToDict(fastaFilePath)
        seqIdToNcbidList = csv.getMapping(mapFilePath, 0, 1, sep='\t', comment='#')

        for seqId, seq in seqIdToSeq.iteritems():
            seqCount += 1
            if seqId in seqIdSet:
                totalIdenticalSeqCount += 1
                continue
            else:
                seqIdSet.add(seqId)

            taxonId = seqIdToNcbidList[seqId][0]

            if taxonId not in taxonIdToOutBuffer:
                outBuffer = csv.OutFileBuffer(os.path.join(outputDir, str(str(taxonId) + '.fna')))
                taxonIdToOutBuffer[taxonId] = outBuffer

            taxonIdToOutBuffer[taxonId].writeText(str('>' + seqId + '\n' + seq + '\n'))
            taxonIdToOutBuffer[taxonId].close()
            storedSeqCount += 1

            if len(string.replace(common.noNewLine(seq),'N','')) == 0:
                print 'zeros', seqId, fastaFilePath, len(common.noNewLine(seq))

        # for buff in taxonIdToOutBuffer.values():
        #     buff.close()

        print 'totalSeq, storedSeq', seqCount, storedSeqCount
        totalSeqCount += seqCount
        totalStoredSeqCount += storedSeqCount


    print 'totalSeqCount, totalStoredSeqCount, totalIdenticalSeqCount', totalSeqCount, totalStoredSeqCount, totalIdenticalSeqCount

    print 'sequences merged'


def getMaxLen(fastaFilePath):
    """
        Gets the length of the sequence that has maximum length in a fasta file.
    """
    maxLen = 0
    for val in fasta.getSequenceToBpDict(fastaFilePath).itervalues():
        if maxLen < int(val):
            maxLen = int(val)
    return maxLen


def sortSeqDesc(usearch5, usearch6, mergedDir, sortedDir, mapFilePathList):
    """
        Sort sequences in the descending order.
    """
    taxonIdSet = getAllTaxonIdSet(mapFilePathList)
    inFilePathList = []
    sortedFilePathList = []
    for taxonId in taxonIdSet:
        inFilePathList.append(os.path.join(mergedDir, str(str(taxonId) + '.fna')))
        sortedFilePathList.append(os.path.join(sortedDir, str(str(taxonId) + '.fna')))

    # sort sequences, longer first
    for inFile, outFile in zip(inFilePathList, sortedFilePathList):
        sortCmd = str(usearch5 + ' -sort ' + inFile + ' -output ' + outFile + ' --maxlen ' + str(int(getMaxLen(inFile) + 1000)) + ' --minlen 10')
        #sortCmd = str(usearch6 + ' -sortbylength ' + inFile + ' -output ' + outFile) # + ' -maxqt ' + str(int(getMaxLen(inFile))))
        print sortCmd
        sortProc = subprocess.Popen(sortCmd, shell=True, cwd=os.path.dirname(usearch6), bufsize=-1)
        sortProc.wait()
        if sortProc.returncode != 0:
            sys.stderr.write(str('Cmd: ' + sortCmd + ' ended with return code: ' + str(sortProc.returncode) + '\n'))
            return
    print 'sequences sorted'


def toCentroids(sortedDir, centroidsDir, mapFilePathList):
    """
        Find centroids for each cluster within each sorted file.
    """
    taxonIdSet = getAllTaxonIdSet(mapFilePathList)

    sortedFilePathList = []
    centroidsFilePathList = []
    for taxonId in taxonIdSet:
        sortedFilePathList.append(os.path.join(sortedDir, str(str(taxonId) + '.fna')))
        centroidsFilePathList.append(os.path.join(centroidsDir, str(str(taxonId) + '.fna')))

    # get centroids
    for inFile, outFile in zip(sortedFilePathList, centroidsFilePathList):
        getSeeds(inFile, outFile)

        # cdHit = None  # '/net/metagenomics/projects/PPSmg/tools/cd-hit-v4.6.1-2012-08-27$ ./cd-hit-est'
        # usearch1 = None  # '/net/programs/Debian-6.0-x86_64/uclust-1.1.579/uclust'
        # identity = 1.0
        #if len(fastaFileToDict(inFile)) == 1:
        #    continue
        #clusterCmd = str(dnaClust + ' -i ' + inFile + ' -s ' +  str(identity) + ' -l > ' + clusterFile)
        #clusterCmd = str(usearch5 + ' -cluster ' + inFile + ' -id ' +  str(identity) + ' -seedsout ' + outFile + ' --maxlen ' + str(int(getMaxLen(inFile) + 10)) + ' --minlen 10')
        #clusterCmd = str(usearch6 + ' -cluster_fast ' + inFile + ' -id ' +  str(identity) + ' -centroids ' + outFile + ' -uc ' + clusterFile + ' -maxqt ' + str(float(getMaxLen(inFile) + 10)) + ' -minqt 10.0 -strand both')
        #clusterCmd = str(usearch1 + ' --input ' + inFile + ' --id ' +  str(identity) +  ' --maxlen ' + str(int(getMaxLen(inFile) + 10)) + ' --minlen 10 --log log_ucl1.txt --rev --uc ' + clusterFile)
        #cdHit
        #clusterProc = subprocess.Popen(clusterCmd, shell=True, cwd = os.path.dirname(usearch6), bufsize=-1)
        #clusterProc.wait()
        #if clusterProc.returncode != 0:
        #    sys.stderr.write(str('Cmd: ' + clusterCmd + ' ended with return code: ' + str(clusterProc.returncode) + '\n'))
        #    return
    print 'centroids computed'


def getSeeds(inSortedFasta, outSeedsFasta):
    """
        @param inSortedFasta: DNA sequences sorted according to the sequence length in the descending order
        @param outSeedsFasta: a fasta file that contains all seeds
    """
    out = csv.OutFileBuffer(outSeedsFasta)
    seedList = []
    seqList = fasta.getSequencesToList(inSortedFasta)  # list of (sequenceName, sequence)

    for seqId, seq in seqList:
        seq = string.upper(seq)

        newSeed = True
        for seedSeq in seedList:

            if len(seedSeq) < len(seq):
                continue

            # if bool(re.search(seq, seedSeq, re.I)) or bool(re.search(str(Seq(seq).reverse_complement()), seedSeq, re.I)):
            if seq in seedSeq or str(Seq(seq).reverse_complement()) in seedSeq:
                newSeed = False
                break

        if newSeed:
            # print 'new seed:', seqId
            seedList.append(seq)
            out.writeText(str('>' + seqId + '\n' + seq + '\n'))
        # else:
        #    print 'no seed:', seqId

    out.close()

    print 'total', len(seqList)
    print 'seed count', len(seedList)
    print 'duplicate', (len(seqList) - len(seedList))


def filterOutSequencesBatch(taxonIdSet, srcDir, dstDir, notAllowedSeqIdSet):
    """
        For each fasta file that is in directory srcDir filters out sequences that are not defined in the allowedSeqIdSet.
    """
    for taxonId in taxonIdSet:
        srcFilePath = os.path.join(srcDir,str(str(taxonId) + '.1.fna'))
        dstFilePath = os.path.join(dstDir,str(str(taxonId) + '.1.fna'))

        seqIdDict = fasta.getSequenceToBpDict(srcFilePath)
        allowedNamesSet = set()
        for id in seqIdDict.iterkeys():
            if id not in notAllowedSeqIdSet:
                allowedNamesSet.add(id)

        fasta.filterOutSequences(srcFilePath, dstFilePath, allowedNamesSet)


def getAllTaxonIdSet(mapFilePathList):
    taxonIdSet = set()
    for mapFile in mapFilePathList:
        taxonIdSet = taxonIdSet.union(set(csv.getColumnAsList(mapFile, colNum=1, sep='\t')))
    return taxonIdSet


def _main():
    mergeS = False
    sortS = False
    clusterS = False
    filterOutSeq = True  # this is optional, e.g. to remove plasmids

    # handle broken pipes
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    #printStatDbk()
    #checkForPlasmids()

    # mapFilePathList = ['/local/johdro/refdata/static/ncbi-genomes-bacteria_20121122/nobackup/dna_acc.nuc.tax',
    #                    '/local/johdro/refdata/static/ncbi-draftgenomes-bacteria_20121122/nobackup/dna-contigs_acc.nuc.tax',
    #                    '/local/johdro/refdata/static/ncbi-draftgenomes-bacteria_20121122/nobackup/dna-scaffolds_acc.nuc.tax',
    #                    '/local/johdro/refdata/static/ncbi-hmp_20121016/nobackup/dna-contigs_acc.nuc.tax',
    #                    '/local/johdro/refdata/static/ncbi-hmp_20121016/nobackup/dna-scaffolds_acc.nuc.tax',
    #                    '/local/johdro/refdata/static/ncbi-refseq-microbial_56/nobackup/dna_acc.nuc.tax'
    #                    ]

    # fastaFilePathList = ['/local/johdro/refdata/static/ncbi-genomes-bacteria_20121122/nobackup/dna_acc.nuc.fna',
    #                      '/local/johdro/refdata/static/ncbi-draftgenomes-bacteria_20121122/nobackup/dna-contigs_acc.nuc.fna',
    #                      '/local/johdro/refdata/static/ncbi-draftgenomes-bacteria_20121122/nobackup/dna-scaffolds_acc.nuc.fna',
    #                      '/local/johdro/refdata/static/ncbi-hmp_20121016/nobackup/dna-contigs_acc.nuc.fna',
    #                      '/local/johdro/refdata/static/ncbi-hmp_20121016/nobackup/dna-scaffolds_acc.nuc.fna',
    #                      '/local/johdro/refdata/static/ncbi-refseq-microbial_56/nobackup/dna_acc.nuc.fna'
    #                      ]

    # input files
    mapFilePathList = ['/net/refdata/static/nonredundant-microbial_20140513/nobackup/ncbi-genomes-bacteria_20140428/dna_acc.tax',
                        '/net/refdata/static/nonredundant-microbial_20140513/nobackup/ncbi-draftgenomes-bacteria_20140508/dna-scaffolds_acc.tax',
                        '/net/refdata/static/nonredundant-microbial_20140513/nobackup/ncbi-hmp_20131125/dna-contigs_acc.tax',
                        '/net/refdata/static/nonredundant-microbial_20140513/nobackup/ncbi-hmp_20131125/dna-scaffolds_acc.tax',
                        '/net/refdata/static/nonredundant-microbial_20140513/nobackup/ncbi-refseq-microbial_64/dna_acc.tax'
                        ]

    fastaFilePathList = ['/net/refdata/static/nonredundant-microbial_20140513/nobackup/ncbi-genomes-bacteria_20140428/dna_acc.fna',
                         '/net/refdata/static/nonredundant-microbial_20140513/nobackup/ncbi-draftgenomes-bacteria_20140508/dna-scaffolds_acc.fna',
                         '/net/refdata/static/nonredundant-microbial_20140513/nobackup/ncbi-hmp_20131125/dna-contigs_acc.fna',
                         '/net/refdata/static/nonredundant-microbial_20140513/nobackup/ncbi-hmp_20131125/dna-scaffolds_acc.fna',
                         '/net/refdata/static/nonredundant-microbial_20140513/nobackup/ncbi-refseq-microbial_64/dna_acc.fna'
                         ]

    # output dirs
    mergedDir = '/net/refdata/static/nonredundant-microbial_20140513/nobackup/merged'  # '/local/igregor/ref_20121122/nobackup/merged'
    sortedDir = '/net/refdata/static/nonredundant-microbial_20140513/nobackup/sorted'  # '/local/igregor/ref_20121122/nobackup/sorted'
    centroidsDir = '/net/refdata/static/nonredundant-microbial_20140513/nobackup/centroids'  # '/local/igregor/ref_20121122/nobackup/centroids_1_0'
    # clustersDir = '/net/refdata/static/nonredundant-microbial_20140513/nobackup/clusters'  # '/local/igregor/ref_20121122/nobackup/clusters_1_0'

    # tools
    usearch5 = '/net/metagenomics/projects/PPSmg/tools/usearch/usearch5.2.32_i86linux32'
    usearch6 = '/net/metagenomics/projects/PPSmg/tools/usearch/usearch6.0.307_i86linux32'
    # dnaClust = '/net/metagenomics/projects/PPSmg/tools/dnaclust_64bit_parallel/dnaclust'

    # merge sequences from multiple files
    if mergeS:
        mergeSequences(mapFilePathList, fastaFilePathList, mergedDir)

    if sortS:
        sortSeqDesc(usearch5, usearch6, mergedDir, sortedDir, mapFilePathList)

    # sort and cluster sequences
    if clusterS:
        toCentroids(sortedDir, centroidsDir, mapFilePathList)
        move(centroidsDir)

    if filterOutSeq:
        taxonIdSet = getAllTaxonIdSet(mapFilePathList)
        srcDir = '/net/refdata/static/nonredundant-microbial_20140513/nobackup/centroids'  # '/local/igregor/ref_20121122/nobackup/centroids_1_0'
        dstDir = '/net/refdata/static/nonredundant-microbial_20140513/nobackup/centroids_noplasmids'  # '/local/igregor/ref_20121122/nobackup/centroids_1_0_no_plasmids'
        notAllowedSet = set(csv.getColumnAsList('/net/refdata/static/nonredundant-microbial_20140513/nobackup/plasmids_accessions.txt', colNum=0))  # /local/igregor/ref_20121122/nobackup/plasmid_accessions2.txt
        filterOutSequencesBatch(taxonIdSet, srcDir, dstDir, notAllowedSet)




def move(dirFasta):
    for f in os.listdir(dirFasta):
        src = os.path.join(dirFasta, f)
        dst = str(src.split('.')[0] + '.1.fna')
        shutil.move(src, dst)




# MAIN
if __name__ == "__main__":
    _main()
    # move('/net/refdata/static/nonredundant-microbial_20140513/nobackup/centroidsmv')

    #getSeeds('/Users/ivan/Documents/nobackup/stringMatchTest/888741.fna',
    #         '/Users/ivan/Documents/nobackup/stringMatchTest/OUT_888741.fna')