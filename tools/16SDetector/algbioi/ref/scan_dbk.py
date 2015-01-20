#!/usr/bin/env python

import os
import sys
import signal
import string
import subprocess
from sets import Set

from Bio import SeqIO
from Bio.Seq import Seq

from algbioi.com.csv import getMapping
from algbioi.com.csv import OutFileBuffer
from algbioi.com.csv import getColumnAsList
from algbioi.com.fasta import fastaFileToDict
from algbioi.com.fasta import getSequenceToBpDict
from algbioi.com.fasta import getSequencesToList
from algbioi.com.fasta import filterOutSequences
from algbioi.com.common import noNewLine


def printStatDbk():
    """
        Print statistics about one DBK file.
    """
    seqIdSet = Set([])
    taxonSet = Set([])
    cumulativeLen = 0
    recordCount = 0
    zeros = 0
    #
    for record in SeqIO.parse(sys.stdin, "genbank"):
        recordCount += 1
        seqId = record.id

        if seqId in seqIdSet:
            print seqId, 'already in set', seqId

        seq = str(record.seq)
        cumulativeLen += len(seq)

        if len(string.replace(noNewLine(seq),'N','')) == 0:
            zeros += 1

        taxonId = None

        for feature in record.features:
            if feature.type == "source":
                for xrefentry in feature.qualifiers["db_xref"]:
                    ( key, val ) = xrefentry.split(":")
                    if key == "taxon":
                        taxonId = int(val)
                        break
            if taxonId != None:
                break


        if taxonId == None:
            print 'could not find taxonId for', seqId
        else:
            taxonSet.add(taxonId)

    print 'record count', recordCount
    print 'seq count', len(seqIdSet)
    print 'taxon id count', len(taxonSet)
    if len(seqIdSet) > 0:
        print 'avg. seq. len', cumulativeLen/len(seqIdSet)
    print 'zeros', zeros


#Find out which sequences are just full of zeros!


def findPlasmids(outPlasmidFilePath):
    """
        Read sequence descriptions from a DBK files (stdin), output sequence ids (record.id) if the corresponding
        description contain "plasmid".
        Plasmids can be also within the sequences!
    """
    #append to a file if it already exists
    if os.path.isfile(outPlasmidFilePath):
        outFileMode='a'
    else:
        outFileMode='w'
    outBuffer = OutFileBuffer(outPlasmidFilePath, bufferText = False, fileOpenMode=outFileMode)

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
    taxonIdToOutBuffer = dict([])
    seqIdSet = Set([])

    totalSeqCount = 0
    totalStoredSeqCount = 0
    totalIdenticalSeqCount = 0

    for mapFilePath, fastaFilePath in zip(mapFilePathList, fastaFilePathList):
        print 'processing', mapFilePath, fastaFilePath
        seqCount = 0
        storedSeqCount = 0

        seqIdToSeq = fastaFileToDict(fastaFilePath)
        seqIdToNcbidList = getMapping(mapFilePath, 0, 1, sep='\t', comment = '#')

        for seqId, seq in seqIdToSeq.iteritems():
            seqCount += 1
            if seqId in seqIdSet:
                totalIdenticalSeqCount += 1
                continue
            else:
                seqIdSet.add(seqId)

            taxonId = seqIdToNcbidList[seqId][0]

            if taxonId not in taxonIdToOutBuffer:
                outBuffer = OutFileBuffer(os.path.join(outputDir,str(str(taxonId) + '.fna')))
                taxonIdToOutBuffer[taxonId] = outBuffer

            taxonIdToOutBuffer[taxonId].writeText(str('>' + seqId + '\n' + seq + '\n'))
            taxonIdToOutBuffer[taxonId].close()
            storedSeqCount += 1

            if len(string.replace(noNewLine(seq),'N','')) == 0:
                print 'zeros', seqId, fastaFilePath, len(noNewLine(seq))


        print 'totalSeq, storedSeq', seqCount, storedSeqCount
        totalSeqCount += seqCount
        totalStoredSeqCount += storedSeqCount


    print 'totalSeqCount, totalStoredSeqCount, totalIdenticalSeqCount', totalSeqCount, totalStoredSeqCount, totalIdenticalSeqCount

    print 'sequences merged'


def getMaxLen(fastaFilePath):
    """
        Gets the length of the sequence that has maximum length in a fasta file.
    """
    max = 0
    for val in getSequenceToBpDict(fastaFilePath).itervalues():
        if max < int(val):
            max = int(val)
    return max


def toCentroids(cdHit, usearch1, usearch5, usearch6, identity, mergedDir, sortedDir, centroidsDir, clustersDir, taxonIdSet, sortS=True, clusterS=True):
    """
        @param identity: percentage identity
    """

    inFilePathList = []
    sortedFilePathList = []
    centroidsFilePathList = []
    clustersFilePathList = []
    for taxonId in taxonIdSet:
        inFilePathList.append(os.path.join(mergedDir,str(str(taxonId) + '.fna')))
        sortedFilePathList.append(os.path.join(sortedDir,str(str(taxonId) + '.fna')))
        centroidsFilePathList.append(os.path.join(centroidsDir,str(str(taxonId) + '.fna')))
        clustersFilePathList.append(os.path.join(clustersDir,str(str(taxonId) + '.uc')))

    #sort sequences, longer first
    if sortS:
        for inFile, outFile in zip(inFilePathList, sortedFilePathList):
            sortCmd = str(usearch5 + ' -sort ' + inFile + ' -output ' + outFile + ' --maxlen ' + str(int(getMaxLen(inFile) + 1000)) + ' --minlen 10')
            #sortCmd = str(usearch6 + ' -sortbylength ' + inFile + ' -output ' + outFile) # + ' -maxqt ' + str(int(getMaxLen(inFile))))
            print sortCmd
            sortProc = subprocess.Popen(sortCmd, shell=True, cwd = os.path.dirname(usearch6) ,bufsize=-1)
            sortProc.wait()
            if sortProc.returncode != 0:
                sys.stderr.write(str('Cmd: ' + sortCmd + ' ended with return code: ' + str(sortProc.returncode) + '\n'))
                return
        print 'sequences sorted'

    #find centroids for each cluster within each sorted file
    if clusterS:
        for inFile, outFile, clusterFile in zip(sortedFilePathList, centroidsFilePathList, clustersFilePathList):

            getSeeds(inFile, outFile)

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
        @param inSortedFasta: DNA sequences sorted according to the sequence length
        @param outSeedsFasta: a fasta file that contains all seeds
    """
    out = OutFileBuffer(outSeedsFasta)
    seedList = []
    seqList = getSequencesToList(inSortedFasta)

    for seqId, seq in seqList:
        seq = string.upper(seq)

        newSeed = True
        for seedSeq in seedList:

            if len(seedSeq) < len(seq):
                continue

            #if bool(re.search(seq, seedSeq, re.I)) or bool(re.search(str(Seq(seq).reverse_complement()), seedSeq, re.I)):
            if seq in seedSeq or str(Seq(seq).reverse_complement()) in seedSeq:
                newSeed = False
                break

        if newSeed:
            #print 'new seed:', seqId
            seedList.append(seq)
            out.writeText(str('>' + seqId + '\n' + seq + '\n'))
        #else:
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
        srcFilePath = os.path.join(srcDir,str(str(taxonId) + '.fna'))
        dstFilePath = os.path.join(dstDir,str(str(taxonId) + '.fna'))

        seqIdDict = getSequenceToBpDict(srcFilePath)
        allowedNamesSet = Set()
        for id in seqIdDict.iterkeys():
            if id not in notAllowedSeqIdSet:
                allowedNamesSet.add(id)

        filterOutSequences(srcFilePath, dstFilePath, allowedNamesSet)



def main():
    mergeS = False
    sortS = False
    clusterS = False
    filterOutSeq = True

    # handle broken pipes
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    #printStatDbk()
    #checkForPlasmids()

    mapFilePathList = ['/local/johdro/refdata/static/ncbi-genomes-bacteria_20121122/nobackup/dna_acc.nuc.tax',
                       '/local/johdro/refdata/static/ncbi-draftgenomes-bacteria_20121122/nobackup/dna-contigs_acc.nuc.tax',
                       '/local/johdro/refdata/static/ncbi-draftgenomes-bacteria_20121122/nobackup/dna-scaffolds_acc.nuc.tax',
                       '/local/johdro/refdata/static/ncbi-hmp_20121016/nobackup/dna-contigs_acc.nuc.tax',
                       '/local/johdro/refdata/static/ncbi-hmp_20121016/nobackup/dna-scaffolds_acc.nuc.tax',
                       '/local/johdro/refdata/static/ncbi-refseq-microbial_56/nobackup/dna_acc.nuc.tax'
                       ]

    fastaFilePathList = ['/local/johdro/refdata/static/ncbi-genomes-bacteria_20121122/nobackup/dna_acc.nuc.fna',
                         '/local/johdro/refdata/static/ncbi-draftgenomes-bacteria_20121122/nobackup/dna-contigs_acc.nuc.fna',
                         '/local/johdro/refdata/static/ncbi-draftgenomes-bacteria_20121122/nobackup/dna-scaffolds_acc.nuc.fna',
                         '/local/johdro/refdata/static/ncbi-hmp_20121016/nobackup/dna-contigs_acc.nuc.fna',
                         '/local/johdro/refdata/static/ncbi-hmp_20121016/nobackup/dna-scaffolds_acc.nuc.fna',
                         '/local/johdro/refdata/static/ncbi-refseq-microbial_56/nobackup/dna_acc.nuc.fna'
                         ]
    mergedDir = '/local/igregor/ref_20121122/nobackup/merged'

    #merge sequences from multiple files
    if mergeS:
        mergeSequences(mapFilePathList, fastaFilePathList, mergedDir)

    #sort and cluster sequences
    if sortS or clusterS or filterOutSeq:
        taxonIdSet = Set([])
        for mapFile in mapFilePathList:
            taxonIdSet = taxonIdSet.union(Set(getColumnAsList(mapFile, colNum=1, sep='\t')))

    if sortS or clusterS:
        identity = 1.0
        sortedDir = '/local/igregor/ref_20121122/nobackup/sorted'
        centroidsDir = '/local/igregor/ref_20121122/nobackup/centroids_1_0'
        clustersDir = '/local/igregor/ref_20121122/nobackup/clusters_1_0'
        usearch1 = '/net/programs/Debian-6.0-x86_64/uclust-1.1.579/uclust'
        usearch5 = '/net/metagenomics/projects/PPSmg/tools/usearch/usearch5.2.32_i86linux32'
        usearch6 = '/net/metagenomics/projects/PPSmg/tools/usearch/usearch6.0.307_i86linux32'
        #dnaClust = '/net/metagenomics/projects/PPSmg/tools/dnaclust_64bit_parallel/dnaclust'
        cdHit = '/net/metagenomics/projects/PPSmg/tools/cd-hit-v4.6.1-2012-08-27$ ./cd-hit-est'

        toCentroids(cdHit, usearch1, usearch5, usearch6, identity, mergedDir, sortedDir, centroidsDir, clustersDir, taxonIdSet, sortS=sortS, clusterS=clusterS)

    if filterOutSeq:
        srcDir = '/local/igregor/ref_20121122/nobackup/centroids_1_0'
        dstDir = '/local/igregor/ref_20121122/nobackup/centroids_1_0_no_plasmids'
        notAllowedSet = Set(getColumnAsList('/local/igregor/ref_20121122/nobackup/plasmid_accessions2.txt', colNum=0))
        filterOutSequencesBatch(taxonIdSet, srcDir, dstDir, notAllowedSet)


# MAIN
if __name__ == "__main__":
    main()

    #getSeeds('/Volumes/Macintosh HD/Users/ivan/Documents/nobackup/stringMatchTest/888741.fna',
    #         '/Volumes/Macintosh HD/Users/ivan/Documents/nobackup/stringMatchTest/OUT_888741.fna')