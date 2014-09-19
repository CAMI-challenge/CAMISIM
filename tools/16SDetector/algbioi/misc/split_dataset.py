#!/usr/bin/env python

import os
import random

from algbioi.com.csv import OutFileBuffer
from algbioi.com.fasta import fastaFileToDict


def splitSet1(scaffFilePath, contigFilePath, outMapFilePath, settings):
    settingsList = [] #touples [coef, i0, i1]
    s = settings.split(',')
    for j in s:
        x, i = j.split(':')
        i0, i1 = i.split('-')
        settingsList.append([x, i0, i1])

    seqNameToSeq = fastaFileToDict(scaffFilePath)
    contigCounter = 0
    mapBuff = OutFileBuffer(outMapFilePath)
    contigBuff = OutFileBuffer(contigFilePath)

    for seqName in seqNameToSeq:
        seq = seqNameToSeq[seqName]
        length = len(seq)
        for tuple in settingsList:
            coef, i0, i1 = tuple
            count = int(float(coef)*(length/1000))
            l = random.randint(int(i0)*1000, int(i1)*1000)
            for i in range(count):
                start = random.randint(0, (length-l))
                end = start + l
                subSeq = seq[start:(end+1)]
                #print len(subSeq)
                contigName = 'contig' + str(contigCounter)
                contigBuff.writeText(str('>' + contigName + '\n' + subSeq + '\n'))
                mapBuff.writeText(str(seqName + '\t' + contigName + '\n'))
                contigCounter += 1

    mapBuff.close()
    contigBuff.close()


def splitSet2(scaffFilePath, contigFilePath, outMapFilePath, subSeqCount, subSeqOverlapp=None):

    contigCounter = 0
    seqNameToSeq = fastaFileToDict(scaffFilePath)
    mapBuff = OutFileBuffer(outMapFilePath)
    contigBuff = OutFileBuffer(contigFilePath)

    for seqName in seqNameToSeq:
        seq = seqNameToSeq[seqName]
        length = len(seq)
        partLen = length/subSeqCount
        for i in range(subSeqCount):
            if i == 0:
                start = 0
            else:
                start = i*partLen
            if i == (subSeqCount - 1):
                end = length
            else:
                end = min(start + partLen, length)

            if subSeqOverlapp != None:
                start -= random.randint(subSeqOverlapp)
                start = max(start, 0)
                end += random.randint(subSeqOverlapp)
                end = min(end, length)

            subSeq = seq[start:end]
            #print length, len(subSeq)
            contigName = 'contig' + str(contigCounter)
            contigBuff.writeText(str('>' + contigName + '\n' + subSeq + '\n'))
            mapBuff.writeText(str(seqName + '\t' + contigName + '\n'))
            contigCounter += 1

    mapBuff.close()
    contigBuff.close()




def test():
    scaffFilePath = '/Users/ivan/Documents/work/binning/data/CowRumen/assembly/cow_rumen_fragmented_velvet_assembly_scaffolds.fas'
    contigFilePath = '/Users/ivan/Documents/work/binning/data/CowRumen/assembly/cow_rumen_artificial_contigs_div2.fas'
    outMapFilePath = '/Users/ivan/Documents/work/binning/data/CowRumen/assembly/cow_rumen_artificial_contigs_map_div2.csv'
    #settings = '0.2:4-6'

    #splitSet1(scaffFilePath, contigFilePath, outMapFilePath, settings)
    splitSet2(scaffFilePath, contigFilePath, outMapFilePath, 2, subSeqOverlapp=None)


def main():
    pass

if __name__ == "__main__":
    if os.getcwd() == '/Users/ivan/Documents/work/python/workspace/pPPS/src':
        test()
    else:
        main()