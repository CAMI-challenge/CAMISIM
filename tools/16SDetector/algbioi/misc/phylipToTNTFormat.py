#!/usr/bin/env python

import re
import sys
import argparse

from algbioi.com.csv import forEachLine
from algbioi.com.csv import OutFileBuffer


def main():

    parser = argparse.ArgumentParser(description='''Transforms a file in the PHYLIP format to a file in the TNT basic format (~Hennig86)''',
                                     epilog='''''')

    parser.add_argument('-i', '--input-file-phylip', nargs=1, required=True,
                        help='Input file in the phylip format', metavar='input.phylip',
                        dest='inputPhylip')

    parser.add_argument('-o', '--output-file-tnt', nargs=1, required=True,
                        help='Output file in the tnt format', metavar='output.tnt',
                        dest='outputTNT')

    parser.add_argument('-s', '--sort-n', action='store_true',
                        help='Sort sequences according to the numbers that are contained in their names (e.g. "Tax12" is sorted according to "12").',
                        dest='s')

    args = parser.parse_args()

    #print str(args.inputPhylip[0]), str(args.outputTNT[0])

    #parse sequence file
    phylipParsedData = PhylipParser()
    forEachLine(args.inputPhylip[0], phylipParsedData)
    if not phylipParsedData.isDataValid:
        return

    if args.s:
        phylipParsedData.sorSeqNameList()

    #store sequences in the TNT basic format
    outBuffer = OutFileBuffer(args.outputTNT[0])
    storeInTNTBasicFormat(phylipParsedData, outBuffer)
    outBuffer.close()


class PhylipParser():
    def __init__(self):
        self.charCount = None
        self.seqCount = None
        self.seqCounter = None
        self.seqNameList = []
        self.seqNameToSeqDict = dict([])

    def parse(self, line):
        if (self.seqCounter == None) and (re.match(r'^[ \t]*[0-9]+[ \t]+[0-9]+', line)):
            self.seqCount = int(re.sub('^[ \t]*([0-9]+)[ \t]+[0-9]+', r'\1' , line))
            self.charCount = int(re.sub('^[ \t]*[0-9]+[ \t]+([0-9]+)', r'\1' , line))
            self.seqCounter = 0
            #print self.seqCount, self.charCount, self.seqCounter
        if (self.seqCounter != None) and (self.seqCounter < self.seqCount) and re.match(r'^[^ \t]+[ \t]+[^ \t\n\r]{' + str(self.charCount) + '}$', line):
            seqName = re.sub(r'^([^ \t]+)[ \t]+[^ \t\n\r]{' + str(self.charCount) + '}$', r'\1', line)
            seq = re.sub(r'^[^ \t]+[ \t]+([^ \t\n\r]{' + str(self.charCount) + '})$', r'\1', line)
            self.seqNameList.append(seqName)
            self.seqNameToSeqDict[seqName] = seq
            self.seqCounter += 1
            #print seqName, seq

    def isDataValid(self):
        if self.seqCounter == None or len(self.seqNameList) == 0:
            sys.stderr.write('No data has been read!\n')
            return False
        if self.seqCounter != None and self.seqCounter != self.seqCount:
            sys.stderr.write('The number of the sequences read differ from the number of the sequences declared in the file:'
                             + str(self.seqCounter) + ' ' + str(self.seqCount) + '\n')
            return False
        if self.seqCounter != None and len(self.seqNameToSeqDict) != self.seqCount:
            sys.stderr.write('The number of the sequences stored in the dictionary differ from the number of the sequences declared in the file:'
                             + str(len(self.seqNameToSeqDict)) + ' ' + self.seqCount)
            return False

        return True

    def sorSeqNameList(self):
        self.seqNameList.sort(key=lambda k: int(re.sub(r'[^0-9]*([0-9]*)[^0-9]*',r'\1',k)))

    def getSeqNameList(self):
        return self.seqNameList

    def getSeq(self,seqName):
        return self.seqNameToSeqDict[seqName]


def storeInTNTBasicFormat(data, outBuffer):

    seqNameList = data.getSeqNameList()
    if len(seqNameList) == 0:
        sys.stderr.write('storeInTNTBasicFormat: data is empty!\n')
        return

    tooLongNames = []
    for name in seqNameList:
        if ';' in name or '\'' in name or '+' in name or '-' in name:
            print 'Taxon', name, 'contains a not allowed character: ";", "\\", "+", or "-"'
        if len(name) > 32:
            tooLongNames.append(name)
    if len(tooLongNames) > 0:
        print 'Some ', str(len(tooLongNames)), 'names are too long, use the "TAXNAME +N" command of TNT'

    statesSet = set([])
    for name in seqNameList:
        seq = data.getSeq(name).upper()
        for i in range(len(seq)):
            statesSet.add(seq[i])

    if len(statesSet) > 16:
        sys.stderr.write('There are too many states' + str(len(statesSet)) + ' in the dataset, maximum number of states is 16')
        return

    idxToChar = ['0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
    stateToChar = dict([])
    stateToChar['A'] = idxToChar[0]
    stateToChar['C'] = idxToChar[2]
    stateToChar['G'] = idxToChar[1]
    stateToChar['T'] = idxToChar[3]
    stateToChar['U'] = idxToChar[3]
    stateToChar['-'] = idxToChar[4]
    i = 5
    for state in statesSet:
        if state not in stateToChar:
            stateToChar[state] = idxToChar[i]
            i += 1

    outBuffer.writeText('xread\r\n')
    outBuffer.writeText('\'TNT basic format (Hennig86)\'\r\n')
    outBuffer.writeText(str(len(data.getSeq(seqNameList[0]))) + ' ' + str(len(seqNameList)) + '\r\n')

    for name in seqNameList:
        seq = data.getSeq(name).upper()
        intSeq = []
        for i in range(len(seq)):
            intSeq.append(stateToChar[seq[i]])

        outBuffer.writeText(name + ' ' + "".join(intSeq) + '\r\n')

    outBuffer.writeText(';')




if __name__ == "__main__":
    main()

