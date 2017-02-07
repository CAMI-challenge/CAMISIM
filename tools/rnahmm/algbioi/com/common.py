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
import re

def removeNonDna(seq):
    """
        Replaces stretches of nonDNA characters with one 'N'.
    """
    return re.sub(r'[^ATGCatgc]+', 'N', seq).upper()


def noNewLine(str):
    """
        Delete all '\n' and '\r' characters in a string.
    """
    return str.replace('\n', '').replace('\r','')


def createTagFilePath(dstDir, fileNameFromPath, tag):
    """
        Returns a path that results from the concatenation of dstDir and a file from filePath where the name of the file
        changes from e.g.: input.n.fas" to "input.n._ids.fas"
    """
    lastDotIdx = fileNameFromPath.rfind('.')
    return os.path.join(os.path.normpath(dstDir),
                        os.path.basename(os.path.normpath(str(fileNameFromPath[0:lastDotIdx] +  fileNameFromPath[lastDotIdx:]
                                                            + '.' + tag))))


def getMothurOutputFilePath(inputFastaFilePath, refTaxonomyFilePath, suffix='.taxonomy'):
    """
        Returns mothur prediction file path that is generated based on the arguments of the Mothur classify command.
    """
    dirName = os.path.dirname(inputFastaFilePath)
    fastaBaseName = os.path.basename(inputFastaFilePath)
    fastaPart = fastaBaseName[0:fastaBaseName.rindex('.')] # without suffix
    taxBaseName = os.path.basename(refTaxonomyFilePath)
    taxPart = taxBaseName[0:taxBaseName.rindex('.')] # without suffix
    if '.' in taxPart:
        taxPart = taxPart[(taxPart.rindex('.') + 1):] # from last comma till the end

    return os.path.join(dirName, str(fastaPart + '.' + taxPart + suffix))


def seqFileCmp(file1, file2):
    """
        Returns true if two files contain the same sequences regardless of their names (and empty spaces).
    """
    seqList1 = seqFileToSeqList(file1)
    seqList2 = seqFileToSeqList(file2)
    if (len(seqList1) != len(seqList2)):
        print "The files contain different number of sequences", len(seqList1), len(seqList2)
        return False
    else:
        print "Number of sequences: ", len(seqList1)

    seqFoundInS2 = set([])
    seqAF = set([])
    idx1 = 0
    for s1 in seqList1:
        idx1 += 1
        idx2 = 0
        for s2 in seqList2:
            idx2 += 1
            if (s1 == s2):
                #print "same:", idx1, idx2
                if s2 not in seqAF:
                    #print idx1, idx2
                #else:
                    seqAF.add(s2)
                if idx2 not in seqFoundInS2:
                    #print "One sequence is in one file more than once", idx1, idx2
                    seqFoundInS2.add(idx2)
                    continue

    s1 = set([])
    s2 = set([])
    for s in seqList1:
        s1.add(s)
    for s in seqList2:
        s2.add(s)
    if len(s1) != len(s2):
        print "The length of unique sequences differ S1:", len(s1), "S2:", len(s2)
    else:
        print "The number of unique sequences is: ", len(s1)


    if len(seqFoundInS2) == len(seqList1):
        print "Both files contain the same sequences"
        return True
    else:
        print "Sequences matches: ", len(seqList1), " Sequences found: ", len(seqFoundInS2)
        return False

def seqFileToSeqList(file):
    seqList = []
    try:
        f = open(os.path.normpath(file),'r')
    except Exception:
        print "Cannot open file:", file
        raise
    else:
        name = ''
        seq = ''
        for line in f:
            line = noNewLine(line)
            if re.match('>', line):
                if seq != '':
                    assert name != ''
                    seqList.append(seq) #store seq
                    seq = ''
                name = line.replace('>','')
            else:
                seq += line
        if seq != '':
            assert name != ''
            seqList.append(seq) #store seq
        return seqList
    finally:
        f.close()



if __name__ == "__main__":
    inputFastaFilePath = '/net/metagenomics/projects/PPSmg/tests/V35/07/working/contigsMappedBlast1000.fna.ids.23S_rRNA.fna'
    refTaxonomyFilePath = '/net/metagenomics/projects/PPSmg/data/silva/lsuparc_silva106_ncbitax.bacteria+archaea.tax'
    print getMothurOutputFilePath(inputFastaFilePath, refTaxonomyFilePath)

    #file = "200643.1.fna"
    #file1="/Users/ivan/Documents/work/binning/database/silva/align/LSURef_106_tax_silva_trunc.fasta"
    #file2="/Users/ivan/Documents/work/binning/database/silva/lsuparc_silva106_ncbitax.bacteria+archaea.fna"
    #print seqFileCmp(file1, file2)
    #seqFileCmp(str("D:\\A_Phylo\\A_Metagenomic\\pPPS\\workspace\\pPPS\\toKaustubh\\variant01\\trainingData\\" + file),
    #           str("D:\\A_Phylo\\A_Metagenomic\\pPPS\\workspace\\pPPS\\toKaustubh\\variant02\\trainingData\\" + file))
    #filePath = "dir/input.n.fas"
    #dstDir = "D:/A_Phylo/A_Metagenomic/pPPS/workspace/pPPS/toKaustubh"
    #tag = "ids"
    #print createTagFilePath(dstDir, filePath, tag)