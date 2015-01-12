#!/usr/bin/env python

import os
import re

from algbioi.com import common


def substCol(tabFile, mapFile, colKeySrc, colValueDst, delColKeys, shuffleList=None, outFile=None):
    """
        Replaces values in a column in a tabFile by the values from the mapFile

        @param tabFile: file with tab separated columns
        @param mapFile: file with two columns key -> value
        @param colKeySrc: column of the first file that contains keys (number started from 0)
        @param colValueDst: column of the output file that contains values
        @param delColKeys: True or False - should the key column be included in the output file
        @param shuffleList: numbers of columns that will be used in the given order from the tabFile (starts from 0; -1 ~ dummy value 0)
    """

    try:
        f = open(os.path.normpath(mapFile),'r')
    except Exception:
            print "Cannot open file:", mapFile
            raise
    else:
        map = dict([])
        for line in f:
            key = common.noNewLine(re.sub(r'([^\t]+)\t[^\t]+', r'\1' , line))
            val = common.noNewLine(re.sub(r'[^\t]+\t([^\t]+)', r'\1' , line))
            map[key] = val

    try:
        f = open(os.path.normpath(tabFile),'r')
    except Exception:
            print "Cannot open file:", tabFile
            raise
    else:
        tabFileLinesOrig = []

        for line in f:
            list = common.noNewLine(line).split("\t")
            if not re.match('#',list[0]):
                tabFileLinesOrig.append(list)
            else:
                print '#line ignored: ', common.noNewLine(line)

        if shuffleList == None:
            tabFileLines = tabFileLinesOrig
        else:
            for i in range(0,len(shuffleList)):
                if shuffleList[i] == colKeySrc:
                    colKeySrc = i
                    break

            tabFileLines = []
            for lineOrig in tabFileLinesOrig:
                line = []
                for i in range(0,len(shuffleList)):
                    idxOrig = shuffleList[i]
                    if idxOrig < 0:
                        line.append(0)
                    else:
                        line.append(lineOrig[idxOrig])
                tabFileLines.append(line)


        for list in tabFileLines:
            key = list[colKeySrc]
            if key not in map:
                print str('#key: "' + key + '" is not contained in the map file')
                continue
            val = map[key]
            if delColKeys:
                for i in range(colKeySrc + 1, len(list)):
                    list[i-1] = list[i]
                del list[len(list)-1]

            list.append(None)
            if colValueDst < len(list):
                for i in reversed(range(colValueDst, len(list)-1)):
                    list[i+1] = list[i]
            list[colValueDst] = val


        strBuf = ''
        for list in tabFileLines:
            for i in range(0,len(list)):
                if i == 0:
                    strBuf += str(list[i])
                else:
                    strBuf += str('\t' + str(list[i]))
            strBuf += '\n'

        if outFile == None:
            print strBuf
        else:
            try:
                f = open(os.path.normpath(outFile),'w')
            except Exception:
                print "Cannot open file for writing:", outFile
                raise
            else:
                f.write(strBuf)




if __name__ == "__main__":
    #tabFile = 'D:/A_Phylo/A_Metagenomic/pPPS/workspace/pPPS/data/contig_tab2.blastn'
    tabFile = 'D:/A_Phylo/A_Metagenomic/pPPS/workspace/pPPS/data/contig_alignments2.blastn10'
    #tabFile ='/AM/metagenomic/work/projects/pPPS/tests/SRM01/johdro/contig_alignments2.blastn'
    mapFile = 'D:/A_Phylo/A_Metagenomic/pPPS/workspace/pPPS/data/gi_taxid_species-refseq.dmp'
    #mapFile = '/AM/metagenomic/work/projects/pPPS/tools/johdro/gi_taxid_species-refseq.dmp'
    #fileOut = 'D:/A_Phylo/A_Metagenomic/pPPS/workspace/pPPS/data/contig_tab.blastn'
    fileOut = None
    #fileOut = '/AM/metagenomic/work/projects/pPPS/tests/SRM01/johdro/contig_tabbed2.blastn'

    colKeySrc = 3
    colValueDst = 12
    delColKeys = False
    shuffleList = [0,3,6,-1,-1,-1,1,2,4,5,8,7]
    #shuffleList = None
    substCol(tabFile, mapFile, colKeySrc, colValueDst, delColKeys, shuffleList, fileOut)