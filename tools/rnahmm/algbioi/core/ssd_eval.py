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

import re
import os
import glob

from algbioi.core.taxonomy import Taxonomy
from algbioi.core.sequences import toScafContigMap
from algbioi.com.config import Config
from algbioi.com import common


def ppsOut2Placements(ppsOutFile, scafContigFile=None):
    """
        Transforms a PPS assignments to a list of pairs <contigName, assigned_ncbid>

        @param ppsOutFile: PPS output file where the first column is the contig/scaffold name and the last column is ncbid
        @param scafContigFile: scaffold contig mapping (tab separated) if None then all sequences are considered as contigs

        @return: list of pairs <contigName, assigned_ncbid>
    """

    #print 'ppsOut2Placements ppsOutFile:', ppsOutFile
    #print 'ppsOut2Placements scafContigFile:', scafContigFile

    if scafContigFile != None:
        scafToContigs = toScafContigMap(scafContigFile)
    else:
        scafToContigs = dict([])

    outList = []
    try:
        f = open(os.path.normpath(ppsOutFile),'r')
    except Exception:
        print "Cannot open file:", ppsOutFile
        raise
    else:
        lineCounter = 0
        for line in f:
            lineCounter += 1
            line = common.noNewLine(line)
            name = re.sub(r'^([^ \t]+)[ \t]+.*[0-9]+[ \t]*$',r'\1' ,line)
            try:
                ncbid = int(re.sub(r'^[^ \t]+.*[ \t]+([0-9]+)[ \t]*$',r'\1' ,line))
            except Exception:
                try:
                    ncbid = abs(int(re.sub(r'^[^ \t]+.*[ \t]+(-1)[ \t]*$',r'\1' ,line)))
                except Exception:
                    print 'ppsOut2Placements: cannot parse placement for line nr:', lineCounter, 'line:', line
                    raise

            if name in scafToContigs:
                contigsList = scafToContigs[name]
                for contig in contigsList:
                    outList.append([contig, ncbid])
                    #print ':',contig,ncbid
            else:
                outList.append([name, ncbid])
                #print '',name,ncbid

    return outList


def ssd2Placements(ssdDir, scafContigFile=None):
    """
        Transforms sample specific data to placements. Sequences` names are not allowed to have gaps ' '

        @param ssdDir: directory that contains sample specific data
        @param scafContigFile: scaffold contig mapping (tab separated) if None then all sequences are considered as contigs

        @return: list of pairs <contigName, assigned_ncbid>
    """

    #collect map: scaffold -> list of contigs
    if scafContigFile != None:
        scafToContigs = toScafContigMap(scafContigFile)
    else:
        scafToContigs = dict([])

    outList = []
    placedContigs = set([])

    for filePath in glob.glob(os.path.join(os.path.normpath(ssdDir),r'*.f[an][sa]')):
        ncbid = int(re.sub(r'^.*[^0-9]([0-9]+)\.[0-9]+\.f[an][sa]$',r'\1' ,filePath)) #int
        try:
            f = open(os.path.normpath(filePath),'r')
        except Exception:
            print "Cannot open file:", filePath
            raise
        else:
            for line in f:
                line = common.noNewLine(line)
                if re.match('>', line):
                    name = re.sub(r'^([^ \t]+)[ \t]*.*$',r'\1',line.replace('>',''))
                    if name in scafToContigs:
                        contigsList = scafToContigs[name]
                    else:
                        contigsList = [name]
                    for contig in contigsList:
                        if contig in placedContigs:
                            print str('contig "' + contig + '" has already been placed')
                        else:
                            placedContigs.add(contig)
                            outList.append([contig, ncbid])
        #count also BP for each contig!!!

    return outList


def cmpPlacements(refPlacement, placement, taxonomy, taxonomicRanks):
    """
        Compare two placements

        @param refPlacement: reference placement is considered to be true, list of pairs <contigName, assigned_ncbid>
        @param placement: list of pairs <contigName, assigned_ncbid>
        @param taxonomy: NCBI taxonomy

        @return: list of n-tuples where each n-touple contains comparison of two placements of one contig
    """

    placementDict = {}  # contig -> ncbid
    for p in placement:
        placementDict[p[0]] = p[1]

    outList = []
    #for all reference placements
    for p in refPlacement:
        contig = p[0]
        ncbid = int(p[1])

        if contig in placementDict:
            pNcbid = int(placementDict[contig])
        else:
            pNcbid = None # placement assignment

        #status = 'N'
        if pNcbid == None:
            pNcbid = 1

        #    outList.append([contig, ncbid, status, None, None, None, None, None])
        #else:
        #get taxonomy path
        pathDictRef = taxonomy.getPathToRoot(ncbid)
        pathDict = taxonomy.getPathToRoot(pNcbid)
        #if  pathDictRef == None:
        #    print 'pathDictRef is None'
        #    pathDictRef = dict([])
        #if pathDict == None:
        #    print 'pathDict is None'
        #    pathDict = dict([])
        #if (len(pathDictRef) == 0 or len(pathDict) == 0):
        #    print "LCA A IS NONE!"
        #    pathDictLCA = dict([])
        #else:
        #    print ncbid, pNcbid
        pathDictLCA = taxonomy.getPathFromLowestCommonAncestorToRoot([ncbid,pNcbid])

        if pathDictLCA is not None:
            pathDictLCALen = len(pathDictLCA)
        else:
            pathDictLCALen = 0

        if pathDictRef is not None:
            pathDictRefLen = len(pathDictRef)
        else:
            pathDictRefLen = 0

        if pathDict is not None:
            pathDictLen = len(pathDict)
        else:
            pathDictLen = 0

        if (pathDictLCALen == pathDictRefLen) or (pathDictLCALen == pathDictLen):
            samePath = True
        else:
            samePath = False
        if pathDictLen < pathDictRefLen:
            status = 'H'
        elif pathDictLen > pathDictRefLen:
            status = 'L'
        else:
            status = 'S'

        #dist = abs((len(pathDictLCA) - len(pathDict))) + abs((len(pathDictLCA) - len(pathDictRef)))
        dist = abs(pathDictLCALen - pathDictLen) + abs(pathDictLCALen - pathDictRefLen)
        if dist == 0:
            status = 'M'
            assert ncbid == pNcbid, str('The ncbi values don`t match %s %s' % (ncbid, pNcbid))

        ncbiName = taxonomy.getNcbidToName(ncbid)
        pNcbiName = taxonomy.getNcbidToName(pNcbid)

        outList.append([contig, ncbid, status, samePath, dist, ncbiName, pNcbiName, pNcbid])

    return outList


def filterCmpList(cmpListR, rankCut, maxDist, taxonomy, outputInverse = False):
    """
        Filters out all placements that are mismatches and all placements to the higher ranks that are
        lower than certain rank (e.g. bacteria or the root) and all placements that are more distant from the
        reference placement than certain threshold.

        @param cmpListR: comparison with the reference placement returned by the method "cmpPlacements"
        @param rankCut: if a contig is assigned on the same path higher but to this rank or higher, it will be filtered out (i.e. 0 for superkingdom)
        @param maxDist: if a contig is assigned on the same path higher, it will be filtered out if the distance is higher than this value
        @param taxonomy
        @param outputInverse: if True, it outputs all entries that would be filtered out if False
    """

    #filter out mismatches (only placements that lie on the same path remain)
    resultList = []
    for entry in cmpListR:
        passFilter = False
        if entry[3]: #same path
            if (entry[2] == 'M') or (entry[2] == 'L'): #match or lower rank
                passFilter = True
            else:
                assert entry[2] == 'H', 'This entry must lie on the same path on a higher rank'  #match higher rank

                #placement ncbid
                pNcbid = entry[7]

                #placement taxonomy path
                pTaxPathDict = taxonomy.getPathToRoot(pNcbid)
                if (pTaxPathDict == None):
                    #sys.stderr.write("filterCmpList: can`t find ncbid: " + str(pNcbid) + "\n")#!!!pNcbid is 1
                    pTaxPathDictLen = 0
                else:
                    pTaxPathDictLen = len(pTaxPathDict)


                #distance from the reference placement to the placement
                dist = entry[4]

                if (dist <= maxDist) and ((pTaxPathDictLen - 1) > rankCut):
                    passFilter = True

        if (passFilter and (not outputInverse)) or ((not passFilter) and outputInverse):
            resultList.append(entry)

    return resultList


def cmp2Summary(cmpPlacementsList, outFile=None):
    """
        Compute statistics and store the results

        @param cmpPlacementsList: list of contigs computed by: cmpPlacements
        @param outFile if None then summary written to the stdout
    """

    ncbidSet = set([])
    idToPlacementList = dict([]) #map ncbid -> list of corresponding placement lists

    for item in cmpPlacementsList:
        ncbid = item[1]
        ncbidSet.add(ncbid)
        if ncbid in idToPlacementList:
            idToPlacementList[ncbid].append(item)
        else:
            tmp = []
            tmp.append(item)
            idToPlacementList[ncbid] = tmp

    #total
    contigsNumT = 0
    matchesNumT = 0
    samePathNumT = 0
    samePathHigherRanksNumT = 0
    samePathHigherRanksDistT = 0
    samePathLowerRanksNumT = 0
    samePathLowerRanksDistT = 0
    distanceSumT = 0
    diffPathRanksDistT = 0

    buf = ''
    for ncbid in ncbidSet:
        placements = idToPlacementList[ncbid]

        contigsNum = 0 #number of contigs
        matchesNum = 0 #number of matches
        samePathNum = 0 #number of contigs placed to the same path
        samePathHigherRanksNum = 0
        samePathHigherRanksDist = 0
        samePathLowerRanksNum = 0
        samePathLowerRanksDist = 0
        distanceSum = 0 #sum of the distances (placement -> )
        diffPathRanksDist = 0
        for p in placements:
            contigsNum += 1
            if p[2] == 'M':
                matchesNum += 1
            elif p[2] == 'H' and p[3]:
                samePathHigherRanksNum += 1
                samePathHigherRanksDist += p[4]
            elif p[2] == 'L' and p[3]:
                samePathLowerRanksNum += 1
                samePathLowerRanksDist += p[4]
            if p[3]:
                samePathNum += 1

            if not p[3]:
                #if p[4] == None:#correct!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                #    # ['Scaffold_343_4082367.lucy.pga.C1965', 538960, 'N', None, None, None, None, None]
                #    print p
                #    continue
                #else:
                diffPathRanksDist += p[4]

            distanceSum += p[4]

        buf += str(placements[0][5] + ' (' +  str(ncbid)+ '):\n')
        buf += str('contigs: ' + str(contigsNum) + '\n')
        buf += str('same path: ' + str((samePathNum*100.0)/contigsNum) + '%\n')
        buf += str('same path matches: ' + str((matchesNum*100.0)/contigsNum) + '%\n')
        if samePathHigherRanksNum > 0:
            buf += str('same path higher ranks: ' + str((samePathHigherRanksNum*100.0)/contigsNum) + '%\n')
        if samePathLowerRanksNum > 0:
            buf += str('same path lower ranks: ' + str((samePathLowerRanksNum*100.0)/contigsNum) + '%\n')
        buf += str('average distance: ' + str((distanceSum*1.0)/contigsNum) + '\n')
        if samePathHigherRanksNum > 0:
            buf += str('same path higher ranks average distance: ' + str((samePathHigherRanksDist*1.0)/samePathHigherRanksNum) + '\n')
        if samePathLowerRanksNum > 0:
            buf += str('same path lower ranks average distance: ' + str((samePathLowerRanksDist*1.0)/samePathLowerRanksNum) + '\n')
        if (contigsNum - samePathNum) > 0:
            buf += str('diff path average distance: ' + str((diffPathRanksDist*1.0)/(contigsNum - samePathNum)) + '\n')
        buf += '\n'

        #total
        contigsNumT += contigsNum
        matchesNumT += matchesNum
        samePathNumT += samePathNum
        samePathHigherRanksNumT += samePathHigherRanksNum
        samePathHigherRanksDistT += samePathHigherRanksDist
        samePathLowerRanksNumT += samePathLowerRanksNum
        samePathLowerRanksDistT += samePathLowerRanksDist
        distanceSumT += distanceSum
        diffPathRanksDistT += diffPathRanksDist


    buf += str('\n')
    buf += str('Summary\n')
    buf += str('contigs: ' + str(contigsNumT) + '\n')
    buf += str('same path: ' + str((samePathNumT*100.0)/contigsNumT) + '%\n')
    buf += str('same path matches: ' + str((matchesNumT*100.0)/contigsNumT) + '%\n')
    if samePathHigherRanksNumT > 0:
        buf += str('same path higher ranks: ' + str((samePathHigherRanksNumT*100.0)/contigsNumT) + '%\n')
    if samePathLowerRanksNumT > 0:
        buf += str('same path lower ranks: ' + str((samePathLowerRanksNumT*100.0)/contigsNumT) + '%\n')
    buf += str('average distance: ' + str((distanceSumT*1.0)/contigsNumT) + '\n')
    if samePathHigherRanksNumT > 0:
        buf += str('same path higher ranks average distance: ' + str((samePathHigherRanksDistT*1.0)/samePathHigherRanksNumT) + '\n')
    if samePathLowerRanksNumT > 0:
        buf += str('same path lower ranks average distance: ' + str((samePathLowerRanksDistT*1.0)/samePathLowerRanksNumT) + '\n')
    if (contigsNumT - samePathNumT) > 0:
        buf += str('diff path average distance: ' + str((diffPathRanksDistT*1.0)/(contigsNumT - samePathNumT)) + '\n')
    buf += '\n'

    try:
        f = open(os.path.normpath(outFile),'w')
        f.write(buf)
        #for itemList in cmpList:
        #    for item in itemList:
        #        f.write(str(item) + '\t')
        #    f.write('\n')
    except Exception:
        print "Cannot open file:", outFile
        raise
    finally:
        f.close()


if __name__ == "__main__":
    config = Config(open(os.path.normpath('D:\\A_Phylo\\A_Metagenomic\\pPPS\\workspace\\pPPS\\config01.cfg')), 'pPPS')
    databaseFile = os.path.normpath(config.get('databaseFile'))
    taxonomicRanks = config.get('taxonomicRanks').split(',')
    t = Taxonomy(databaseFile, taxonomicRanks)
    ppsOutFile = 'D:/A_Phylo/A_Metagenomic/reindeer/predictions/pps05/contigsOut/SRM_Large_Contigs_namesOnly.fna.out'
    scafContigFile = 'D:/A_Phylo/A_Metagenomic/reindeer/data/scaffolds-contigs.tab'
    placement = ppsOut2Placements(ppsOutFile, scafContigFile)
    #for pair in placement:
    #    print pair[0], pair[1]
    ssdDir = 'D:/A_Phylo/A_Metagenomic/reindeer/predictions/pps05/sampleSpecificData/'
    refPlacement = ssd2Placements(ssdDir, scafContigFile)
    #for pair in refPlacement:
    #    print pair[0], pair[1]
    cmpList = cmpPlacements(refPlacement, placement, t, taxonomicRanks)
    #for itemList in cmpList:
    #    for item in itemList:
    #        print ' ', item,
    #    print ''
    outFile = 'D:/A_Phylo/A_Metagenomic/reindeer/predictions/pps05/cmp_summary.txt'
    cmp2Summary(cmpList, outFile)
