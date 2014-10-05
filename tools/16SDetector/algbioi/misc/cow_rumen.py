#!/usr/bin/env python

"""
    Contains methods that handle the Cow Rumen data set.
"""

import os
import glob
import re

from algbioi.com.csv import forEachLine
from algbioi.com.csv import OutFileBuffer
from algbioi.com.taxonomy_ncbi import TaxonomyNcbi


def getBinToOrderDict():
    """
        Maps each WGS bin to an ncbid at the rank order.
    """
    return {'AC2a':171549, 'ADa':29, 'AFa':136, 'AGa':171549, 'AH': 171549, 'AIa':186802,
            'AJ':171549, 'AMa':136, 'AN':186802, 'APb':186802, 'AQ':171549, 'AS1a':186802,
            'ATa':186802, 'AWa':186802, 'BOa':186802}


class BinParser():
    """
        Parse each bin and write sample output files to the outBuffers.
    """
    def __init__(self, outBufferOrderIds, outBufferWgsIds, ncbid, wgsName):
        self.outBufferOrderIds = outBufferOrderIds
        self.outBufferWgsIds = outBufferWgsIds
        self.ncbid = ncbid
        self.wgsName = wgsName

    def parse(self, line):
        if re.match(r'>.*', line):
            name = re.sub(r'^>([^ ]+)$',r'\1', line)
            self.outBufferOrderIds.writeText(str(name + '\t' + str(self.ncbid) + '\n'))
            self.outBufferWgsIds.writeText(str(name + '\t' + self.wgsName + '\n'))


def binsToPredictions():

    outFileOrderIds = '/Users/ivan/Documents/work/binning/data/CowRumen/cowRumenOrderNcbids.txt'
    outFileWgsIds = '/Users/ivan/Documents/work/binning/data/CowRumen/cowRumenWgsIds.txt'

    wgsBinsDir = '/Users/ivan/Documents/work/binning/data/CowRumen/cow_rumen_genome_bins'

    binToOrderDict = getBinToOrderDict()

    outBufferOrderIds = OutFileBuffer(outFileOrderIds)
    outBufferWgsIds = OutFileBuffer(outFileWgsIds)

    #for each file name in directory
    for filePath in glob.glob(os.path.join(os.path.normpath(wgsBinsDir),'*.fas')):
        wgsName = re.sub(r'(^[a-zA-Z0-9]+)-.*',r'\1', os.path.basename(filePath))
        forEachLine(filePath, BinParser(outBufferOrderIds, outBufferWgsIds, binToOrderDict[wgsName], wgsName))
        print wgsName #, filePath

    outBufferOrderIds.close()
    outBufferWgsIds.close()


class DictParser():
    def __init__(self):
        self._dict = dict([])

    def parse(self, line):
        if re.match(r'^[^\t]+\t[^\t]+$', line):
            key = re.sub(r'^([^\t]+)\t[^\t]+$',r'\1',line)
            val = re.sub(r'^[^\t]+\t([^\t]+)$',r'\1',line)
            self._dict[key] = val

    def getDict(self):
        return self._dict


class _BufferedTaxonomy():
    def __init__(self, rank, taxonomy):
        self._rank = rank
        self._taxonomy = taxonomy
        self._dict = dict([])

    def getNcbidAtSpecRank(self, ncbid):
        if ncbid in self._dict:
            return self._dict[ncbid]
        else:
            ncbidRank = self._getNcbidAtRank(ncbid, self._rank, self._taxonomy)
            self._dict[ncbid] = ncbidRank
            return ncbidRank

    def _getNcbidAtRank(self, ncbid, rank, taxonomy):

        if taxonomy.getRank(ncbid) == rank:
            return ncbid
        else:
            parentNcbid = taxonomy.getParentNcbid(ncbid)
            if parentNcbid == None:
                return None
            else:
                return self._getNcbidAtRank(parentNcbid, rank, taxonomy)



def evaluatePPSClusters():

    binToOrderDict = {'AC2a':171549, 'ADa':29, 'AFa':136, 'AGa':171549, 'AH': 171549, 'AIa':186802, 'AJ':171549,
                  'AMa':136, 'AN':186802, 'APb':186802, 'AQ':171549, 'AS1a':186802, 'ATa':186802,
                  'AWa':186802, 'BOa':186802}

    #pred = '/Users/ivan/Documents/work/binning/tests/CowRumen/01/output/cow_rumen_fragmented_velvet_assembly_scaffolds.fas.pOUT'
    pred = '/Users/ivan/Documents/work/binning/tests/CowRumen/02/output/cow_rumen_fragmented_velvet_assembly_scaffolds.fas.pOUT'
    #pred = '/Users/ivan/Documents/work/binning/data/CowRumen/binning/cow_rumen_fragmented_velvet_assembly_scaffolds.fas.nox.fna._just_pred.out'

    ref = '/Users/ivan/Documents/work/binning/data/CowRumen/cowRumenWgsIds.txt'

    #seqToBp = getSequenceToBpDict('/Users/ivan/Documents/work/binning/data/CowRumen/assembly/cow_rumen_fragmented_velvet_assembly_scaffolds.fas')

    databaseFile = "/Users/ivan/Documents/work/binning/taxonomy/ncbi_taxonomy_20110629/ncbitax_sqlite.db"

    taxonomy = TaxonomyNcbi(databaseFile)

    rank = 'genus'
    rankTaxonomy = _BufferedTaxonomy(rank, taxonomy)

    predDict = forEachLine(pred, DictParser()).getDict()
    refDict = forEachLine(ref, DictParser()).getDict()

    refWgsToScaffCount = dict([])
    for s in refDict:
        wgs = refDict[s]
        if wgs in refWgsToScaffCount:
            refWgsToScaffCount[wgs] += 1
        else:
            refWgsToScaffCount[wgs] = 1

    ncbidSet = set([])

    #get bins at rank X
    for key in predDict:
        ncbid = rankTaxonomy.getNcbidAtSpecRank(predDict[key])
        if ncbid == None:
            continue
        ncbidSet.add(int(ncbid))

    wgsNameSet = set(refDict.values())

    print 'Pipeline Bins'
    for ncbid in ncbidSet:
        bufHead = str(taxonomy.getScientificName(ncbid) + '(' + str(ncbid) + ')' + '-----------------------------------------')
        bufBody = ''

        for wgsName in wgsNameSet:
            wgsAssignedNameList = []
            for seqName in predDict:
                if rankTaxonomy.getNcbidAtSpecRank(int(predDict[seqName])) == ncbid:
                    if seqName in refDict and refDict[seqName] == wgsName:
                        wgsAssignedNameList.append(seqName)
            if len(wgsAssignedNameList) > 0:
                bufBody += str(wgsName + ' ' + str(len(wgsAssignedNameList)) + '/' + str(refWgsToScaffCount[wgsName]) + ' (' +
                               str(round((len(wgsAssignedNameList)*100.0)/(refWgsToScaffCount[wgsName]),1)) +
                               ')\n')

        if bufBody != '':
            print bufHead
            print bufBody

    print 'Consistency -----------------------------------'
    print '-----------------------------------'

    for wgsName in wgsNameSet:
        #wgsCounter = 0
        #for seqName in refDict:
        #    if refDict[seqName] == wgsName:
        #        wgsCounter+=1

        bufHead = str(wgsName + ' ' + ' (' + str(refWgsToScaffCount[wgsName]) + ') ----------------------')
        bufBody = ''
        assignedCount = 0

        for ncbid in ncbidSet:
            wgsAssignedNameList = []
            for seqName in predDict:
                if rankTaxonomy.getNcbidAtSpecRank(int(predDict[seqName])) == ncbid:
                    if seqName in refDict and refDict[seqName] == wgsName:
                        wgsAssignedNameList.append(seqName)

            if len(wgsAssignedNameList) > 0:
                bufBody += str(taxonomy.getScientificName(ncbid) + ' (' + str(ncbid) + ') ' + str(round(((100.0*len(wgsAssignedNameList))/refWgsToScaffCount[wgsName]),1)) + '\n')

            assignedCount += len(wgsAssignedNameList)

        if assignedCount < refWgsToScaffCount[wgs]:
            bufBody += str('unassigned ' + str(round(((refWgsToScaffCount[wgsName] - assignedCount)*100.0)/refWgsToScaffCount[wgsName],1)) + '' )
        print bufHead
        print bufBody


    taxonomy.close()



if __name__ == "__main__":
    evaluatePPSClusters()
    #binsToPredictions()