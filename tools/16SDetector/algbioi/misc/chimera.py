#!/usr/bin/env python

"""To detect Chimeras in the Mercier dataset - is NOT complete."""

from algbioi.com.csv import getColumnAsList
from algbioi.com.csv import getMapping
from algbioi.com.csv import filterOutLines
from algbioi.com.csv import OutFileBuffer
from algbioi.com.fasta import getSequenceToBpDict


def filterAlignments():
    """
        Read in all entries from the alignmentsFile. Copy entries to the filteredAlignmentsFile that were aligned to
        one of the reference genomes.
    """
    genomeNcbidsFile = '/net/metagenomics/projects/PPSmg/data/V35/genome_ncbids.txt'
    giToNcbidFile = '/net/metagenomics/projects/PPSmg/data/V35/chimera/gi_taxid_dna_leaves.dmp'
    alignmentsFile = '/net/metagenomics/projects/PPSmg/data/V35/chimera/nobackup/alignments.csv'
    filteredAlignmentsFile = '/net/metagenomics/projects/PPSmg/data/V35/chimera/nobackup/alignmentsFiltered.csv'

    #read in genome ncbids
    genomeNcbidSet = set(getColumnAsList(genomeNcbidsFile, colNum=0, comment='#'))
    giToNcbidDict = getMapping(giToNcbidFile, 0, 1, sep='\t', comment = '#')
    ncbidToGiDict = getMapping(giToNcbidFile, 1, 0, sep='\t', comment = '#')

    genomeGiIdSet = set([])
    for ncbid in genomeNcbidSet:
        assert ncbid in ncbidToGiDict, str('There is no mapping for ncbid: ' +  str(ncbid))
        #if len(ncbidToGiDict[ncbid]) > 1:
        #    print str('There are more gi ids' + str(ncbidToGiDict[ncbid]) + ' for ncbid: ' +  str(ncbid))
        assert len(ncbidToGiDict[ncbid]) != 0, str('For ncbid:' + str(ncbid) + ' there is no gi id')
        for gi in ncbidToGiDict[ncbid]:
            genomeGiIdSet.add(gi)

    filterOutLines(alignmentsFile, filteredAlignmentsFile, genomeGiIdSet, entryModifyFunction=None, colNum=4, sep='\t', comment='#')


class AlignEntry():
    def __init__(self, line, giToNcbidDict):
        list = line.split('\t')
        self.queryId = int(list[0])
        self.queryBegin = int(list[1])
        self.queryEnd = int(list[2])
        self.queryLen = int(list[3])
        self.refId = int(list[4])
        self.refNcbid = int(giToNcbidDict[str(self.refId)][0])
        self.score = float(list[7])
    def __str__(self):
        return str(str(self.queryId) + ' ' + str(self.queryBegin) + ' ' + str(self.queryEnd) + ' ' + str(self.queryLen) + ' ' + str(self.refId)
                   + ' ' + str(self.score) + ' ' + str(self.refNcbid))
    def getAlignLen(self):
        return self.queryEnd - self.queryBegin + 1


def getStat():
    alignAcceptParam = 0.001
    filteredAlignmentsFile = '/Users/ivan/Documents/work/binning/data/V35/chimera/nobackup/alignmentsFiltered.csv'
    giToNcbidFile = '/Users/ivan/Documents/work/binning/data/V35/chimera/gi_taxid_dna_leaves.dmp'
    giToNcbidDict = getMapping(giToNcbidFile, 0, 1, sep='\t', comment = '#')
    lines = getColumnAsList(filteredAlignmentsFile, colNum=0, comment='#',sep='\n')
    queryToAlign = dict([])
    for i in lines:
        queryId = int(i.split()[0])
        if queryId in queryToAlign:
            queryToAlign[queryId].append(AlignEntry(i,giToNcbidDict))
        else:
            queryToAlign[queryId] = [AlignEntry(i,giToNcbidDict)]

    for i in queryToAlign:
        queryToAlign[i].sort(key=lambda x: x.score, reverse=True)

    inspectList = []
    for queryId in queryToAlign:
        #does the best alignment spans the whole contig
        bestAlign = queryToAlign[queryId][0]
        queryLen = bestAlign.queryLen
        alignLen = bestAlign.queryEnd - bestAlign.queryBegin + 1
        if queryLen != alignLen:
            inspectList.append(queryId)

    print '#inspectList', len(inspectList)
    count = 0
    #
    queryIdToIntervalsList = dict([])
    for queryId in inspectList:
        bestAlign = queryToAlign[queryId][0]
        #begin, end, whoIdx
        intervalsList = [(bestAlign.queryBegin, bestAlign.queryEnd, bestAlign)] #

        for align in queryToAlign[queryId][1:len(queryToAlign[queryId])]:
            #is it`s ncbid different?

            skipAlign = False
            for interval in intervalsList:
                if interval[2].refNcbid == align.refNcbid:
                    skipAlign = True
            if skipAlign:
                continue

            alignBegin = align.queryBegin
            alignEnd = align.queryEnd
            considerAlign = True

            for interval in intervalsList:
                if interval[0] <= alignBegin and alignBegin <= interval[1]:
                    alignBegin = interval[1] + 1 #end

                if interval[0] <= alignEnd and alignEnd <= interval[1]:
                    alignEnd = interval[0] - 1 #begin

                if alignBegin > alignEnd or alignBegin < 1 or alignEnd > bestAlign.queryLen:
                    considerAlign = False
                    break

            if considerAlign and ((alignEnd - alignBegin + 1) >= bestAlign.getAlignLen()*alignAcceptParam):
                intervalsList.append((alignBegin, alignEnd, align)) # it is long enough

        queryIdToIntervalsList[queryId] = intervalsList

    count = 0
    for queryId in inspectList:
        intervalsList = queryIdToIntervalsList[queryId]
        if len(intervalsList) == 1:
            continue
        count += 1
        print queryId, '------------------------------------'
        for interval in intervalsList:
            #intervalsList.sort(key=lambda x: x[0], reverse=False)
            print interval[0], interval[1], interval[2].refNcbid, interval[2].score, interval[2].queryLen

    print "count",count




#    for queryId in inspectList:
#        bestAlign = queryToAlign[queryId][0]
#
#        secondBestAlign = None
#        for align in queryToAlign[queryId][1:len(queryToAlign[queryId])]:
#            if align.refNcbid != bestAlign.refNcbid:
#                secondBestAlign = align
#                break
#
#        if secondBestAlign == None:
#            continue
#
#        if bestAlign.score*0.9 < secondBestAlign.score:
#            #if queryId in notConsider:
#            #    continue
#            for i in queryToAlign[queryId]:
#                print i
#            print '-----------------------------------------'
#            count += 1
#    print 'count', count






    #list = getColumnAsList(filteredAlignmentsFile, colNum=0, sep='\n',comment='#')
    #for i in list:
    #    print i

#transforms a fasta file
def fastaToIntLength():
    #fastaFilePath = '/Users/ivan/Documents/work/binning/data/V35/contigs_1000.txt'
    #outFile = '/Users/ivan/Documents/work/binning/data/V35/distr/contigs_1000_lengths.txt'
    fastaFilePath = ''
    outFile = ''
    seqIdToBp = getSequenceToBpDict(fastaFilePath)
    print 'seq read'
    out = OutFileBuffer(outFile)
    for id in seqIdToBp:
        out.writeText(str(str(seqIdToBp[id]) + '\n'))
    out.close()
    print 'done'



if __name__ == "__main__":
    getStat()
    #filterAlignments()
    #main()