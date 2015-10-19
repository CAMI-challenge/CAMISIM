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
import sys
import re
import glob

from algbioi.core.sequences import seqWeightThenLenCmp
from algbioi.core.ref_seq import RefSequences
from algbioi.com import common


class PPSInput():
    """
    Attributes
    ...

    """
    def __init__(self, sequences, taxonomicRanks, forbiddenDict, exSSDContigNameList, summaryAllFile=None):

        self.sequences = sequences
        self.taxonomicRanks = taxonomicRanks

        self.ncbidToSequences = {}      # map: ncbid -> list of sequences
        self.ncbidToBp = {}             # map: ncbid -> # of BP that map to this ncbi node and its children
        self.ncbidToRankId = {}         # map: ncbid -> rankId
        self.ncbidToName = {}           # map: ncbid -> name
        self.ncbidToParentNcbid = {}    # map: ncbid -> parent_ncbid or None
        self.ncbidToJumpedRanksNum = {} # map: ncbid -> number of ranks that one has to be jumped to get to the next
                                        # defined rank
        #self.rankToNodesDict = {}
        #self.ncbidToSequences = {}
        #map: rank -> list of Nodes (ncbids)
        #map: rank -> set of ncbids

        # transform forbiddenDict to forbidden sets
        forbiddenDictSet = {}
        if forbiddenDict is not None:
            for ncbid in forbiddenDict:
                forbiddenDictSet[ncbid] = set(forbiddenDict[ncbid])

        # transform the list into a set
        if exSSDContigNameList is None:
            exSSDContigNameSet = set()
        else:
            exSSDContigNameSet = set(exSSDContigNameList)

        for sequence in self.sequences.sequences:

            taxonomyPath = sequences.getTaxonomyPath(sequence.id)
            if taxonomyPath is None:
                continue  # the sequence was not assigned

            if sequence.name in exSSDContigNameSet:
                print str('Exclude according exSSDContigNameList:' + sequence.name)
                continue  # don't consider this sequence

            for rankIdx in range(0, len(taxonomyPath.keys())):

                ncbid = int(taxonomyPath[self.taxonomicRanks[rankIdx]].ncbid)

                # if ncbid and seq.id are together from the forbidden dict then continue!!!
                if (ncbid in forbiddenDictSet) and (sequence.id in forbiddenDictSet[ncbid]):
                    ignoreSeq = True
                else:
                    ignoreSeq = False

                if taxonomyPath[self.taxonomicRanks[rankIdx]].isCopy():  # if this rank is not specified in the taxonomy
                    if ncbid not in self.ncbidToJumpedRanksNum:
                        # get down to a specified real node
                        idx = rankIdx
                        while True:
                            idx += 1
                            if not taxonomyPath[self.taxonomicRanks[idx]].isCopy():
                                break
                        # go up and count fake nodes (where there is not corresponding node in the taxonomy)
                        fn = 0
                        idx -=1
                        while idx >=0:
                            if taxonomyPath[self.taxonomicRanks[idx]].isCopy():
                                idx -= 1
                                fn += 1
                            else:
                                break
                        self.ncbidToJumpedRanksNum[ncbid] = fn
                    continue

                if not ignoreSeq:
                    if ncbid not in self.ncbidToSequences:
                        self.ncbidToSequences[ncbid] = []
                    self.ncbidToSequences[ncbid].append(sequence)  # ncbid -> list of sequences

                if ncbid not in self.ncbidToRankId:
                    self.ncbidToRankId[ncbid] = rankIdx  # one ncbid -> more ranks ids
                else:
                    assert self.ncbidToRankId[ncbid] == rankIdx

                name = taxonomyPath[self.taxonomicRanks[rankIdx]].name
                if ncbid not in self.ncbidToName:
                    self.ncbidToName[ncbid] = name
                else:
                    assert self.ncbidToName[ncbid] == name

                if rankIdx > 0:
                    idx = rankIdx
                    parentNcbid = None
                    while idx > 0:
                        if not taxonomyPath[self.taxonomicRanks[int(idx - 1)]].isCopy():
                            parentNcbid = int(taxonomyPath[self.taxonomicRanks[int(idx - 1)]].ncbid)  # parent !!!
                            break
                        idx -= 1
                else:
                    parentNcbid = None

                if ncbid not in self.ncbidToParentNcbid:
                    self.ncbidToParentNcbid[ncbid] = parentNcbid
                else:
                    assert self.ncbidToParentNcbid[ncbid] == parentNcbid

        for ncbid in self.ncbidToSequences:  # sums up all sequences that can be used to model this clade
            seqList = self.ncbidToSequences[ncbid]
            sum = 0
            for seq in seqList:
                sum += seq.seqBp
            self.ncbidToBp[ncbid] = sum

        self.toSummary(self.ncbidToSequences, summaryAllFile)  # alsoPrint=True


    def toSummary(self, ncbidList, toFile=None, alsoPrint=False):
        """
            Summary (e.g. which clades can be modeled)

            @param ncbidList: list of all ncbids that are considered in the summary
                (e.g. all clades that can be modeled)
            @param toFile: the summary will be stored to this file (if != None)
            @param alsoPrint: whether the summary should be printed to the stdout
        """
        sumEntry = []
        for ncbid in ncbidList:
            seqNum = len(self.ncbidToSequences[ncbid])
            bp = self.ncbidToBp[ncbid]
            rankId = self.ncbidToRankId[ncbid]
            entry = []
            ncbidCurrent = ncbid

            while rankId >= 0:
                name = self.ncbidToName[ncbidCurrent]
                entry.append(str(name + ' (' + str(ncbidCurrent) + ')'))
                rankId -= 1
                if ncbidCurrent in self.ncbidToJumpedRanksNum:
                    for i in range(0,self.ncbidToJumpedRanksNum[ncbidCurrent]):
                        entry.append(str('not spec. "' + str(self.taxonomicRanks[rankId]) + '"'))
                        rankId -= 1

                ncbidParent = self.ncbidToParentNcbid[ncbidCurrent]
                ncbidCurrent = ncbidParent

            entry.append(seqNum)
            entry.append(bp)
            entry.reverse()
            sumEntry.append(entry)

        sumEntry.sort(reverse=True)

        # write it to a file
        if toFile is not None:
            try:
                f = open(os.path.normpath(toFile),'w')
                f.write('# collective sequence length, number of sequences, lineage: scientific name (ncbi taxon id)\n')
                for e in sumEntry:
                    f.write(str(e) + '\n')
            except Exception:
                print "Cannot create a file or write to it:", toFile
                raise

        #print it
        if alsoPrint:
            print 'Summary-------------------------------------------------'
            for e in sumEntry:
                print e
            print '-------------------------------------------------'


    def getLeafs(self, ncbidList):
        """
            Returns a list of leaf nodes (clades) from a list of ncbids (representing internal or leaf nodes)

            @param ncbidList: list of ncbids

            @return: list of leaf ncbids
        """
        #for each ncbid: store a path from its parent to the root
        parentsSet = set()
        for ncbid in ncbidList:
            id = ncbid
            while id in self.ncbidToParentNcbid:
                id = self.ncbidToParentNcbid[id]
                parentsSet.add(id)

        #store all leaf nodes
        leafs = []
        for ncbid in ncbidList:
            if ncbid not in parentsSet:
                leafs.append(ncbid)

        return leafs


    def createPPSInputFiles(self, outFilePath, outTrainDataDir, rankIdAll, rankIdCut, rankIdCutMinBp, minPercentInLeaf,
                            maxLeafClades, minBpToModel, minGenomesWgs, minBpPerSpecies, wgsGenomesDir, forbiddenDict, dbFile, taxonomicRanks,
                            fastaLineMaxChar, minSSDfileSize, maxSSDfileSize, weightStayAll, summaryTrainFile=None):
        """
            Create the input for PhyloPythiaS: a list of clades (ncbids) that will be modeled
            and corresponding sample specific data.

            @param outFilePath: a file path where the list of the clades will be stored
            @param outTrainDataDir: a directory where the sample specific data will be stored
            @param rankIdAll: up to this rank all ncbids will be included (if there is enough data to model them)
            @param rankIdCut: the ncbids will be considered only up to this rank if there is enough data assigned to them
            @param rankIdCutMinBp: the nodes (from rankIdAll to rankIdCut) will be considered if there is at least this amount of data assigned to them
            @param minPercentInLeaf: (as rankIdCutMinBp) but furthermore this percentage of data have to be assigned to the clades
            @param maxLeafClades: this is the maximum number of leaf clades that will be considered
            @param minBpToModel: a clade can be modeled if there is at least this amount of the sample specific data assigned to it
            @param minGenomesWgs: a clade can be modeled if there is at least this number of genomes/draft genomes in the public DB
            @param minBpPerSpecies: all reference sequences must sum up to at least this number, otherwise the species won't be counted (or partly?)
            @param wgsGenomesDir: a directory where the genomes/draft genomes are stored (ncbid.1.fas/fna)
            @param forbiddenDict: a dict (ncbid -> list of seq.ids) of data that can`t be used as the sample specific data
        """

        # get all ncbids whose rank is at least rankIdAll or whose rank is at least rankIdCut and there is at least
        # rankIdCutMinBp of the sample specific data
        outNcbids = []

        # here also add all ncbids that were assigned with a very high weight, they will stay in the set if possible
        outNcbidsStaySet = set()
        for seq in self.sequences.sequences:
            weight = seq.getTaxonomyPathWeight()
            if (weight is not None) and (weight >= weightStayAll):
                taxPathDict = seq.getTaxonomyPath()
                for rankId in range(min(len(taxPathDict), (rankIdCut+1))):
                    ncbid = taxPathDict[taxonomicRanks[rankId]].ncbid
                    outNcbidsStaySet.add(ncbid)

        # print outNcbidsStaySet

        for ncbid in self.ncbidToSequences:
            rankId = self.ncbidToRankId[ncbid]
            if rankId <= rankIdAll:
                outNcbids.append(ncbid)
                outNcbidsStaySet.add(ncbid)
                continue
            if (rankId <= rankIdCut) and (self.ncbidToBp[ncbid] >= rankIdCutMinBp):  # currently doesn't have sense
                outNcbids.append(ncbid)
                continue
            if ncbid in outNcbidsStaySet:
                outNcbids.append(ncbid)

        # delete all ncbids for which there is not enough training data,
        # i.e. sample specific data or genomes/draft genomes
        temp = []
        # dbData = DBData(wgsGenomesDir, dbFile)
        refData = RefSequences(wgsGenomesDir, dbFile)
        for ncbid in outNcbids:
            if self.ncbidToBp[ncbid] >= minBpToModel:
                temp.append(ncbid)
                continue
            #if dbData.getGenomeWgsCount(ncbid, minGenomesWgs) >= minGenomesWgs:
            if refData.isRefSufficient(ncbid, minGenomesWgs, minBpPerSpecies):
                temp.append(ncbid)
                continue
            print('There is not enough data to model: %s' % ncbid)
        outNcbids = temp
        refData.close()

        # delete all ncbids that are not leafs
        outNcbids = self.getLeafs(outNcbids)

        # get leaf ncbids to which at least X% of the sample specific data (that is assigned to the leafs) is assigned
        # (except for ncbids that are in the "stay set" - outNcbidsStaySet)
        while True:
            # sums up the data in leafs
            bpInLeafs = 0
            for ncbid in outNcbids:
                if ncbid in self.ncbidToBp:
                    bpInLeafs += min(self.ncbidToBp[ncbid], maxSSDfileSize)
                else:
                    print 'no sample specific data for:', ncbid, ' ', bpInLeafs
            minBpLeaf = (bpInLeafs/100.0)*minPercentInLeaf

            # for leaf ncbids that don't contain enough data, put their parents to the candidate list
            candidateSet = set()
            breakLoop = True
            for ncbid in outNcbids:
                if self.ncbidToBp[ncbid] < minBpLeaf:
                    if ncbid not in outNcbidsStaySet:
                        candidateSet.add(self.ncbidToParentNcbid[ncbid])
                        breakLoop = False
                    else:
                        candidateSet.add(ncbid)
                else:
                    candidateSet.add(ncbid)

            # transform the candidate set to the ncbids list
            outNcbids = []
            for ncbid in candidateSet:
                if ncbid is not None:
                    outNcbids.append(ncbid)
            outNcbids = self.getLeafs(outNcbids)

            if breakLoop:
                break

        # take only the first X clades (defined by the maxLeafClades) for which we have the most sample specific data
        while len(outNcbids) > maxLeafClades:
            minBp = sys.maxint
            minBpNcbid = -1
            minIdx = -1
            idx = 0
            for ncbid in outNcbids:
                bp = self.ncbidToBp[ncbid]
                if bp < minBp and (ncbid not in outNcbidsStaySet):
                    minBp = bp
                    minBpNcbid = ncbid
                    minIdx = idx
                idx += 1
            if minBpNcbid < 0:  # this means that all ncbids in the outNcbids are also in outNcbidsStaySet
                idx = 0
                for ncbid in outNcbids:
                    bp = self.ncbidToBp[ncbid]
                    if bp < minBp:
                        minBp = bp
                        minBpNcbid = ncbid
                        minIdx = idx
                    idx += 1
            if minBpNcbid < 0:
                print 'Did not find an ncbid with min BP'
                break
            assert ncbid in self.ncbidToParentNcbid, str('cannot find parent for ncbid: ' + ncbid)

            # replace ncbid with its parent and try again
            outNcbids[minIdx] = self.ncbidToParentNcbid[ncbid] #eats up one leaf in one loop
            outNcbids = self.getLeafs(outNcbids)

        # Store the Clades to Model
        print('Nodes--------------------------------------------')
        for ncbid in outNcbids:
            print ncbid
        f = None
        try:
            f = open(os.path.normpath(outFilePath),'w')
            i = 0
            for ncbid in outNcbids:
                if i == 0:
                    f.write(str(ncbid))
                    i += 1
                else:
                    f.write('\n' + str(ncbid))
        except Exception:
            print "Cannot create a file or write to it:", outFilePath
            raise
        finally:
            if f is not None:
                f.close()

        print('Sample specific data------------------------------------')

        # Store the sample specific data
        self.storeSSD(outTrainDataDir, outNcbids, self.ncbidToSequences, fastaLineMaxChar, forbiddenDict, minSSDfileSize, maxSSDfileSize)

        # Store the summary of the clades that will be modeled
        self.toSummary(outNcbids, summaryTrainFile, alsoPrint=True)


    def storeSSD(self, SSDDir, ncbidList, ncbidToSeqDict, fastaLineMaxChar, forbiddenDict, minSSDfileSize, maxSSDfileSize):
        """
            Store the sample specific data to a directory. Creates files of form: "ncbid.1.fna"

            @param SSDDir: directory to store the sample specific data
            @param ncbidList: list of all ncbids for which the SSD should be stored
            @param ncbidToSeqDict: ncbid -> list of sequences
            @param forbiddenDict: is a dict of lists (map: ncbid -> list of sequence.id) that contain sequences that
            @param maxSSDfileSize: in bp maximum size of one file, short sequences are eliminated
            @param minSSDfileSize: in bp min size of a file, if there is less data the file won't be created
                will not be used as sample specific data for a respective ncbid (if None ~ not considered)
        """
        try:
            os.mkdir(os.path.normpath(SSDDir))
        except OSError:
            for filePath in glob.glob(os.path.join(os.path.normpath(SSDDir), r'*.1.fna')):
                os.remove(filePath)

        for ncbid in ncbidList:

            seqList = ncbidToSeqDict[ncbid]

            # filter out sequences that will not be used as the SSD
            if (forbiddenDict is not None) and (ncbid in forbiddenDict):  # this entries are already excluded in the init function !!!
                forbiddenSet = set(forbiddenDict[ncbid])
                temp = []
                for seq in seqList:
                    if seq.id not in forbiddenSet:
                        temp.append(seq)
                    else:
                        print('skip seq_id: %s for ncbid: %s' % (seq.id, ncbid))
                if len(temp) == 0:
                    continue  # no SSD will be stored for this ncbid
                else:
                    seqList = temp

            # do not create this file
            countBp = 0
            for seq in seqList:
                countBp += seq.seqBp
            if countBp < minSSDfileSize:
                continue

            # filter out short sequences
            if countBp > maxSSDfileSize:
                seqList.sort(cmp=seqWeightThenLenCmp, reverse=True)
                seqBpCount = 0
                tmp = []
                # this filters out only short sequences
                for seq in seqList:
                    if (seqBpCount + seq.seqBp) > maxSSDfileSize:
                        continue  # was break
                    else:
                        seqBpCount += seq.seqBp
                        tmp.append(seq)

                seqList = tmp

                #this removes sequences shorter than 1000 and then by random until I get less sequences than maxSSDfileSize
                #while countBp > maxSSDfileSize and len(seqList) > 1:
                #    if seqList[len(seqList) - 1].seqBp <= 1000:
                #        seqRem = seqList.pop(len(seqList) - 1)
                #    else:
                #        seqRem = seqList.pop(random.randint(0, len(seqList) - 1))
                #    countBp -= seqRem.seqBp


            filePath = os.path.join(os.path.normpath(SSDDir), str(str(ncbid) + '.1.fna'))
            range = fastaLineMaxChar
            f = None
            try:
                f = open(os.path.normpath(filePath),'w')
                k = 0
                counter = 0
                for seq in seqList:
                    #if seqCount != None and counter == seqCount:
                    #    break
                    #else:
                    #    counter += 1

                    entry = str('>' + str(seq.scaffold.id) + '_' + str(seq.id))
                    if k == 0:
                        f.write(entry)
                        k += 1
                    else:
                        f.write('\n' + entry)

                    s = seq.getSeq()
                    l = seq.seqBp
                    i = 0
                    while i < l:
                        j = i + range
                        f.write('\n' + s[i:j])
                        i += range
            except Exception as e:
                print "Cannot create a file or write to it:", filePath
                raise e
            finally:
                if f is not None:
                    f.close()


def storeDictToAFile(filePath, dictOfLists):
    """
        Stores a dictionary of lists (i.e. map: key -> list of items) to a file in format: (key tab item)

        @param filePath: a file in which the dictionary will be stored in format: (key tab item)
        @param dictOfLists: dict to be stored that represents mapping: (key -> list of items)
    """
    f = None
    try:
        f = open(os.path.normpath(filePath), 'w')
        k = 0
        for key in dictOfLists:
            list = dictOfLists[key]
            for item in list:
                if k == 0:
                    f.write(str(key) + '\t' + str(item))
                    k = 1
                else:
                    f.write('\n' + str(key) + '\t' + str(item))
    except Exception as e:
        print "Cannot create a file or write to it:", filePath
        raise e
    finally:
        if f is not None:
            f.close()


def loadDictFromAFile(filePath):
    """
        Returns a dictionary that is stored in a file.

        @param filePath: a file in which a dictionary is stored in format: (key tab item)

        @return: dict that represents mapping: (key -> list of items)
    """
    try:
        dictOfLists = dict([])
        f = open(os.path.normpath(filePath), 'r')

        for line in f:
            pair = re.findall('[^\t]+', common.noNewLine(line))
            assert len(pair) == 2, str('There are not two values separated by \t at line: ' + line)
            key = int(pair[0])
            val = int(pair[1])
            if key in dictOfLists:
                dictOfLists[key].append(val)
            else:
                list = []
                list.append(val)
                dictOfLists[key] = list

        return dictOfLists
    except Exception:
        print "Cannot create a file or write to it:", filePath
        raise
    finally:
        f.close()


def updateForbiddenList(forbiddenList, filePath):
    """
        Appends (key tab item) pairs that are given as a list to a file. If the file doesn`t exists it creates it.
        Values that are already contained in the list are not appended. The updated file then represents a dict of lists
        (i.e. map: key -> list of items)

        @param forbiddenList: a list of several fields where the first entry is in format ([0-9]+_[0-9]+)
            and the second field is a number. The second field represent keys. The second number in the first field
            represent corresponding items
        @param filePath: a file that will be updated
    """
    #if file exists then read values from it
    if os.path.isfile(os.path.normpath(filePath)):
        dictOfLists = loadDictFromAFile(filePath)
    else:
        dictOfLists = dict([])

    #add values from the forbidden list to the dictionary
    for entry in forbiddenList: #entry format: ['114_290', 186803, 'H', True, 4, u'Lachnospiraceae', u'Bacteria', 2]
        ncbid = int(entry[1])
        try:
            seqId = int(re.sub('[0-9]+_([0-9]+)',r'\1', entry[0]))
        except Exception:
            #sys.stderr.write('updateForbiddenList: entry has wrong format and we skip it: ' + str(entry[0]) + '\n')#!!!
            #e.g. sequences from the expert sample specific data
            continue
        if ncbid in dictOfLists:
            dictOfLists[ncbid].append(seqId)
        else:
            list = []
            list.append(seqId)
            dictOfLists[ncbid] = list

    #remove possible duplicates in the individual lists
    for ncbid in dictOfLists:
        list = dictOfLists[ncbid]
        set = set(list)
        if len(list) != len(set):
            print str('There are duplicates in the forbidden list: ' + str(ncbid) + ' ' + str(len(list)) + ' ' + str(len(set)))
            temp = []
            set = set()
            for i in list:
                if i not in set:
                    temp.append(i)
                    set.add(i)
            dictOfLists[ncbid] = temp

    #store the values to the file
    storeDictToAFile(filePath, dictOfLists)


def test():
    """
        @deprecated: old test
    """
    pass

    #config = Config(open(os.path.normpath('D:\\A_Phylo\\A_Metagenomic\\pPPS\\workspace\\pPPS\\config01.cfg')), 'pPPS')

    #read sequences
    #sequences = Sequences(config) #!!!

    #write ids file
    #sequences.writeSequences(config.get('inputIdsFastaFile'))

    #taxonomy
    #taxonomy = Taxonomy(config.get('databaseFile'), config.get('taxonomicRanks').split(','))

    #placeSequences(sequences, taxonomy, config.get('final_RAxML_outputs'))

    #if eval(config.get('placeContigsFromTheSameScaffold')):
    #    sequences.placeContigsFromTheSameScaffold(taxonomy)

    #taxonomy.close()

    #pps = PPSInput(sequences, config.get('taxonomicRanks').split(','), config.get('summaryAllFile'))

    #pps.createPPSInputFiles(config.get('nodesFile'), config.get('trainingDataDir'),
    #                    int(config.get('rankIdAll')), int(config.get('rankIdCut')), int(config.get('rankIdCutMinBp')),
    #                    int(config.get('minTrainingBp')),
    #                    int(config.get('fastaLineMaxChar')), config.get('summaryTrainFile'))


#(self, outFilePath, outTrainDataDir, rankIdAll, rankIdCut, rankIdCutMinBp, fastaLineMaxChar):

    #print placements for each sequence that is assigned/placed
    #print '--------------------------------'
    #for s in sequences.sequences:
    #    dict = sequences.getTaxonomyPath(s.id)
    #    if dict == None:
    #        continue
    #    print s.name, '...'
    #    for rank in Config.taxonomicRanks:
    #        if rank not in dict:
    #            break
    #        print rank, dict[rank].ncbid, dict[rank].name
    #    print ''




if __name__ == "__main__":
  #test()
  pass