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

import sys
import os
import re
import subprocess
import glob

from algbioi.com.config import Config
from algbioi.com.csv import forEachLine
from algbioi.com.csv import OutFileBuffer
from algbioi.com.csv import getMapping
from algbioi.com.fasta import fastaFileToDict


class MGCluster():
    """ Main class """

    def __init__(self, config, mgWorkingDir, s16Prefix, sequences, taxonomy, sampleName):
        self._config = config
        self._mgWorkingDir = mgWorkingDir
        self._taxonomicRanks = config.get('taxonomicRanks').split(',')
        self._minBpToModel = int(config.get('minBpToModel'))
        self._clustDir = os.path.normpath(os.path.join(mgWorkingDir,'clust'))
        self._s16Prefix = s16Prefix
        self._sequences = sequences
        self._mgToDM = None
        self._mgToCluster = None
        self._taxonomy = taxonomy
        self._sampleName = sampleName
        self._mgList = None
        self._mgToMaxThreshold = None
        self._seqIdToTaxPathDict = None
        self._seqIdToWeight = None
        self._initDone = False

        if not os.path.exists(self._clustDir):
            try:
                os.mkdir(self._clustDir)
            except OSError:
                print 'Can`t create directory:', self._clust
                raise

#---------------------------------------------------------------------------------

    def _init(self, align=True, dm=True, cluster=True):
        """
            Init data, compute: alignment, distance matrix, clusters.
        """
        if self._initDone:
            return
        self._initDone = True

        fastaPathList = [] # fasta files containing regions that correspond to particular marker genes
        self._mgList = [] # list of names of marker genes
        mgToFastaPath = dict([]) # marker gene name -> fasta file path

        #collect regions from Amphora mg
        for fastaFile in glob.glob(os.path.join(os.path.normpath(self._mgWorkingDir),'*.gff')):
            fastaPathList.append(fastaFile)
        for path in fastaPathList:
            name = re.sub('([^\.]+)\..*$', r'\1' , os.path.basename(path))
            mg = re.sub(r'([^_]+)_dna', r'\1',name)
            dir = os.path.dirname(path)
            self._mgList.append(mg)
            mgToFastaPath[mg] = path

        #add 16S
        s16List = ['5S_rRNA', '16S_rRNA', '23S_rRNA']
        for mg in s16List:
            mgToFastaPath[mg] = str(self._s16Prefix + '.' + mg + '.fna')
            self._mgList.append(mg)

        #For each marker gene create filtered fasta file that contains for each mg and sequence at most one region.
        mgToFilteredFastaPath = dict([])
        mgToSeqNameToTaxPathDict = dict([]) #mg -> seqName (~region name) -> pred
        for mg in self._mgList:
            mgToSeqNameToTaxPathDict[mg] = dict([])

        for seq in self._sequences.sequences:
            id = str(str(seq.scaffold.id) + '_' + str(seq.id))
            for mg,tag,pred in zip(seq.getCandidateTaxPathSourceList(), seq.getCandidateTaxPathTagList(),
                                    seq.getCandidateTaxPathDictList()):
                mgToSeqNameToTaxPathDict[mg][str(id + '_' + tag)] = pred

        #for each marker gene: choose only one sequence region for each mg and sequence
        #all sequences are predicted at least at superkingdom
        for mg in self._mgList:
            seqNameToPred = mgToSeqNameToTaxPathDict[mg] #sequence region predictions for this mg
            seqNameToSeq = fastaFileToDict(mgToFastaPath[mg]) #read the fasta file
            outPath = os.path.normpath(os.path.join(self._clustDir, str(mg + '.filter.fna')))
            mgToFilteredFastaPath[mg] = outPath
            out = OutFileBuffer(outPath)
            seqBaseToSeqName = dict([]) # sequence base (scaffId_seqId) -> region name
            for seqName in seqNameToSeq:
                seqBase = re.sub(r'^([0-9]+_[0-9]+)[^0-9].*',r'\1', seqName)
                if seqBase not in seqBaseToSeqName:
                    seqBaseToSeqName[seqBase] = []
                seqBaseToSeqName[seqBase].append(seqName)
            for seqBase in seqBaseToSeqName:
                seqId = int(re.sub(r'^[0-9]+_([0-9]+)',r'\1', seqBase))
                seqBaseTaxPathDict = self._sequences.getSequence(seqId).getTaxonomyPath()
                list = seqBaseToSeqName[seqBase]
                candidateSeq = [] # sequence region is predicted at least at rank superkingdom
                for seqName in list:
                    if seqName not in seqNameToPred:
                        taxPathDict = None
                    else:
                        taxPathDict = seqNameToPred[seqName]
                    if taxPathDict != None:
                         candidateSeq.append(seqName)
                if len(candidateSeq) == 0:
                    continue
                candidateSeq2 = [] # sequence regions predicted at least at the same rank as the whole sequence
                for seqName in candidateSeq:
                    taxPathDict = seqNameToPred[seqName]
                    if ((seqBaseTaxPathDict == None)
                        or (len(taxPathDict) >= len(seqBaseTaxPathDict))): #predict at least at the same level
                        candidateSeq2.append(seqName)
                if len(candidateSeq2) > 0: #take the longest sequence
                    sMax = candidateSeq2[0]
                    for s in candidateSeq2[1:]:
                        if len(seqNameToSeq[s]) > len(seqNameToSeq[sMax]):
                            sMax = s
                else: #all sequence regions are predicted higher than the sequence
                    sMax = candidateSeq[0] #sequence region with the most specific prediction
                    for s in candidateSeq[1:]:
                        taxPathDictMax = seqNameToPred[sMax]
                        taxPathDictS = seqNameToPred[s]
                        if taxPathDictS == None:
                            continue
                        if taxPathDictMax == None:
                            sMax = s
                            continue
                        if len(taxPathDictMax) < len(taxPathDictS):
                            sMax = s

                    candidateSeq3 = [] #get all sequence regions with the most specific prediction
                    taxPathDictMax = seqNameToPred[sMax]
                    for s in candidateSeq:
                        taxPathDictS = seqNameToPred[s]
                        if taxPathDictMax == None:
                            candidateSeq3.append(s)
                        elif len(taxPathDictS) == len(taxPathDictMax):
                            candidateSeq3.append(s)
                    sMax = candidateSeq3[0]
                    for s in candidateSeq3[1:]: #take the longest sequence
                        if len(seqNameToSeq[sMax]) < len(seqNameToSeq[s]):
                            sMax = s

                out.writeText(str('>' + str(sMax) + '\n' + str(seqNameToSeq[sMax]) + '\n'))

            out.close()

        mgToAlignPath = dict([])
        for mg in self._mgList:
            mgToAlignPath[mg] = os.path.normpath(os.path.join(self._clustDir, str(mg + '.align.fna')))

        #build alignment
        if align:
            for mg in self._mgList:
                alignCmd = str(self._config.get('aligner') + ' -in ' + mgToFilteredFastaPath[mg]
                + ' -out ' + mgToAlignPath[mg] + ' -quiet')
                assert os.name == 'posix'
                predictProc = subprocess.Popen(alignCmd, cwd=self._mgWorkingDir, shell=True, bufsize=-1) #stdout=subprocess.STDOUT, stderr=subprocess.STDOUT)
                predictProc.wait()
                print 'Muscle return code for', mg, ':', predictProc.returncode
                if predictProc.returncode != 0:
                    sys.stderr.write(str(alignCmd + ' \n'))

        #compute DM
        if dm:
            for mg in self._mgList:
                mothur = os.path.join(os.path.normpath(self._configRRNA16S.get('mothurInstallDir')), 'mothur')
                mothurCmd = str('time ' + mothur + ' "#dist.seqs(fasta=' + mgToAlignPath[mg]
                                + ', processors=2, countends=F, calc=nogaps, cutoff=0.3, output=lt)"')
                assert os.name == 'posix'
                mothurProc = subprocess.Popen(mothurCmd, shell=True, bufsize=-1, cwd=self._mgWorkingDir)
                mothurProc.wait()
                print 'Mothur return code dist:', mg, mothurProc.returncode
                #distFilePath = os.path.join(os.path.dirname(mgToAlignPath[mg]), str(mg + '.align.phylip.dist'))
                #self._mgToDM[mg] = forEachLine(distFilePath, DM())
                #self._mgToDM[mg].printDM()

        #cluster
        if cluster:
            for mg in self._mgList:
                distFilePath = os.path.join(os.path.dirname(mgToAlignPath[mg]), str(mg + '.align.phylip.dist'))
                mothur = os.path.join(os.path.normpath(self._configRRNA16S.get('mothurInstallDir')), 'mothur')
                mothurCmd = str('time ' + mothur + ' "#cluster(phylip=' + distFilePath
                                + ', method=furthest, hard=t, precision=1000)"')
                assert os.name == 'posix'
                mothurProc = subprocess.Popen(mothurCmd, shell=True, bufsize=-1, cwd=self._mgWorkingDir)
                mothurProc.wait()
                print 'Mothur return code cluster:', mg, mothurProc.returncode

        #read DM and clusters

        #sequence predictions
        self._seqIdToTaxPathDict = dict([])
        self._seqIdToWeight = dict([])
        for seq in self._sequences.sequences:
            id = int(seq.id)
            self._seqIdToTaxPathDict[id] = seq.getTaxonomyPath()
            self._seqIdToWeight[id] = seq.getTaxonomyPathWeight()

        #similarity thresholds
        thresholds = self._configMG.get('mgSimilarityThresholds')
        self._mgToMaxThreshold = dict([])
        tmpDict = getMapping(self._configMG.get('mgSimilarityThresholds'), 0, 1, sep='\t', comment = '#')
        for k in tmpDict:
            self._mgToMaxThreshold[k] = float(tmpDict[k][0])

        self._mgToDM = dict([])
        self._mgToCluster = dict([])
        for mg in self._mgList:
            file = os.path.join(os.path.dirname(mgToAlignPath[mg]), str(mg + '.align.phylip.dist'))
            self._mgToDM[mg] = forEachLine(file, DM())
            file = os.path.join(os.path.dirname(mgToAlignPath[mg]), str(mg + '.align.phylip.fn.list'))
            self._mgToCluster[mg] = forEachLine(file, MCluster(self._seqIdToTaxPathDict, self._mgToMaxThreshold[mg]))



#---------------------------------------------------------------------------------

    def refineSpecificPred(self):
        self._init(align=False, dm=False, cluster=False)
        seqToCandidatePred = dict([])
        for mg in self._mgList:
            mCluster = self._mgToCluster[mg]
            tCluster = mCluster.getLastNoConflictClustering()
            threshold = tCluster.getThreshold()
            if threshold > self._mgToMaxThreshold[mg]: #just in the case the first clustering was already conflicting
                continue
            #seqNameToGroupId = tCluster.getSeqNameToGroupId()
            groupIdToSeqNameSet = tCluster.getGroupIdToSeqNameSet()
            for groupId in groupIdToSeqNameSet:
                group = groupIdToSeqNameSet[groupId]
                if len(group) < 2:
                    continue
                seqNameMaxTaxPathDict = None #the lowest prediction within the group (all lie on the common path to the root)
                seqNameToTaxPathDict = dict([])
                #weightList = []
                for seqName in group:
                    seqId = int(re.sub(r'^[0-9]+_([0-9]+)$', r'\1', seqName))
                    taxPathDict = self._seqIdToTaxPathDict[seqId]
                    #weightList.append(self._seqIdToWeight[seqId])
                    seqNameToTaxPathDict[seqName] = taxPathDict
                    if (seqNameMaxTaxPathDict == None) or (len(seqNameToTaxPathDict[seqNameMaxTaxPathDict]) < len(taxPathDict)):
                        seqNameMaxTaxPathDict = seqName
                maxTaxPathDict = seqNameToTaxPathDict[seqNameMaxTaxPathDict]

                weightList = []
                for seqName in group:
                    if len(seqNameToTaxPathDict[seqName]) >= len(maxTaxPathDict):
                        weightList.append(self._seqIdToWeight[int(re.sub(r'^[0-9]+_([0-9]+)$', r'\1', seqName))])

                for seqName in group:
                    if seqName == seqNameMaxTaxPathDict:
                        continue
                    if len(seqNameToTaxPathDict[seqName]) < len(maxTaxPathDict):
                        if seqName not in seqToCandidatePred:
                            seqToCandidatePred[seqName] = []
                        seqToCandidatePred[seqName].append((mg, maxTaxPathDict, min(weightList)))

        #resolve candidate predictions
        for seqName in seqToCandidatePred:
            list = seqToCandidatePred[seqName]
            if len(list) == 1:
                taxPathDict = self._taxonomy.replicateTaxPathDict(list[0][1])
                weight = list[0][2]
            else:
                #get lowest common ancestor
                taxPathDictList = []
                weightList = []
                for t in list:
                    taxPathDictList.append(t[1])
                    weightList.append(t[2])
                taxPathDict = self._taxonomy.getLongestCommonPathFromMultipleAssignments(taxPathDictList)
                weight = min(weightList)
            seqId = int(re.sub(r'^[0-9]+_([0-9]+)$', r'\1', seqName))
            scaffoldId = int(re.sub(r'^([0-9]+)_[0-9]+$', r'\1', seqName))
            currentTaxPathDict = self._sequences.getSequence(seqId).getTaxonomyPath()
            if (taxPathDict != None) and (len(currentTaxPathDict) < len(taxPathDict)):
                print 'Spec. pred override:', seqId, currentTaxPathDict, '->', taxPathDict
                self._sequences.setTaxonomyPathOverride(seqId, scaffoldId, taxPathDict, weight)


#----------------------------------------------------------------

    def reconstructOTU(self, mgToConsider = ['16S_rRNA', '23S_rRNA', 'rpoB']):
        self._init(align=False, dm=False, cluster=False)
        #collect all ncbids
        ncbids = set([])
        for seq in self._sequences.sequences:
            taxPathDict = seq.getTaxonomyPath()
            if taxPathDict != None:
                for t in taxPathDict:
                    ncbids.add(taxPathDict[t].ncbid)

        #remove elements that correspond to the superkingdom
        ncbids.discard(1)     #Root
        ncbids.discard(2)     #Bacteria
        ncbids.discard(2157)  #Archaea
        ncbids.discard(2759)  #Eukaryota

        #for each ncbid collect a list of sequences their lowest assignment is to this ncbid
        ncbidToSeqList = dict([])
        innerNcbidSet = set([]) #set of ncbids that are not leafs
        speciesNcbidSet = set([])
        for seq in self._sequences.sequences:
            taxPathDict = seq.getTaxonomyPath()
            if taxPathDict == None:
                continue
            ncbid = taxPathDict[self._taxonomicRanks[len(taxPathDict) - 1]]

            for rankIdx in range(len(taxPathDict)-1):
                tmp = taxPathDict[self._taxonomicRanks[rankIdx]].ncbid
                innerNcbidSet.add(tmp)
                if self._taxonomicRanks[rankIdx] == 'species':
                    speciesNcbidSet.append(tmp)

            if ncbid not in ncbidToSeqList:
                ncbidToSeqList[ncbid] = []
            ncbidToSeqList[ncbid].append(str(str(seq.scaffold.id) + '_' + str(seq.id)))

        #try to resolve subclades of all ncbids
        for ncbid in ncbids:

            #skip ncbids that are at the rank species !!!
            if ncbid in speciesNcbidSet:
                continue

            if ncbid not in ncbidToSeqList:
                continue

            seqNameList = ncbidToSeqList[ncbid]
            for mg in mgToConsider:
                tCluster = self._mgToCluster[mg].getLastNoConflictClustering()
                seqNameToGroupId = tCluster.getSeqNameToGroupId()
                groupIdToSeqNameSet = tCluster.getGroupIdToSeqNameSet()

                groupIdToSeqNameMgList = dict([])

                for seqName in seqNameList:
                    if seqName in seqNameToGroupId:
                        groupId = seqNameToGroupId[seqName]
                        if groupId not in groupIdToSeqNameMgList:
                            groupIdToSeqNameMgList[groupId] = []
                        groupIdToSeqNameMgList[groupId].append(seqName)

                groupIdToBp = dict([]) #for each group store its size
                for groupId in groupIdToSeqNameMgList:
                    seqNameList = groupIdToSeqNameMgList[groupId]
                    bp = 0
                    for seqName in seqNameList:
                        seqId = int(re.sub(r'^[0-9]+_([0-9]+)$', r'\1', seqName))
                        bp += self._sequences.getSequence(seqId).seqBp
                    groupIdToBp[groupId] = bp

                candidateGroupIdOtuSet = set([])
                for groupId in groupIdToBp:
                    if groupIdToBp[groupId] >= self._minBpToModel:
                        candidateGroupIdOtuSet.add(groupId)

                if ((ncbid in innerNcbidSet) and (len(candidateGroupIdOtuSet) >= 1) or
                    (ncbid not in innerNcbidSet) and (len(candidateGroupIdOtuSet) >= 2)):
                    #create new OTUs
                    newNcbid = self._taxonomy.createNewOtuDBEntry(ncbid, self._sampleName, rank='species')
                    taxPathDict = self._taxonomy.getPathToRoot(newNcbid)
                    for groupId in candidateGroupIdOtuSet:
                        seqNameList = groupIdToSeqNameMgList[groupId]
                        weightList = []
                        for seqName in seqNameList:
                            seqId = int(re.sub(r'^[0-9]+_([0-9]+)$', r'\1', seqName))
                            weightList.append(self._sequences.getSequence(seqId).getTaxonomyPathWeight())

                        for seqName in seqNameList:
                            seqId = int(re.sub(r'^[0-9]+_([0-9]+)$', r'\1', seqName))
                            scaffoldId = int(re.sub(r'^([0-9]+)_[0-9]+$', r'\1', seqName))
                            weight = None
                            print 'NewOtu:', ncbid, seqId, self._sequences.getSequence(seqId).getTaxonomyPath(), '->', taxPathDict
                            self._sequences.setTaxonomyPathOverride(seqId, scaffoldId,
                                                                    self._taxonomy.replicateTaxPathDict(taxPathDict), min(weightList))

                    break # don`t try to infer OTUs for



                #get all clusters at this node
                #compute the size of all clusters at this node
                #try to extend the clusters using other marker genes
                #store suggested OTUs





            #get sequences that are predicted to this ncbid but not lower



#---------------------------------------------------------------------------------


    def reconstructOTU_OLD(self, lowestRank='genus'):
        self._init(align=False, dm=False, cluster=False)

        allowedRanks=['root','superkingdom','phylum','class','order','family','genus','species']
        #for each clade, write which OTUs could be reconstructed
        #collect clades
        ncbids = set([])
        for seq in self._sequences.sequences:
            taxPathDict = seq.getCandidateTaxPathDictList()
            for t in taxPathDict:
                print taxPathDict, t
                ncbids.add(taxPathDict[t].ncbid)
        seqNodeListD = dict([])
        seqLowerListD = dict([])
        seqUpperListD = dict([])
        for ncbid in ncbids:
            #if ncbid != 171549:
            #    continue
            seqNodeListD[ncbid] =  set([])
            seqLowerListD[ncbid] = set([])
            seqUpperListD[ncbid] = set([])

            for seq in self._sequences.sequences:
                path = seq.getTaxonomyPath()
                if path == None or len(path) <= 1:
                    continue
                found = False
                for rank in allowedRanks:
                    if path[rank].node == ncbid:
                        found=True
                        seqNodeListD[ncbid].add(seq.id)
                    elif not found:
                        seqUpperListD[ncbid].add(seq.id)
                    else:
                        seqLowerListD[ncbid].add(seq.id)

        for ncbid in ncbids:
            if ncbid != 171549: #!remove then!
                continue
            s16Clust = self._mgToCluster['16S_rRNA']

            lastThreshold = 0.0
            for threshold in s16Clust.thresholdsList:
                if threshold < 0.3:
                    lastThreshold = threshold
                    continue
                if threshold > 0.5:
                    lastThreshold = threshold

                cluster = self.thresholdToTCluster[threshold]
                #look at each cluster
                relClustIdSet = set([])
                for i in range(cluster.clusterIdCount):
                    #is some sequence from the cluster at my node? filter out clusters
                    for seqId in seqNodeListD[ncbid]:
                        if seqId in cluster.clusterIdToSeqSet[i]:
                            relClustIdSet.add(i)
                            break

                #inspect relevant clusters "relClustIdSet"
                filteredClustId = set([])
                for i in relClustIdSet:
                    seqs = cluster.clusterIdToSeqSet[i]
                    relSeq = 0
                    wrongSeq = 0
                    for s in seqs:
                        if s in seqNodeListD[ncbid] or s in seqUpperListD[ncbid] or s in seqLowerListD[ncbid]:
                            relSeq += 1
                        else:
                            wrongSeq += 1
                    if wrongSeq == 0: # corrected wrongSet
                        filteredClustId.add(i)
                #have clusters to consider:
                bp = 0
                for i in filteredClustId:
                    s = cluster.clusterIdToSeqSet[i]
                    sAtNode = set([])
                    for seq in s:
                        if seq in seqNodeListD[ncbid]:
                            sAtNode.add(seq) # take just sequences that were at the node
                            bp += self._sequences.getSequence(seq).seqBp

                    print ncbid, 'clust', sAtNode, bp
                    if bp > 100000:
                        #create new OTU and assign all sequences in the cluster to it
                        newNcbid = self._taxonomy.createNewOtuDBEntry(ncbid, self._sampleName, rank='species')
                        for seq in sAtNode:
                            taxPathDictOTU = self._taxonomy.getPathToRoot(newNcbid)
                            self._sequences.setTaxonomyPathOverride(seq.id, seq.scaffold.id, taxPathDictOTU, 100.0)

                break

#---------------------------------------------------------------------------------

class DM():
    """
        Phylip Distance matrix for a marker gene (and line parser).
    """
    def __init__(self):
        self._counter = -1
        self._seqNum = None
        self._seqNames = []
        self._matrix = []
        self._nameToIndex = dict([])

    def parse(self, line):
        if self._counter == -1:
            self._seqNum = int(line)
        else:
            tokens = line.split()
            name = re.sub(r'^([0-9]+_[0-9]+)_.*',r'\1', tokens[0])
            self._seqNames.append(name)
            self._nameToIndex[name] = self._counter
            list = []
            for e in tokens[1:]:
                list.append(float(e))
            self._matrix.append(list)
        self._counter += 1

    def printDM(self):
        for name1 in self._seqNames:
            for name2 in self._seqNames:
                if name1 != name2:
                    print self.getDist(name1,name2)

    def getDist(self, seqName1, seqName2):
        if seqName1 not in self._nameToIndex or seqName2 not in self._nameToIndex:
            return None
        idx1 = self._nameToIndex[seqName1]
        idx2 = self._nameToIndex[seqName2]
        #print idx1, idx2
        return self._matrix[max(idx1,idx2)][min(idx1,idx2)]

    def getSeqNameList(self):
        return self._seqNames


#---------------------------------------------------------------------------------

class MCluster():
    """
        Clusters for different thresholds
    """
    def __init__(self, seqIdToTaxPathDict, maxSimilarityThreshold):
        self._thresholdsList = []
        self._thresholdIdxToTCluster = dict([])
        self._lineCounter = -1
        self._seqIdToTaxPathDict = seqIdToTaxPathDict
        self._lastNoConflictThresholdIdx = 0
        self._maxSimilarityThreshold = maxSimilarityThreshold

    def parse(self, line):
        if self._lineCounter == -1: #skip unique clusters
            pass
        else:
            c = TCluster(line)
            self._thresholdsList.append(c.getThreshold())
            self._thresholdIdxToTCluster[self._lineCounter] = c
        self._lineCounter += 1

    def finalize(self):
        #finds largest threshold at which all groups within a clustering are consistent
        for i in range(len(self._thresholdsList)):
            if self._thresholdsList[i] > self._maxSimilarityThreshold:
                break
            tCluster = self._thresholdIdxToTCluster[i]
            if self._isConsistent(tCluster):
                self._lastNoConflictThresholdIdx = i
            else:
                break

    def getLastNoConflictClustering(self):
        return self._thresholdIdxToTCluster[self._lastNoConflictThresholdIdx]

    def _isConsistent(self, tCluster):
        groupIdToSeqNameSet = tCluster.getGroupIdToSeqNameSet()
        for groupId in groupIdToSeqNameSet:
            taxPathDictList = []
            maxLenTaxPathDict = None

            for seqName in groupIdToSeqNameSet[groupId]:
                t = self._seqIdToTaxPathDict[int(re.sub(r'^[0-9]+_([0-9]+)$', r'\1', seqName))]
                taxPathDictList.append(t)
                if (maxLenTaxPathDict == None) or (len(maxLenTaxPathDict) < len(t)):
                    maxLenTaxPathDict = t

            assert len(taxPathDictList) > 0
            if len(taxPathDictList) == 1:
                continue

            allowedNcbids = set([])
            for node in maxLenTaxPathDict:
                allowedNcbids.add(node.ncbid)

            for t in taxPathDictList:
                for node in t:
                    if node.ncbid not in allowedNcbids:
                        return False

        return True
        #for each mg, get the highest threshold at which there is no conflicting cluster with known labels


        #            for mg in mgList:
        #                mCluster = self._mgToCluster[mg]
        #                mCluster.setNoConflictThreshold()

    def getThresholdsList(self):
        """
            List of thresholds of different clusterings.
        """
        return self._thresholdsList


    def getClusterAtThreshold(self, thresholdIdx):
        """
            Gets cluster where the thresholdIdx correspond to the index in the array returned from getThresholdsList.
        """
        return self._thresholdIdxToTCluster[thresholdIdx]

#---------------------------------------------------------------------------------

class TCluster():
    """
        One clustering of marker genes at a specific threshold.
    """
    def __init__(self, line):
        tokens = line.split(',')
        self._threshold = float(re.sub(r'^([^\t]+)\t[^\t]+\t.*', r'\1', tokens[0]))
        tokens[0] = re.sub(r'^[^\t]+\t[^\t]+\t(.*)', r'\1', tokens[0])
        self.groupIdCount = 0
        self.seqNameToGroupId = dict([])
        self.groupIdToSeqNameSet = dict([])
        for token in tokens:
            names = token.split('\t')
            self.groupIdToSeqNameSet[self.groupIdCount] = set([])
            for name in names:
                #print name
                if re.match(r'^[0-9]+_.*$', name):
                    seqName = re.sub(r'^([0-9]+_[0-9]+)_.*$',r'\1', name)
                    self.seqNameToGroupId[seqName] = self.groupIdCount
                    self.groupIdToSeqNameSet[self.groupIdCount].add(seqName)
            self.groupIdCount += 1

    def getThreshold(self):
        return self._threshold

    def getSeqNameToGroupId(self):
        return self.seqNameToGroupId

    def getGroupIdToSeqNameSet(self):
        return self.groupIdToSeqNameSet

#---------------------------------------------------------------------------------

def test():
    config = Config(open('/Users/ivan/Documents/work/binning/tests/CowRumen/03/config.cfg'), 'pPPS')
    mgWorkingDir = '/Users/ivan/Documents/work/binning/tests/CowRumen/03/working/mgWorking'
    s16Prefix = '/Users/ivan/Documents/work/binning/tests/CowRumen/03/working/cow_rumen_fragmented_velvet_assembly_scaffolds.fas.ids'
    clust = MGCluster(config, mgWorkingDir, s16Prefix)
    clust.preprocess(align=False, dm=False, cluster=False, readData=True)
    #clust.buildSpecificPred()

    clust.reconstructOTU()


if __name__ == "__main__":
    test()