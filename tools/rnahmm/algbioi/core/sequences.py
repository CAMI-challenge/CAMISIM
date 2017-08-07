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

from algbioi.core.sequence import Sequence
from algbioi.core.scaffold import Scaffold
from algbioi.com.config import Config
from algbioi.com.common import noNewLine


class Sequences():

    def __init__(self, inputFastaFile, inputFastaScaffoldsFile, scaffoldsToContigsMapFile, taxonomicRanks, minSeqLen):
        """
            Represents a set of contigs/sequences and scaffolds.

            @param config: configuration file.
            @param inputFastaFile: fasta file with contigs/sequences
            @param inputFastaScaffoldsFile: fasta file with scaffolds
            @param scaffoldsToContigsMapFile: mapping file
            @raise:
        """
        self.sequences = []  # list of contigs/sequences
        self.scaffolds = []  # list of corresponding scaffolds
        self.scaffoldsDict = {}  # map: sequence.id -> sequence
        self.sequencesDict = {}  # map: sequence.id -> sequence
        self.placedSeqSet = set()  # set of sequence.id to which a taxonomy path was assigned

        self.scaffNameToScaff = {}

        self._scaffoldPattern = '^(.*)'  # config.get('scaffoldPattern')
        self._taxonomicRanks = taxonomicRanks  # config.get('taxonomicRanks').split(',')
        self._fastaLineMaxChar = 80  # int(config.get('fastaLineMaxChar'))
        self._minSeqLen = minSeqLen  # int(config.get('minSeqLen'))
        self._minScaffLen = self._minSeqLen  # int(config.get('minScaffLen'))
        self._inputFastaFile = inputFastaFile
        self._inputFastaScaffoldsFile = inputFastaScaffoldsFile
        self._scaffoldsToContigsMapFile = scaffoldsToContigsMapFile

        #reads a list of sequences/contigs
        self._seqCounter = -1
        self._scaffCounter = -1

        #read sequences/contigs
        self._readContigsScaffolds(self._inputFastaFile, True) # True ~ read contigs

        #read scaffolds if available
        if self._inputFastaScaffoldsFile is not None:
            self._readContigsScaffolds(self._inputFastaScaffoldsFile, False) #False ~ read scaffolds
            scaffRead = True
        else:
            scaffRead = False

        #read scaffold-contigs mapping if available (scaff_name tab contig_name)
        #@return: (map: scaff_name -> list of contig names)
        if self._scaffoldsToContigsMapFile is not None:
            scaffToContigListDict = toScafContigMap(self._scaffoldsToContigsMapFile)
            scaffContigMapRead = True
        else:
            scaffContigMapRead = False

        #bind contigs with the corresponding scaffolds
        if (not scaffRead) and (not scaffContigMapRead):

            #assigns scaffolds
            scaffoldDict = dict()
            for s in self.sequences:
                sn = re.findall(self._scaffoldPattern, s.name);
                if len(sn) != 1:
                    raise Exception('The pattern to match scaffold name is ambiguous or wrong, variants:', sn)
                scaffoldName = sn[0]
                if scaffoldName in scaffoldDict: #sequence belongs to an existing scaffold
                    scaffold = scaffoldDict[scaffoldName]
                    s.setScaffold(scaffold)
                    scaffold.addNextContig(s)
                else:
                    #create a new scaffold
                    scaffold = self._addScaff(scaffoldName, s, None)
                    scaffoldDict[scaffoldName] = scaffold
                    s.setScaffold(scaffold)

        elif scaffContigMapRead:

            #assign contigs to the scaffolds (according: scaffToContigListDict)
            #the scaffolds either already exists or not

            #get the map: contig -> scaffold
            contigToScaffDict = toContigScafMap(scaffToContigListDict)

            #for all sequences/contigs
            #    get the name of the corresponding scaffold
            #    if the scaffold with this name doesn`t exist then create it
            #    bind the scaffold and the contig


            assSet = set()  # temp

            count = 0
            for s in self.sequences:

                #got the scaffold name that is not in the mapping!
                if s.name not in contigToScaffDict:
                    count += 1  # str('There is no mapping for contig: ' + s.name + ' to the corresponding scaffold')
                    scaffold = self._addScaff(s.name, None, None)  # name the corresponding scaffold as the contig name !!!
                else:
                    scaffoldName = contigToScaffDict[s.name]
                    assSet.add(s.name)
                    if scaffoldName in self.scaffNameToScaff:
                        scaffold = self.scaffNameToScaff[scaffoldName]
                    else:
                        scaffold = self._addScaff(scaffoldName, None, None)
                #bind the sequence/contig with the corresponding scaffold
                s.setScaffold(scaffold)
                scaffold.addNextContig(s)

            if count > 0:
                print str('There were ' + str(count) + ' contigs for which there was no defined scaffold in the mapping file')
                for n in contigToScaffDict:
                    if n not in assSet:
                        print n, #temp

        else:
            assert False


        #build index on top of all sequences and scaffolds
        for s in self.scaffolds:
            self.scaffoldsDict[s.id] = s
        for s in self.sequences:
            self.sequencesDict[s.id] = s

        #check for duplicate sequences #need to check all ids not only the last one !!!
        hashToSeqIdDict = dict([]) #map: hash -> list of seq ids
        for s in self.sequences:
            h = s.getHash()
            if h in hashToSeqIdDict:
                list = hashToSeqIdDict[h]
                for id in list:
                    sh = self.sequencesDict[id]
                    if s.getSeq() == sh.getSeq():
                        print str('Sequences "' + s.name + '" (' + str(s.id) + ') and "' + sh.name + '" (' + str(sh.id) + ') are identical!')
                list.append(s.id)
            else:
                temp = []
                temp.append(s.id)
                hashToSeqIdDict[h] = temp

        #check for duplicate scaffolds
        if scaffRead:
            hashToSeqIdDict = dict()  # map: hash -> list of scaff ids
            for s in self.scaffolds:
                if not s.getScaffSeqDef():  # the sequence for this scaffold is not defined
                    continue
                h = s.getHash()
                if h in hashToSeqIdDict:
                    list = hashToSeqIdDict[h]
                    for id in list:
                        sh = self.scaffoldsDict[id]
                        s1 = s.getScaff()
                        s2 = sh.getScaff()
                        assert ((s1 is not None) and (s2 is not None)), 'The scaffold sequence must be defined'
                        if s.getSeq() == sh.getSeq():
                            print str('Sequences "' + s.name + '" (' + str(s.id) + ') and "' + sh.name + '" (' + str(sh.id) + ') are identical!')
                    list.append(s.id)
                else:
                    temp = []
                    temp.append(s.id)
                    hashToSeqIdDict[h] = temp


    def _readContigsScaffolds(self, filePath, readContigs = True):
        """
            Read contigs or scaffolds from a file.
        """
        try:
            f = open(os.path.normpath(filePath),'r')
        except Exception:
            print "Cannot open file:", filePath
            raise
        else:
            name = ''
            seq = ''
            for line in f:
                line = noNewLine(line)
                if re.match('>', line):
                    if seq != '':
                        assert name != ''
                        if readContigs:
                            self._addSeq(name, seq)  # store seq
                        else:
                            self._addScaff(name, None, seq)
                        seq = ''
                    name = line.replace('>','')
                else:
                    seq += line
            if seq != '':
                assert name != ''
                if readContigs:
                    self._addSeq(name, seq) #store seq
                else:
                    self._addScaff(name, None, seq)
        finally:
            f.close()


    def _addSeq(self, name, seq):
        """
            Adds one sequence/contig.
        """
        s = Sequence((self._seqCounter + 1), name, seq)
        if s.seqBp >= self._minSeqLen:
            self.sequences.append(s)
            self._seqCounter += 1

    def get_sequence_count(self):
        return self._seqCounter

    def _addScaff(self, name, contig, scaffSeq):
        """
            Adds one scaffold.
        """
        self._scaffCounter += 1
        scaff = Scaffold(self._scaffCounter, name, contig, scaffSeq)
        assert name not in self.scaffNameToScaff, str('The scaffold name ' + name + ' already exists')
        self.scaffNameToScaff[name] = scaff
        self.scaffolds.append(scaff)
        if scaff.getScaffSeqDef() and (scaff.seqBp < self._minScaffLen):
            scaff.removeScaffSeq()
        return scaff


    def getSequence(self, seqId):
        assert seqId in self.sequencesDict, str('seqId: "' + seqId + '" is not in the Sequences')
        return self.sequencesDict[seqId]


    def setRemoveNonDna(self, removeNonDnaChars):
        """
            non-DNA characters won`t be considered in the..
        """
        for seq in self.sequences:
            seq.setRemoveNonDna(removeNonDnaChars)

        for scaff in self.scaffolds:
            scaff.setRemoveNonDna(removeNonDnaChars)


    def setCandidateTaxonomyPath(self, seqId, scaffoldId, taxPathDict, weight, source=None, tag=None):
        assert seqId in self.sequencesDict, str('seqId: ' + str(seqId))
        assert seqId not in self.placedSeqSet
        assert scaffoldId in self.scaffoldsDict
        sequence = self.sequencesDict[seqId]
        assert sequence.scaffold.id == scaffoldId
        assert len(taxPathDict.keys()) >= 1
        sequence.setCandidateTaxonomyPath(taxPathDict, weight, source, tag)


    def setScaffCandidateTaxonomyPath(self, scaffoldId, taxPathDict, taxonomy, weight):
        assert scaffoldId in self.scaffoldsDict
        assert len(taxPathDict.keys()) >= 1
        scaffold = self.scaffoldsDict[scaffoldId]
        assert len(scaffold.contigs) > 0
        for sequence in scaffold.contigs:
            tp = taxonomy.replicateTaxPathDict(taxPathDict)
            sequence.setCandidateTaxonomyPath(tp, weight)
        return len(scaffold.contigs)


    def _filterPlacements(self, weightsList, placementList, topPercentThreshold):

        assert len(placementList) == len(weightsList)
        cWeights = []
        assignments = []
        weight = None

        if len(placementList) >= 1:

            for e in weightsList:
                if e is None:
                    cWeights.append(50.0)  # if a weight is not known, it is considered to be 50.0
                else:
                    cWeights.append(e)
            #assignments with at least this weight will be considered
            assert len(cWeights) >= 1
            threshold = float(max(cWeights))*(1.0 - topPercentThreshold)
            weights = []
            for i in range(len(placementList)):
                if cWeights[i] >= threshold:
                    assignments.append(placementList[i])
                    weights.append(cWeights[i])
            weight = min(weights)

        return [assignments, weight]


    def setTaxonomyPathsFromCandidatePaths(self, taxonomy, topPercentThreshold):
        for seq in self.sequences:
            candidateTaxPathDictList = seq.getCandidateTaxPathDictList()
            candidateTaxPathDictListWeights = seq.getCandidateTaxPathDictWeightsList()
            assert len(candidateTaxPathDictList) == len(candidateTaxPathDictListWeights)
            result = self._filterPlacements(candidateTaxPathDictListWeights, candidateTaxPathDictList, topPercentThreshold)
            assignments = result[0]
            weight = result[1]
            taxPathDict = taxonomy.getLongestCommonPathFromMultipleAssignments(assignments)
            if taxPathDict is not None:
                self.setTaxonomyPath(seq.id, seq.scaffold.id, taxPathDict, weight)


    def setTaxonomyPath(self, seqId, scaffoldId, taxPathDict, weight):
        assert seqId in self.sequencesDict
        assert seqId not in self.placedSeqSet
        assert scaffoldId in self.scaffoldsDict

        sequence = self.sequencesDict[seqId]
        assert sequence.scaffold.id == scaffoldId
        assert len(taxPathDict.keys()) >= 1
        sequence.setTaxonomyPath(taxPathDict, weight)
        self.placedSeqSet.add(seqId)

    def delTaxonomyPath(self, sequence):
        assert sequence.id in self.sequencesDict
        assert sequence.id in self.placedSeqSet
        # the sequence is not placed anymore
        sequence.delTaxonomyPath()
        self.placedSeqSet.remove(sequence.id)


    def setTaxonomyPathOverride(self, seqId, scaffoldId, taxPathDict, weight):
        assert seqId in self.sequencesDict
        assert seqId in self.placedSeqSet
        assert scaffoldId in self.scaffoldsDict

        sequence = self.sequencesDict[seqId]
        assert sequence.scaffold.id == scaffoldId
        assert len(taxPathDict.keys()) >= 1
        sequence.setTaxonomyPath(taxPathDict, weight)


    def getTaxonomyPath(self, seqId):
        assert seqId in self.sequencesDict, 'Ask for a wrong id of a sequence!'
        if seqId in self.placedSeqSet:
            return self.sequencesDict[seqId].getTaxonomyPath()
        else:
            return None


    def placeContigsFromTheSameScaffold(self, taxonomy, agThreshold, assignedPartThreshold, topPercentThreshold):
        """
            Place all contigs from the same scaffold.

            @param agThreshold: agreement threshold - all assignments must lie on the same path and below the line given
                                by this threshold
            @param assignedPartThreshold:
        """
        for scaffold in self.scaffolds:
            placedSeqList = []
            notPlacedSeqList = []
            placedBp = 0
            notPlacedBp = 0
            for seq in scaffold.contigs:
                if seq.id in self.placedSeqSet:
                    placedSeqList.append(seq)
                    placedBp += seq.seqBp
                else:
                    notPlacedSeqList.append(seq)
                    notPlacedBp += seq.seqBp

            # assigned part of a scaffold in <0,1>
            if (placedBp + notPlacedBp) == 0:
                assingnedPart = 0
                bf = "Zero length sequences:\n"
                for zls in scaffold.contigs:
                     bf += str(zls.name + str(zls.seqBp) + '\n')
            else:
                assignedPart = (float(placedBp)/(float(placedBp) + float(notPlacedBp)))
            if assignedPart >= assignedPartThreshold:
                assignNotAssigned = True
            else:
                assignNotAssigned = False

            # place all not placed sequences based on one placed sequence
            if (len(placedSeqList) == 1) and (len(notPlacedSeqList) > 0) and assignNotAssigned:
                taxPathDict = placedSeqList[0].getTaxonomyPath()
                assert (taxPathDict is not None) and (len(taxPathDict) > 0)
                for seq in notPlacedSeqList:
                    t = taxonomy.replicateTaxPathDict(taxPathDict)
                    self.setTaxonomyPath(seq.id, scaffold.id, t, placedSeqList[0].getTaxonomyPathWeight())

            # place all sequences to the agreement based on the placed sequences
            elif len(placedSeqList) > 1:

                taxPathDictList = []
                weightsList = []
                for seq in placedSeqList:
                    tp = seq.getTaxonomyPath()
                    assert (tp is not None) and (len(tp) > 0)
                    taxPathDictList.append(tp)
                    weightsList.append(seq.getTaxonomyPathWeight())

                # get only placements with top perCent weights
                resultA = self._filterPlacements(weightsList, taxPathDictList, topPercentThreshold)
                assignmentsA = resultA[0]
                weightA = resultA[1]

                taxPathDict = taxonomy.getLongestCommonPathFromMultipleAssignments2(assignmentsA, agThreshold)

                if taxPathDict is not None:

                    for seq in placedSeqList:
                        t = taxonomy.replicateTaxPathDict(taxPathDict)
                        self.setTaxonomyPathOverride(seq.id, scaffold.id, t, weightA)

                    if assignNotAssigned:
                        for seq in notPlacedSeqList:
                            t = taxonomy.replicateTaxPathDict(taxPathDict)
                            self.setTaxonomyPath(seq.id, scaffold.id, t, weightA)
                else:
                    for seq in placedSeqList:
                        self.delTaxonomyPath(seq)




    def writeScaffolds(self, outFile, writeIds=True, outputFileContigSubPattern=None):
        self._writeSequencesOrScaffolds(outFile, False, writeIds, outputFileContigSubPattern) #False ~ write scaffolds


    def writeSequences(self, outFile, writeIds=True, outputFileContigSubPattern=None):
        self._writeSequencesOrScaffolds(outFile, True, writeIds, outputFileContigSubPattern) #True ~ write sequences


    def _writeSequencesOrScaffolds(self, outFile, writeSeq=True, writeIds=True, outputFileContigSubPattern=None):
        try:
            f = open(os.path.normpath(outFile), 'w')
            k = 0
            if writeSeq:
                list = self.sequences  # write sequences
            else:
                list = self.scaffolds  # write scaffolds

            for seq in list:
                if seq.seqBp == 0:
                    continue

                if writeIds:
                    if writeSeq:
                        name = str(str(seq.scaffold.id) + '_' + str(seq.id))  # sequence
                    else:
                        name = str(seq.id)  # scaffold
                else:
                    if outputFileContigSubPattern is None:
                        name = seq.name
                    else:
                        name = re.sub(outputFileContigSubPattern, r'\1' , seq.name)
                if k == 0:
                    f.write('>' + name)
                    k += 1
                else:
                    f.write('\n>' + name)
                s = seq.getSeq()
                l = seq.seqBp
                i = 0
                range = self._fastaLineMaxChar
                while i < l:
                    j = i + range
                    f.write('\n' + s[i:j])
                    i += range
        except Exception:
            print "Cannot create a file or write to it:", outFile
            raise
        finally:
            f.close()


    def writeSeqNameSeqId(self, mapFilePath):
        """
            Creates a file that contains mapping: sequence name -> sequence id.
        """
        self._writeSeqOrScaffNameTabId(mapFilePath, True) #True ~ do it for sequences


    def writeScaffNameScaffId(self, mapFilePath):
        """
            Creates a file that contains mapping: scaffold name -> scaffold id.
        """
        self._writeSeqOrScaffNameTabId(mapFilePath, False) #False ~ do it for scaffolds


    def _writeSeqOrScaffNameTabId(self, mapFilePath, forSequences = True):
        """
            Creates a file that represents mapping: sequence/contig or scaffold name -> corresponding_id (name tab id).
            @param forSequences: True ~ do it for sequences, False ~ do it for scaffolds.
        """
        try:
            f = open(os.path.normpath(mapFilePath), 'w')
            k = 0
            if forSequences:
                list = self.sequences  #for sequences
            else:
                list = self.scaffolds  #for scaffolds

            for seq in list:
                if k == 0:
                    k += 1
                    f.write(seq.name + '\t' + str(seq.id))
                else:
                    f.write('\n' + seq.name + '\t' + str(seq.id))
        except Exception:
            print "Cannot create a file or write to it:", mapFilePath
            raise
        finally:
            f.close()


    def writeScaffoldContigMap(self, outFilePath):
        """
            Creates file: scaffold_id tab contig_id.
        """
        try:
            f = open(os.path.normpath(outFilePath), 'w')
            k = 0
            for scaffold in self.scaffolds:
                for contig in scaffold.contigs:
                    entry = str(str(scaffold.id) + '\t' + str(scaffold.id) + '_' + str(contig.id))
                    if k == 0:
                        k += 1
                        f.write(entry)
                    else:
                        f.write('\n' + entry)
        except Exception:
            print "Cannot create a file or write to it:", outFilePath
            raise
        finally:
            f.close()


    def writePlacementsPPOut(self, outFile, taxaRanks, outputFileContigSubPattern):

        try:
            f = open(os.path.normpath(outFile), 'w')

            f.write('#Output of pPPS\n#\n'),
            header = str('#ID' + '\t' + 'root')
            for rank in taxaRanks:
                header += str('\t' + rank)
            f.write(header)

            for seq in self.sequences:
                entry = str('\n' + re.sub(outputFileContigSubPattern, r'\1' , seq.name))
                taxPathDict = seq.getTaxonomyPath()
                if taxPathDict is None:
                    entry += str('\t')
                else:
                    entry += str('\t' + 'root')
                for rank in taxaRanks:
                    if (taxPathDict is not None) and (rank in taxPathDict) and (not taxPathDict[rank].isCopy()):
                        entry += str('\t' + taxPathDict[rank].name)
                    else:
                        entry += '\t'
                f.write(entry)
        except Exception:
            print "Cannot create a file or write to it:", outFile
            raise
        finally:
            f.close()


    def writePlacementsOut(self, outFile, taxaRanks, outputFileContigSubPattern):

        try:
            f = open(os.path.normpath(outFile), 'w')
            f.write('# SEQUENCEID	TAXID')
            # k = 0

            for seq in self.sequences:

                taxPathDict = seq.getTaxonomyPath()
                ncbid = 1
                for rank in taxaRanks:
                    if ((taxPathDict is not None) and (rank in taxPathDict)):
                        ncbid = taxPathDict[rank].ncbid
                    else:
                        break

                if ncbid == 1:
                    continue

                entry = (noNewLine(re.sub(outputFileContigSubPattern, r'\1' , seq.name)) + '\t' + str(ncbid))

                # if k == 0:
                #     f.write(entry)
                #     k += 1
                # else:
                f.write('\n' + entry)

        except Exception:
            print "Cannot create a file or write to it:", outFile
            raise
        finally:
            f.close()


    def discardInconsistentPlacements(self, discardContigThreshold, discardScaffThreshold, taxonomicRanks):
        discardedContigsCount = 0
        print 'discard placements: "discardContigThreshold:"', discardContigThreshold, ' discardScaffThreshold:', discardScaffThreshold
        totalBp = 0
        cumulativeCons = 0.0
        #for all scaffolds
        for scaff in self.scaffolds:
            #what is the deapest assignment of a contig
            depth = 0
            for contig in scaff.contigs:
                if contig.getTaxonomyPath() is not None:
                    d = len(contig.getTaxonomyPath())
                    if d > depth:
                        depth = d
            #get deapest assignment with the longest sequence
            length = 0
            deapest = None
            for contig in scaff.contigs:
                if contig.getTaxonomyPath() is not None:
                    if len(contig.getTaxonomyPath()) == depth:
                        if contig.seqBp > length:
                            deapest = contig
                            length = contig.seqBp
            if deapest == None:
                print 'discardInconsistentPlacements: No contig of scaffold assigned, scaff id:', scaff.id
                continue

            consistentList = []
            inconsistentList = []
            consistentBp = 0
            inconsistentBp = 0
            td = deapest.getTaxonomyPath()
            for contig in scaff.contigs:
                c = contig.getTaxonomyPath()
                if c != None:
                    rank = taxonomicRanks[len(c) - 1]
                if (c == None) or (td[rank].ncbid == c[rank].ncbid):
                    consistentList.append(contig)
                    consistentBp += contig.seqBp
                else:
                    inconsistentList.append(contig)
                    inconsistentBp += contig.seqBp

            ciBp = consistentBp + inconsistentBp
            consistency = (float(consistentBp)/float(ciBp))
            totalBp += ciBp
            cumulativeCons += consistency*float(ciBp)

            #discard inconsistent contigs
            if consistency < discardContigThreshold:
                for contig in inconsistentList:
                    self.delTaxonomyPath(contig)
                    discardedContigsCount += 1
            #discard inconsistent contigs
            if consistency < discardScaffThreshold:
                for contig in consistentList:
                    self.delTaxonomyPath(contig)
                    discardedContigsCount += 1

        print "Total consistency of PPS:", cumulativeCons/float(totalBp)

        return discardedContigsCount


def toScafContigMap(scafContigFile):
    """
        Reads scaffold contig mapping.

        @param scafContigFile: scaffold-contig mapping (tab separated)

        @return: map: scaffold -> list of contigs
    """
    scafToContigs = dict()
    try:
        f = open(os.path.normpath(scafContigFile),'r')
    except Exception:
        print "Cannot open file:", scafContigFile
        raise
    else:
        for line in f:
            line = noNewLine(line)
            scaffold = re.sub(r'^[ ]*([^\t]+)\t[^\t]*',r'\1', line)# gap deleted !!!
            contig = re.sub(r'^[ ]*[^\t]+\t([^\t]*)',r'\1', line)
            if scaffold in scafToContigs:
                scafToContigs[scaffold].append(contig)
            else:
                temp = []
                temp.append(contig)
                scafToContigs[scaffold] = temp

    return scafToContigs


def toContigScafMap(scaffToContigListDict):
    """
        Transforms the input dictionary (key -> list of items) into the inverse dictionary (item -> key).
        To each item (contig) only one key (scaffold) exists.

        @param scaffToContigListDict: map: scaffold -> list of contigs

        @return: map: contig name -> scaffold name
    """
    contigToScaffDict = dict()
    for scaff in scaffToContigListDict:
        for contig in scaffToContigListDict[scaff]:
            assert contig not in contigToScaffDict, str('Contig: ' + contig + ' is not allowed to be twice in the directory!')
            contigToScaffDict[contig] = scaff
    return contigToScaffDict


def replaceIdsWithNames(outputFileContigSubPattern, nameToIDsFile, targetFile, outFile):
    """
        @deprecated: NOT IMPLEMENTED YET!!!
        replace ids with names
        @param nameToIdsFile: file that contains lines: contigName tab contigID
        @param targetFile: file that contain in the first column scaffoldID_contigID which will be replaced by its name
        @param outFile: file that contain the first column in the form scaffoldID_contigID with the name
        (that can be modified by substitution defined in the config file .. according to outputFileContigSubPattern)
    """
    idToName = dir([])
    assert False, 'NOT IMPLEMENTED YET'
    #try:
    #    f = open(os.path.normpath(nameToIDsFile), 'r')
    #    for line in f:
    #        if re.match('^#', line):
    #            continue
    #        name = re.sub(outputFileContigSubPattern, r'\1' , noNewLine(re.sub(r'^([^ \t]+)\t[0-9]+$',r'\1', line)))
    #        id = int(noNewLine(re.sub(r'^[^ \t]+\t([0-9]+)$',r'\1', line)))
    #        idToName[id] = name
    #except Exception:
    #    print "Cannot create a file or write to it:", outFile
    #    raise
    #finally:
    #    f.close()

    #now: go through the targetFile and for each line do:
    #    extract contigID and the rest of the line ^[0-9]+_[0-9]+([^0-9].*)$
    #    write name + rest of the line + \n to the outFile !!!!!!!!!!



#compare two sequences according to their length
#def seqLenCmp(seq1, seq2):
#    return seq1.seqBp - seq2.seqBp


def seqWeightThenLenCmp(seq1, seq2):
    """
        Compares two sequences according to the assignment weight and if it`s the same or not defined,
        it compares them according to their length.
    """
    weight1 = seq1.getTaxonomyPathWeight()
    weight2 = seq2.getTaxonomyPathWeight()
    if (weight1 is not None) and (weight2 is not None):
        if (weight1 - weight2) < 0:
            return -1
        else:
            return 1
    else:
        return seq1.seqBp - seq2.seqBp



#---------------------------------------------------
def _test():
#s = Sequence(1, 'seqName','AATTGGCCC\n\rAAA\n')
#print 'sequence name: ', s.seqName
#print 'sequence:', s.getSeq()
#print 'seqBp:', s.seqBp, '({0})'.format(len(s.seqCompressed))

    config = Config(open(os.path.normpath('D:\\A_Phylo\\A_Metagenomic\\pPPS\\workspace\\pPPS\\config01.cfg')), 'pPPS')
    outputFileContigSubPattern = config.get('outputFileContigSubPattern')

    nameToIDsFile = os.path.normpath('D:\\A_Phylo\\A_Metagenomic\\pPPS\\workspace\\pPPS\\wdir02\\inputTW.fas.cToIds')
    targetFile = os.path.normpath('D:\\A_Phylo\\A_Metagenomic\\pPPS\\workspace\\pPPS\\wdir02\\inputTW.fas.ids.out')
    outFile = os.path.normpath('D:\\A_Phylo\\A_Metagenomic\\pPPS\\workspace\\pPPS\\wdir02\\inputTW.fas.pOUT')

    replaceIdsWithNames(outputFileContigSubPattern, nameToIDsFile, targetFile, outFile)



    #s = Sequences(configpPPS)
    ##s = Sequences('D:\\A_Phylo\\A_Metagenomic\\pPPS\\workspace\\pPPS\\input.fas_sub')
    #s.writeSequences(configpPPS.get('inputIdsFastaFile'))

    #for seq in s.sequences:
    #    print seq.id, seq.name, seq.scaffold.id, seq.getSeq()
    #for scaffold in s.scaffolds:
    #    print scaffold.id, scaffold.name, ":",
    #    for contig in scaffold.contigs:
    #        print contig.name, contig.id, ",",
    #    print ""


if __name__ == "__main__":
    _test()