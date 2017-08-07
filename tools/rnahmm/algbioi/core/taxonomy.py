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
import sqlite3
import datetime

from algbioi.com.config import Config


class Node():
    def __init__(self, ncbid, rank, name):
        self.ncbid = ncbid
        self.rank = rank
        self.name = name
        self._copy = False

    def __str__(self):
        return str(self.rank + ' ' + str(self.ncbid) + ' ' + self.name + ' ' + str(self._copy))

    def copy(self, newRank):
        node = Node(self.ncbid, newRank, self.name)
        node._copy = True
        return node

    def isCopy(self):
        return self._copy

    def clone(self):
        node =  Node(self.ncbid, self.rank, self.name)
        node._copy = self._copy
        return node

#---------------------------------------------------
class Taxonomy():
    def __init__(self, databaseFile, taxonomicRanks):#taxonomicRanks is a list of taxonomic ranks
        self.taxonomicRanks = taxonomicRanks
        self.taxonomicRanksSet = set(taxonomicRanks)
        self.ncbidToPathToRoot = dict([]) #buffer: ncbid -> path to root
        self.newOtuCounter = 0
        try:
            self.conn = sqlite3.connect(os.path.normpath(databaseFile))
            self.cursor = self.conn.cursor()
        except Exception:
            print "Failed to create connection to a database:", databaseFile
            raise

    def createNewOtuDBEntry(self, parentNcbid, sampleName, rank='species'):
        #get max taxon_id
        self.cursor.execute(str('SELECT max(taxon_id) FROM taxon'))
        result = self.cursor.fetchall()
        maxTaxonId = int(result[0][0])
        self.cursor.execute(str('SELECT max(ncbi_taxon_id) FROM taxon'))
        result = self.cursor.fetchall()
        maxNcbiTaxonId = int(result[0][0])
        self.cursor.execute(str('SELECT max(taxon_id) FROM taxon_name'))
        result = self.cursor.fetchall()
        maxTaxonId2 = int(result[0][0])
        #print maxTaxonId, maxTaxonId2, maxNcbiTaxonId
        assert maxTaxonId == maxTaxonId2

        #insert new entries
        self.newOtuCounter += 1
        #table taxon:
        taxon_id = maxTaxonId + 1
        ncbi_taxon_id = maxNcbiTaxonId + 1
        parent_taxon_id = parentNcbid
        node_rank = rank
        now = datetime.datetime.now()
        timeStamp = str(str(now.year) + str(now.month) + str(now.day) + '_' + str(now.hour) + str(now.minute) + str(now.second))
        name = str('OTU' + str(self.newOtuCounter) + '_' + timeStamp + '_' + sampleName)
        name_class = 'scientific name'
        #insert into tables
        self.cursor.execute(str('INSERT INTO taxon (taxon_id, ncbi_taxon_id, parent_taxon_id, node_rank)' +
        ' VALUES (?,?,?,?)'),(taxon_id, ncbi_taxon_id, parent_taxon_id, node_rank))
        self.cursor.execute(str('INSERT INTO taxon_name (taxon_id, name, name_class) VALUES (?,?,?)'),
                            (taxon_id, name, name_class))
        self.conn.commit()

        #print taxon_id, ncbi_taxon_id
        return ncbi_taxon_id
        #785593 1050858
        #delete from taxon where taxon_id = 785592;

    def close(self):
        self.cursor.close()
        self.conn.close()

    #hang the new ncbid below the parent
    #return ncbid
    def getNewNCBID(self, parentNcbid):
        pass


    def getNcbidToName(self, ncbid):
        taxonNcbid = ncbid
        self.cursor.execute(str('SELECT TN.name FROM taxon_name TN, taxon T WHERE T.ncbi_taxon_id=?' +
                                ' AND T.taxon_id = TN.taxon_id AND TN.name_class="scientific name"'),(taxonNcbid,))
        result = self.cursor.fetchall()
        assert len(result) == 1, str('Cannot find name for ncbi: ' + str(taxonNcbid))
        name = result[0][0]
        return name



    def getPathToRoot(self, ncbid):
        assert ncbid != None
        pathDict = dict([]) #map: rank -> Node
        #if we are in the root, empty dict is returned
        if ncbid == 1:
            return None #!!!!!!!!!!

        #check buffer
        if ncbid in self.ncbidToPathToRoot:
            return self.replicateTaxPathDict(self.ncbidToPathToRoot[ncbid])

        taxonNcbid = ncbid
        taxonRank = ''
        history = [] #temp

        while True:
            history.append(taxonNcbid) #temp
            if taxonNcbid == 1: #the root of the taxonomy reached
                break
            self.cursor.execute('SELECT taxon_id FROM taxon T WHERE T.ncbi_taxon_id=?',(taxonNcbid,))
            result = self.cursor.fetchall()
            if len(result) != 1:
                sys.stderr.write('Taxonomy: Cannot find taxon_id for ncbi:' + str(taxonNcbid) + ' result:' + str(result) + ' \n')
                return None # !!!!!!!!!
            taxonId = result[0][0]

            self.cursor.execute('SELECT node_rank FROM taxon T WHERE T.taxon_id=?', (taxonId,))
            result = self.cursor.fetchall()
            assert len(result) == 1, str('Cannot find rank for taxon_id: ' + str(taxonId))
            taxonRank = result[0][0]

            self.cursor.execute('SELECT name FROM taxon_name TN WHERE TN.taxon_id=? AND name_class="scientific name"', (taxonId,))
            result = self.cursor.fetchall()
            assert len(result) == 1, str('Cannot find scientific name for taxon_id: ' + str(taxonId))
            taxonName = result[0][0]

            if taxonRank in self.taxonomicRanksSet:
                pathDict[taxonRank] = Node(taxonNcbid, taxonRank, taxonName)
            #else:
            #    if taxonNcbid != 1:
            #        if taxonRank == 'no rank':
            #            print 'Rank ignored:', taxonRank, 'for:', taxonName, 'for ncbid:', taxonNcbid
            #        else:
            #            print 'Unknown rank:', taxonRank, 'for:', taxonName, 'for ncbid:', taxonNcbid

            #get parent
            self.cursor.execute('SELECT parent_taxon_id FROM taxon T WHERE T.taxon_id=?', (taxonId,))
            result = self.cursor.fetchall()
            assert len(result) == 1, str('Cannot find parent for taxon_id', taxonId)
            taxonNcbid = result[0][0]

        #test if the path is correct and if a rank is missing along the path, it adds it as a copy of an existing lower node
        revRanks = list(self.taxonomicRanks)
        revRanks.reverse()
        deapestRank = None
        for rank in revRanks:
            if rank in pathDict:
                deapestRank = rank
                break
        assert deapestRank != None, '%s' % (ncbid)

        lastNode = None
        for rank in revRanks:

            if rank == deapestRank:
                lastNode = pathDict[rank]
                continue

            if lastNode == None:
                continue

            if rank not in pathDict:
                pathDict[rank] = lastNode.copy(rank)
                #print str('Set rank "' + rank + '" for "' + pathDict[deapestRank].rank + ' ' + pathDict[deapestRank].name +
                #    ' (' + str(pathDict[deapestRank].ncbid) + ')" as: "' + lastNode.rank + ' ' + lastNode.name + ' (' + str(lastNode.ncbid) + ')"')
            else:
                lastNode = pathDict[rank]

        #printing
        #print '------------------------'
        #for r in self.taxonomicRanks:
        #    if r not in pathDict:
        #        break
        #    print pathDict[r].rank, pathDict[r].name, pathDict[r].ncbid
        ########

        #store the result into a buffer:
        self.ncbidToPathToRoot[ncbid] = self.replicateTaxPathDict(pathDict)

        return pathDict


    def getPathToRootSemicolonSeparated(self, ncbid):
        """
            Return the path from a clade to the root as semicolon separated ncbids.
        """
        pathDict = self.getPathToRoot(int(ncbid))
        outBuffer = ''
        for rank in self.taxonomicRanks:
            if rank in pathDict:
                node = pathDict[rank]
                assert node.rank == rank
                outBuffer += str(node.ncbid) + ';'
                if int(node.ncbid) == int(ncbid):
                    break

        if outBuffer == '':
            return None
        else:
            return outBuffer


    def getLongestCommonPathFromMultipleAssignments2(self, taxPathDictList, threshold):
        """
        !!!!!!!!!!!

        @param taxPathDictList: list of taxPathDicts
        @param threshold: percent that says which fraction of the assigned sequences must lie at least at this rank (or lower)

        @return: - longest path s.t. all placements up to this rank lie on this path, moreover at least X%
        of the sequences that are assigned are assigned at least to this rank or lower
        """

        assert len(taxPathDictList) > 0
        if len(taxPathDictList) == 1:
            return taxPathDictList[0]

        #longest common path
        lpTaxPathDict = self.getLongestCommonPathFromMultipleAssignments(taxPathDictList)
        if lpTaxPathDict == None:
            print 'getLongestCommonPathFromMultipleAssignments2: Longest common path is 0'
            return None
        #min number of sequences assigned to the final level:

        minCount = float(len(taxPathDictList))*float(threshold)
        depth = len(lpTaxPathDict)
        while depth >= 0:
            count = 0
            for taxPathDict in taxPathDictList:
                if len(taxPathDict) >= depth:
                    count += 1
            if count >= minCount:
                break
            else:
                depth -= 1

        if depth == 0:
            sys.stderr.write('Taxonomy:getLongestCommonPathFromMultipleAssignments2: The depth must be at least 1'
                             + str(taxPathDictList) + '\n')
            return None

        taxPathDict = dict([])
        for rank in self.taxonomicRanks:
            taxPathDict[rank] = lpTaxPathDict[rank].clone()
            depth -= 1
            if depth == 0:
                break
        assert depth == 0, 'The number of elements in the reference dictionary is too small'
        return taxPathDict


    def replicateTaxPathDict(self, taxPathDict):
        """
            Replicate the dictionary. !!!!!!!!!!!
        """
        taxPathDict2 = dict([])
        for rank in taxPathDict:
            taxPathDict2[rank] = taxPathDict[rank].clone()
        return taxPathDict2


    def getLongestCommonPathFromMultipleAssignments(self, taxPathDictList):
        """
            input: list of taxPathDicts
            output: taxPathDict - longest path s.t. all placements up to this rank lie on this path
        """
        if len(taxPathDictList) == 0:
            return None

        ncbiLowest = 1
        longestPath = taxPathDictList[0]
        for t in taxPathDictList[1:]:
            if len(longestPath) < len(t):
                longestPath = t

        for r in self.taxonomicRanks:
            if r not in longestPath:
                break
            ncbi1 = longestPath[r].ncbid
            ncbiN = ncbi1
            for entry in taxPathDictList:
                if r in entry:
                    ncbiI = entry[r].ncbid
                    if ncbiI != ncbi1:
                        ncbiN = ncbiI
                        break
            if ncbiN != ncbi1:
                break
            else:
                ncbiLowest = ncbi1

        return self.getPathToRoot(ncbiLowest)


    def getPathFromLowestCommonAncestorToRoot(self, listOfNcbid):

        if len(listOfNcbid) == 1:
            return self.getPathToRoot(listOfNcbid[0])

        listOfPlacements = []
        for ncbid in listOfNcbid:
            path = self.getPathToRoot(ncbid)
            if path != None:
                listOfPlacements.append(path)

        if len(listOfPlacements) == 1:
            return listOfPlacements[0]
        elif len(listOfPlacements) == 0:
            return None

        lowestCommonNcbid = 1
        for rank in self.taxonomicRanks:
            if rank not in listOfPlacements[0]:
                break
            ncbid1 = listOfPlacements[0][rank].ncbid
            ncbidN = ncbid1
            for placementI in listOfPlacements[1:]:
                if rank not in placementI:
                    ncbidN = 1
                    break
                ncbidI = placementI[rank].ncbid
                if ncbidI != ncbid1:
                    ncbidN = ncbidI
                    break

            if ncbid1 != ncbidN:
                break
            else:
                lowestCommonNcbid = ncbid1

        return self.getPathToRoot(lowestCommonNcbid)


    def getNcbidFromScientificName(self, scientificName):
        """
            Gets NCBID of a clade for which we know its scientific name.
            @return ncbid or None
        """
        self.cursor.execute(str('SELECT T.ncbi_taxon_id FROM taxon_name TN, taxon T ' +
                                'WHERE TN.name_class="scientific name" AND TN.name=? AND TN.taxon_id=T.taxon_id'),
                                (scientificName,))
        result = self.cursor.fetchall()
        if len(result) == 1:
            return int(result[0][0])
        else:
            assert len(result) == 0, str('The scientific name ' + scientificName + 'is ambiguous!')
            return None


#---------------------------------------------------
def test1():

    config = Config(open(os.path.normpath('D:\\A_Phylo\\A_Metagenomic\\pPPS\\workspace\\pPPS\\config01.cfg')), 'pPPS')

    databaseFile = os.path.normpath(config.get('databaseFile'))
    taxonomicRanks = config.get('taxonomicRanks').split(',')
    t = Taxonomy(databaseFile, taxonomicRanks)

    s = 'Assignment of query to the lowest common ancestor of Bacteroides thetaiotaomicron (226186), Porphyromonas gingivalis (242619) and Parabacteroides distasonis (435591).'
    #s = 'Assignment of query to the lowest common ancestor of Bacteroides thetaiotaomicron , Porphyromonas gingivalis (242619) and Parabacteroides distasonis (435591).'
    #s = 'Assignment of query to the lowest common ancestor of Bacteroides thetaiotaomicron , Porphyromonas gingivalis  and Parabacteroides distasonis (435591).'
    #s = 'Assignment of query to the lowest common ancestor of Bacteroides thetaiotaomicron (226186) and Parabacteroides distasonis (435591).'
    #s = 'Assignment of query to the lowest common ancestor of Halobacterium sp. (64091), Thermococcus kodakarensis (69014), Pyrococcus horikoshii (70601), Methanothermobacter thermautotrophicus (187420), Methanopyrus kandleri AV19 (190192), Methanosarcina mazei (192952), Archaeoglobus fulgidus (224325), Methanocaldococcus jannaschii (243232), Methanococcoides burtonii (259564), Methanococcus maripaludis S2 (267377), Haloarcula marismortui (272569), Methanospirillum hungatei (323259), Methanosphaera stadtmanae (339860), Natronomonas pharaonis (348780), Methanosaeta thermophila PT (349307), Candidatus methanoarchaeon RC1 (351160), Haloquadratum walsbyi (362976), Methanoculleus marisnigri JR1 (368407), Methanocorpusculum labreanum (410358) and Methanobrevibacter smithii (420247).'
    #s = 'Assignment of query to the lowest common ancestor of Thermococcus kodakarensis (69014) and Pyrococcus horikoshii (70601).'
    list = re.findall(r'\([0-9]+\)', s)
    list2 = []
    for i in list:
        str = re.sub(r'\(([0-9]+)\)', r'\1', i)
        list2.append(str)
        print str
    print '-------------------------'


    pathDict = t.getPathFromLowestCommonAncestorToRoot(list2)
    #pathDict = t.getPathToRoot(170187)

    for k in taxonomicRanks:
        if k not in pathDict:
            break
        n = pathDict[k]
        print n.ncbid, n.rank, n.name

def test2():
    config = Config(open(os.path.normpath('D:\\A_Phylo\\A_Metagenomic\\pPPS\\workspace\\pPPS\\config01.cfg')), 'pPPS')

    databaseFile = os.path.normpath(config.get('databaseFile'))
    taxonomicRanks = config.get('taxonomicRanks').split(',')
    t = Taxonomy(databaseFile, taxonomicRanks)

    taxPathDictList = []
    taxPathDictList.append(t.getPathToRoot(33958))#Lactobacillaceae
    taxPathDictList.append(t.getPathToRoot(91061))#Bacilli
    taxPathDictList.append(t.getPathToRoot(2))#Bacteria
    #taxPathDictList.append(t.getPathToRoot(1385))#Bacillales
    taxPathDictList.append(t.getPathToRoot(1578))#Lactobacilus
    #taxPathDictList.append(t.getPathToRoot(31979))#Clostridiaceae
    taxPathDictList.append(t.getPathToRoot(2))#Bacteria

    taxPathDict = t.getLongestCommonPathFromMultipleAssignments(taxPathDictList)

    for key in taxPathDict:
        print key, taxPathDict[key]

def test3():
    config = Config(open(os.path.normpath('/Users/ivan/Documents/work/binning/tests/CowRumen/03/config.cfg')), 'pPPS')
    databaseFile = os.path.normpath(config.get('databaseFile'))
    taxonomicRanks = config.get('taxonomicRanks').split(',')
    t = Taxonomy(databaseFile, taxonomicRanks)
    parentNcbid = 1239 #Firmicutes
    sampleName = 'test_sample'
    rank = 'species'
    t.createNewOtuDBEntry(parentNcbid, sampleName, rank)
    t.close()


if __name__ == "__main__":
  test3()