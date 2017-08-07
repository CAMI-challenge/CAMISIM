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
import glob
import sqlite3

from algbioi.com.config import Config
from algbioi.com import taxonomy_ncbi


class RefSequences():
    def __init__(self, refDir, databaseFilePath):
        """
            Provides information about NCBI reference sequences (genomes or draft genomes).

            @param refDir: directory that contains reference sequences,
                each file has format ncbi_taxon_id.[0-9]+.fna(fas), for instance 382638.1.fna or 2110.1.fas
            @param databaseFilePath: ncbi taxonomy file in sqlite3 format
        """
        assert os.path.isdir(refDir)
        assert os.path.isfile(databaseFilePath)
        self._taxonIdSet = set()  # taxonIds in the reference
        self._taxonIdToSize = {}  # taxonId -> cumulative file size
        for fileName in os.listdir(refDir):
            if fileName.endswith(('.fna', '.fas')):
                taxonId = int(fileName[0:fileName.index('.')])
                self._taxonIdSet.add(taxonId)
                fileSize = int(os.path.getsize(os.path.join(refDir, fileName)))
                if taxonId in self._taxonIdToSize:
                    self._taxonIdToSize[taxonId] += fileSize
                else:
                    self._taxonIdToSize[taxonId] = fileSize
        self._taxonomy = taxonomy_ncbi.TaxonomyNcbi(databaseFilePath, considerNoRank=True)
        self._childrenBuffer = {}  # taxonId -> set of children taxon Ids
        self._rankBuffer = {}  # taxonId -> rank


    def _getChildren(self, taxonId):
        """
            Returns a set of taxonIds of child clades.
            @rtype: set
        """
        if taxonId in self._childrenBuffer:
            return self._childrenBuffer[taxonId]
        else:
            childrenList = self._taxonomy.getChildrenNcbids(taxonId)
            if childrenList is None:
                children = None
            else:
                children = set(childrenList)
            self._childrenBuffer[taxonId] = children
            return children

    def _getRank(self, taxonId):
        """
            Returns taxonomic rank of given taxonId.
            @rtype: str
        """
        if taxonId in self._rankBuffer:
            return self._rankBuffer[taxonId]
        else:
            rank = self._taxonomy.getRank(taxonId)
            self._rankBuffer[taxonId] = rank
            return rank

    def isRefSufficient(self, taxonId, minTotalCount, minBpPerSpeciesCount):
        """
            Returns true if the reference contains sufficient number of sequences to model a clade.

            @param taxonId: search reference for this clade
            @param minTotalCount: minimum number of reference sequences for this clade
            @param minBpPerSpeciesCount: a sequence file will be counted only proportionally if it contains less
                than this number of bp
            @rtype: bool
        """
        assert isinstance(taxonId, int)
        assert isinstance(minTotalCount, int)
        assert isinstance(minBpPerSpeciesCount, int)
        info = self._Info(taxonId, minTotalCount, minBpPerSpeciesCount)
        self._collectSpecies(taxonId, info)
        #print info.inspectedSpecTaxonId, info.collected  # debug !!!
        return info.isSufficient()

    class _Info():
        def __init__(self, taxonId, minTotalCount, minBpPerSpeciesCount):
            """
                Helper class that contains the current state of the search.
            """
            self.taxonId = taxonId
            self.minTotalCount = minTotalCount
            self.minBpPerSpeciesCount = minBpPerSpeciesCount
            self.collected = 0.0001  # number of sequences collected so far
            self.specTaxonIdToCount = {}
            self.inspectedSpecTaxonId = set()  # species taxonIds for which there is enough reference data

        def isSufficient(self):
            """"  Returns true if enough data has been found so far. """
            if self.collected > self.minTotalCount:
                return True
            else:
                return False

        def isSpecSufficient(self, specTaxonId):
            """ Returns true if enough data has been found for this species. """
            return specTaxonId in self.inspectedSpecTaxonId

        def update(self, specTaxonId, size):
            """"
                Updates the current state.

                @param specTaxonId: species taxon id for which some data was found
                @param size: file size of a subspecies
            """
            if specTaxonId in self.inspectedSpecTaxonId:  # this species has already been inspected
                return

            if size < self.minBpPerSpeciesCount:
                count = float(size) / float(self.minBpPerSpeciesCount)
            else:
                count = 1.0001

            if specTaxonId in self.specTaxonIdToCount:
                if self.specTaxonIdToCount[specTaxonId] < 1.0:
                    self.specTaxonIdToCount[specTaxonId] += count
            else:
                self.specTaxonIdToCount[specTaxonId] = count

            if self.specTaxonIdToCount[specTaxonId] > 1.0:
                self.collected += 1.0
                self.inspectedSpecTaxonId.add(specTaxonId)

    def _collectSpecies(self, taxonId, info):
        """
            @param taxonId: taxon id of species or higher ranks (or any if species not defined)
            @type info: _Info
        """
        if info.isSufficient():
            return

        if self._getRank(taxonId) == "species":
            self._collectStrains(taxonId, taxonId, info)
        else:
            children = self._getChildren(taxonId)
            if children is None:
                self._collectStrains(taxonId, taxonId, info)
            else:
                for id in children:
                    self._collectSpecies(id, info)
                    if info.isSufficient():
                        break

    def _collectStrains(self, specTaxonId, currentTaxonId, info):
        """
            @param specTaxonId: species taxonId (or lower if species not defined)
            @param currentTaxonId: current taxon id is at species rank or lower
            @type info: _Info
        """
        if info.isSpecSufficient(specTaxonId) or info.isSufficient():
            return

        if currentTaxonId in self._taxonIdSet:  # there is a reference
            size = self._taxonIdToSize[currentTaxonId]
            # print currentTaxonId, size  # debug
            info.update(specTaxonId, size)

        if info.isSpecSufficient(specTaxonId):
            return

        children = self._getChildren(currentTaxonId)
        if children is not None:
            for id in children:
                self._collectStrains(specTaxonId, id, info)

    def getNonBacterialNonArchaeal(self):
        """ Returns a list of non Bacterial or non Archaeal ncbi taxon ids that are in the reference. """
        returnList = []
        for taxonId in self._taxonIdSet:
            parentSet = self._taxonomy.getParentsNcbidSet(taxonId)
            if (2 not in parentSet) and (2157 not in parentSet):
                returnList.append(taxonId)
        return 'count:', len(returnList), 'list of ncbi taxon ids: ', returnList

    def close(self):
        self._taxonomy.close()


def testBacteriaArchaea():
    refDir = '/Users/ivan/Documents/nobackup/sequences'
    databaseFilePath = '/Users/ivan/Documents/work/binning/taxonomy/20121122/ncbitax_sqlite.db'
    r = RefSequences(refDir, databaseFilePath)
    print r.getNonBacterialNonArchaeal()
    r.close()
    #('count:', 7, 'list: ', [117575, 1051631, 396359, 683735, 977801, 1126885, 658056])


def testRefSequences():
    refDir = '/Users/ivan/Documents/nobackup/sequences'
    databaseFilePath = '/Users/ivan/Documents/work/binning/taxonomy/20121122/ncbitax_sqlite.db'
    taxonId = 119066  # 2130  # 186803
    minTotalCount = 3  # 43
    minBpPerSpeciesCount = 280000
    r = RefSequences(refDir, databaseFilePath)
    print r. isRefSufficient(taxonId, minTotalCount, minBpPerSpeciesCount)
    r.close()

    # next to check: 658089
    # set([237576, 658088, 796942, 33039, 712976, 796945, 936594, 33043, 45851, 658081, 658082, 658083, 658085, 658086, 658087, 652712, 658089, 1322, 29360, 177972, 712414, 658655, 39486, 168384, 742723, 40520, 552395, 652706, 665937, 33038, 166486, 410072, 43995, 665950, 665951, 796943, 796944, 53443, 360807, 88431, 105841, 301301, 43305])
    # True


if __name__ == "__main__":
    # testBacteriaArchaea()
    testRefSequences()


class DBData():
    """
        @deprecated: use class RefSequences instead
    """

    def __init__(self, ncbiProcessDir, databaseFile):
        self._databaseFile = databaseFile
        # buffer genome/WGS info
        count = 0
        self._ncbiBuffer = set()  # set of ncbids for which there is a genome/WGS in the database
        for filePath in glob.glob(os.path.join(os.path.normpath(ncbiProcessDir),'*.[0-9]*.f[an][sa]')):
            count += 1
            ncbid = int(re.sub(os.path.join(r'.*[^0-9]([0-9]+).[0-9]+.f[an][sa]$'),r'\1' ,filePath))
            self._ncbiBuffer.add(ncbid)
            #print ncbid
        print('%s ref. sequences found in directory: %s' % (count, os.path.normpath(ncbiProcessDir)))

    def getGenomeWgsCount(self, ncbid, threshold):
        """
            Gets the number of genomes/wgs available in the directory.

            @param threshold: it returns max. threshold genomes/wgs (after this threshold is reached, it returns)
            @param dir: directory that contains genomes/wgs in the form: "ncbid.[0-9]*.f[an][sa]"

            @return: the number of genomes/wgs from different species that are subclades of the input ncbid
        """
        try:
            conn = sqlite3.connect(os.path.normpath(self._databaseFile))
            cursor = conn.cursor()

            speciesIdsList = []
            self._collectSpecies(speciesIdsList, cursor, ncbid, threshold)
            return len(speciesIdsList)
        except Exception:
            print "Failed to create connection to a database:", self._databaseFile
            raise
        finally:
            cursor.close()
            conn.close()

    def _genomeExists(self, listOfNcbids):
        """
            Gets the first ncbid of species/subspecies from the input list that is contained in the directory.

            @param dir: directory that contain genomes/wgs files

            @return None or the first ncbid for which there is a genome or a draft genome in the directory
        """
        #for ncbid in listOfNcbids:
        #    if glob.glob(os.path.join(os.path.normpath(dir), str(str(ncbid) + '.[0-9]*.f[an][sa]'))):
        #        return ncbid
        for ncbid in listOfNcbids:
            if ncbid in self._ncbiBuffer:
                return ncbid
        return None

    def _collectSpecies(self, speciesIds, cursor, root, threshold):
        """
            @param speciesIds: output list of ncbids (species or subspecies) for which there are genomes/wgs in
                the directory
        """
        if len(speciesIds) >= threshold:
            return
        cursor.execute('SELECT node_rank FROM taxon T WHERE T.ncbi_taxon_id=?',(root,))
        result = cursor.fetchall()
        if len(result) == 0:
            return
        assert len(result) == 1
        if result[0][0] == 'species':
            list = []
            list.append(root)
            self._collectSubSpecies(list, cursor, root)
            ncbid = self._genomeExists(list)  # now check if at least one genome exists from the list
            if ncbid is not None:
                speciesIds.append(ncbid)
        else:
            cursor.execute('SELECT ncbi_taxon_id FROM taxon T WHERE T.parent_taxon_id=?',(root,))
            result = cursor.fetchall()
            for item in result:
                self._collectSpecies(speciesIds, cursor, int(item[0]), threshold)


    def _collectSubSpecies(self, list, cursor, root):
        """
            Get ncbids of all subspecies of the root
            @list: output list with all ncbids of the subspecies
        """
        cursor.execute('SELECT ncbi_taxon_id FROM taxon T WHERE T.parent_taxon_id=?',(root,))
        result = cursor.fetchall()
        if len(result) > 0:
            for item in result:
                id = int(item[0])
                list.append(id)
                self._collectSubSpecies(list, cursor, id)


def test(ncbid):
    config = Config(open(os.path.normpath('D://A_Phylo//A_Metagenomic//pPPS//workspace//pPPS//config01.cfg')), 'pPPS')
    databaseFile = os.path.normpath(config.get('databaseFile'))
    ncbiProcessDir = os.path.normpath('D://A_Phylo//A_Metagenomic//pPPS//workspace//pPPS//wdir02//ncbiProcDir')
    dbData = DBData(ncbiProcessDir, databaseFile)
    threshold = 3
    print dbData.getGenomeWgsCount(ncbid, threshold)

    #config = Config(open(os.path.normpath('//AM//metagenomic//work//projects//pPPS//tests//TW//TW01//config.cfg')), 'pPPS')

    #threshold = 3
    #dir = 'D://A_Phylo//A_Metagenomic//pPPS//workspace//pPPS//genomes'
    #dir = '//AM//metagenomic//work//projects//pPPS//tests//TW//TW01//ncbiProcDir'

    #databaseFile = os.path.normpath(config.get('databaseFile'))
    #taxonomicRanks = config.get('taxonomicRanks').split(',')

    #count = getGenomeWgsCount(ncbid, threshold, dir, databaseFile, taxonomicRanks)

    #print count, 'genomes/wgs for ncbid:', ncbid


#if __name__ == "__main__":
  #test(122)
  #haveData(126)
  #haveData(84999) #Coriobacteriales
  #haveData(171549) #Bacteroidales
  #haveData(815) #Bacteriodaceae
  #haveData(171551) #Porphyromonadaceae
  #haveData(171552) #Prevotellaceae
  #haveData(171550) #Rikenellaceae
  ###test(976) #Bacteroidetes
  #haveData(200666) #Sphingobacteriales
  #haveData(768503) #Cytophagia
  #haveData(117743) #Flavobacteria
  #haveData(475963) #Caldilineales
  #haveData(292625) #Anaerolineae
  #haveData(200795) #Chloroflexi
  #haveData(204431) #Fibrobacteraceae (59374, 834)
  #haveData(186803) #Lachnospiraceae
  #haveData(541000) #Ruminococcaceae
  #haveData(186802) #Clostridiales
  #haveData(31979) #Clostridiaceae
  #haveData(186806) #Eubacteriaceae
  #haveData(186807) #Peptococcaceae
  #haveData(31977) #Veillonellaceae
  #haveData(186801) #Clostridia
  #haveData(128827) #Erysipelotrichaceae
  #haveData(1239) #Firmicutes
  #haveData(91061) #Bacilli
  #haveData(255528) #Victivallaceae (340101)
  #haveData(126) #Planctomycetaceae
  #haveData(481) #Neisseriaceae
  #haveData(213121) #Desulfobulbaceae (577650, 177439, 589865)
  #haveData(213421) #Desulfuromonaceae
  #haveData(69541) #Desulfuromonadales
  #haveData(72294) #Campylobacteraceae
  #haveData(1224) #Proteobacteria
  #haveData(28211) #Alphaproteobacteria
  #haveData(1236) #Gammaproteobacteria
  #haveData(137) #Spirochaetaceae
  #haveData(186333) #Anaeroplasmataceae
  #haveData(186332) #Anaeroplasmatales
  #haveData(31969) #Mollicutes

  #haveData(278082) #Victivallales
  #haveData(256845) #Lentisphaerae