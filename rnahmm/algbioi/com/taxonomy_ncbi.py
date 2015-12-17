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
import sqlite3

TAXONOMIC_RANKS = ['root','superkingdom','phylum','class','order','family','genus','species']


class TaxonomyNcbi():
    """
        Represents an interface to the sqlite3 database in which the NCBI taxonomy is stored.
        (NOTE that methods and variables starting with "_" are local and shouldn`t be used from the outside)

        @author: Ivan
    """
    def __init__(self, databaseFile, allowedRanks=TAXONOMIC_RANKS, considerNoRank=False):
        """
            @param databaseFile: usually file named "ncbitax_sqlite.db"
            @param allowedRanks: taxonomic ranks that will be considered (where 'root' is the root of the taxonomy)
            @param considerNoRank: consider ranks 'no rank' if true
        """
        self._allowedRanks = set(allowedRanks)
        if considerNoRank:
            self._allowedRanks.add('no rank')
        try:
            self.conn = sqlite3.connect(os.path.normpath(databaseFile))
            self.cursor = self.conn.cursor()
        except Exception:
            sys.stderr.write(str('TaxonomyNcbi: Failed to create connection to database: ' + databaseFile))
            raise

    def getScientificName(self, ncbid, checkRank=False):
        """
            @return: scientific name or None
            @rtype: str
        """
        if ncbid == -1:
            ncbid = 1
            sys.stderr.write('ncbid -1 converted to 1')

        if checkRank and (not self.isRankNcbidAllowed(ncbid)):
            return None

        self.cursor.execute(str('SELECT TN.name FROM taxon_name TN, taxon T WHERE T.ncbi_taxon_id=?' +
                                ' AND T.taxon_id = TN.taxon_id AND TN.name_class="scientific name"'),(ncbid,))
        result = self.cursor.fetchall()
        if len(result) == 1:
            return result[0][0]
        else:
            sys.stderr.write(str('TaxonomyNcbi: Cannot find name for ncbi: ' + str(ncbid)))
            return None

    def getNcbid(self, scientificName, checkRank = False):
        """
            @return: ncbid or None
            @rtype: int
        """
        self.cursor.execute(str('SELECT T.ncbi_taxon_id FROM taxon_name TN, taxon T ' +
                                'WHERE TN.name_class="scientific name" AND TN.name=? AND TN.taxon_id=T.taxon_id'),
                                (scientificName,))
        result = self.cursor.fetchall()
        if len(result) == 1:
            ncbid = int(result[0][0])
            if checkRank and (not self.isRankNcbidAllowed(ncbid)):
                return None
            return ncbid
        else:
            if len(result) == 0:
                sys.stderr.write(str('TaxonomyNcbi: Cannot find scientific name "' + scientificName + '" in the database.\n'))
            else:
                sys.stderr.write(str('TaxonomyNcbi: scientific name "' + scientificName + '" is ambiguous!\n'))
            return None

    def getNcbid2(self, name, checkRank = False):
        """
            @param name: doesn`t have to be a scientific name
            @return: ncbid or None
            @rtype: int
        """
        self.cursor.execute(str('SELECT T.ncbi_taxon_id FROM taxon_name TN, taxon T ' +
                                'WHERE TN.name=? AND TN.taxon_id=T.taxon_id'),
                                (name,))
        result = self.cursor.fetchall()
        if len(result) == 1:
            ncbid = int(result[0][0])
            if checkRank and (not self.isRankNcbidAllowed(ncbid)):
                return None
            return ncbid
        else:
            if len(result) == 0:
                sys.stderr.write(str('TaxonomyNcbi: Cannot find scientific name "' + name + '" in the database.\n'))
            else:
                sys.stderr.write(str('TaxonomyNcbi: scientific name "' + name + '" is ambiguous!\n'))
            return None

    def getChildrenNcbids(self, ncbid):  # SELECT T1.ncbi_taxon_id from taxon T1 where T1.parent_taxon_id=818;
        self.cursor.execute(str('SELECT T1.ncbi_taxon_id from taxon T1 where T1.parent_taxon_id=?'),(ncbid,))
        result = self.cursor.fetchall()
        if len(result) == 0:
            return None
        else:
            resultList = []
            for i in result:
                resultList.append(i[0])
            return resultList

    def getParentNcbid(self, ncbid):
        """
            @return: ncbid or None
            @rtype: int
        """
        if ncbid == 1:
            return None
        taxonId = self._getTaxonId(ncbid)

        while True:
            if ncbid == 1 or ncbid == None or taxonId == None: #the root of the taxonomy reached
                return None
            ncbid = self._getParentNcbid(taxonId)

            taxonId = self._getTaxonId(ncbid)
            rank = self._getRank(taxonId)
            if (rank in self._allowedRanks) or (ncbid == 1 and 'root' in self._allowedRanks):
                return ncbid

    def getParentsNcbidSet(self, ncbid):
        """
            @return: set of parent ncbi taxon ids.
            @rtype: set
        """
        s = set()
        currentId = ncbid
        while True:
            currentId = self.getParentNcbid(currentId)
            if currentId is None:
                break
            else:
                s.add(currentId)
        return s

    def getRank(self, ncbid, checkRank=False):
        """
            @return: rank or None
            @rtype: str
        """
        if checkRank and (not self.isRankNcbidAllowed(ncbid)):
            return None
        return self._getRank(self._getTaxonId(ncbid))

    def isRankNcbidAllowed(self, ncbid):
        """
            @rtype: bool
        """
        taxonId = self._getTaxonId(ncbid)
        rank = self._getRank(taxonId)
        if rank in self._allowedRanks:
            return True
        else:
            return False

    def isRankAllowed(self, rank):
        """
            @rtype: bool
        """
        if rank in self._allowedRanks:
            return True
        else:
            return False

    def exists(self, ncbid):
        if self._getTaxonId(ncbid) is None:
            return False
        else:
            return True

    def close(self):
        """
            Close the database after you stop using it.
        """
        self.cursor.close()
        self.conn.close()

    def _getTaxonId(self, ncbid):
        if ncbid is None:
            return None
        if ncbid == -1:
            ncbid = 1
            #sys.stderr.write('ncbid=(-1) converted to ncbid=(1)\n')
        self.cursor.execute('SELECT taxon_id FROM taxon T WHERE T.ncbi_taxon_id=?',(ncbid,))
        result = self.cursor.fetchall()
        if len(result) != 1:
            #sys.stderr.write('TaxonomyNcbi: Cannot find taxon_id for ncbi:' + str(ncbid) + ' result:' + str(result) + ' \n')
            return None
        return int(result[0][0])

    def _getParentNcbid(self, taxonId):
        if taxonId is None:
            return None
        self.cursor.execute('SELECT parent_taxon_id FROM taxon T WHERE T.taxon_id=?', (taxonId,))
        result = self.cursor.fetchall()
        if len(result) != 1:
            sys.stderr.write(str('TaxonomyNcbi: Cannot find parent for taxon_id' + str(taxonId)))
            return None
        return int(result[0][0])

    def _getRank(self, taxonId):
        if taxonId is None:
            return None
        self.cursor.execute('SELECT node_rank FROM taxon T WHERE T.taxon_id=?', (taxonId,))
        result = self.cursor.fetchall()
        if len(result) != 1:
            sys.stderr.write(str('TaxonomyNcbi: Cannot find rank for taxon_id: ' + str(taxonId)))
            return None
        return str(result[0][0])


def test():
    databaseFile = "/Users/ivan/Documents/work/binning/taxonomy/ncbi_taxonomy_20110629/ncbitax_sqlite.db"
    taxonomy = TaxonomyNcbi(databaseFile)

    print 'Scientific name for ncbid 286730 (Alkaliflexus imshenetskii) is:', taxonomy.getScientificName(286730)

    print 'Ncbid of Lachnospiraceae (186303) is:', str(taxonomy.getNcbid('Lachnospiraceae'))

    print 'Get Parent of 167965 which is 74152:', str(taxonomy.getParentNcbid(167965))

    print 'Get rank of Lachnosipraceae (family):', str(taxonomy.getRank(186803))

    taxonomy.close()


def test2():
    databaseFile = "/Users/ivan/Documents/work/binning/taxonomy/ncbi_taxonomy_20110629/ncbitax_sqlite.db"
    taxonomy = TaxonomyNcbi(databaseFile)
    ncbid = 1
    parent = taxonomy.getParentNcbid(ncbid)
    print ncbid, taxonomy.getScientificName(ncbid)
    print parent, taxonomy.getScientificName(parent)
    print parent, taxonomy.getRank(parent)

    taxonomy.close()

def test3():
    databaseFile = "/Users/ivan/Documents/work/binning/taxonomy/ncbi_taxonomy_20110629/ncbitax_sqlite.db"
    taxonomy = TaxonomyNcbi(databaseFile)
    print taxonomy.getParentsNcbidSet(186803)

def test4():
    taxonomy = TaxonomyNcbi('/Users/ivan/Documents/work/binning/taxonomy/20121122/ncbitax_sqlite.db')
    print taxonomy.childrenNcbids(83763)


if __name__ == "__main__":
    pass
    #test()
    #test2()
    #test3()
    #test4()