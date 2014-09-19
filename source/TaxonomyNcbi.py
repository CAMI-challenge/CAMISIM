#!/usr/bin/env python

import os
import sys
import re
import sqlite3
from sets import Set


#Represents an interface to the sqlite3 database in which the NCBI database is stored.
#(NOTE that methods and variables starting with "_" are local and shouldn`t be used from the outside)
class TaxonomyNcbi():

	#Constructor
	#@param databaseFile: usually file named "ncbitax_sqlite.db"
	#@param allowedRanks: taxonomic ranks that will be considered (where 'root' is the root of the taxonomy)
	#@param considerNoRank: consider ranks 'no rank' if true
	def __init__(self, databaseFile, allowedRanks=['root','superkingdom','phylum','class','order','family','genus','species'],
				 considerNoRank=False):

		self._allowedRanks = Set(allowedRanks)
		if considerNoRank:
			self._allowedRanks.add('no rank')
		try:
			self.conn = sqlite3.connect(os.path.normpath(databaseFile))
			self.cursor = self.conn.cursor()
		except Exception:
			sys.stderr.write(str('TaxonomyNcbi: Failed to create connection to database: ' + databaseFile))
			raise


	#@return scientific name or None
	def getScientificName(self, ncbid, checkRank = False):

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

	def get_ncbi_id_of_ranks(self, ncbi_id, ranks=['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']):
		result = [''] * len(ranks)
		if ncbi_id == 1:
			return result
		taxon_id = self._getTaxonId(ncbi_id)
		rank = self._getRank(taxon_id)
		if rank in ranks:
			result[ranks.index(rank)] = ncbi_id
		while True:
			if ncbi_id == 1 or ncbi_id is None or taxon_id is None: #the root of the taxonomy reached
				return result
			ncbi_id = self._getParentNcbid(taxon_id)
			taxon_id = self._getTaxonId(ncbi_id)
			rank = self._getRank(taxon_id)
			if rank in ranks:
				result[ranks.index(rank)] = ncbi_id
		return result

	#@return ncbid or None
	def getNcbid(self, scientificName, checkRank = False, verbose=True):
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
			if len(result) == 0 and verbose:
				sys.stderr.write(str('TaxonomyNcbi: Cannot find scientific name "' + scientificName + '" in the database.\n'))
			else:
				sys.stderr.write(str('TaxonomyNcbi: scientific name "' + scientificName + '" is ambiguous!\n'))
			return None

	#@return ncbid, rank or None
	def get_ncbi_id_and_rank(self, scientificName, checkRank = False, verbose=True):
		self.cursor.execute(str('SELECT T.ncbi_taxon_id, T.node_rank FROM taxon_name TN, taxon T ' +
								'WHERE TN.name_class="scientific name" AND TN.name=? AND TN.taxon_id=T.taxon_id'),
								(scientificName,))
		result = self.cursor.fetchall()
		if len(result) == 1:
			ncbid = int(result[0][0])
			rank = str(result[0][1])
			if checkRank and rank not in self._allowedRanks:
				return None
			return ncbid, rank
		else:
			if len(result) == 0 and verbose:
				sys.stderr.write(str('TaxonomyNcbi: Cannot find scientific name "' + scientificName + '" in the database.\n'))
			else:
				if len(result) == 2:
					ncbi_1 = int(result[0][0])
					rank_1 = str(result[0][1])
					ncbi_2 = int(result[1][0])
					rank_2 = str(result[1][1])
					if ncbi_1 == self.getParentNcbid(ncbi_2) or ncbi_2 == self.getParentNcbid(ncbi_1):
						ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
						if ranks.index(rank_1) > ranks.index(rank_2):
							return ncbi_1, rank_1
						return ncbi_2, rank_2
				sys.stderr.write(str('TaxonomyNcbi: scientific name "' + scientificName + '" is ambiguous ({})!\n'.format(len(result))))
			return None

	#@return ncbid or None
	def getParentNcbid(self, ncbid):
		if ncbid == 1:
			return None
		taxonId = self._getTaxonId(ncbid)

		while True:
			if ncbid == 1 or ncbid is None or taxonId is None: #the root of the taxonomy reached
				return None
			ncbid = self._getParentNcbid(taxonId)
			taxonId = self._getTaxonId(ncbid)
			rank = self._getRank(taxonId)
			if rank in self._allowedRanks:
				return ncbid


	#@return rank or None
	def getRank(self, ncbid, checkRank = False):
		if checkRank and (not self.isRankNcbidAllowed(ncbid)):
			return None
		return self._getRank(self._getTaxonId(ncbid))


	#@return: True or False
	def isRankNcbidAllowed(self, ncbid):
		taxonId = self._getTaxonId(ncbid)
		rank = self._getRank(taxonId)
		if rank in self._allowedRanks:
			return True
		else:
			return False


	#@return True or False
	def isRankAllowed(self, rank):
		if rank in self._allowedRanks:
			return True
		else:
			return False


	def get_ncbi_of_rank(self, rank, ncbi):
		#if ncbi <= 2:
		#	return ncbi
		current_rank = self.getRank(ncbi)
		#print current_rank
		while current_rank != rank and current_rank != "root" and current_rank is not None:
			ncbi = self.getParentNcbid(ncbi)
			current_rank = self.getRank(ncbi)
			if ncbi == 1:
				return None
		#if current_rank == "root" or current_rank is None:
		#	return None
		return ncbi


	#close the database after you stop using it
	def close(self):
		self.cursor.close()
		self.conn.close()


	def _getTaxonId(self, ncbid):
		if ncbid == None:
			return None
		self.cursor.execute('SELECT taxon_id FROM taxon T WHERE T.ncbi_taxon_id=?',(ncbid,))
		result = self.cursor.fetchall()
		if len(result) != 1:
			sys.stderr.write('TaxonomyNcbi: Cannot find taxon_id for ncbi:' + str(ncbid) + ' result:' + str(result) + ' \n')
			return None
		return int(result[0][0])


	def _getParentNcbid(self, taxonId):
		if taxonId == None:
			return None
		self.cursor.execute('SELECT parent_taxon_id FROM taxon T WHERE T.taxon_id=?', (taxonId,))
		result = self.cursor.fetchall()
		if len(result) != 1:
			sys.stderr.write(str('TaxonomyNcbi: Cannot find parent for taxon_id' + str(taxonId)))
			return None
		return int(result[0][0])


	def _getRank(self, taxonId):
		if taxonId == None:
			return None
		self.cursor.execute('SELECT node_rank FROM taxon T WHERE T.taxon_id=?', (taxonId,))
		result = self.cursor.fetchall()
		if len(result) != 1:
			sys.stderr.write(str('TaxonomyNcbi: Cannot find rank for taxon_id: ' + str(taxonId)))
			return None
		return str(result[0][0])


def test():
	databaseFile = "D:/A_Phylo/A_Metagenomic/pPPS/workspace/pPPS/ncbi/ncbi_taxonomy_20110629/ncbitax_sqlite.db"
	taxonomy = TaxonomyNcbi(databaseFile)

	print 'Scientific name for ncbid 286730 (Alkaliflexus imshenetskii) is:', taxonomy.getScientificName(286730)

	print 'Ncbid of Lachnospiraceae (186303) is:', str(taxonomy.getNcbid('Lachnospiraceae'))

	print 'Get Parent of 167965 which is 74152:', str(taxonomy.getParentNcbid(167965))

	print 'Get rank of Lachnosipraceae (family):', str(taxonomy.getRank(186803))


	taxonomy.close()



if __name__ == "__main__":
  test()