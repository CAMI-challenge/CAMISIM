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


    Computes the confusion (comparison) matrices for different taxonomic ranks.
    A confusion matrix presents comparison of two prediction methods,
    usually taxonomic predictions of a particular method and true (or reference) assignments.
    Rows represent true (reference) assignments, columns represent predictions by a particular method.
"""

import os
import argparse
from algbioi.com import taxonomy_ncbi
from algbioi.com import csv
from algbioi.com import fasta as fas


class _TaxonomyWrapCM():
    def __init__(self, databaseFile):
        """
            Taxonomy wrapper that buffers frequently used operations for this module.
            @param databaseFile: database in the sqlite3 format
        """
        self._taxonomy = taxonomy_ncbi.TaxonomyNcbi(databaseFile)
        # buffers
        self._rankToRankId = {}
        self._rankIdToRank = {}
        self._taxonIdToParentTaxonId = {}
        self._taxonIdToRankId = {}
        self._taxonIdToScientificName = {}
        # map: rank <-> rankId
        rankId = 0
        for rank in taxonomy_ncbi.TAXONOMIC_RANKS:
            self._rankToRankId[rank] = rankId
            self._rankIdToRank[rankId] = rank
            rankId += 1

    def getParent(self, taxonId):
        """
            @return: parent taxonId
            @rtype int
        """
        parentId = self._taxonIdToParentTaxonId.get(taxonId, None)
        if parentId is None:
            parentId = self._taxonomy.getParentNcbid(taxonId)
            self._taxonIdToParentTaxonId[taxonId] = parentId
        return parentId

    def getRankIdOfTaxonId(self, taxonId):
        """
            @return: rankId of the taxonId (ids correspond to method getRankId)
            @rtype: int
        """
        rankId = self._taxonIdToRankId.get(taxonId, None)
        if rankId is None:
            rank = self._taxonomy.getRank(taxonId)
            rankId = self._rankToRankId.get(rank, None)
            self._taxonIdToRankId[taxonId] = rankId
        return rankId

    def getRankId(self, rank):
        """
            @type rank: str
            @return: id of the given taxonomic rank
            @rtype: int
        """
        return self._rankToRankId.get(rank, None)

    def getSortedScientificNames(self, taxonIdSet):
        """
            Gets a list of sorted scientific names that correspond to the input taxonIds.
            If there is a problem with the retrieval of the scientific names, taxonIds are used instead.

            @param taxonIdSet: set of taxonIds
            @type taxonIdSet: set of int
            @return: list of sorted scientific names or taxonIds, mapping scientific name to taxonId
            @rtype: (list of str, dict)
        """
        names = []
        nameToTaxonId = {}
        for id in taxonIdSet:
            name = self._getScientificName(id)
            if name is None:
                name = str(id)  # use taxonId if the scientific name cannot be found
            names.append(name)
            nameToTaxonId[name] = id

        if len(taxonIdSet) != len(nameToTaxonId):  # the scientific names were not unique, taxonIds will be used instead
            names = []
            nameToTaxonId = {}
            for id in taxonIdSet:
                names.append(str(id))
                nameToTaxonId[str(id)] = id

        names.sort()
        return names, nameToTaxonId

    def _getScientificName(self, taxonId):
        """
            @return: scientific name
            @rtype: str
        """
        name = self._taxonIdToScientificName.get(taxonId, None)
        if name is None:
            name = self._taxonomy.getScientificName(taxonId)
            self._taxonIdToScientificName[taxonId] = name
        return name

    def close(self):
        self._taxonomy.close()


class ConfusionMatrix():
    def __init__(self, seqNameToBp, seqNameToPred, seqNameToRefPred, taxonomy, ranksList=None):
        """
            Initializes the main class that computes the confusion matrices.

            @param seqNameToBp: contains mapping, sequence name to bp (as int); or a fasta file
                @type seqNameToBp: dict; or a fasta file
            @param seqNameToPred: contains mapping, sequence name to taxonId; or a tab separated prediction file
                @type seqNameToPred: dict; or a tab separated file, first column ~ sequence name, last column taxonId
            @param seqNameToRefPred: contains mapping, sequence name to taxon Id; or a tab separated reference file
                @type seqNameToRefPred: dict; or a tab separated file, first column ~ sequence name, last column taxonId
            @param ranksList: list of ranks for which the confusion matrices will be computed (None ~ all default ranks)
                @type ranksList: list of str
            @param taxonomy: database file in the sqlite3 format; or taxonomy returned by function "getTaxonomy"
        """
        # Check input options and read in the data (if appropriate)
        self._initFailed = False  # replace this with exceptions!
        if isinstance(seqNameToBp, dict):
            self._seqNameToBp = seqNameToBp
        elif isinstance(seqNameToBp, str) and os.path.isfile(seqNameToBp):
            self._seqNameToBp = fas.getSequenceToBpDict(seqNameToBp)
        else:
            print("Can't get sequence info from:", seqNameToBp)
            self._initFailed = True
            return
        if isinstance(seqNameToPred, dict):
            self._seqNameToPred = seqNameToPred
        elif isinstance(seqNameToPred, str) and os.path.isfile(seqNameToPred):
            self._seqNameToPred = csv.predToDict(seqNameToPred)
        else:
            print("Can't get prediction info from:", seqNameToPred)
            self._initFailed = True
            return
        if isinstance(seqNameToRefPred, dict):
            self._seqNameToRefPred = seqNameToRefPred
        elif isinstance(seqNameToRefPred, str) and os.path.isfile(seqNameToRefPred):
            self._seqNameToRefPred = csv.predToDict(seqNameToRefPred)
        else:
            print("Can't get reference prediction info from:", seqNameToRefPred)
            self._initFailed = True
            return
        if isinstance(taxonomy, str) and os.path.isfile(taxonomy):
            self._taxonomy = _TaxonomyWrapCM(taxonomy)
        elif isinstance(taxonomy, _TaxonomyWrapCM):
            self._taxonomy = taxonomy
        else:
            print("Can't use taxonomy: ", taxonomy)
        if ranksList is None:
            ranksList = taxonomy_ncbi.TAXONOMIC_RANKS[1:]  # default ranks
        else:
            allowedRanksSet = set(taxonomy_ncbi.TAXONOMIC_RANKS[1:])  # custom ranks
            for rank in ranksList:
                if rank not in allowedRanksSet:
                    print('Rank: "' + str(rank) + '" is not allowed!')
                    self._initFailed = True
                    return
        rankIdsList = []  # rankIds that will be considered
        for rank in ranksList:
            rankIdsList.append(self._taxonomy.getRankId(rank))
        self._allowedRankIdsSet = set(rankIdsList)

        # get predictions at different taxonomic ranks
        # rankId -> (seqId -> taxonIdAtRank)
        self._rankIdToPredMap = {}
        self._rankIdToRefMap = {}
        for rankId in rankIdsList:
            self._rankIdToPredMap[rankId] = {}
            self._rankIdToRefMap[rankId] = {}

        # get predictions at given ranks
        for seqId, taxonId in self._seqNameToPred.iteritems():
            while (taxonId is not None) and (taxonId != 1):
                rankId = self._taxonomy.getRankIdOfTaxonId(taxonId)
                if rankId in self._allowedRankIdsSet:
                    self._rankIdToPredMap[rankId][seqId] = taxonId
                taxonId = self._taxonomy.getParent(taxonId)

        # get reference predictions at given ranks
        for seqId, taxonId in self._seqNameToRefPred.iteritems():
            while (taxonId is not None) and (taxonId != 1):
                rankId = self._taxonomy.getRankIdOfTaxonId(taxonId)
                if rankId in self._allowedRankIdsSet:
                    self._rankIdToRefMap[rankId][seqId] = taxonId
                taxonId = self._taxonomy.getParent(taxonId)

    def generateConfusionMatrix(self, rank, prefixOutputPath):
        """
            Generates confusion matrix at given rank.
            The object must have been initialized considering this rank.

            @param prefixOutputPath: prefix of the output file path
        """
        if self._initFailed:
            return
        rankId = self._taxonomy.getRankId(rank)
        if rankId not in self._allowedRankIdsSet:
            print("Can't consider rank: " + rank)
            return
        if not os.path.isdir(os.path.dirname(prefixOutputPath)):
            print("Output prefix is wrong, the corresponding directory doesn't exist: " +
                  os.path.dirname(prefixOutputPath))
            return

        # entries of the confusion matrix
        tableCountMap = {}  # (taxonId_ref, taxonId_pred) -> count
        tableBpMap = {}  # (taxonId_ref, taxonId_pred) -> bp

        # predictions (and reference) at the given rank
        seqNameToPred = self._rankIdToPredMap[rankId]
        seqNameToRef = self._rankIdToRefMap[rankId]
        predTaxonIdSet = set()
        refTaxonIdSet = set()

        # fill in entries of the confusion matrix
        for seqId, bp in self._seqNameToBp.iteritems():
            predId = seqNameToPred.get(seqId, None)
            refId = seqNameToRef.get(seqId, None)
            if predId is not None:
                predTaxonIdSet.add(predId)  # stores it only if it's predicted at this rank
            if refId is not None:
                refTaxonIdSet.add(refId)
            key = (refId, predId)
            if key not in tableCountMap:
                tableCountMap[key] = 1
                tableBpMap[key] = bp
            else:
                tableCountMap[key] += 1
                tableBpMap[key] += bp

        # get taxonIds contained in prediction and reference prediction, common for both, unique for pred. and ref.
        commonTaxonIdSet = predTaxonIdSet.intersection(refTaxonIdSet)
        uniquePredIdSet = predTaxonIdSet.difference(commonTaxonIdSet)
        uniqueRefIdSet = refTaxonIdSet.difference(commonTaxonIdSet)

        # get taxonIds contained in predictions and reference predictions as lists of scientific names
        commonNames, commonMap = self._taxonomy.getSortedScientificNames(commonTaxonIdSet)
        uniquePredNames, uniquePredMap = self._taxonomy.getSortedScientificNames(uniquePredIdSet)
        uniqueRefNames, uniqueRefMap = self._taxonomy.getSortedScientificNames(uniqueRefIdSet)

        # headers
        predHeader = commonNames + uniquePredNames + ['unassigned']  # predictions
        refHeader = commonNames + uniqueRefNames + ['unassigned']  # reference
        predHeaderTaxonIds = []
        refHeaderTaxonIds = []
        for name in commonNames:
            id = commonMap[name]
            predHeaderTaxonIds.append(id)
            refHeaderTaxonIds.append(id)
        for name in uniquePredNames:
            predHeaderTaxonIds.append(uniquePredMap[name])
        for name in uniqueRefNames:
            refHeaderTaxonIds.append(uniqueRefMap[name])
        predHeaderTaxonIds.append(None)  # predicted as unassigned
        refHeaderTaxonIds.append(None)  # unassigned in reference

        # count matches
        matchCount = 0
        matchBp = 0
        for taxonId in commonTaxonIdSet:
            count = tableCountMap.get((taxonId, taxonId), None)
            if count is not None:
                bp = tableBpMap.get((taxonId, taxonId), None)
                assert bp is not None
                matchCount += count
                matchBp += bp

        # count mismatches
        mismatchCount = 0
        mismatchBp = 0
        for predTaxonId in predHeaderTaxonIds[:-1]:
            for refTaxonId in refHeaderTaxonIds[:-1]:
                if predTaxonId == refTaxonId:
                    continue
                assert (predTaxonId is not None) and (refTaxonId is not None)
                count = tableCountMap.get((refTaxonId, predTaxonId), None)
                if count is not None:
                    bp = tableBpMap.get((refTaxonId, predTaxonId), None)
                    assert bp is not None
                    mismatchCount += count
                    mismatchBp += bp

        # count pred total, ref total
        predTotalCount = 0
        predTotalBp = 0
        refTotalCount = 0
        refTotalBp = 0
        for predTaxonId in predHeaderTaxonIds:
            for refTaxonId in refHeaderTaxonIds:
                count = tableCountMap.get((refTaxonId, predTaxonId), None)
                if count is None:
                    continue
                bp = tableBpMap.get((refTaxonId, predTaxonId), None)
                assert bp is not None
                if predTaxonId is not None:
                    predTotalCount += count
                    predTotalBp += bp
                if refTaxonId is not None:
                    refTotalCount += count
                    refTotalBp += bp

        # total
        totalCount = 0
        totalBp = 0
        for bp in self._seqNameToBp.values():
            totalCount += 1
            totalBp += bp

        # write the confusion matrix to a file
        out = csv.OutFileBuffer(os.path.normpath(prefixOutputPath + '.' + str(rank) + '_cmp.csv'))

        header = 'ref/pred'
        for e in predHeader:
            header += ', ' + e
        out.writeText(header + '\n')

        for i in range(len(refHeaderTaxonIds)):
            line = refHeader[i]
            refTaxonId = refHeaderTaxonIds[i]
            for j in range(len(predHeaderTaxonIds)):
                predTaxonId = predHeaderTaxonIds[j]
                count = tableCountMap.get((refTaxonId, predTaxonId), None)
                line += ', '
                if count is not None:
                    bp = tableBpMap.get((refTaxonId, predTaxonId), None)
                    assert bp is not None
                    line += str(int(round(float(bp) / 1000.0))) + 'k (' + str(count) + ')'
            out.writeText(line + '\n')

        out.writeText(',\n')
        out.writeText('Matches, ' + str(int(round(float(matchBp) / 1000.0))) + 'k, ' + str(matchCount) + ', ' +
                      self._div(matchBp, matchBp + mismatchBp, 1) + ' %k' + ', ' +
                      self._div(matchCount, matchCount + mismatchCount, 1) + ' %\n')

        out.writeText('Mismatches, ' + str(int(round(float(mismatchBp) / 1000.0))) + 'k, ' + str(mismatchCount) + ', ' +
                      self._div(mismatchBp, matchBp + mismatchBp, 1) + ' %k' + ', ' +
                      self._div(mismatchCount, matchCount + mismatchCount, 1) + ' %\n')

        out.writeText('Pred. assigned, ' + str(int(round(float(predTotalBp) / 1000.0))) + 'k, ' + str(predTotalCount) + ', ' +
                      self._div(predTotalBp, totalBp, 1) + ' %k, ' +
                      self._div(predTotalCount, totalCount, 1) + ' %\n')

        out.writeText('Ref. assigned, ' + str(int(round(float(refTotalBp) / 1000.0))) + 'k, ' + str(refTotalCount) + ', ' +
                      self._div(refTotalBp, totalBp, 1) + ' %k, ' +
                      self._div(refTotalCount, totalCount, 1) + ' %\n')

        out.writeText('Total fasta, ' + str(int(round(float(totalBp) / 1000.0))) + 'k, ' + str(totalCount) + '\n')
        out.close()

    def _div(self, dividend, divisor, roundNDigits):
        if abs(divisor) < 0.000001:
            return 'NaN'
        else:
            return str(round(100.0 * (float(dividend) / float(divisor)), roundNDigits))

    def getTaxonomy(self):
        return self._taxonomy

    def close(self, closeTaxonomy=True):
        if closeTaxonomy:
            self._taxonomy.close()


def _main():
    """ Main method of the script, see the module documentation and argument description. """
    parser = argparse.ArgumentParser(description='Computes the confusion matrix, detailed comparison of two placements '
                                                 '(e.g. true/reference assignments vs. predictions by a particular '
                                                 'method).', epilog='Module description: ' + __doc__)

    parser.add_argument('-f', '--fasta', nargs=1, type=file, required=True, help='Fasta file.', metavar='sequences.fna',
                        dest='f')

    parser.add_argument('-p', '--predictions', nargs=1, type=file, required=True,
                        help='Tab separated prediction file (first column sequence name, last column predicted '
                             'ncbi taxon id).', metavar='predictions.csv', dest='p')

    parser.add_argument('-t', '--true-ref-assignments', nargs=1, type=file, required=True,
                        help='Tab separated true/reference assignments file (first column sequence name, '
                             'last column predicted ncbi taxon id.', metavar='true_ref_assignments.csv', dest='t')

    parser.add_argument('-d', '--database', nargs=1, type=file, required=True,
                        help='Database file containing the NCBI taxonomy in the sqlite3 format.',
                        metavar='ncbitax_sqlite.db', dest='d')

    parser.add_argument('-o', '--output', nargs=1, required=True,
                        help='Prefix of the output path.', metavar='sample0', dest='o')

    parser.add_argument('-r', '--ranks', nargs=1, help='Compute the measures only for these ranks (given as comma '
                                                       'separated strings) Default ~ consider all ranks.',
                        metavar='order,family,genus', dest='r')

    args = parser.parse_args()

    assert len(args.f) == 1 and len(args.p) == 1 and len(args.t) == 1 and len(args.d) == 1 and len(args.o) == 1

    if args.r and len(args.r) == 1:
        ranks = str(args.r[0]).strip("'").strip('"').split(',')
    else:
        ranks = taxonomy_ncbi.TAXONOMIC_RANKS[1:]

    confusionMatrix = ConfusionMatrix(args.f[0].name, args.p[0].name, args.t[0].name, args.d[0].name, ranks)

    for rank in ranks:
        confusionMatrix.generateConfusionMatrix(rank, args.o[0])

    confusionMatrix.close()


if __name__ == "__main__":
    _main()
