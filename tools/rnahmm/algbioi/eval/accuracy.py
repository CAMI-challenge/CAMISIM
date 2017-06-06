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
import argparse

from algbioi.com import csv
from algbioi.com import fasta
from algbioi.com import taxonomy_ncbi


class _TaxonomyWrapperA():
    """
        Wraps the functionality of the database.
    """

    def __init__(self, databaseFile):
        self._taxonomy = taxonomy_ncbi.TaxonomyNcbi(databaseFile)
        self._rankToId = {}
        self._ncbidToRankId = {}
        self._predAtRankId = {}  # rankId -> ncbid -> ncbid at given rank
        self._noDefAtRankId = {}  # rankId -> set of ncbids for which the ncbid at given rank is not defined
        self._ncbidToNcbidParent = {}  # ncbid -> parent ncbid

        id = 0
        for rank in taxonomy_ncbi.TAXONOMIC_RANKS:
            self._rankToId[rank] = id
            self._predAtRankId[id] = {}
            self._noDefAtRankId[id] = set()
            id += 1

    def _getRankId(self, ncbid):
        """
            Gets a rankId given an ncbi taxon id
            @rtype: int
        """
        rankId = self._ncbidToRankId.get(ncbid, None)
        if rankId is not None:
            return rankId
        else:
            rank = self._taxonomy.getRank(ncbid)
            if rank is None:
                return None
            else:
                rankId = self._rankToId.get(rank, None)
                self._ncbidToRankId[ncbid] = rankId
                return rankId

    def _getParent(self, ncbid):
        """
            Gets direct parent ncbi taxon id.
        """
        parent = self._ncbidToNcbidParent.get(ncbid, None)
        if parent is None:
            parent = self._taxonomy.getParentNcbid(ncbid)
            self._ncbidToNcbidParent[ncbid] = parent
        return parent

    def getPredDictAtRank(self, seqToNcbid, rank):
        """
            Gets predictions at the given rank as a dictionary.

            @param seqToNcbid: contain mapping, sequence name -> ncbi taxon id
            @type seqToNcbid: dict
            @param rank: the resulting dictionary will contain predictions a this rank
            @type rank: str
            @return: mapping, sequence name -> ncbi taxon id at given rank
            @rtype: dict
        """
        rankId = self._rankToId[rank]
        retDict = {}
        predAtRankBuff = self._predAtRankId[rankId]
        noDefAtRankBuff = self._noDefAtRankId[rankId]

        for seq, ncbid in seqToNcbid.iteritems():

            pred = predAtRankBuff.get(ncbid, None)
            if pred is not None:
                retDict[seq] = pred  # we already know the ncbid at given rank
                continue

            if ncbid in noDefAtRankBuff:
                continue  # the ncbid is not defined at this rank (we already know)

            ncbidRankId = self._getRankId(ncbid)
            #if ncbidRankId is None:
            #    noDefAtRankBuff.add(ncbid)
            #    continue  # we have just found out that the ncbid is not defined at this rank

            if ncbidRankId == rankId:  # the ncbid is defined already at the right rank
                predAtRankBuff[ncbid] = ncbid
                retDict[seq] = ncbid
                continue

            if (ncbidRankId is None) or (ncbidRankId > rankId):  # the right ncbid may be defined at a higher rank
                current = self._getParent(ncbid)
                while current is not None:
                    currentRankId = self._getRankId(current)
                    if currentRankId == rankId:
                        retDict[seq] = current
                        predAtRankBuff[ncbid] = current
                        break
                    current = self._getParent(current)
                if current is None:
                    noDefAtRankBuff.add(ncbid)  # we have just found out that the ncbid is not defined at this rank

        return retDict

    def close(self):
        self._taxonomy.close()


class Accuracy():
    """
        Implements computation of the "precision" and "recall" according to different definitions.
    """

    def __init__(self, seqIdToBp, seqIdToPred, seqIdToTruePred, taxonomy, correctLabelThreshold=None):
        """
            Initializes the accuracy object.
            @param seqIdToBp: dictionary or a fasta file
            @param seqIdToPred: dictionary or a prediction file
            @param seqIdToTruePred: dictionary or a true prediction file
            @param taxonomy: database file in the sqlite3 format, or taxonomy object retrieved from not closed Accuracy
        """
        if isinstance(seqIdToBp, dict):
            self._seqToBp = seqIdToBp
        else:
            assert os.path.isfile(seqIdToBp)
            self._seqToBp = fasta.getSequenceToBpDict(seqIdToBp)

        if isinstance(seqIdToPred, dict):
            self._seqToPred = seqIdToPred
        else:
            assert os.path.isfile(seqIdToPred)
            self._seqToPred = csv.predToDict(seqIdToPred)

        if isinstance(seqIdToTruePred, dict):
            self._seqToTrue = seqIdToTruePred
        else:
            assert os.path.isfile(seqIdToTruePred)
            self._seqToTrue = csv.predToDict(seqIdToTruePred)

        if isinstance(taxonomy, _TaxonomyWrapperA):
            self._taxonomy = taxonomy
        else:
            assert os.path.isfile(taxonomy)
            self._taxonomy = _TaxonomyWrapperA(taxonomy)

        # correct the predictions self._seqToPred
        if correctLabelThreshold is not None:
            self._seqToPred = self._correctPredictions(
                self._seqToBp, self._seqToPred, self._seqToTrue, self._taxonomy, correctLabelThreshold)


    def _correctPredictions(self, seqIdToBp, seqIdToPred, seqIdToTruePred, taxonomy, correctLabelThreshold):
        """

        """
        newPred = {}

        ranks = taxonomy_ncbi.TAXONOMIC_RANKS[1:]
        ranks.reverse()
        for rank in ranks:

            # get true clades at given rank
            seqIdToLabelRank = taxonomy.getPredDictAtRank(seqIdToTruePred, rank)

            # get pred clades at given rank
            seqIdToPredRank = taxonomy.getPredDictAtRank(seqIdToPred, rank)

            # map: true taxonId -> seqId
            labelToSeqIdList = {}
            for seqId, taxonId in seqIdToLabelRank.iteritems():
                if taxonId in labelToSeqIdList:
                    labelToSeqIdList[taxonId].append(seqId)
                else:
                    labelToSeqIdList[taxonId] = [seqId]

            for taxonId, seqIdList in labelToSeqIdList.iteritems():

                idToBp = {}
                sumBp = 0
                for seqId in seqIdList:
                    id = seqIdToPredRank.get(seqId, None)
                    if id is None:
                        continue
                    bp = seqIdToBp[seqId]
                    if id in idToBp:
                        idToBp[id] += bp
                    else:
                        idToBp[id] = bp
                    sumBp += bp

                entryList = []
                for id, bp in idToBp.iteritems():
                    entryList.append((id, bp))
                if len(entryList) == 0:
                    continue
                entryList.sort(key=lambda x: x[1], reverse=True)

                id, bp = entryList[0]
                if id == taxonId:
                    continue
                percent = float(bp) / float(sumBp)
                if percent >= correctLabelThreshold:
                    for seqId in seqIdList:
                        if seqIdToPred.get(seqId, None) == id:
                            newPred[seqId] = taxonId

        for seqId, taxonId in seqIdToPred.iteritems():
            if seqId not in newPred:
                newPred[seqId] = taxonId

        return newPred


    def getAccuracy(self, rank, minFracClade=None, minFracPred=None, asBp=True, weightAccordingBinSize=True):
        """
            Precision (specificity) and Recall (sensitivity) according to PhyloPythiaS and PhyloPythia papers.

            The number of classes correspond to the number of classes in the true reference and param "minFracClades".

            @param rank: on which taxonomic rank the predictions should be considered
            @param minFracClade: a clade is considered only if the dataset (true labels) contain at least this
                          fraction of sequences that belong to the clade
            @param minFracPred: a clade is considered only if the corresponding predicted bins contain at least this
                         fraction of the overall sequences (None ~ this criteria is not considered and only
                         true "reference" bins are used for the comparison).
            @param asBp: count it according to the sequence lengths
            @param weightAccordingBinSize: weight individual bins according to their bin size

            @return: [precision, recall, classPrecisionNum, classRecallNum]
        """
        predAtRankDict = self._taxonomy.getPredDictAtRank(self._seqToPred, rank)
        trueAtRankDict = self._taxonomy.getPredDictAtRank(self._seqToTrue, rank)
        tp = {}  # class label -> count of sequences correctly assigned to clade i
        t = {}  # class label -> true count of sequences of clade i
        p = {}  # class label -> count of sequences assigned to clade i
        tpOther = 0  # count of sequences correctly unassigned
        tOther = 0  # true count of sequences that are unassigned at given rank

        # iterate over all sequences
        for seq, seqLen in self._seqToBp.iteritems():
            # bp
            if asBp:
                bp = seqLen
            else:
                bp = 1

            # true
            i = None
            if seq in trueAtRankDict:
                i = trueAtRankDict[seq]
                if i not in t:
                    t[i] = bp
                else:
                    t[i] += bp
            else:
                tOther += bp
                if seq not in predAtRankDict:
                    tpOther += bp

            # pred
            j = None
            if seq in predAtRankDict:
                j = predAtRankDict[seq]
                if j not in p:
                    p[j] = bp
                else:
                    p[j] += bp

            # match
            if i == j and i is not None:
                if i not in tp:
                    tp[i] = bp
                else:
                    tp[i] += bp

        classesP = p.keys()  # classes for precision
        classesR = t.keys()  # classes for recall

        # filter out least abundant TRUE clades
        if minFracClade is not None:
            sumT = tOther  # true bin containing all sequences undefined at this rank
            for i in classesR:
                sumT += t[i]
            rmList = []
            for i in classesR:
                if (sumT == 0) or (float(t[i]) / float(sumT) < minFracClade):
                    rmList.append(i)
            for i in rmList:
                classesR.remove(i)
            if (sumT == 0) or (float(tOther) / float(sumT) < minFracClade):
                tOther = 0

        # filter out least abundant PREDICTED clades
        if minFracPred is not None:
            sumT = 0
            for i in classesP:
                sumT += p[i]
            rmList = []
            for i in classesP:
                if (sumT == 0) or (float(p[i]) / float(sumT) < minFracPred):
                    rmList.append(i)
            for i in rmList:
                classesP.remove(i)

        # zero missing entries
        for i in classesR:
            if i not in tp:
                tp[i] = 0
            if i not in p:
                p[i] = 0
        for i in classesP:
            if i not in tp:
                tp[i] = 0
            if i not in t:
                t[i] = 0

        wp = {}  # weights for precision
        wr = {}  # weights for recall
        if weightAccordingBinSize:
            # compute weights of individual bins that correspond to the number of bp/sequences
            # assigned to individual bins
            sumP = 0.0
            sumR = 0.0
            for i in classesP:
                sumP += p[i]
            for i in classesR:
                sumR += t[i]
            sumR += tOther

            for i in classesP:
                wp[i] = float(p[i]) / sumP
            for i in classesR:
                wr[i] = float(t[i]) / sumR
            if tOther > 0:
                wrOther = float(tOther) / sumR
        else:
            # all bins are equally important
            for i in classesP:
                wp[i] = 1.0 / float(len(classesP))

            for i in classesR:
                if tOther > 0:
                    w = 1.0 / float(len(classesR) + 1)
                    wr[i] = w
                    wrOther = w
                else:
                    wr[i] = 1.0 / float(len(classesR))
        if len(classesR) == 0 and tOther > 0:
            wrOther = 1.0

        # precision
        precision = 0.0
        for i in classesP:
            if p[i] > 0:
                precision += (float(tp[i]) / float(p[i])) * wp[i]

        # recall
        recall = 0.0
        classesRCount = len(classesR)
        for i in classesR:
            recall += (float(tp[i]) / float(t[i])) * wr[i]
        if tOther > 0:
            recall += (float(tpOther) / float(tOther)) * wrOther
            classesRCount += 1
            #
        return [precision, recall, len(classesP), classesRCount]

    def getAccuracyPrint(self, ranks, minFracClade, minFracPred, overview=True, asBp=True, weightAccordingBinSize=True):
        """
            Gets the precision and recall values printed as a string

            @param ranks: compute the precision and recall at these ranks

            @rtype: str
        """
        buff = '# precision, recall, #classes precision, #classes recall, seq. count/bp, weighted bins\n'
        for rank in ranks:
            if overview:  # overview
                buff += str(rank + ',--,--,--,----------,----------\n')
                buff += self.getAccuracyPrintEntry(rank, minFracClade, minFracPred, False, False)  # asBp, weighted
                buff += self.getAccuracyPrintEntry(rank, minFracClade, minFracPred, True, False)
                buff += self.getAccuracyPrintEntry(rank, minFracClade, minFracPred, False, True)
                buff += self.getAccuracyPrintEntry(rank, minFracClade, minFracPred, True, True)
            else:  # custom
                buff += self.getAccuracyPrintEntry(rank, minFracClade, minFracPred,
                                                   asBp=asBp, weightAccordingBinSize=weightAccordingBinSize)
        return buff

    def getAccuracyPrintEntry(self, rank, minFracClade, minFracPred, asBp=True, weightAccordingBinSize=True):
        p, r, cp, cr = self.getAccuracy(rank, minFracClade, minFracPred, asBp, weightAccordingBinSize)
        if asBp:
            c = 'bp'
        else:
            c = 'count'
        if weightAccordingBinSize:
            w = 'weighted'
        else:
            w = 'not weighted'
        return str('%s, %s, %s, %s, "%s", "%s"\n' % (round(p * 100.0, 1), round(r * 100.0, 1), cp, cr, c, w))

    def getTaxonomy(self):
        return self._taxonomy

    def close(self, closeTaxonomy=True):
        if closeTaxonomy:
            self._taxonomy.close()


def _main():
    parser = argparse.ArgumentParser(
        description='Computes precision and recall measures according to different definitions.', epilog='')

    parser.add_argument('-f', '--fasta', nargs=1, type=file, required=True, help='Fasta file.', metavar='contigs.fna',
                        dest='f')

    parser.add_argument('-p', '--predictions', nargs=1, type=file, required=True,
                        help='Tab separated prediction file (first column sequence name, last column predicted ncbid).',
                        metavar='pred.csv', dest='p')

    parser.add_argument('-t', '--true-assignments', nargs=1, type=file, required=True,
                        help='Tab separated true assignments file (first column sequence name, '
                             'last column predicted ncbid.', metavar='true_assignments.csv', dest='t')

    parser.add_argument('-d', '--database', nargs=1, type=file, required=True,
                        help='Database file containing the NCBI taxonomy in the sqlite3 format.',
                        metavar='ncbitax_sqlite.db', dest='d')

    parser.add_argument('-r', '--ranks', nargs=1, help='Compute the measures only for these ranks (given as comma '
                                                       'separated strings) Default ~ consider all ranks.',
                        metavar='order,family,genus', dest='r')

    parser.add_argument('-c', '--min-frac-clade', nargs=1,
                        help='A clade is considered in the computation of "Recall" only if the reference (the true '
                             'assignments) contain at least this fraction of sequences that belong to the '
                             'corresponding clade at the corresponding rank. (e.g. value 0.01 means that all clades '
                             'that are considered in the computations of recall represent at least 1%% of the overall '
                             'dataset) Default ~ 0.01', metavar='0.01', dest='c')

    parser.add_argument('-b', '--min-frac-bin', nargs=1,
                        help='In the computation of "Precision". A clade is considered only if the corresponding '
                             'predicted bins contain at least this fraction of the overall predicted sequences at the '
                             'corresponding rank (Default ~ 0.01)', metavar='0.01', dest='b')

    parser.add_argument('-s', '--consider-seq-len', action='store_true',
                        help='Compute the measures based on the sequence lengths (in bp). '
                             '(Default ~ based on sequence counts)', dest='s')

    parser.add_argument('-w', '--weight-bins', action='store_true',
                        help='The measures are computed using weighted averages over bin sizes. Size of true bins is '
                             'used to compute "Recall". Size of predicted bins is used to compute "Precision". '
                             '(Default ~ not weighted)', dest='w')

    parser.add_argument('-o', '--overview', action='store_true',
                        help='Compute the measures according to several default settings. '
                             'You can still set the (-c) and (-b) options.', dest='o')
    args = parser.parse_args()

    if args.r:
        ranks = str(args.r[0].name).strip("'").strip('"').split(',')
    else:
        ranks = taxonomy_ncbi.TAXONOMIC_RANKS[1:]

    if args.c:
        minFracClade = float(args.c[0])
    else:
        minFracClade = 0.01

    if args.b:
        minFracPred = float(args.b[0])
    else:
        minFracPred = 0.01

    acc = Accuracy(args.f[0].name, args.p[0].name, args.t[0].name, args.d[0].name)

    print(acc.getAccuracyPrint(ranks, minFracClade, minFracPred,
                               overview=bool(args.o), asBp=bool(args.s), weightAccordingBinSize=bool(args.w)))
    acc.close()


def _test():
    # -f /Users/ivan/Documents/work/binning/data/simMC/AMGN_AMD.Arachne.contigs.fna
    # -p /Users/ivan/Documents/work/binning/tests/simMC/AMD05/output/AMGN_AMD.Arachne.contigs.fna.pOUT
    # -t /Users/ivan/Documents/work/binning/data/simMC/AMD.Arachne.genus
    # -d /Users/ivan/Documents/work/binning/taxonomy/ncbi_taxonomy_20110629/ncbitax_sqlite.db
    fastaFilePath = '/Users/ivan/Documents/work/binning/data/simMC/AMGN_AMD.Arachne.contigs.fna'
    # just marker genes
    predFilePath = '/Users/ivan/Documents/work/binning/tests/simMC/AMD05/output/AMGN_AMD.Arachne.contigs.fna.pOUT'
    # with taxator
    #predFilePath = '/Users/ivan/Documents/work/binning/tests/simMC/AMD06/output/AMGN_AMD.Arachne.contigs.fna.pOUT'
    trueFilePath = '/Users/ivan/Documents/work/binning/data/simMC/AMD.Arachne.genus'
    databaseFile = '/Users/ivan/Documents/work/binning/taxonomy/ncbi_taxonomy_20110629/ncbitax_sqlite.db'
    ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    acc = Accuracy(fastaFilePath, predFilePath, trueFilePath, databaseFile)
    print(acc.getAccuracyPrint(ranks, minFracClade=0.01, minFracPred=0.01, overview=True))
    acc.close()


if __name__ == "__main__":
    _main()
    #_test()