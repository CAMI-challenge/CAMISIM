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


    Computes the scaffold-contig consistency based on new definitions (2012).

    @author: Ivan
    @version: 26.4.2013 (internal subversion: 17)
"""

import os
import sys
import argparse

from algbioi.com.csv import predToDict
from algbioi.com.csv import getMapping
from algbioi.com.fasta import getSequenceToBpDict
from algbioi.com.taxonomy_ncbi import TaxonomyNcbi


class _TaxonomyWrapper():
    """
        Wraps the taxonomy to buffer (speed up) taxonomy calls.
    """
    def __init__(self, databaseFile):
        self._taxonomy = TaxonomyNcbi(databaseFile)
        self._ncbidToNcbidParent = dict()
        self._closed = False

    def getParent(self, ncbid):
        """
            Gets direct parent ncbi taxon id.
        """
        if ncbid in self._ncbidToNcbidParent:
            return self._ncbidToNcbidParent[ncbid]
        else:
            parent = self._taxonomy.getParentNcbid(ncbid)
            self._ncbidToNcbidParent[ncbid] = parent
            return parent

    def getDist(self, ncbid, ncbidSet):
        """
            Gets the distance between ncbid and the closest clade in the ncbidSet which
            represents a path from the root to a clade.
        """
        current = ncbid
        dist = 0
        while (current not in ncbidSet):
            if current is None:
                sys.stderr.write('Consistency:_TaxonomyWrapper:getDist: current is "None" '
                                 + str(ncbid) + ' ' + str(ncbidSet) + '\n')
                break
            current = self.getParent(current)
            dist += 1
        return dist

    def getDistantParent(self, ncbid, intDist):
        """
            Gets parent of "ncbid" that is in distance "intDist".
        """
        parentNcbid = ncbid
        for i in range(int(intDist)):
            parentNcbid = self.getParent(parentNcbid)
        return parentNcbid

    def getDistTowardsRoot(self, ncbid, parentNcbid):
        """
            Gets the distance from ncbid to parentNcbid that is its parent and lies on the same path to the root.
        """
        current = ncbid
        dist = 0
        while (current != parentNcbid) and (current != 1):
            if current is None:
                sys.stderr.write('Consistency:_TaxonomyWrapper:getDistTowardsRoot: current is "None" '
                                 + str(ncbid) + ' ' + str(parentNcbid) + '\n')
                break
            current = self.getParent(current)
            dist += 1
        return float(dist)

    def getScientificName(self, ncbid):
        return self._taxonomy.getScientificName(ncbid)

    def close(self):
        """ To free resources. """
        self._taxonomy.close()
        self._closed = True

    def isClosed(self):
        return self._closed


class ScScaffold():
    """
        Represents one scaffold.
    """
    def __init__(self, name, ncbid, pathSet, weightedDistToPath, contigsNameList, contigNameToNcbid,
                 ncbidToBp, ncbidToWeight, ncbidToDist, ncbidToLeafDist):
        self._name = name
        self._ncbid = ncbid
        self._pathSet = pathSet
        self._weightedDistToPath = weightedDistToPath
        self._contigsNameList = contigsNameList
        self._contigNameToNcbid = contigNameToNcbid
        self._ncbidToBp = ncbidToBp
        self._ncbidToWeight = ncbidToWeight
        self._ncbidToDist = ncbidToDist
        self._ncbidToLeafDist = ncbidToLeafDist

    def getName(self):
        return self._name

    def getNcbid(self):
        """
            Returns ncbid according to which the consistency of this Scaffold was computed.
        """
        return self._ncbid

    def getPathSet(self):
        """
            Get the path from the root to the ncbid according to which this scaffold was assigned,
            the path is represented as a set.
        """
        return self._pathSet

    def getToLeafDist(self, ncbid):
        """
            Get the distance from the argument ncbid to the ncbid to which the scaffold was assigned.
        """
        if ncbid in self._ncbidToLeafDist:
            return self._ncbidToLeafDist[ncbid]
        else:
            return None

    def getToPathDist(self, ncbid):
        """
            Get the distance from the argument ncbid to the path according to which this scaffold was assigned.
        """
        if ncbid in self._ncbidToDist:
            return self._ncbidToDist[ncbid]
        else:
            return None

    def getContigsNameList(self):
        """
            Gets all contigs' names of this scaffold as a list.
        """
        return self._contigsNameList

    def getCollectiveLength(self):
        """
            Get the sum of the lengths (in bp) of all constituent contigs.
        """
        length = 0.0
        for ncbid in self._ncbidToBp:
            length += self._ncbidToBp[ncbid]
        return length

    def getConsistencyTotal(self, asCount=False):
        """
            Consistency: def 1.
        """
        consistentCount = 0.0
        for contigName in self._contigsNameList:
            ncbid = self._contigNameToNcbid[contigName]
            if ncbid in self._pathSet:
                consistentCount += 1
        if asCount:
            return float(consistentCount)
        else:
            return float(consistentCount) / len(self._contigsNameList)

    def getConsistencyTotalBp(self, asBpCount=False):
        """
            Consistency: def 2.
        """
        consistentBp = 0.0
        sumBp = 0.0
        for ncbid in self._ncbidToBp:
            bp = self._ncbidToBp[ncbid]
            if ncbid in self._pathSet:
                consistentBp += bp
            sumBp += bp
        if asBpCount:
            return float(consistentBp)
        else:
            return float(consistentBp) / sumBp

    def getConsistencyAvgDist(self, asTotalCount=False):
        """
            Consistency: def3.
        """
        sumDist = 0.0
        for contigName in self._contigsNameList:
            ncbid = self._contigNameToNcbid[contigName]
            sumDist += self._ncbidToDist[ncbid]
        if asTotalCount:
            return float(sumDist)
        else:
            return float(sumDist) / len(self._contigsNameList)

    def getConsistencyWeightedAvgDist(self):
        """
            Consistency: def4.
        """
        return self._weightedDistToPath

    def getConsistencyAvgDistLeaf(self, asTotalCount=False):
        """
            Consistency: def5.
        """
        sumDist = 0.0
        for contigName in self._contigsNameList:
            ncbid = self._contigNameToNcbid[contigName]
            sumDist += self._ncbidToLeafDist[ncbid]
        if asTotalCount:
            return float(sumDist)
        else:
            return float(sumDist) / len(self._contigsNameList)

    def getConsistencyAvgWeightedDistLeaf(self):
        """
            Consistency: def6.
        """
        sumDist = 0.0
        for ncbid in self._ncbidToWeight:
            sumDist += float(self._ncbidToWeight[ncbid]) * self._ncbidToLeafDist[ncbid]
        return float(sumDist)


class Consistency():
    def __init__(self, contigNameToBp, contigNameToNcbid, scaffToContigList, taxonomy,
                 minScaffContigCount=None, minScaffBpLen=None, cladesSet=None, considerContigWithNoScaff=True,
                 ignoreScaffPredToRoot=True):
        """
            Initializes the main Consistency class.

            @param contigNameToBp: dictionary that maps contig names to bp (int);
                or a fasta file that contain contigs
            @param contigNameToNcbid: dictionary that maps contig names to ncbids (int);
                or a prediction file - first column contig name, last column ncbid
            @param scaffToContigList: dictionary that maps scaffold names to list of contig names;
                or a file - first column scaffold name, second column contig name
            @param minScaffContigCount: consider only scaffolds that contain at least this number of contigs
            @param minScaffBpLen: consider only scaffolds with at least this collective length (in bp)
            @param cladesSet: consider only scaffolds that contain at least one contig from this set
            @param considerContigWithNoScaff: consider also contigs that are not assigned to scaffolds
                (as artificial scaffolds)
            @param ignoreScaffPredToRoot: ignore scaffolds that are assigned based on the root (uninformative)
        """
        # check input options
        assert minScaffContigCount is None or isinstance(minScaffContigCount, int)
        assert minScaffBpLen is None or isinstance(minScaffBpLen, int)
        assert cladesSet is None or isinstance(cladesSet, set)
        assert isinstance(considerContigWithNoScaff, bool)
        assert isinstance(ignoreScaffPredToRoot, bool)

        if isinstance(contigNameToBp, dict):
            self._contigNameToBp = contigNameToBp
        elif isinstance(contigNameToBp, str) and os.path.isfile(contigNameToBp):
            self._contigNameToBp = getSequenceToBpDict(contigNameToBp)
        else:
            print("Can't get contig info from: ", contigNameToBp)
            return
        if isinstance(contigNameToNcbid, dict):
            self._contigToPred = contigNameToNcbid
        elif isinstance(contigNameToNcbid, str) and os.path.isfile(contigNameToNcbid):
            self._contigToPred = predToDict(contigNameToNcbid)
        else:
            print("Can't get prediction info from: ", contigNameToNcbid)
            return
        if isinstance(scaffToContigList, dict):
            self._scaffToContigsList = scaffToContigList
        elif isinstance(scaffToContigList, str) and os.path.isfile(scaffToContigList):
            self._scaffToContigsList = getMapping(scaffToContigList, 0, 1, '\t')
        else:
            print("Can't get scaffold config mapping from: ", scaffToContigList)
            return

        if isinstance(taxonomy, _TaxonomyWrapper) and (not taxonomy.isClosed()):
            self._taxonomy = taxonomy
        elif isinstance(taxonomy, str) and os.path.isfile(taxonomy):
            self._taxonomy = _TaxonomyWrapper(taxonomy)
        else:
            print("Can't use taxonomy:", taxonomy)
            return

        # check the consistency of the data!

        # if a contig that is defined in the mapping doesn't exist (in the fasta file) we remove it
        for scaff, contigsList in self._scaffToContigsList.iteritems():
            removeList = []
            for contig in contigsList:
                if contig not in self._contigNameToBp:
                    removeList.append(contig)

            for contig in removeList:
                contigsList.remove(contig)

        # if a contig was predicted but there is no scaffold assigned to it then this
        # contig is assigned to an "artificial scaffold"
        if considerContigWithNoScaff:
            scaffContigSet = set()
            for s, l in self._scaffToContigsList.iteritems():
                for c in l:
                    scaffContigSet.add(c)
            aloneContigSet = set()
            for c in self._contigToPred:
                if c not in scaffContigSet:
                    aloneContigSet.add(c)

            for c in aloneContigSet:
                scaffName = str('scaffold_' + c)  # make up a scaffold name
                assert scaffName not in self._scaffToContigsList, 'The names of contigs are ambiguous!'
                self._scaffToContigsList[scaffName] = [c]

        # filter out scaffolds according to the input constrains
        self._scaffolds = dict()
        for scaffName, contigsList in self._scaffToContigsList.iteritems():
            if minScaffContigCount is not None:
                if len(contigsList) < minScaffContigCount:
                    continue

            if minScaffBpLen is not None:
                sum = 0
                for contig in contigsList:
                    sum += self._contigNameToBp[contig]
                if sum < minScaffBpLen:
                    continue

            if cladesSet is not None:
                passScaff = False
                for contig in contigsList:
                    if (contig in self._contigToPred) and (self._contigToPred[contig] in cladesSet):
                        passScaff = True
                        break
                if not passScaff:
                    continue

            # process the scaffold, but if everything in the scaffold was assigned to the root, then ignore it!
            s = self._processScaffold(scaffName)
            if not ((s.getNcbid() == 1) and ignoreScaffPredToRoot):
                self._scaffolds[scaffName] = s

    def _getPred(self, contigName):
        """
            Gets prediction for contig "contigName" or 1 if it wasn't assigned.
        """
        if contigName in self._contigToPred:
            return self._contigToPred[contigName]
        else:
            return 1

    def _processScaffold(self, scaffName):
        """
            Gets a "Scaffold" - an object that contain pre-processed scaffold contig consistency information.
            @rtype: ScScaffold
        """
        allNcbidSet = set()
        for contigName in self._scaffToContigsList[scaffName]:
            allNcbidSet.add(self._getPred(contigName))
        parentNcbidSet = set([1])
        leafNcbidSet = set()
        for ncbid in allNcbidSet:
            if ncbid == 1:
                continue
            current = self._taxonomy.getParent(ncbid)
            while current not in parentNcbidSet:
                if current is None:
                    sys.stderr.write('Consistency:Consistency:_processScaffold: current is "None" ' + str(ncbid) + ' '
                                     + str(current) + ' ' + str(parentNcbidSet) + '\n')
                    break
                parentNcbidSet.add(current)
                current = self._taxonomy.getParent(current)
        for ncbid in allNcbidSet:
            if ncbid not in parentNcbidSet:
                leafNcbidSet.add(ncbid)

        ncbidToBp = dict()
        sumBp = 0
        for contigName in self._scaffToContigsList[scaffName]:
            ncbid = self._getPred(contigName)
            bp = self._contigNameToBp[contigName]
            sumBp += int(bp)
            if ncbid not in ncbidToBp:
                ncbidToBp[ncbid] = int(bp)
            else:
                ncbidToBp[ncbid] += int(bp)

        ncbidToWeight = dict()
        for ncbid in allNcbidSet:
            ncbidToWeight[ncbid] = float(ncbidToBp[ncbid])/float(sumBp)

        # for all paths defined by leaf ncbids compute the weighted distance to all other ncbids
        minDistW = sys.float_info.max
        # minDist = sys.maxint
        minDistNcbid = 1
        minPath = set([1])
        for ncbid in leafNcbidSet:
            path = set([1])
            current = ncbid
            while current != 1:
                if current is None:
                    sys.stderr.write('Consistency:Consistency:_processScaffold: '
                                     'current is "None" (while current != 1) ' + str(current) + '\n')
                    break
                path.add(current)
                current = self._taxonomy.getParent(current)

            distW = 0.0
            dist = 0.0
            for ncbidA in allNcbidSet:
                d = float(self._taxonomy.getDist(ncbidA, path))
                distW += d * ncbidToWeight[ncbidA]
                dist += d
            if distW < minDistW:
                minDistW = distW
                # minDist = dist
                minDistNcbid = ncbid
                minPath = path

        # if everything is assigned to the root, then the distance is 0
        if len(leafNcbidSet) == 0:
            minDistW = 0.0
            # minDist = 0.0
            assert minDistNcbid == 1

        # for each ncbid compute the distance to the path
        ncbidToDist = dict()
        for ncbid in allNcbidSet:
            ncbidToDist[ncbid] = float(self._taxonomy.getDist(ncbid, minPath))

        contigNameToNcbid = dict()
        for contigName in self._scaffToContigsList[scaffName]:
            contigNameToNcbid[contigName] = self._getPred(contigName)

        # for each ncbid compute the distance to the leaf (minDistNcbid)
        ncbidToLeafDist = dict()
        for ncbid in allNcbidSet:
            d = int(ncbidToDist[ncbid])
            lcaNcbid = self._taxonomy.getDistantParent(ncbid, d)
            ncbidToLeafDist[ncbid] = float(d + self._taxonomy.getDistTowardsRoot(minDistNcbid, lcaNcbid))

        return ScScaffold(scaffName, minDistNcbid, minPath, minDistW, self._scaffToContigsList[scaffName],
                          contigNameToNcbid, ncbidToBp, ncbidToWeight, ncbidToDist, ncbidToLeafDist)

    def getScaffoldsDict(self):
        """
            Gets all scaffolds (as a dict) that were filtered out according to the parameters.
        """
        return self._scaffolds

    def getScaffoldsPrint(self):
        """
            Gets a list of scaffolds to be printed out.
            @rtype: str
        """
        buff = ''
        scaffList = []
        for scaffName in self._scaffolds:
            scaffList.append(scaffName)
        scaffList.sort()

        for scaffName in scaffList:
            scaff = self._scaffolds[scaffName]
            contigCount = len(scaff.getContigsNameList())
            pathSet = scaff.getPathSet()
            scaffNcbid = scaff.getNcbid()
            buff += str(scaff.getName() + ', ' + str(scaffNcbid) + ', ' +
                        str(round(float(scaff.getCollectiveLength()) / 1000.0, 3)) + 'kbp,  (' +
                        str(int(scaff.getConsistencyTotal(asCount=True))) + '/' + str(contigCount) + ')')
            if abs(scaff.getConsistencyTotal() - 1) > 0.0001:
                buff += str(', ' + str(round(scaff.getConsistencyTotal() * 100, 0)) + '%, ' +
                            str(round(scaff.getConsistencyTotalBp() * 100, 0)) + '%bp')
            if scaff.getConsistencyAvgDist() > 0.0001:
                buff += ',  pathD:, ' + str(round(scaff.getConsistencyAvgDist(), 2)) + ', ' + \
                        str(round(scaff.getConsistencyWeightedAvgDist(), 2)) + 'w'
            if scaff.getConsistencyAvgDistLeaf() > 0.0001:
                buff += str(',  leafD:,' + str(round(scaff.getConsistencyAvgDistLeaf(), 2)) + ', ' +
                            str(round(scaff.getConsistencyAvgWeightedDistLeaf(), 2)) + 'w')

            buff += ',  ('
            i = 0
            contigList = scaff.getContigsNameList()
            contigList.sort()
            for contig in contigList:
                contigNcbid = self._getPred(contig)
                bp = 0
                if contig in self._contigNameToBp:
                    bp = self._contigNameToBp[contig]
                buff += str(contig + ' ' + str(int(bp)) + 'bp ' + str(contigNcbid))
                if contigNcbid == scaffNcbid:
                    buff += '*'
                elif contigNcbid in pathSet:
                    buff += str('+' + str(int(self._taxonomy.getDistTowardsRoot(scaffNcbid, contigNcbid))))
                else:
                    buff += str('-' + str(int(scaff.getToLeafDist(contigNcbid))) + '-' +
                                str(int(scaff.getToPathDist(contigNcbid))))
                if (i + 1) == contigCount:
                    buff += ')'
                else:
                    buff += '; '
                i += 1
            buff += '\n'

        return buff

    def getGroupedScaffoldsPrint(self):
        """
            Gets scaffolds grouped according to their ncbid, to be printed out.
        """
        ncbidToScaffList = dict()  # ncbid -> list of scaffolds
        for scaffName, scaff in self._scaffolds.iteritems():
            scaffNcbid = scaff.getNcbid()
            if scaffNcbid not in ncbidToScaffList:
                ncbidToScaffList[scaffNcbid] = [scaff]
            else:
                ncbidToScaffList[scaffNcbid].append(scaff)

        # ncbids
        ncbidSet = set(ncbidToScaffList.keys())

        # scientific name -> ncbid list (there can be more than one ncbids for the same scientific name)
        scientificNameToNcbidList = dict()
        nameList = []
        for ncbid in ncbidSet:
            name = self._taxonomy.getScientificName(ncbid)
            if name in scientificNameToNcbidList:
                scientificNameToNcbidList[name].append(ncbid)
            else:
                scientificNameToNcbidList[name] = [ncbid]
                nameList.append(name)
        nameList.sort()

        buff = ''
        nameList.append('Summary')
        scientificNameToNcbidList['Summary'] = [1]

        for name in nameList:
            for ncbid in scientificNameToNcbidList[name]:
                if name != 'Summary':
                    scaffolds = ncbidToScaffList[ncbid]
                else:
                    ncbid = str('all clades (' + str(len(set(ncbidToScaffList.keys()))) + ')')
                    scaffolds = list(self._scaffolds.values())

                if len(scaffolds) == 0:  # there are no scaffolds for this ncbid
                    continue

                totalContigCount = 0.0  # count the number of all contigs in all scaffolds
                totalBpLen = 0.0  # sum up lengths of all scaffolds' lengths
                totalConsistentContigCount = 0.0  # number of contigs that are consistent
                totalConsistentBpLen = 0.0  # number of Bp that are consistent
                totalPathDist = 0.0  # distance from all contigs to the path
                totalPathDistWeighted = 0.0  # weighted distance from all contigs to the path
                totalLeafDist = 0.0  # distance from all contigs to the respective leafs (of the path)
                totalLeafDistWeighted = 0.0  # weighted distance from all contigs to the respective leaf (of the path)

                for scaff in scaffolds:
                    contigCount = float(len(scaff.getContigsNameList()))
                    collectiveLength = float(scaff.getCollectiveLength())
                    totalContigCount += contigCount
                    totalConsistentContigCount += float(scaff.getConsistencyTotal(asCount=True))
                    totalBpLen += collectiveLength
                    totalConsistentBpLen += float(scaff.getConsistencyTotalBp(asBpCount=True))
                    totalPathDist += float(scaff.getConsistencyAvgDist(asTotalCount=True))
                    totalPathDistWeighted += collectiveLength * float(scaff.getConsistencyWeightedAvgDist())
                    totalLeafDist += float(scaff.getConsistencyAvgDistLeaf(asTotalCount=True))
                    totalLeafDistWeighted += collectiveLength * scaff.getConsistencyAvgWeightedDistLeaf()

                buff += str(name + ', (' + str(ncbid) + '), scaffolds: ' + str(len(scaffolds)) + ', contigs: (' +
                            str(int(totalConsistentContigCount)) + '/' + str(int(totalContigCount)) + '), ' +
                            str(round(((totalConsistentContigCount / totalContigCount) * 100.0), 2)) + '%, (' +
                            str(round(totalConsistentBpLen / 1000.0, 1)) + '/' + str(round(totalBpLen / 1000.0, 1)) +
                            ' kb), ' + str(round(((totalConsistentBpLen / totalBpLen) * 100.0), 2)) +
                            '% bp,  pathDist:, ' + str(round(totalPathDist / totalContigCount, 2)) + ', ' +
                            str(round(totalPathDistWeighted / totalBpLen, 2)) + 'w, leafDist:, ' +
                            str(round(totalLeafDist / totalContigCount, 2)) + ', ' +
                            str(round(totalLeafDistWeighted / totalBpLen, 2)) + 'w')
                buff += '\n'
        return buff

    def getTaxonomy(self):
        return self._taxonomy

    def close(self, closeTaxonomy=True):
        """
            Release resources.
        """
        if closeTaxonomy and (not self._taxonomy.isClosed()):
            self._taxonomy.close()


def _main():
    """
        Main function.
    """
    parser = argparse.ArgumentParser(description='Computes the scaffold-contig consistency based on '
                                                 'the "maximum support path".',
                                     epilog=__doc__)

    parser.add_argument('-f', '--fasta', nargs=1, type=file, required=True,
                        help='Fasta file containing contigs.', metavar='contigs.fna',
                        dest='f')

    parser.add_argument('-p', '--predictions', nargs=1, type=file, required=True,
                        help='Tab separated prediction file (first column contig name, '
                             'last column predicted ncbi taxon id.',
                        metavar='pred.csv',
                        dest='p')

    parser.add_argument('-m', '--mapping', nargs=1, type=file, required=True,
                        help='Tab separated scaffold-contig mapping file (first column scaffold name, '
                             'second column contig name.',
                        metavar='scaffoldToContig.csv',
                        dest='m')

    parser.add_argument('-d', '--database', nargs=1, type=file, required=True,
                        help='Database file in the sqlite3 format.', metavar='ncbitax_sqlite.db',
                        dest='d')

    parser.add_argument('-c', '--min-scaff-contig-count', nargs=1,
                        help='Consider only scaffolds that contain at least this number of contigs.', metavar='X',
                        dest='c')

    parser.add_argument('-b', '--min-scaff-bp', nargs=1,
                        help='Consider only scaffolds that has at least this collective length '
                             '(sum of the lengths of all constituent contigs).',
                        metavar='X',
                        dest='b')

    parser.add_argument('-n', '--allowed-ncbids-list', nargs=1,
                        help='Consider only scaffolds that contain at least one contig that was assigned to one '
                             'of this clades (given as comma separated NCBIDs)',
                        metavar='186803,83763,128827',
                        dest='n')

    parser.add_argument('-a', '--consider-artificial-scaff', action='store_true',
                        help='Consider contigs that are not assigned to scaffolds. '
                             '(For each such contig an artificial scaffold is created.)',
                        dest='a')

    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print scaffold-contig consistency info for each scaffold.',
                        dest='v')

    args = parser.parse_args()

    if args.c:
        assert len(args.c) == 1
        minScaffContigCount = int(args.c[0])
    else:
        minScaffContigCount = None

    if args.b:
        assert len(args.b) == 1
        minScaffBpLen = int(args.b[0])
    else:
        minScaffBpLen = None

    if args.n:
        cladesSet = set()
        assert len(args.n) == 1
        for i in str(args.n[0]).split(','):
            cladesSet.add(int(i))
    else:
        cladesSet = None

    assert len(args.f) == 1 and len(args.p) == 1 and len(args.m) == 1 and len(args.d) == 1  # make this nicer

    cons = Consistency(args.f[0].name, args.p[0].name, args.m[0].name, args.d[0].name, minScaffContigCount,
                       minScaffBpLen, cladesSet, args.a)
    if args.v:
        print cons.getScaffoldsPrint()
    print cons.getGroupedScaffoldsPrint()
    cons.close()


def _test1():
    """
        For debugging purposes.
    """
    db = '/Users/ivan/Documents/work/binning/taxonomy/ncbi_taxonomy_20110629/ncbitax_sqlite.db'

    minScaffContigCount = 26
    minScaffBpLen = None
    cladesSet = None  # set([133925])
    considerContigWithNoScaff = False

    #SRM
    fasta = '/Users/ivan/Documents/work/binning/data/Reindeer/contigs.fna'
    pred = '/Users/ivan/Documents/work/binning/tests/SRM1/07/output/contigs.fna.pOUT'
    map = '/Users/ivan/Documents/work/binning/data/Reindeer/scaffolds-contigs.tab'

    #HG
    #fasta = '/Users/ivan/Documents/work/binning/data/HumanGut/working/contigs.fna'
    #pred = '/Users/ivan/Documents/work/binning/tests/HumanGut/02/output/contigs.fna.pOUT'
    #pred = '/Users/ivan/Documents/work/binning/data/HumanGut/working/PPS_contigs_no_minus_one.txt'
    #map = '/Users/ivan/Documents/work/binning/data/HumanGut/working/scafftocontig.txt'

    #TW
    #fasta = '/Users/ivan/Documents/work/binning/data/TW/assembly/contigs_tw.fna'
    #pred = '/Volumes/Macintosh HD/Users/ivan/Documents/work/binning/tests/TW/TW13/output/contigs_tw.fna.pOUT'
    #map = '/Users/ivan/Documents/work/binning/data/TW/scaff_contig.tab'

    cons = Consistency(fasta, pred, map, db, minScaffContigCount, minScaffBpLen, cladesSet, considerContigWithNoScaff)
    print cons.getScaffoldsPrint()
    print cons.getGroupedScaffoldsPrint()
    cons.close()


def _test2():
    db = '/Users/ivan/Documents/work/binning/taxonomy/ncbi_taxonomy_20110629/ncbitax_sqlite.db'
    contigNameToBp = {'contig1': 120, 'contig2': 150, 'contig3': 200}
    contigNameToNcbid = {'contig1': 186803, 'contig2': 186802, 'contig3': 200795}
    scaffToContigList = {'scaff1': ['contig1', 'contig2', 'contig3']}
    cons = Consistency(contigNameToBp, contigNameToNcbid, scaffToContigList, db)
    scaffolds = cons.getScaffoldsDict()  # get all scaffolds as a list
    for scaffoldName, scaffold in scaffolds.iteritems():
        print scaffoldName, scaffold.getConsistencyTotalBp()
    print cons.getScaffoldsPrint()
    taxonomy = cons.getTaxonomy()
    cons.close(closeTaxonomy=False)

    cons = Consistency(contigNameToBp, contigNameToNcbid, scaffToContigList, taxonomy)
    print cons.getScaffoldsPrint()
    print cons.getGroupedScaffoldsPrint()
    cons.close()  # this will also close the taxonomy


if __name__ == "__main__":
    _main()
    #_test1()
    #_test2()