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


    This module enables to mask parts of the reference sequences at different ranks, which is needed for
    the simulated data set evaluation.

    Input:
        (1) A rank that should be masked.
        (2) A list of clades that should be removed at different ranks.
        (3) Directory with reference sequences (as needed for PPS).
        (4) Directory with Silva or Amphora mg database

    Output:
        (1) A list of clades that should be masked from the reference sequences at given ranks
            according to the taxonomy and reference sequences.
        (2) A directory that does not contain masked reference sequences (as symlinks to the original sequences)
        (3) A directory with marker gene databases with masked sequences
"""

import os
import re
import glob
import argparse

from algbioi.com import csv
from algbioi.com import fasta as fas
from algbioi.com import taxonomy_ncbi as tax


_STRAIN = 'strain'
_RANKS = ['phylum', 'class', 'order', 'family', 'genus', 'species', _STRAIN]


class _TaxonomyWrapMD():
    def __init__(self, taxonomy):
        """ A taxonomy wrapper.
            @param taxonomy: file path to the ncbi taxonomy database in the sqlite3
            @type taxonomy: str
        """
        self._taxonomy = tax.TaxonomyNcbi(taxonomy, considerNoRank=True)
        self._existsTaxonIdSet = set()
        self._taxonIdToDirectChildrenSet = {}

    def exists(self, taxonId):
        """ @type taxonId: int
            @rtype: bool
        """
        assert isinstance(taxonId, int)
        if taxonId in self._existsTaxonIdSet:
            return True
        elif self._taxonomy.exists(taxonId):
            self._existsTaxonIdSet.add(taxonId)
            return True
        else:
            return False

    def isLeaf(self, taxonId):
        """ @type taxonId: int
            @rtype: bool
        """
        assert isinstance(taxonId, int)
        if self._taxonomy.getChildrenNcbids(taxonId) is None:
            return True
        else:
            return False

    def getTaxonIdAtRank(self, taxonId, rank):
        """ @type taxonId: int
            @type rank: str
            @rtype: int
        """
        assert isinstance(taxonId, int)
        assert (rank in _RANKS), str('Rank "%s" is not supported!' % rank)
        if rank == _STRAIN:
            assert self.isLeaf(taxonId)
            return taxonId
        current = taxonId
        while current is not None:
            if rank == self._taxonomy.getRank(current):
                return current
            else:
                current = self._taxonomy.getParentNcbid(current)
        return None

    def getAllChildren(self, taxonId):
        """ @type taxonId: int
            @rtype: set of int
        """
        assert isinstance(taxonId, int)
        return self._getAllChildren(taxonId, set())

    def close(self):
        self._taxonomy.close()

    def _getAllChildren(self, taxonId, childrenSet):
        """ @type childrenSet: set of int"""
        children = self._getDirectChildren(taxonId)
        if children is not None:
            childrenSet.update(children)
            for i in children:
                self._getAllChildren(i, childrenSet)
        return childrenSet

    def _getDirectChildren(self, taxonId):
        if taxonId in self._taxonIdToDirectChildrenSet:
            return self._taxonIdToDirectChildrenSet[taxonId]
        else:
            lst = self._taxonomy.getChildrenNcbids(taxonId)
            if lst is None:
                return None
            else:
                resultSet = set(map(int, lst))
                self._taxonIdToDirectChildrenSet[taxonId] = resultSet
                return resultSet


def _refFilePathToTaxonId(refFilePath):
    """
        Extracts a taxonId from a filePath that represents a fasta file (*.fna or *.fas) that contains reference
        sequences for given taxonId as used for PPS.

        @param refFilePath: a path to a file in format PATH/"taxonId"\.[0-9]+\.fna (or fas);
        @return: taxonId
        @rtype: int
    """
    assert refFilePath.endswith(('.fna', '.fas')), str(
        'The fasta files are allowed to end either with .fna or .fas ' + refFilePath)
    return int(os.path.basename(refFilePath).rsplit('.', 2)[0])


def _mgSeqIdToTaxonId(seqId):
    """
        Extracts a taxonId from sequence id used in the Amphora or Silva mg databases (ends with '|ncbid:taxonId")

        @param seqId: sequence id used in mg databases
        @return: taxonId
        @rtype: int
        """
    return int(seqId.rsplit('|', 1)[1].rsplit(':', 1)[1])


def _main():
    """ Main function (script interface), see module description."""
    parser = argparse.ArgumentParser(description=__doc__, epilog='')

    parser.add_argument('-a', '--action', nargs=1, required=True, choices=['cl', 'mr', 'mg'],
                        help='''To determine ONE action that will be performed. (cl) generates a list of leaf level
                        clades that should be masked. (mr) masked reference sequences as sim-links.
                        (mg) mask Amphora or Silva mg db.''', dest='action')

    parser.add_argument('-c', '--input-clades-list', nargs=1, type=file, required=True,
                        help='A list of clades that will be masked from the databases.', metavar='clades.txt',
                        dest='clades')

    parser.add_argument('-d', '--input-dir', action='store', nargs=1, required=True,
                        help='Directory that contains input files (reference sequences or mg reference sequences).',
                        metavar='in_dir', dest='inDir')

    parser.add_argument('-r', '--rank', action='store', nargs=1, required=True, choices=_RANKS,
                        help='Rank to be masked.', metavar='genus', dest='rank')

    parser.add_argument('-o', '--output-dir', action='store', nargs=1, required=True,
                        help='Directory that contains the output files (masked (mg), reference sequences, '
                             'or a list of masked clades)', metavar='out_dir', dest='outDir')

    parser.add_argument('-t', '--taxonomy-file', nargs=1, type=file, required=True,
                        help='NCBI taxonomy database file in the sqlite3 format.', metavar='ncbitax_sqlite.db',
                        dest='taxonomy')

    # parse and read arguments
    args = parser.parse_args()
    action = args.action[0]
    inDir = args.inDir[0]
    outDir = args.outDir[0]
    rank = args.rank[0]
    cladesFilePath = args.clades[0].name
    taxonomyFilePath = args.taxonomy[0].name
    maskDb(action, inDir, outDir, rank, cladesFilePath, taxonomyFilePath, verbose=True)


def maskDb(action, inDir, outDir, rank, clades, taxonomyFilePath, verbose=False):
    """
        Main function (function interface), see module description.

        @param action: one action that will be performed [cl, mr, mg] ~ (generate list, mask seq, mask mg)
        @type action str
        @param inDir: directory containing input files
        @type inDir: str
        @param outDir: directory containing output files
        @type: outDir: str
        @param rank: the data will be excluded at this rank
        @type rank: str
        @param clades: a file containing clades that will be masked (one ncbi taxon id at a line),
            or a set of ncbi taxon ids that will be masked
        @type clades: file or set of int
        @param taxonomyFilePath: taxonomy database file in the sqlite3 format
        @type taxonomyFilePath: str
    """
    # check input parameters
    assert action in ['cl', 'mr', 'mg'], str('Given action is not supported: ' + action)
    if action == 'mr':
        assert os.name == 'posix', 'Symbolic links can be created only on posix systems, action "mr" is not valid!'
    for dir in [inDir, outDir]:
        assert os.path.isdir(dir), str("Directory doesn't exists: " + dir)
    assert rank in _RANKS, str('Not supported rank: ' + rank)
    assert os.path.isfile(taxonomyFilePath), str("Taxonomy database file doesn't exist: " + taxonomyFilePath)
    assert isinstance(clades, set) or (isinstance(clades, str) and os.path.isfile(clades)), str(
        "Parameter 'clades' can be either a file or a set of ncbi taxonIds to be excluded.")

    # maps a rank to a lower rank
    toLowerRank = {}
    for i in range(1, len(_RANKS)):
        toLowerRank[_RANKS[i-1]] = _RANKS[i]

    taxonomy = _TaxonomyWrapMD(taxonomyFilePath)

    # leaf clades to mask
    if isinstance(clades, set):
        inCladesSet = set(map(int, clades))
    else:
        inCladesSet = set(map(int, csv.getColumnAsList(clades)))

    # clades in the reference
    refCladesSet = set()
    if action in ['cl', 'mr']:
        # get the list of all taxon ids that appear in the directory (as PPS reference)
        for fastaFilePath in glob.glob(os.path.join(os.path.normpath(inDir), r'*.f[na][as]')):  # *.fas or *.fna
            refCladesSet.add(_refFilePathToTaxonId(fastaFilePath))  # taxonId.1.fna or taxonId.1.fas
    elif action in ['mg']:
        # get the list of all taxon ids that appear in any file in the input directory as taxonomy ".tax"
        for mapFilePath in glob.glob(os.path.join(os.path.normpath(inDir), r'*.tax')):  # *.tax
            refCladesSet.update(set(map(_mgSeqIdToTaxonId, csv.getColumnAsList(mapFilePath, sep='\t'))))
    else:
        assert False, str('Not supported action: ' + action)

    # checks whether taxonIds are in the taxonomy
    for taxonId in inCladesSet:
        assert taxonomy.exists(taxonId), str(
            'taxonId: %s from clades list is not contained in the taxonomy!' % taxonId)
    for taxonId in refCladesSet:
        assert taxonomy.exists(taxonId), str(
            'taxonId: %s from the reference is not contained in the taxonomy!' % taxonId)

    # checks whether the taxonIds are leafs (doesn't have to be (unless you want to mask at the strain level))
    for taxonId in inCladesSet:
        if not taxonomy.isLeaf(taxonId):
            print('Taxon id %s does not represent a leaf clade in the taxonomy.' % taxonId)

    if verbose:
        print('Initial checks done.')

    # taxonIds that should be excluded
    toExcludeSet = set()
    for taxonId in inCladesSet:
        taxonIdAtRank = taxonomy.getTaxonIdAtRank(taxonId, rank)
        if taxonIdAtRank is None:  # the lineage is not defined at this rank ! try a lower rank !
            print('Taxon id: "%s" is not defined at rank: "%s"' % (taxonId, rank))
            currentRank = rank  # find a lower rank at which it's defined
            while currentRank in toLowerRank:
                currentRank = toLowerRank[currentRank]
                taxonIdAtRank = taxonomy.getTaxonIdAtRank(taxonId, currentRank)
                if taxonIdAtRank is not None:
                    break
            if taxonIdAtRank is None:
                taxonIdAtRank = taxonId
                currentRank = _STRAIN
            print('Taxon id: %s will be masked at rank: %s' % (taxonId, currentRank))

        # all child clades (and itself)
        toExcludeSet.add(int(taxonIdAtRank))
        toExcludeSet.update(set(map(int, taxonomy.getAllChildren(taxonIdAtRank))))

    # all clades that should be excluded (there is at least one sequence for each taxonId in the reference)
    toExcludeSet.intersection_update(refCladesSet)
    if verbose:
        print('Data to mask collected done.')

    print('To exclude: ', len(toExcludeSet))

    # exclude data from the reference
    if action == 'cl':
        # generates a list of taxonIds
        out = csv.OutFileBuffer(os.path.join(outDir, 'exclude_list.txt'))
        for taxonId in toExcludeSet:
            out.writeText(str(taxonId) + '\n')
        out.close()
    elif action == 'mr':
        # masked reference sequences (create sim links to files that were not excluded)
        for fastaFilePath in glob.glob(os.path.join(os.path.normpath(inDir), r'*.f[na][as]')):  # *.fas or *.fna
            taxonId = _refFilePathToTaxonId(fastaFilePath)  # taxonId.1.fna or taxonId.1.fas
            if taxonId not in toExcludeSet:
                # assert os.name == 'posix'
                os.symlink(fastaFilePath, os.path.join(outDir, os.path.basename(fastaFilePath)))
    elif action == 'mg':
        # exclude sequences from the marker gene databases
        for mapFilePath in glob.glob(os.path.join(os.path.normpath(inDir), r'*.tax')):

            # get entries that can stay in the mapping and fasta files
            allowedEntriesSet = set(map(_mgSeqIdToTaxonId, csv.getColumnAsList(mapFilePath, sep='\t')))
            allowedEntriesSet.difference_update(toExcludeSet)

            # filter out entries from the mapping file
            csv.filterOutLines(mapFilePath, os.path.join(outDir, os.path.basename(mapFilePath)),
                               allowedEntriesSet, entryModifyFunction=_mgSeqIdToTaxonId, colNum=0, sep='\t')

            # filter out entries from the fasta file
            fastaFilePath = str(mapFilePath.rsplit('.', 1)[0] + '.fna')
            fas.filterOutSequences(fastaFilePath, os.path.join(outDir, os.path.basename(fastaFilePath)),
                                   allowedEntriesSet, seqNameModifyFunction=_mgSeqIdToTaxonId)
    else:
        assert False, 'Not supported action!'

    taxonomy.close()
    if verbose:
        print('Data masked done.')


def isRankAllowed(rank):
    """
        Returns true if the taxonomic rank is allowed to be excluded using this module

        @param rank: taxonomic rank
        @type rank: str
        @rtype: bool
    """
    return bool(rank in _RANKS)


###------------------------------------------------------------------------


def _getLabelsCreateFasta():
    """
        To process the original mercier dataset with 59 strains. Take only contigs that were mapped to the reference
        genomes. Output a fasta file and a mapping file.
    :rtype : None
    """
    # input fasta file
    fastaFilePath = '/net/metagenomics/projects/PPSmg/data/V35/contigs_1000.txt' #contigs_1000.txt
    seqIdToSeq = fas.fastaFileToDict(fastaFilePath)

    # contigs mapped to genome names
    nameLabelsFilePath = '/net/metagenomics/projects/PPSmg/data/V35/contigs_1000_blast_labels.txt' #contigs_1000_blast_labels.txt
    seqIdToNameLabels = csv.getMapping(nameLabelsFilePath, 0, 1, sep='\t', comment = '#')

    # mapping: genome name -> taxon id
    genomeListFilePath = '/net/metagenomics/projects/PPSmg/data/V35/genome_list2.txt' #genome_list.txt
    nameLabelToNcbid = csv.getMapping(genomeListFilePath, 0, 2, sep=';', comment = '#')

    # to store mapped sequences
    outFastaFilePath = '/net/metagenomics/projects/PPSmg/data/V35/contigsMappedBlast1000.fna' #contigsMappedBlast1000.fna
    outFasta = csv.OutFileBuffer(outFastaFilePath)
    # to stored taxonomic mapping of mapped sequences
    outLabelsFilePath = '/net/metagenomics/projects/PPSmg/data/V35/contigsMappedBlast1000Labels.txt' #contigsMappedBlast1000Labels.txt
    outLabels = csv.OutFileBuffer(outLabelsFilePath)

    for seqId in seqIdToSeq:
        if seqId in seqIdToNameLabels:
            outFasta.writeText(str('>' + str(seqId) + '\n' + seqIdToSeq[seqId] + '\n'))

    outFasta.close()
    print 'fasta created'

    for seqId in seqIdToSeq:
        if seqId in seqIdToNameLabels:
            nameLabel = seqIdToNameLabels[seqId][0]
            ncbid = nameLabelToNcbid[nameLabel][0]
            outLabels.writeText(str(str(seqId) + '\t' + str(ncbid) + '\n'))

    outLabels.close()
    print 'labels created'


def getFirstLabelAtAllowedRank():
    rank='species' # !!!!!!!

    predFile1 = '/net/metagenomics/projects/PPSmg/data/V35/contigsMappedBlast1000Labels.txt'
    predFile2 = '/net/metagenomics/projects/PPSmg/data/V35/contigsMappedBlast1000LabelsSpecies.txt'
    seqIdToLabel = csv.getMapping(predFile1, 0, 1, sep='\t', comment = '#')
    outPred = csv.OutFileBuffer(predFile2)

    taxonomy = tax.TaxonomyNcbi('/net/metagenomics/projects/PPSmg/data/nobackup/NCBI20120828/ncbiTax/ncbitax_sqlite.db')

    for seqId in seqIdToLabel:
        ncbid = int(seqIdToLabel[seqId][0])
        while not taxonomy.isRankNcbidAllowed(ncbid):
            ncbid = taxonomy.getParentNcbid(ncbid)
        outPred.writeText(str(seqId + '\t' + str(ncbid) + '\n'))

    taxonomy.close()
    outPred.close()


def removeLines(mg):
    removeListFilePath = '/net/metagenomics/projects/PPSmg/data/V35/genome_ncbids.txt'
    #removeListFilePath = '/net/metagenomics/projects/PPSmg/data/V35/genome_accession_silva.txt'
    srcFilePath = str('/net/metagenomics/projects/PPSmg/data/markerGenes/db/' + mg + '_bact+arch_dnaV.tax')
    dstFilePath = str('/net/metagenomics/projects/PPSmg/data/V35/genomesRemoved/markerGenes/db/' + mg + '_bact+arch_dnaV.tax')
    #srcFilePath = str('/net/metagenomics/projects/PPSmg/data/silva/' + mg + '_silva106_ncbitax.bacteria+archaea.tax' )
    #dstFilePath = str('/net/metagenomics/projects/PPSmg/data/V35/genomesRemoved/silva/' + mg + '_silva106_ncbitax.bacteria+archaea.tax' )
    pattern = r'.*ncbid:([0-9]+)$'
    #pattern = r'^([^\-]+)\-.*$'

    removeSet = set(csv.getColumnAsList(removeListFilePath, colNum=0, comment='#'))
    col0 = csv.getColumnAsList(srcFilePath, colNum=0, sep='\t', comment='#')
    col1 = csv.getColumnAsList(srcFilePath, colNum=1, sep='\t', comment='#')
    out = csv.OutFileBuffer(dstFilePath)
    removed = 0
    for col0,col1 in zip(col0,col1):
        if re.sub(pattern, r'\1', col0) not in removeSet:
            out.writeText(str(col0 + '\t' + col1 + '\n'))
        else:
            removed += 1

    out.close()
    print mg, 'removeLines', removed


def removeEntries(mg):
    """
        Removes sequences from the marker gene files at the level from species, genus, family etc.
    """
    removeListPath = '/net/metagenomics/projects/PPSmg/data/V35/genome_ncbids_species.txt'
    srcFilePath = str('/net/metagenomics/projects/PPSmg/data/markerGenes2/db/' + mg + '_bact+arch_dnaV.tax')
    dstFilePath = str('/net/metagenomics/projects/PPSmg/data/V35/mgScenarios/speciesRemoved/db/' + mg + '_bact+arch_dnaV.tax')
    out = csv.OutFileBuffer(dstFilePath)
    removeSet = set(csv.getColumnAsList(removeListPath, colNum=0, comment='#'))
    removeSetInt = set()
    removeSetIds = set()
    removed = 0
    for s in removeSet:
        if s != '':
            removeSetInt.add(int(s))
    col0 = csv.getColumnAsList(srcFilePath, colNum=0, sep='\t', comment='#')
    col1 = csv.getColumnAsList(srcFilePath, colNum=1, sep='\t', comment='#')
    for col0,col1 in zip(col0,col1):
        lineSetInt = set()
        for s in col1.split(';'):
            if s != '':
                lineSetInt.add(int(s))
        if len(removeSetInt.intersection(lineSetInt)) > 0: #the intersection is not empty
            removed += 1
            removeSetIds.add(col0)
        else:
            out.writeText(str(col0 + '\t' + col1 + '\n'))
    out.close()

    print mg, 'removedEntries', removed

    srcFilePath = str('/net/metagenomics/projects/PPSmg/data/markerGenes2/db/' + mg + '_bact+arch_dnaV.noalign.fna')
    dstFilePath = str('/net/metagenomics/projects/PPSmg/data/V35/mgScenarios/speciesRemoved/db/' + mg + '_bact+arch_dnaV.noalign.fna')
    out = csv.OutFileBuffer(dstFilePath)
    seqIdToSeq = fas.fastaFileToDict(srcFilePath)
    removed=0
    for seqId in seqIdToSeq:
        if seqId in removeSetIds:
            removed += 1
        else:
            out.writeText(str('>' + str(seqId) + '\n' + str(seqIdToSeq[seqId]) + '\n'))

    out.close()

    print mg, 'removedSeq', removed


def removeSequences(mg):
    removeListFilePath = '/net/metagenomics/projects/PPSmg/data/V35/genome_ncbids.txt'
    #removeListFilePath = '/net/metagenomics/projects/PPSmg/data/V35/genome_accession_silva.txt'
    srcFilePath = str('/net/metagenomics/projects/PPSmg/data/markerGenes/db/' + mg + '_bact+arch_dnaV.noalign.fna')
    dstFilePath = str('/net/metagenomics/projects/PPSmg/data/V35/genomesRemoved/markerGenes/db/' + mg + '_bact+arch_dnaV.noalign.fna')
    #srcFilePath = str('/net/metagenomics/projects/PPSmg/data/silva/' + mg + '_silva106_ncbitax.bacteria+archaea.fna' )
    #dstFilePath = str('/net/metagenomics/projects/PPSmg/data/V35/genomesRemoved/silva/' + mg + '_silva106_ncbitax.bacteria+archaea.fna' )
    pattern = r'.*ncbid:([0-9]+)$'
    #pattern = r'^([^\-]+)\-.*$'

    removeSet = set(csv.getColumnAsList(removeListFilePath, colNum=0, comment='#'))
    seqIdToSeq = fas.fastaFileToDict(srcFilePath)
    out = csv.OutFileBuffer(dstFilePath)
    removed = 0
    for seqId in seqIdToSeq:
        if re.sub(pattern, r'\1', str(seqId)) not in removeSet:
            out.writeText(str('>' + str(seqId) + '\n' + str(seqIdToSeq[seqId]) + '\n'))
        else:
            removed += 1

    out.close()
    print mg, 'removeSequences', removed


def filterSequences():
    """
        To filter sequences with a specific label.
    """
    inFileName = '/net/metagenomics/projects/PPSmg/data/V35/contigsMappedBlast1000.fna'
    outFileName = '/net/metagenomics/projects/PPSmg/data/V35/nostocRemoved/contigsMappedBlast1000NostocRm.fna'
    mapFileName = '/net/metagenomics/projects/PPSmg/data/V35/contigsMappedBlast1000Labels.txt'
    labelRemove = 103690
    #seq id -> label
    labelToIdsDict = csv.getMapping(mapFileName, 1, 0, sep='\t', comment = '#')
    allowedNamesSet = set()
    for i in labelToIdsDict:
        if int(i) != int(labelRemove):
            for j in labelToIdsDict[i]:
                allowedNamesSet.add(j)

    fas.filterOutSequences(inFileName, outFileName, allowedNamesSet)


def _test1():
    pass


def _test2():

    #getFirstLabelAtAllowedRank()
    list = ['dnaG', 'infC', 'pgk', 'rplA', 'rplC', 'rplE', 'rplK', 'rplM', 'rplP', 'rplT', 'rpoB', 'rpsC', 'rpsI', 'rpsK', 'rpsS', 'tsf', 'frr', 'nusA', 'pyrG', 'rplB',
               'rplD', 'rplF', 'rplL', 'rplN', 'rplS', 'rpmA', 'rpsB', 'rpsE', 'rpsJ', 'rpsM', 'smpB', '5S']
    #list = ['ssuref','lsuparc']
    for mg in list:
        removeLines(mg)
        removeSequences(mg)

def _test3():
    t = _TaxonomyWrapMD('/Users/ivan/Documents/work/binning/taxonomy/20121122/ncbitax_sqlite.db')
    id = 83763
    print 'exists:', t.exists(id)
    print 'isLeaf:', t.isLeaf(id)
    print 'all children:', t.getAllChildren(id)


if __name__ == "__main__":
    _main()
    #_test3()

