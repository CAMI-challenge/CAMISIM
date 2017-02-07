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


    Creates the marker gene database for the 16S and 23S genes as well as for the Amphora marker genes.

    The input 16S and 23S marker genes are generated using the "arb" software from the files downloaded from the
    SILVA database (http://www.arb-silva.de). The database with the "Amphora" marker genes was generated using script
    "download_mg.py"
    that downloads sequences from NCBI.
"""
import os
import glob
import argparse

from Bio import Entrez
from Bio import SeqIO

from algbioi.com import csv
from algbioi.com import fasta as fas
from algbioi.com import taxonomy_ncbi as tax


class _TaxonomyWrap():
    def __init__(self, taxonomyFile,
                 allowedRanks=['root','superkingdom','phylum','class','order','family','genus','species']):
        """
            Represents the ncbi taxonomy, buffers entries to efficiently compute the path to the root.
            @param taxonomyFile: database in sqlite3 format
            @param allowedRanks:
        """
        self._taxonomy = tax.TaxonomyNcbi(taxonomyFile, allowedRanks)
        self._allowedRanks = allowedRanks
        self._parentDict = {}  # map: taxonId -> parent taxonId
        self._rankDict = {}  # map: taxonId -> rank

    def _getParent(self, taxonId):
        if taxonId in self._parentDict:
            return self._parentDict[taxonId]
        else:
            parent = self._taxonomy.getParentNcbid(taxonId)
            self._parentDict[taxonId] = parent
            return parent

    def _getRank(self, taxonId):
        if taxonId in self._rankDict:
            return self._rankDict[taxonId]
        else:
            rank = self._taxonomy.getRank(taxonId)
            self._rankDict[taxonId] = rank
            return rank

    def close(self):
        self._taxonomy.close()

    def getPathToRoot(self, taxonId):
        """
            Returns a path to the root in the taxonomy. Considers only allowed ranks. If e.g. a species is not defined
            at a specified rank, the taxon id from a lower rank is used at this rank instead.
            @param taxonId: taxon id of a node where the path starts
            @type taxonId: int
            @return: semicolon separated taxon ids, from second allowed rank till taxonId, ends with a semicolon
            @rtype: str
        """
        assert isinstance(taxonId, int)
        rankToTaxonId = {}
        current = taxonId
        while current != 1:
            rank = self._getRank(current)
            if rank in self._allowedRanks:
                rankToTaxonId[rank] = current
            current = self._getParent(current)
            if current is None:
                return None
        ranks  = self._allowedRanks[1:]
        ranks.reverse()
        pathList = []
        current = taxonId
        for rank in ranks:
            if rank in rankToTaxonId:
                current = rankToTaxonId[rank]
            pathList.append(current)

        pathList.reverse()
        return str(";".join(map(str, pathList)) + ';')


def _main():
    """ See the module description."""
    parser = argparse.ArgumentParser(description=__doc__, epilog="""""")

    parser.add_argument('-i', '--input-data-dir', action='store', nargs=1, required=True,
        help="""Directory that contains fasta files and corresponding mapping files, for each "*.tax" (or "*.csv")
                 file there must be a "*.fna" file with the same name. All files with suffix "tax" (or "*.csv")
                 will be considered. (Takes only Bacteria and Archaea)""",
        metavar='input_dir',
        dest='inDir')

    parser.add_argument('-o', '--output-dir', action='store', nargs=1, required=True,
        help='Directory that contains the output files.',
        metavar='out_dir',
        dest='outDir')

    parser.add_argument('-s', '--source-type', required=True, nargs=1, choices=["s","a"],
        help='To determine the source, use "s" for the Silva database and "a" for the Amphora database.',
        dest='srcType')

    parser.add_argument('-t', '--taxonomy-file', nargs=1, type=file, required=True,
        help='NCBI taxonomy database file in the sqlite3 format.', metavar='ncbitax_sqlite.db',
        dest='taxonomy')

    parser.add_argument('-n', '--not-considered-taxonIds', action='store', nargs=1,
        help='Comma separated leaf level or top level taxonIds (as a string) what fill be filtered out. (optional)',
        metavar='"2759,10239,77133,155900,408172,32644, 408170,433727,749907,556182,702656,410661,652676,410659,797283'\
                ',408171,703336,256318,32630,433724,766747,488339,942017,1076179,717931,455559,527640,904678,552539,'\
                '54395,198431,358574,415540,511564,369433,380357,81726,198834,271928,311313,2759,749906,1077529,'\
                '1077529,361146,511563,361147"',
        dest='filterOut')

    # parse arguments
    args = parser.parse_args()
    inDir = args.inDir[0]
    outDir =  args.outDir[0]
    srcType = args.srcType[0]
    filterOutTaxonIdsSet = set()
    try:
        if args.filterOut:
            filterOutTaxonIdsSet.update(set(map(int, str(args.filterOut[0]).split(','))))
    except:
        print('Taxon ids that are to be filtered out are in a wrong format! Comma separated integers are needed!')
        raise

    taxonomy = _TaxonomyWrap(args.taxonomy[0].name)
    for dir in [inDir, outDir]:
        assert os.path.isdir(dir), 'Path: "' + dir + '" does not exists!'

    # create db for each gene
    mapDict = {}  # map: seqId -> ncbid
    for mapFilePath in glob.glob(os.path.join(os.path.normpath(inDir), r'*.[ct][sa][vx]')):  # *.csv or *.tax

        assert mapFilePath.endswith(('.csv', '.tax')), \
            'The mapping files can either end with .csv or .tax ' + mapFilePath

        base = os.path.basename(mapFilePath).rsplit('.', 1)[0]  # cut out dir path and suffix
        fastaDict = fas.fastaFileToDict(os.path.join(os.path.dirname(mapFilePath), (base + '.fna'))) # map: seqId -> seq
        print("Processing: %s seq count: %s" % (base, str(len(fastaDict))))

        if 'a' in srcType:  # Amphora
            mapDict = {}
            for k in csv.getColumnAsList(mapFilePath, colNum=0, sep='\t'):
                v =  int(k.rsplit('|', 1)[1].split(':')[1]) # get ncbid
                assert ((k not in mapDict) or (mapDict[k] == v)), str(
                    'There are at least two different values for key: ' + str(k) + ' in ' + mapFilePath)
                mapDict[k] = v
        elif 's' in srcType:  # Silva
            mapTmp = csv.getMapping(mapFilePath, 0, 2, '\t')
            mapDict = {}
            for k, v in mapTmp.iteritems():
                mapDict[k] = int(v[0])
        else:
            assert False, 'Unsupported source type!'

        # same number of entries in both files (fasta and mapping) ?
        if len(mapDict) != len(fastaDict):
            print(str('%s: The mapping file and the corresponding fasta file have different number of entries: ' +
                      '"%s" "%s" these files will be skipped!') % (base, str(len(mapDict)), str(len(fastaDict))))
            continue

        # are duplicates in the mapping file ?
        count = len(csv.getColumnAsList(mapFilePath))
        if len(mapDict) != count:
            print('%s: The mapping file contained duplicates! unique: %s non-unique: %s' % (
                base, str(len(mapDict)), str(count)))

        # store data to the output directory
        outDna = csv.OutFileBuffer(os.path.join(outDir, str(base + '.fna')))
        outTax = csv.OutFileBuffer(os.path.join(outDir, str(base + '.tax')))
        count = 0
        filteredLeaf = 0
        filteredSup = 0
        notMapped = 0
        noBacArch = 0
        for seqId, taxonId in mapDict.iteritems():
            if taxonId in filterOutTaxonIdsSet:
                filteredLeaf += 1
                continue
            path = taxonomy.getPathToRoot(taxonId)
            if path is None:
                print('Could not find: %s for seqId: %s record skipped!' % (str(taxonId), seqId))
                notMapped += 1
                continue
            topLevel = int(path.split(';', 1)[0])
            if topLevel in filterOutTaxonIdsSet:
                filteredSup += 1
                continue
            if topLevel not in [2, 2157]:  # Bacteria, Archaea
                noBacArch += 1
                print('NoBactArch: ', topLevel)

            seq = fastaDict[seqId]
            if 'a' in srcType:  # Amphora
                id = seqId
            elif 's' in srcType:  # Silva
                id = str(seqId + '|ncbid:' + str(taxonId))

            outTax.writeText(str(id + '\t' + path + '\n'))
            outDna.writeText(str('>' + id + '\n' + seq + '\n'))
            count += 1

        outDna.close()
        outTax.close()
        print('Stored entries: %s filtered out: %s leaf, %s top level, not mapped: %s' %
              (count, filteredLeaf, filteredSup, notMapped))
        if noBacArch > 0:
            print('WARN: stored %s of non Bacterial and non Archaeal sequences: ' % (noBacArch))

        # Silva:
        #-i /Users/ivan/Documents/work/binning/database/silva111/arbGenerated -s s -t /Users/ivan/Documents/work/binning/taxonomy/20121122/ncbitax_sqlite.db
        # -o /Users/ivan/Documents/work/binning/database/silva111/db -n ...

        # Amphora
        # -i /Users/ivan/Documents/work/binning/database/markerGenes3/mGenesExtracted -s a -t /Users/ivan/Documents/work/binning/taxonomy/20121122/ncbitax_sqlite.db
        # -o /Users/ivan/Documents/work/binning/database/markerGenes3/db

    taxonomy.close()
    print 'done'

def _testTaxonomyWrap():
    taxonomy = _TaxonomyWrap('/Users/ivan/Documents/work/binning/taxonomy/20121122/ncbitax_sqlite.db')
    print(str(taxonomy.getPathToRoot(870603)))

if __name__ == "__main__":
    _main()
    #_testTaxonomyWrap()