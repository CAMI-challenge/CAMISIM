#!/usr/bin/env python

"""
    Compares two versions of the NCBI databases in terms of the ncbids listed in a file.
"""

import argparse

from algbioi.com.csv import getColumnAsList
from algbioi.com.taxonomy_ncbi import TaxonomyNcbi


def _main():

    parser = argparse.ArgumentParser(description=__doc__, epilog='''For each ncbi taxon id from the taxon-ids-file
        verifies whether the ncbi taxon id is in both databases and whether the paths to the root in both databases
        contain the same nodes (in terms of ncbi taxon ids).''')

    parser.add_argument('-f', '--first-database', nargs=1, type=file, required=True,
        help='Database file in the sqlite3 format.', metavar='first_ncbitax_sqlite.db',
        dest='db1')

    parser.add_argument('-s', '--second-database', nargs=1, type=file, required=True,
        help='Database file in the sqlite3 format.', metavar='second_ncbitax_sqlite.db',
        dest='db2')

    parser.add_argument('-t', '--taxon-ids-file', nargs=1, type=file, required=True,
        help='A file that contains one ncbi taxon id per line.', metavar='ncbi_taxon_ids.csv',
        dest='ids')

    args = parser.parse_args()

    tax1 = TaxonomyNcbi(args.db1[0].name)
    tax2 = TaxonomyNcbi(args.db2[0].name)
    ncbiTaxonList = getColumnAsList(args.ids[0].name, colNum=0, comment='#')

    for id in ncbiTaxonList:
        name1 = tax1.getScientificName(id)
        name2 = tax2.getScientificName(id)
        #if name1 is None:
        #    print("Id: " + str(id) + " not found in the first db.")
        if name2 is None:
            print("Id: " + str(id) + " not found in the second db.")
        if name1 is None and name2 is None:
            print("Id: " + str(id) + " found in neither of the db.")

        #if name1 != name2:
        #    print('For id:' + str(id) + ' scientific names differ,  first: ' +  str(name1) + ' second: ' +  str(name2))

        #rank1 = tax1.getRank(id)
        #rank2 = tax2.getRank(id)
        #if rank1 != rank2:
        #    print('For id:' + str(id) + ' ranks differ,  first: ' +  str(rank1) + ' second: ' +  str(rank2))

        parentSet1 = tax1.getParentsNcbidSet(id)
        parentSet2 = tax2.getParentsNcbidSet(id)
        if parentSet1 !=  parentSet2 and name1 != name2 and name1 is not None:
            print('For id:' + str(id) +
                  ' parent sets differ,  first: ' +  str(parentSet1) + ' second: ' +  str(parentSet2))
            print('For id:' + str(id) + ' scientific names differ,  first: ' +  str(name1) + ' second: ' +  str(name2))




# Silva and MG: /Users/ivan/Documents/work/binning/taxonomy/20120828/ncbitax_sqlite.db

# Ref.: /Users/ivan/Documents/work/binning/taxonomy/20121122/ncbitax_sqlite.db

# Ref. ncbids: /Users/ivan/Documents/work/binning/ref/NCBI20121122/ncbi_taxon_ids.csv



def _test():
    #-f /Users/ivan/Documents/work/binning/taxonomy/20120828/ncbitax_sqlite.db -s /Users/ivan/Documents/work/binning/taxonomy/20121122/ncbitax_sqlite.db -t /Users/ivan/Documents/work/binning/ref/NCBI20121122/ncbi_taxon_ids.csv
    pass


if __name__ == "__main__":
    #_test()
    _main()
