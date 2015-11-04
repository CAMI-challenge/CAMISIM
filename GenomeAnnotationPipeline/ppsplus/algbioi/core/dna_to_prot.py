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
from Bio import SeqIO


def dnaToProt(dnaFile, protFile):

    try:
        f = open(os.path.normpath(protFile), 'w')
        fr = open(dnaFile, "r")

        for record in SeqIO.parse(fr, "fasta"):
            name = record.description
            dna0 = record.seq
            dna1 = dna0.reverse_complement()

            prot0 = dna0.translate()
            f.write('>' + name + '_p0\n')
            f.write(str(prot0) + '\n')

            prot1 = dna1.translate()
            f.write('>' + name + '_pr0\n')
            f.write(str(prot1) + '\n')

            #frameshift
            prot2 = dna0[1:].translate()
            f.write('>' + name + '_p1\n')
            f.write(str(prot2) + '\n')

            prot3 = dna0[2:].translate()
            f.write('>' + name + '_p2\n')
            f.write(str(prot3) + '\n')

            prot4 = dna1[1:].translate()
            f.write('>' + name + '_pr1\n')
            f.write(str(prot4) + '\n')

            prot5 = dna1[2:].translate()
            f.write('>' + name + '_pr2\n')
            f.write(str(prot5) + '\n')
        f.close()
    except Exception:
        print "Cannot create a file or write to it:", protFile
        raise

if __name__ == "__main__":

    from sys import argv

    dnaFile = argv[1]
    protFile = argv[2]

    #dnaFile="D://A_Phylo//A_Metagenomic//pPPS//workspace//pPPS//wdir02//inputTW.fas.ids"
    #protFile="D://A_Phylo//A_Metagenomic//pPPS//workspace//pPPS//wdir02//inputTW.fas.ids.prot"

    #dnaFile=os.path.normpath("C://Documents and Settings//Administrator//Desktop//temp//SRM_Large_Contigs.fna.ids")
    #protFile=os.path.normpath("C://Documents and Settings//Administrator//Desktop//temp//SRM_Large_Contigs.fna.ids.prot")

    #dnaFile=os.path.normpath("C://Documents and Settings//Administrator//Desktop//temp//HGcontigs.fna.ids")
    #protFile=os.path.normpath("C://Documents and Settings//Administrator//Desktop//temp//HGcontigs.fna.ids.prot")

    #dnaFile=os.path.normpath("C://Documents and Settings//Administrator//Desktop//temp//biofilmcontigs.fna.ids")
    #protFile=os.path.normpath("C://Documents and Settings//Administrator//Desktop//temp//biofilmcontigs.fna.ids.prot")

    #dnaFile=os.path.normpath("C://Documents and Settings//Administrator//Desktop//temp//TM7contigs.fna.ids")
    #protFile=os.path.normpath("C://Documents and Settings//Administrator//Desktop//temp//TM7contigs.fna.ids.prot")

    dnaToProt(dnaFile, protFile)