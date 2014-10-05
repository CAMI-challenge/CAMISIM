#!/usr/bin/env python


import os
import sys
from Bio import SeqIO

def toStockholmFormat(fastaFile, stockholmFile):

    try:
        fw = open(os.path.normpath(stockholmFile), 'w')
        fr = open(os.path.normpath(fastaFile), "r")
        fw.write('# STOCKHOLM 1.0')
        for record in SeqIO.parse(fr, "fasta"):
            name = record.description
            name = name.replace(' ','_')
            seq = str(record.seq)
            seq = seq.replace('-','.')
            fw.write('\n' + name + '\t' + seq)

        fw.write('\n//')
        fw.close()
        fr.close()
    except Exception:
        print sys.stderr.write(str("Cannot read from file:" + fastaFile + " or write to file:" + stockholmFile))
        raise

if __name__ == "__main__":

    from sys import argv

    #fastaFile = argv[1]
    #stockholmFile = argv[2]

    fastaFile = os.path.normpath('C://Documents and Settings//Administrator//Desktop//temp//amphoraMGen//dnaG.aln')
    stockholmFile = os.path.normpath('C://Documents and Settings//Administrator//Desktop//temp//amphoraMGenS//dnaG.aln')

    toStockholmFormat(fastaFile, stockholmFile)