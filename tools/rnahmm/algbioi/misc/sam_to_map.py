"""
    Gets labels for a simulated dataset from a SAM mapping file and accession to ncbi taxon id file.
"""
import os
import argparse
from algbioi.com import csv


def parseSam(samFile):
    """
        @param samFile: sam file from an assembler
        @return: mapping: seqId -> accession
    """
    d = {}
    for line in open(samFile):
        tokens = line.split()
        if line.startswith('@') or len(tokens) < 3:
            continue
        name = tokens[0]
        labelTokens = tokens[2].split('|')
        if len(labelTokens) < 4:
            print('Label token "%s" is of a unknown type' % tokens[2])
            continue
        label = labelTokens[3]
        if (name in d) and (d[name] != label):
            print('Contig "%s" is present with two different labels "%s" and "%s", the first will be kept.' %
                  (name, d[name], label))
        else:
            d[name] = label
    return d


def samToMap(samFile, accToNcbiFile, outMapFile):
    """

        @param samFile: sam file from an assembler
        @param accToNcbiFile: mapping: accessions -> ncbi taxon ids
        @param outMapFile: output file or directory
    """
    accToNcbi = csv.getMapping(accToNcbiFile, 0, 1, sep='\t')
    contigToAcc = parseSam(samFile)
    out = csv.OutFileBuffer(outMapFile)
    for contigId, acc in contigToAcc.iteritems():
        taxonId = accToNcbi.get(acc, None)
        if taxonId is None:
            print("No mapping for %s %s" % (contigId, acc))
        else:
            out.writeText(contigId + '\t' + taxonId[0] + '\n')
    out.close()


def _main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog='')

    parser.add_argument('-s', '--sam-file', nargs=1, type=file, required=True, help='A SAM file.',
                        metavar='config.sam', dest='s')

    parser.add_argument('-m', '--accession-to-taxon-id-map', nargs=1, type=file, required=True,
                        help='Accession id to NCBI taxon id mapping file, there should be a mapping for each accession '
                             'number contained in the SAM file. ',
                        dest='m')

    parser.add_argument('-o', '--output-file', required=True, action='store', nargs=1,
                        help='Output file name. ',
                        dest='o')

    args = parser.parse_args()

    samFile = args.s[0].name
    accToTaxonIdFile = args.m[0].name
    outArg = args.o[0]
    if os.path.isdir(outArg):
        outFile = os.path.join(outArg, 'map.csv')
    else:
        outFile = os.path.join(os.getcwd(), outArg)

    samToMap(samFile, accToTaxonIdFile, outFile)


if __name__ == "__main__":
    _main()