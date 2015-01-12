import re
import os
import sys

from Bio.Seq import Seq

from algbioi.com import fasta as fas
from algbioi.com import csv
from algbioi.com import taxonomy_ncbi as tax
import traceback
import traceback

def sayHello2(a, b, d=4, e=5, c=3):
    """
    General description.

    @param t:
    @param n:
    @param uu: uu description
        @type uu: int
    @param ss: ssu param
        @type ss: str
    @param bb: bb something
        type bb: float
    """

    s = set()

    s.add('a')
    print s
    s.add('b')
    s = Seq('ATGC')
    print a

    #b='jkgh'

    print("%s %s %s %s %s" % (a, b, c, d, e))

    return int(2345)



def test2():
    print 'new test'
    # s = 'lsuparc_silva106_ncbitax.bacteria+archaea.tax'
    # print s[(s.rindex('.', 0, s.rindex('.')) + 1):s.rindex('.')]
    # print 'done'
    sayHello2('a', 'b', 44, e=33, c=33)
    # print 'more'
    sayHello2('a','a','a')



def stat():
    """


    """
    #file = '/Volumes/hhu-hera/PPSmg/data/nobackup/mercier51Strains/contigs_soapdenovo-20121119.fna'
    #file = '/Volumes/hhu-hera/PPSmg/data/mercier50/nobackup/seed0/soap_.contig'
    #file = '/Volumes/hhu-hera/PPSmg/data/mercier50/nobackup/seed5/soap_.contig'
    #file = '/Volumes/hhu-hera/PPSmg/data/mercier50/nobackup/seed5_63mer/soap_.contig'
    #file = '/Volumes/hhu-hera/PPSmg/data/mercier50/nobackup/seed5_127/soap_seed5_127.contig'
    #fas.fastaFileToDict(file)
    file = '/Volumes/hhu-hera/PPSmg/data/mercier50/nobackup/seed0/ReadsR_seed0.fa'
    c = 0
    bp = 0
    minLen = 1000
    maxLen = 0
    totalBp = 0
    totalCount = 0
    for k, v in fas.fastaFileToDict(file).iteritems():
        l = len(v)
        totalCount += 1
        if l > 1000:
            c += 1
            bp += l
        totalBp += l
        if l < minLen:
            minLen = l
        elif l > maxLen:
            maxLen = l

    print('Bigger than 1000bp (contigs, bp):', c, bp)
    print('maxLen, minLen, avgLen:', minLen, maxLen, (totalBp / totalCount))
    print('total:', totalCount)





def toPercent(costList, costIdxHundredPercent=2):
    """ For PTree, relative cost comparison. """
    percent = costList[costIdxHundredPercent] / 100.0
    return map(lambda x: round(x / percent, 3), costList)


def toPercent2(timeStrList, timeIdxHundredPercent=2):
    """ For PTree, relative time comparison. """
    timeList = timeStrList.split(';')
    secList = []
    for t in timeList:
        m = 0
        s = 0
        h = 0
        for i in t.split(' '):
            if 's' in i:
                s = float(re.sub(r'([0-9\.]+)s', r'\1', i))
            elif 'm' in i:
                m = float(re.sub(r'([0-9]+)m', r'\1', i))
            elif 'h' in i:
                h = float(re.sub(r'([0-9]+)h', r'\1', i))
            else:
                assert False, 'Unknown entry'
        sec = s + (m * 60) + (h * 60 * 60)
        secList.append(sec)
    percent = secList[timeIdxHundredPercent] / 100.0
    return map(lambda x: round(x / percent, 3), secList)


def filterOutSequences(fastaFile, predFile, outFastaFile, outPredFile, minBp=1000):
    seqIdToSeq = fas.fastaFileToDict(fastaFile)
    seqIdToBp = fas.getSequenceToBpDict(fastaFile)
    seqIdToPred = csv.predToDict(predFile)
    outFasta = csv.OutFileBuffer(outFastaFile)
    outPred = csv.OutFileBuffer(outPredFile)

    totalBp = 0
    taken = 0
    takenBp = 0
    for seqId, seq in seqIdToSeq.iteritems():
        bp = seqIdToBp[seqId]
        totalBp += bp
        if bp >= minBp:
            taxonId = seqIdToPred[seqId]
            outFasta.writeText('>' + str(seqId) + '\n' + str(seq) + '\n')
            outPred.writeText(str(seqId) + '\t' + str(taxonId) + '\n')
            taken += 1
            takenBp += bp
    outFasta.close()
    outPred.close()

    print('Total sequences: ', len(seqIdToSeq))
    print('Total size: ', totalBp)
    print('Taken sequences: ', taken)
    print('Taken size:', takenBp)


def fastaBySeqNameList(inSeqIdList, inFastaFile, outFastaFile):
    """
        Generates outFastaFile that contains sequences that are in the inSeqIdList and inFastaFile.
    """
    seqIdList = csv.getColumnAsList(inSeqIdList)
    seqIdToSeq = fas.fastaFileToDict(inFastaFile)
    out = csv.OutFileBuffer(outFastaFile)
    for seqId in seqIdList:
        seq = seqIdToSeq.get(seqId, None)
        if seq is None:
            print("Can't find sequence for seqId: %s" % seqId)
        else:
            out.writeText('>' + str(seqId) + '\n' + str(seq) + '\n')
    out.close()


def testExeption():
    d = {}
    try:
        print (str(1/0))
    except Exception as ex:
        #print(ex.message)
        #print(ex.args)
        print traceback.print_exc(file=sys.stdout)
    print('done')


def sortColumn(inCsvFilePath, outCsvFilePath=None, colNum=0, separator=',', sortReverse=False):
    """
        Sorts a csv file according to the given column.

        @param inCsvFilePath: input csv file path
        @param outCsvFilePath: output csv file path sorted according to the given column (if None then prints to stdout)
        @param colNum: the column according to which the file will be sorted
        @param separator: column separator
        @param reverse: whether the file should be sorted in the reverse order
    """
    pass




# """
# 	Simple script to estimate inconsistent predictions at each rank
#
# 	@author: Johannes (based on Ivan's code)
# """

from algbioi.eval.consistency import Consistency
import argparse
from sys import stdin, stdout, stderr, exit
from os import fdopen
import signal

def countConsistentPerRank( cons ):
    consistent = {}
    unconsistent = {}
    weights = cons._contigNameToBp
    taxonomy = cons.getTaxonomy()._taxonomy
    for scaffold in cons.getScaffoldsDict().values():
        consistent_taxids = scaffold.getPathSet()
        print consistent_taxids
        for contig in scaffold.getContigsNameList():
            taxid = cons._getPred( contig )
            rank = taxonomy.getRank( taxid )

            if rank == None:
                print taxid

            if taxid in consistent_taxids:
                try:
                    consistent[rank] += weights[contig]
                except KeyError:
                    consistent[rank] = weights[contig]
            else:
                try:
                    unconsistent[rank] += weights[contig]
                except KeyError:
                    unconsistent[rank] = weights[contig]
    return consistent, unconsistent



def _main2():
    """
        Main function.
    """
    parser = argparse.ArgumentParser(description='Computes the scaffold-contig consistency based on '
    											'the "maximum support path".',
    								epilog=__doc__)

    parser.add_argument('-w', '--weights', nargs=1, type=str, required=True,
    					help='Tab separated weights per sequence (likely its length).',
    					metavar='contigs.fna.seqlen',
    					dest='w')

    parser.add_argument('-p', '--predictions', nargs=1, type=str, required=False,
    					default=[ stdin ],
    					help='Tab separated prediction files (first column contig name, last column predicted ncbi taxon id. If not specified it will read from stdin',
    					metavar='pred.tax',
    					dest='p')

    parser.add_argument('-m', '--mapping', nargs=1, type=str, required=True,
    					help='Tab separated scaffold-contig mapping file (first column scaffold name, second column contig name.',
    					metavar='group_to_seqname.tsv',
    					dest='m')

    parser.add_argument('-d', '--database', nargs=1, type=str, required=True,
    					help='Database file in the sqlite3 format.', metavar='ncbi-taxonomy.sqlite',
    					dest='d')

    signal.signal(signal.SIGPIPE, signal.SIG_DFL) #handle broken pipes

    args = parser.parse_args()

    cons = Consistency( args.w[0], args.p[0], args.m[0], args.d[0], None, None, None, False)

    consistent, unconsistent = countConsistentPerRank( cons )
    for rank in frozenset( consistent.keys() + unconsistent.keys() ):
        stdout.write( "%s\t%i\t%i\n" % (rank, consistent.get(rank,0), unconsistent.get(rank,0)) )


def indexToDna(index, length, mapping={0: 'A', 1: 'T', 2: 'G', 3: 'C'}):
    s = str(bin(index))[2:]  # get binary number, trim starting 0b
    s = ((length * 2) - len(s)) * '0' + s  # add starting zeros
    dna = ''
    for i in range(length):  # for each dna character
        dna += mapping[int(s[i*2:(i+1)*2], 2)]  # map binary to dna
    return dna

def indexLineToDnaLine(indexLine, kmerLen):
    line = ''
    for entry in indexLine.split('\t'):
        index, count = entry.split(':')
        if line != '':
            line += '\t'
        line += indexToDna(int(index), kmerLen) + ':' + str(count)
    return line

#{'AG': 3, 'TT': 2, 'GA': 4, 'TG': 3, 'TA': 1, 'TC': 2, 'CT': 2}

def DecimalToBinary(inTxtFilePath, n, outTxtFilePath):
    """





    @param outTxtFilePath:
    @param n:
    @param inTxtFilePath:
    """
    inpath = open(inTxtFilePath, 'r')
    lines = inpath.readlines()

    for line in lines:
        # d1 = dict (((lambda i: (int(i[0]), int(i[1])))(element.split(':')) for element in line.split('\t')))
        d1 = dict(((lambda i: (((bin(int(i[0])))[2:].zfill(2 * n)), int(i[1])))(element.split(':')) for element in
                   line.split('\t')))
        # for i in d1:

        # d2 = dict (int(i[0]),int(i[1]))

    print d1.keys()

    lists = []

    values = []
    for BinaryData in d1.keys():
        # lists = [BinaryData[i:i + 2] for i in range(0, len(BinaryData), 2)]
        lists.append([BinaryData[i:i + 2] for i in range(0, len(BinaryData), 2)])
        values.append(d1.get(BinaryData))

    print lists

    knownData = {'00': 'A', '01': 'T', '10': 'G', '11': 'C'}

    # Yao has difficulties from here and yao does not know how to write dict into file.

    results = []

    for data in lists:

        #value = [knownData.get(y) for y in data if y in knownData.keys()]
        value = []
        for y in data:
            if y in knownData.keys():
                value.append(knownData.get(y))


        results.append(value)

    print results

    Kmer = []

    for Nucleotide in results:

        Kmer.append("".join(Nucleotide))

    print Kmer

    new_dic = dict(zip(Kmer,values))

    print new_dic

    inpath.close()

    outpath = open(outTxtFilePath, 'w')

    for i in Kmer:

        outpath.write(i+'\t')



        # for k, v in new_dic.iteritems():
        #     outpath.write(str(k) + ':' + str(v) + '\t')

def toMostAbundant(d='AG:3	TA:1	TG:3	TC:2	GA:4	TT:2	CT:2'):
    l = []
    for entry in d.split('\t'):
        k, v = entry.split(':')
        l.append((k,v))
    l.sort(key=lambda x: x[1], reverse=True)
    print l


def testException():
    try:
        int(None)
    except Exception as e:
        # pass
        traceback.print_exc(sys.stdout)


def checkTrainData(slFileDir):
    ti = 0
    tf = 0
    for f in os.listdir(slFileDir):
        print str(f)
        lines = 0
        for line in open(os.path.join(slFileDir, f), 'r'):
            lineList = line.split('\t')
            for entry in lineList:
                e = entry.split(':')
                if len(e) > 2:
                    print str(entry), str(lines)
                if len(e) == 2:
                    try:
                        ti += int(e[0])
                        tf += float(e[1])
                    except Exception as e:
                        print 'Exception: ', e.message, str(entry), str(lines)
            lines += 1
        print str(f), ':',str(lines), '--------------'


if __name__ == "__main__":
    checkTrainData('/Users/ivan/Documents/nobackup/vm_hg/train_data')
    #testException()

    #t()

    # DecimalToBinary(inTxtFilePath='/Users/ivan/Documents/nobackup/1.txt', n=2, outTxtFilePath='/Users/ivan/Documents/nobackup/1_out.txt')
    # print indexLineToDnaLine('2:3	4:1	6:3	7:2	8:4	5:2	13:2', 2)
    # print indexToDna(93, 6)
    # print re.sub(r'^>.*\n', 'N', '>name\nATGC')
    # fastaBySeqNameList('/Volumes/hhu-hera/data/CowRumen/chunked070513/chunks2000.scaffold_names',
    #                    '/Users/ivan/Documents/work/binning/data/CowRumen/assembly/cow_rumen_fragmented_velvet_assembly_scaffolds.fas',
    #                    '/Volumes/hhu-hera/data/CowRumen/chunked070513/chunks2000.scaffolds')
    #
    #refToClades('/Volumes/hera - net/metagenomics/projects/PPSmg/data/nobackup/NCBI20121122/sequences',
    #          '/Users/ivan/Documents/nobackup/species_list.txt',
    #          '/Users/ivan/Documents/work/binning/taxonomy/20121122/ncbitax_sqlite.db',
    #          rank='species')
    #testExeption()
    #stat()
    #test2()
    #print toPercent([12,24,50,209,3], 2)
    #filterOutSequences('/Users/ivan/Documents/work/binning/data/mercier51Strains/contigs_soapdenovo-20121119.fna',
    #                   '/Users/ivan/Documents/work/binning/data/mercier51Strains/binning_soapdenovo-20121119.tax',
    #                   '/Users/ivan/Documents/work/binning/data/mercier51Strains/contigs_soapdenovo-20121119_1000bp.fna',
    #                   '/Users/ivan/Documents/work/binning/data/mercier51Strains/binning_soapdenovo-20121119_1000bp.tax')
    pass