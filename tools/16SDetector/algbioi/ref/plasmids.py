#!/usr/bin/env python

import sys
import signal
import string
from Bio import SeqIO


#from FastaFileFunctions import fastaFileToDict
#from FastaFileFunctions import getSequenceToBpDict
#from TabSepFileFunctions import getMapping
from algbioi.com.csv import OutFileBuffer
#from TabSepFileFunctions import getColumnAsList
#from Common import noNewLine


def scanForPlasmids():
    """
        Reads in a genbank file from stdin, store an accession number of a sequence that contain word
        'plasmid' (ignorecase) in the record.description Store also for each sequence a list of locations
        that contain 'plasmid' (ignorecase) in feature values

        USAGE:
        nohup time zcat /local/johdro/refdata/static/ncbi-genomes-bacteria_20121122/dna.gbff.gz \
        /local/johdro/refdata/static/ncbi-draftgenomes-bacteria_20121122/dna-contigs.gbff.gz \
        /local/johdro/refdata/static/ncbi-draftgenomes-bacteria_20121122/dna-scaffolds.gbff.gz \
        /local/johdro/refdata/static/ncbi-hmp_20121016/dna-contigs.gbff.gz \
        /local/johdro/refdata/static/ncbi-hmp_20121016/dna-scaffolds.gbff.gz \
        /local/johdro/refdata/static/ncbi-refseq-microbial_56/dna.gbff.gz \
        | python /net/metagenomics/projects/PPSmg/scripts/scripts25/plasmids.py &
    """

    #plasmidAccessionFile = '/Users/ivan/Documents/nobackup/refseq/bacterialGenomes/accession_test.txt'
    plasmidAccessionFile = '/local/igregor/ref_20121122/nobackup/plasmid_accessions.txt'
#    plasmidRegionsFile = ''
    seqCount = 0
    uniqueSeqCount = 0
    plasmidSeqCount = 0

    accessionSet = set([])
#    accessionToLocationList = dict([])

    outPlasmidAccessions = OutFileBuffer(plasmidAccessionFile)
#    outPlasmidLocations = OutFileBuffer(plasmidRegionsFile)
    for record in SeqIO.parse(sys.stdin, "genbank"):
        seqCount += 1

        if str(record.id) in accessionSet:
            continue
        else:
            accessionSet.add(str(record.id))

        uniqueSeqCount += 1

        #is the whole sequence a plasmid
        if 'PLASMID' in str(record.description).upper():
            outPlasmidAccessions.writeText(str(str(record.id) + '\n'))
            plasmidSeqCount += 1
#        else:
#            #get plasmid regions
#            for feature in record.features:
#
#                for key in feature.qualifiers:
#                    if str(key) not in ['misc_feature']: #check this
#                        continue
#                    val = feature.qualifiers[key]
#                    if 'PLASMID' in string.upper(str(val)):
#                        if str(record.id) in accessionToLocationList:
#                            accessionToLocationList[str(record.id)].append(str(feature.location))
#                        else:
#                            accessionToLocationList[str(record.id)] = str(feature.location)

#    for key, val in accessionToLocationList.iteritems():
#        outPlasmidLocations.writeText(str(key + '\t' + str(val).replace("'",'').replace(']','').replace('[','')   + '\n'))

    outPlasmidAccessions.close()
#    outPlasmidLocations.close()

    print 'seqCount', seqCount
    print 'uniqueSeqCount', uniqueSeqCount
    print 'plasmidSeqCount', plasmidSeqCount


def scan():
    recordCount = 0
    recordPlasmid = 0
    breakCount = 50
    featureType = set([])
    keySet = set([])

    for record in SeqIO.parse(sys.stdin, "genbank"):
        recordCount += 1

        seqId = record.id

        if string.find(record.description, 'plasmid') != -1:
            #print record.description
            recordPlasmid += 1

        for feature in record.features:
            #print feature.location, 'location'

            for key in feature.qualifiers:
                val = feature.qualifiers[key]
                if 'PLASMID' in string.upper(str(val)):
                    print feature
                    #keySet.add((feature.type, key, str(val)))
                    #keySet.add((feature.type, key))
            #print 'location', feature.location

            #print 'type', feature.type
            #if feature.type == "CDS":
            #    print feature

            #    for x in feature.qualifiers:
            #        print x, feature.qualifiers[x]

            #for xrefentry in feature.qualifiers:
            #    print xrefentry
                #    ( key, val ) = xrefentry.split(":")
                #    if key == "taxon":
                #        taxonId = int(val)
                #        break

            #for qualifier in feature.qualifiers:
            #    print qualifier

            #featureType.add(feature.type)
            #print feature

            #for xrefentry in feature.qualifiers["db_xref"]:
            #    print featurexrefentry



            #if feature.type == "source":
                #for xrefentry in feature.qualifiers["db_xref"]:
                #    ( key, val ) = xrefentry.split(":")
                #    if key == "taxon":
                #        taxonId = int(val)
                #        break


        if recordCount == breakCount:
            break

    print 'recordCount', recordCount
    #print 'recordPlasmid', recordPlasmid
    #print 'featureType', featureType
    #for key in keySet:
    #    print key




       # if seqId in seqIdSet:
       #     print seqId, 'already in set', seqId

       # seq = str(record.seq)
       # cumulativeLen += len(seq)

      #  if len(string.replace(noNewLine(seq),'N','')) == 0:
 #           zeros += 1

       # taxonId = None

       # for feature in record.features:
       #     if feature.type == "source":
       #         for xrefentry in feature.qualifiers["db_xref"]:
       #             ( key, val ) = xrefentry.split(":")
       #             if key == "taxon":
       #                 taxonId = int(val)
       #                 break
       #     if taxonId != None:
       #         break


       # if taxonId == None:
       #     print 'could not find taxonId for', seqId
       # else:
       #     taxonSet.add(taxonId)

#    print 'record count', recordCount
#    print 'seq count', len(seqIdSet)
#    print 'taxon id count', len(taxonSet)
#    if len(seqIdSet) > 0:
#        print 'avg. seq. len', cumulativeLen/len(seqIdSet)
#    print 'zeros', zeros



# MAIN
if __name__ == "__main__":
    # handle broken pipes
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    scanForPlasmids()
    #scan()