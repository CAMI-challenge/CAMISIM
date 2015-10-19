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
import sys
import time
import operator
import argparse
import random
from xml.dom import minidom
from xml.dom import Node
from Bio import Entrez
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
#Entrez.email = 'A.N.Other@example.com'
from algbioi.core.taxonomy import Taxonomy


def test():
    markerGeneName =  'rpsC' #'rpsI' #'rpsS'  # 'rpsK'
    annotationDir = os.path.normpath('D:/A_Phylo/A_Metagenomic/data/markerGenes/annotation')
    outDir = os.path.normpath('D:/A_Phylo/A_Metagenomic/data/markerGenes/mGenesExtracted')
    taxonomy = Taxonomy(os.path.normpath('D:/A_Phylo/A_Metagenomic/data/ncbiTaxonomy20111007/ncbitax_sqlite.db'), ['superkingdom','phylum','class','order','family','genus','species'])
    relaxGeneNames = False
    recSkipCount = 0
    firstErrorStop = False
    createGeneDb(markerGeneName, annotationDir, outDir, taxonomy, relaxGeneNames, recSkipCount, firstErrorStop)


def protNameEntryToGid(protEntry):
    return protEntry.split('|')[2].split(':')[1]

def dnaNameEntryToGid(dnaEntry):
    return dnaEntry.split('|')[2].split(':')[1]


def createGeneDb(markerGeneName, annotationDir, outDir, taxonomy, relaxGeneNames = False, recSkipCount=0, firstErrorStop=False):
    """
        First X records in the annotation file are skipped.

        DOCUMENTATION
        http://www.ncbi.nlm.nih.gov/gene?term=rpsK
        http://docs.python.org/library/xml.dom.minidom.html
        http://eutils.ncbi.nlm.nih.gov/corehtml/query/static/efetchseq_help.html
        http://biopython.org/wiki/Seq

        @param relaxGeneNames: if True then don`t check if the gene names in the
    """
    annotationXMLfile = os.path.join(annotationDir, str(markerGeneName + '.xml'))

    #get the new taxonomy file !!!
    #taxonomy = Taxonomy(os.path.normpath('D:/A_Phylo/A_Metagenomic/data/ncbiTaxonomy20111007/ncbitax_sqlite.db'), ['superkingdom','phylum','class','order','family','genus','species'])

    outBuffDnaFasta = OutputBuffer(os.path.join(outDir, str(markerGeneName + '_bact+arch_dna.fna')), 'dna', 'fasta', taxonomy)
    outBuffDnaTax = OutputBuffer(os.path.join(outDir, str(markerGeneName + '_bact+arch_dna.tax')), 'dna', 'taxonomy', taxonomy)

    outBuffProtFasta = OutputBuffer(os.path.join(outDir, str(markerGeneName + '_bact+arch_prot.fna')), 'prot', 'fasta', taxonomy)
    outBuffProtTax = OutputBuffer(os.path.join(outDir, str(markerGeneName + '_bact+arch_prot.tax')), 'prot', 'taxonomy', taxonomy)

    storedDNAEntries = 0
    storedProtEntries = 0

    geneNameSynSet = set([])
    geneNameSynSet.add(markerGeneName.lower())

    #fileXML = os.path.normpath('D:/A_Phylo/A_Metagenomic/data/markerGenes/annotation/18S.xml')
    #markerGeneName = '18S'

    #annotationXML = os.path.normpath('D:/A_Phylo/A_Metagenomic/data/markerGenes/annotation/rpsK.xml')
    #annotationXML = os.path.normpath('D:/A_Phylo/A_Metagenomic/data/markerGenes/testAnnotation.xml')
    #markerGeneName = 'rpsK'

    #annotationXML = os.path.normpath('D:/A_Phylo/A_Metagenomic/data/markerGenes/annotation/frr.xml')
    #markerGeneName = 'frr'

#    annotationXML = os.path.normpath('D:/A_Phylo/A_Metagenomic/data/markerGenes/annotation/rplP.xml')
#    markerGeneName = 'rplP'

#    annotationXML = os.path.normpath('D:/A_Phylo/A_Metagenomic/data/markerGenes/annotation/rplP.xml')
#    markerGeneName = 'rplP'

#    annotationXML = os.path.normpath('D:/A_Phylo/A_Metagenomic/data/markerGenes/annotation/rplS.xml')
#    markerGeneName = 'rplS'

    #annotationXML = os.path.normpath('D:/A_Phylo/A_Metagenomic/data/markerGenes/annotation/rpsB.xml')
    #markerGeneName = 'rpsB'

    xmlDoc = minidom.parse(annotationXMLfile)
    xmlRoot = xmlDoc.documentElement
    geneAnnotationList = xmlRoot.getElementsByTagName('Entrezgene')

    recordCount = 0


    for geneAnnotation in geneAnnotationList: #for geneAnnotation in geneAnnotationList:

        #print record info at the end
        dump = False
        recordCount += 1
        if recordCount <= recSkipCount:
            continue

        gid = getLeafNodeValueByPath(geneAnnotation,
                        'Entrezgene_track-info/Gene-track/Gene-track_geneid')

        if gid == None:
            sys.stderr.write('gid is missing at record number: ' + str(recordCount) + 'record skipped')
            continue

#        print 'taxName: ', getLeafNodeValueByPath(geneAnnotation,
#                        'Entrezgene_source/BioSource/BioSource_org/Org-ref/Org-ref_taxname')

        lineage = getLeafNodeValueByPath(geneAnnotation,
                        'Entrezgene_source/BioSource/BioSource_org/Org-ref/Org-ref_orgname/OrgName/OrgName_lineage')

        if lineage == None:
            sys.stderr.write(str(gid) + ': lineage entry is missing, record skipped')
            continue

        if lineage.find('Bacteria') == -1 and lineage.find('Archaea') == -1:
            print str(str(gid) + ': lineage skipped (not Bacteria or Archeae): ' + lineage )
            continue

#        accesionNumber = getLeafNodeValueByPath(geneAnnotation,
#                        'Entrezgene_locus/Gene-commentary/Gene-commentary_accession')
        intervalFrom2 = getLeafNodeValueByPath(geneAnnotation,
                        'Entrezgene_locus/Gene-commentary/Gene-commentary_products/Gene-commentary/Gene-commentary_genomic-coords/Seq-loc/Seq-loc_int/Seq-interval/Seq-interval_from')

        intervalTo2 = getLeafNodeValueByPath(geneAnnotation,
                        'Entrezgene_locus/Gene-commentary/Gene-commentary_products/Gene-commentary/Gene-commentary_genomic-coords/Seq-loc/Seq-loc_int/Seq-interval/Seq-interval_to')

        strandStr = getLeafNodeValueByPath(geneAnnotation,
                        'Entrezgene_locus/Gene-commentary/Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval/Seq-interval_strand/Na-strand','value')
        strandStr2 = getLeafNodeValueByPath(geneAnnotation,
                        'Entrezgene_locus/Gene-commentary/Gene-commentary_products/Gene-commentary/Gene-commentary_genomic-coords/Seq-loc/Seq-loc_int/Seq-interval/Seq-interval_strand/Na-strand','value')

        if operator.xor(strandStr == None, strandStr2 == None) or (strandStr == strandStr2):
            if strandStr == None:
                strandStr = strandStr2
            if strandStr == 'minus':
                strandStr = str(2)
            else:
                strandStr = str(1)#default strand is plus
        else:
            sys.stderr.write(str(gid) + ': Got different values for strand, record skipped: ' + strandStr + ' ' + strandStr2 + '\n')
            continue

        #always take the first one
        intervalGiList = dict([])
        intervalGiList[0] = getLeafNodeValueByPath(geneAnnotation,
                        'Entrezgene_locus/Gene-commentary/Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval/Seq-interval_id/Seq-id/Seq-id_gi')
        intervalGiList[1] = getLeafNodeValueByPath(geneAnnotation,
                        'Entrezgene_locus/Gene-commentary/Gene-commentary_products/Gene-commentary/Gene-commentary_genomic-coords/Seq-loc/Seq-loc_int/Seq-interval/Seq-interval_id/Seq-id/Seq-id_gi')
        intervalGiList[2] = getLeafNodeValueByPath(geneAnnotation,
                        'Entrezgene_gene-source/Gene-source/Gene-source_src-int')
        intervalGiList[3] = getLeafNodeValueByPath(geneAnnotation,
                        'Entrezgene_non-unique-keys/Dbtag/Dbtag_tag/Object-id/Object-id_id')

        intervalGi = None
        for i in intervalGiList:
            if intervalGiList[i] != None:
                if intervalGi == None:
                    intervalGi = intervalGiList[i]
                else:
                    assert intervalGi == intervalGiList[i], str(str(intervalGiList[0]) + ' ' + str(intervalGiList[1]) + ' '
                                                                + str(intervalGiList[2]) + ' ' + str(intervalGiList[3]))

        seqLocWhole = getLeafNodeValueByPath(geneAnnotation,
                        'Entrezgene_locus/Gene-commentary/Gene-commentary_products/Gene-commentary/Gene-commentary_seqs/Seq-loc/Seq-loc_whole/Seq-id/Seq-id_gi')

        #TESTS
        useDnaSeq = True
        useProtSeq = True

        #GET DNA SEQUENCE
        if intervalGi != None and intervalFrom2!= None and intervalTo2 != None and strandStr != None:
            intervalFrom = int(intervalFrom2) + 1
            intervalTo = int(intervalTo2) + 1
            try:
                #http://eutils.ncbi.nlm.nih.gov/corehtml/query/static/efetchseq_help.html
                handle = Entrez.efetch(db="nucleotide", id=str(intervalGi), seq_start=str(intervalFrom),
                                   seq_stop=str(intervalTo), strand=strandStr, rettype="gb", retmode="xml")
                dnaXml = handle.read()
                dnaXmlRoot = minidom.parseString(dnaXml).documentElement
            except Exception:
                sys.stderr.write(str(gid) + ': unable to parse Entrez result for DNA: Entrez.efetch(db="nucleotide", id=str(' + str(intervalGi)
                           + '), seq_start=str(' + str(intervalFrom2) + '+1), seq_stop=str(' + str(intervalTo2)
                           + '+1), strand=' + str(strandStr) + ', rettype="gb", retmode="xml")\n')
                dna = None
                dnaTranslation = None
                dnaTranslationTable = None
                dnaNcbid = None
                dnaGeneId = None
                dnaIntervalGi = None
                dnaGeneName = None
                dump = True
                useDnaSeq = False
            else:
                dna = getLeafNodeValueByPath(dnaXmlRoot,'GBSeq/GBSeq_sequence')

                dnaGeneName = getValueByKey2(dnaXmlRoot, 'GBSeq/GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier',
                                               'GBQualifier_name', 'gene', 'GBQualifier_value', None, 'gene', markerGeneName)
                #print 'h1: dnaGeneName', dnaGeneName

                if dnaGeneName == None:
                    dnaGeneName = getValueByKey(dnaXmlRoot, 'GBSeq/GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier',
                                               'GBQualifier_name', 'gene', 'GBQualifier_value')
                #print 'h2: dnaGeneName', dnaGeneName

                dnaGeneNameSyn = None
                if dnaGeneName != None and dnaGeneName.lower() != markerGeneName.lower():
                    dnaGeneNameSyn = getValueByKey2(dnaXmlRoot, 'GBSeq/GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier',
                                               'GBQualifier_name', 'gene_synonym', 'GBQualifier_value', None, 'gene', dnaGeneName)
                #print 'h3 dnaGeneNameSyn: ', dnaGeneNameSyn

                if dnaGeneName == None:
                    dnaSectionGeneName = markerGeneName
                else:
                    dnaSectionGeneName = dnaGeneName
                #print 'h4 dnaSectionGeneName', dnaSectionGeneName

                if dnaGeneNameSyn != None and dnaGeneNameSyn.lower() == markerGeneName.lower():
                    dnaGeneName = str(markerGeneName)
                    geneNameSynSet.add(dnaSectionGeneName.lower())
                #print 'h5 dnaGeneName:', dnaGeneName

                dnaTranslation = getValueByKey2(dnaXmlRoot, 'GBSeq/GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier',
                                               'GBQualifier_name', 'translation', 'GBQualifier_value', None, 'gene', dnaSectionGeneName)
                dnaTranslationTable = getValueByKey2(dnaXmlRoot, 'GBSeq/GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier',
                                               'GBQualifier_name', 'transl_table', 'GBQualifier_value', None, 'gene', dnaSectionGeneName)
                dnaNcbid = getValueByKey(dnaXmlRoot, 'GBSeq/GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier',
                                               'GBQualifier_name', 'db_xref', 'GBQualifier_value', 'taxon:')
                dnaGeneId = getValueByKey2(dnaXmlRoot, 'GBSeq/GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier',
                                               'GBQualifier_name', 'db_xref', 'GBQualifier_value', 'GeneID:', 'gene', dnaSectionGeneName)
                dnaIntervalGi = getValueByKey2(dnaXmlRoot, 'GBSeq/GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier',
                                               'GBQualifier_name', 'db_xref', 'GBQualifier_value', 'GI:', 'gene', dnaSectionGeneName)


                if (dna == None) or (dnaTranslationTable == None) or (dnaNcbid == None):
                    print(str(gid) + ': DNA: some important entries are missing\n')
                    dump = True
        else:
            sys.stderr.write(str(gid) + ': unable to retrieve DNA: Entrez.efetch(db="nucleotide", id=str(' + str(intervalGi)
                       + '), seq_start=str(' + str(intervalFrom2) + '+1), seq_stop=str(' + str(intervalTo2)
                       + '+1), strand=' + str(strandStr) + ', rettype="gb", retmode="xml")\n')
            dna = None
            dnaTranslation = None
            dnaTranslationTable = None
            dnaNcbid = None
            dnaGeneId = None
            dnaIntervalGi = None
            dnaGeneName = None
            dump = True
            useDnaSeq = False


        #GET PROTEIN SEQUENCE
        if seqLocWhole != None:
            try:
                #http://eutils.ncbi.nlm.nih.gov/corehtml/query/static/efetchseq_help.html
                handle = Entrez.efetch(db="nucleotide", id=str(seqLocWhole), rettype="gb", retmode="xml")
                protXml = handle.read()
                protXmlRoot = minidom.parseString(protXml).documentElement
            except Exception:
                sys.stderr.write(str(gid) + ': unable to parse Entrez result for PROT: Entrez.efetch(db="nucleotide", id=str(' + str(seqLocWhole) + '), rettype="gb", retmode="xml")\n')
                prot = None
                protTranslationTable = None
                protNcbid = None
                protGeneId = None
                protGeneName = None
                useProtSeq = False
                dump = True
            else:
                prot = str(getLeafNodeValueByPath(protXmlRoot,'GBSeq/GBSeq_sequence')).upper()
                protTranslationTable = getValueByKey(protXmlRoot, 'GBSeq/GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier',
                                               'GBQualifier_name', 'transl_table', 'GBQualifier_value')
                protNcbid = getValueByKey(protXmlRoot, 'GBSeq/GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier',
                                               'GBQualifier_name', 'db_xref', 'GBQualifier_value', 'taxon:')
                protGeneId = getValueByKey(protXmlRoot, 'GBSeq/GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier',
                                               'GBQualifier_name', 'db_xref', 'GBQualifier_value', 'GeneID:')
                protGeneName = getValueByKey(protXmlRoot, 'GBSeq/GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier',
                                               'GBQualifier_name', 'gene', 'GBQualifier_value')
                protGeneNameSyn = getValueByKey(protXmlRoot, 'GBSeq/GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier',
                                             'GBQualifier_name', 'gene_synonym', 'GBQualifier_value')
                if (protGeneName != None and protGeneName.lower() != markerGeneName.lower()
                    and protGeneNameSyn != None and protGeneNameSyn.lower() == markerGeneName.lower()):
                    geneNameSynSet.add(protGeneNameSyn.lower())
                    protGeneName = markerGeneName

                if (prot == None) or (protTranslationTable == None) or (protNcbid == None):
                    print(str(gid) + ': PROT: some important entries are missing\n')
                    dump = True
        else:
            sys.stderr.write(str(gid) + ': unable to retrieve PROT: Entrez.efetch(db="nucleotide", id=str(' + str(seqLocWhole) + '), rettype="gb", retmode="xml")\n')
            prot = None
            protTranslationTable = None
            protNcbid = None
            protGeneId = None
            protGeneName = None
            useProtSeq = False
            dump = True

        #MATCHING DNA/PROT SEQUENCES

        #translate dna seq to prot seq
        dnaTranslated = None
        if dnaTranslationTable == None and protTranslationTable == None:
            print str(str(gid) + ': no record for translation tables')
            dump = True
        else:
            if dnaTranslationTable != None and protTranslationTable != None:
                if int(dnaTranslationTable) != int(protTranslationTable):
                    sys.stderr.write(str(gid) + ': different translation tables: ' + str(dnaTranslationTable) + ' '+ str(protTranslationTable) + '\n')
                    useDnaSeq = False
                    useProtSeq = False
                    dump = True
            if dnaTranslationTable != None:
                translationTable = int(dnaTranslationTable)
            else:
                translationTable = int(protTranslationTable)
            try:
                #translate dna to protein
                codingDna = Seq(str(dna), generic_dna)
                dnaTranslated = codingDna.translate(table=translationTable, to_stop=True, cds=True)
            except Exception, err:
                sys.stderr.write(str(gid) + ': Exception in DNA translate: ' + str(err) + '\n')
                dnaTranslated = None
                useDnaSeq = False
                dump = True
        if dnaTranslated == None:
            print str(str(gid) + ' could not translate dna sequence')
            dump = True

        #is the dna record consistent
        if dnaTranslated != None and dnaTranslation != None:
            if str(dnaTranslated) != str(dnaTranslation):
                sys.stderr.write(str(gid) + ': wrong translation within dna record\n')
                useDnaSeq = False
                dump = True

        #does the prot record match the dna record
        if prot != None:
            if (dnaTranslation != None) and (str(prot) != str(dnaTranslation)):
                sys.stderr.write(str(gid) + ': prot sequence doesn`t match to the dna "translation" sequence\n')
                useDnaSeq = False
                useProtSeq = False
                dump = True
            if (dnaTranslated != None) and (str(prot) != str(dnaTranslated)):
                sys.stderr.write(str(gid) + ': prot sequence doesn`t match to the dna "translated" sequence\n')
                useDnaSeq = False
                useProtSeq = False
                dump = True


        #NCBIDs dna record vs prot record
        if dnaNcbid != None and protNcbid != None:
            if int(dnaNcbid) != int(protNcbid):
                sys.stderr.write(str(gid) + ': the ncbids are different dnaNcbid:' + str(dnaNcbid) + ' protNcbid' + str(protNcbid) + '\n')
                useDnaSeq = False
                useProtSeq = False
                dump = True
        elif dnaNcbid == None and protNcbid == None:
                sys.stderr.write(str(gid) + ': both ncbids are missing\n')
                useDnaSeq = False
                useProtSeq = False
                dump = True
        else:
            print str(str(gid) + ': one of the ncbids is missing dnaNcbid:' + str(dnaNcbid) + ' protNcbid:' + str(protNcbid))
            dump = True

        #DNA and GENE IDs
        if dnaGeneId != None and protGeneId != None and int(dnaGeneId) != int(protGeneId):
            sys.stderr.write(str(gid) + ': different gene ids,  dnaGeneId:' + str(dnaGeneId) + ' protGeneId:' + str(protGeneId) + '\n')
            useDnaSeq = False
            useProtSeq = False
            dump = True


        if ((dnaGeneName != None and dnaGeneName.lower() not in geneNameSynSet)
            or (protGeneName != None and protGeneName.lower() not in geneNameSynSet)
            or (dnaGeneName != None and protGeneName != None and dnaGeneName != protGeneName)):
            if not relaxGeneNames:
                sys.stderr.write(str(gid) + ': gene names are wrong, dnaGeneName:' + str(dnaGeneName)
                             + ' protGeneName:' + str(protGeneName) + '\n')
                useDnaSeq = False
                useProtSeq = False
                dump = True
            else:
                print str(str(gid) + ': gene names are not correct (no action is taken), dnaGeneName:' + str(dnaGeneName)
                             + ' protGeneName:' + str(protGeneName) + '\n')

        #RECORD DUMP
        if dump:
            print 'recordCount:', recordCount
            print 'DUMP gid:', gid

            print 'DNA'
            print 'dna:           ', dna
            print 'dnaTranslation:', dnaTranslation
            print 'dnaTranslated: ', dnaTranslated
            print 'dnaTranslationTable:', dnaTranslationTable
            print 'dnaNcbid', dnaNcbid
            print 'dnaGeneId', dnaGeneId
            print str('dnaSeq: ' + str(intervalGi) + ':' + str(intervalFrom2) + '+1 - ' + str(intervalTo2) + '+1 s' + str(strandStr))
            print dnaGeneName

            print 'PROTEIN'
            print '               ', prot
            print 'protTranslationTable', protTranslationTable
            print 'protNcbid', protNcbid
            print 'protGeneId', protGeneId
            print 'seqLocWhole', seqLocWhole
            print protGeneName
            print '---------------------------'
            if firstErrorStop:
                try:
                    print dnaXml
                except Exception:
                    print 'dnaXml not defined'
                print '---------------------------'
                try:
                    print protXml
                except Exception:
                    print 'protXml not defined'
                break

        if dnaNcbid != None:
            ncbid = dnaNcbid
        else:
            ncbid = protNcbid

        #STORE DNA SEQUENCE
        if (useDnaSeq and dna != None and ncbid != None and intervalGi != None
            and intervalFrom2!= None and intervalTo2 != None and strandStr != None and gid != None):
            outBuffDnaFasta.writeDnaRecord(dna, ncbid, intervalGi, str(int(intervalFrom2)+1), str(int(intervalTo2)+1), strandStr, gid)
            outBuffDnaTax.writeDnaRecord(  dna, ncbid, intervalGi, str(int(intervalFrom2)+1), str(int(intervalTo2)+1), strandStr, gid)
            storedDNAEntries += 1
        else:
            sys.stderr.write(str(gid) + ': cannot use this record to get DNA sequence\n')

        #STORE PROT SEQUENCE
        if useProtSeq and prot != None and ncbid != None and seqLocWhole != None and gid != None:
            outBuffProtFasta.writeProtRecord(prot, ncbid, seqLocWhole, gid)
            outBuffProtTax.writeProtRecord(  prot, ncbid, seqLocWhole, gid)
            storedProtEntries += 1
        else:
            sys.stderr.write(str(gid) + ': cannot use this record to get PROT sequence\n')


        #store the metainfo about the prot or dna sequences
        #if useDnaSeq or useProtSeq:
        #    pass

        #sleep for a while - don`t overload the server
        time.sleep(float(random.random()/10))

    outBuffDnaFasta.close()
    outBuffDnaTax.close()
    outBuffProtFasta.close()
    outBuffProtTax.close()

    print str('Stored DNA entries: ' + str(storedDNAEntries) + ' Stored Prot entries: ' + str(storedProtEntries)
              + ' Record count: ' + str(recordCount))


def _getLeafNodeValueByPath(node, pathList, resultList, attrName, collectNodes = False):
    """
        Parse a path in a xml document and return all corresponding entries.
    """
    try:
        if len(pathList) == 0:
            if collectNodes:
                resultList.append(node)
            elif attrName != None:
                if node.hasAttribute(attrName):
                    resultList.append(node.getAttribute(attrName))
            else:
                assert len(node.childNodes) == 1
                resultList.append(node.childNodes[0].nodeValue)
        else:
            head = pathList[0]
            tail = pathList[1:]

            children = node.getElementsByTagName(head)
            candidateChildList = []
            for child in children:
                if child.parentNode == node:
                    candidateChildList.append(child)
            for child in candidateChildList:
                _getLeafNodeValueByPath(child, list(tail), resultList, attrName, collectNodes)

    except Exception:
        sys.stderr.write('Cannot parse path: ' + str(pathList) + ' at node: ' + str(node) + '\n')
        raise


def getValueByKey(node, path, keyTag, key, valueTag, valuePrefix = None):
    """
        #xmlRoot,'GBSeq/GBSeq_feature-table/GBFeature/GBFeature_quals/GBQualifier','GBQualifier_name','translation','GBQualifier_value'
    """
    pathList = path.rsplit('/')
    resultList = []
    attrName = None
    collectNodes = True
    resultValueDict = dict([])
    _getLeafNodeValueByPath(node, pathList, resultList, attrName, collectNodes)

    count = 0
    for tag in resultList:
        count += 1
        value = None
        if tag.nodeType == Node.ELEMENT_NODE:
            nameList = tag.getElementsByTagName(keyTag)
            assert len(nameList) == 1
            name = nameList[0]
            if name.childNodes[0].nodeValue == key:
                #print name.childNodes[0].nodeValue
                valueList = tag.getElementsByTagName(valueTag)
                assert len(valueList) <= 1, str(valueList)
                if len(valueList) == 0:
                    continue
                value = valueList[0].childNodes[0].nodeValue
                if valuePrefix != None:
                    if value.count(valuePrefix) == 1:
                        value = value.replace(valuePrefix, '', 1)
                        #break
                    else:
                        value = None
                #else:
                #    break
        if value != None:
            resultValueDict[count] = value

    list = []
    for i in resultValueDict:
        list.append(resultValueDict[i])
    resultValueSet = set(list)
    if len(list) == 0:
        return None
    if len(resultValueSet) != 1:
        sys.stderr.write('getValueByKey ambiguous result: ' + str(resultValueSet))
        return None
    else:
        return list[0]


def getValueByKey2(node, path, keyTag, key, valueTag, valuePrefix, refKey, refValue):
    pathList = path.rsplit('/')
    groupTag = pathList[len(pathList)-1]
    pathList = pathList[0:len(pathList)-1]
    resultList = []
    attrName = None
    collectNodes = True
    _getLeafNodeValueByPath(node, pathList, resultList, attrName, collectNodes)

    confirmedValues = []
    notContradictedValues = []

    for group in resultList:
        value = None
        refFound = False
        contraRefFount = False
        elements = group.getElementsByTagName(groupTag)
        for element in elements:
            if element.nodeType != Node.ELEMENT_NODE:
                continue
            keyTagList = element.getElementsByTagName(keyTag)
            valueTagList = element.getElementsByTagName(valueTag)
            if len(keyTagList) == 0 or len(valueTagList) == 0:
                continue
            assert len(keyTagList) == 1 and len(valueTagList) == 1, str(str(keyTagList) + str(valueTagList))
            assert len(keyTagList[0].childNodes) == len(valueTagList[0].childNodes) == 1, str(str(keyTagList[0].childNodes) + str(valueTagList[0].childNodes))
            k = keyTagList[0].childNodes[0].nodeValue
            v = valueTagList[0].childNodes[0].nodeValue
            if (k == key) and ((valuePrefix == None) or (v.count(valuePrefix) == 1)):
                if valuePrefix != None:
                    v = v.replace(valuePrefix, '', 1)
                if value == None:
                    value = v
                else:
                    assert value == v, str(str(value) + ' ' + str(v))

            if k == refKey:
                if v.lower() == refValue.lower():
                    refFound = True
                else:
                    contraRefFount = True

        if value != None:
            if refFound:
                confirmedValues.append(value)
            elif not contraRefFount:
                notContradictedValues.append(value)

    if len(confirmedValues) > 0:
        if  len(set(confirmedValues)) == 1:
            return confirmedValues[0]
        else:
            sys.stderr.write('getValueByKey2 ambiguous result None returned: ' + str(confirmedValues))
            return None
    elif len(notContradictedValues) > 0:
        if len(set(notContradictedValues)) == 1:#, str(notContradictedValues)
            return notContradictedValues[0]
        else:
            sys.stderr.write('multiple possible values in getValueByKey2, None returned: ' + str(notContradictedValues) + '\n')
            return None
    else:
        return None


def getLeafNodeValueByPath(node, path, attrName=None):
    """
        Gets the element value or an attribute value or None if there is not one occurance of this element or attribute.
    """
    pathList = path.rsplit('/')
    resultList = []
    _getLeafNodeValueByPath(node, pathList, resultList, attrName)
    resultSet = set(resultList)
    if len(resultSet) == 1:
        return resultList[0]
    elif len(resultSet) == 0:
        #sys.stderr.write('Cannot find path: ' + str(path) + '\n')
        return None
    else:
        sys.stderr.write('Ambiguous result, possible values: ' + str(resultList) + ' for path:' + path + '\n')
        return None


class OutputBuffer():
    """
        To output a sequence or a taxonomy record to a file.
    """
    def __init__(self, outFilePath, seqType, fileType, taxonomy):
        try:
            self.outFile = open(outFilePath, 'w')
        except Exception:
            sys.stderr.write('Cannot open a file for writing: ' + str(outFilePath) + '\n')
            raise
        if seqType == 'dna':
            self.seqType = 'dna'
        elif seqType == 'prot':
            self.seqType = 'prot'
        else:
            raise Exception('The supported "seqType" is either "dna" or "prot", but not' + str(seqType))
        if fileType == 'stockholm':
            self.fileType = 'sto'
            try:
                self.outFile.write('# STOCKHOLM 1.0\n')
            except Exception:
                sys.stderr.write('Cannot write to file: ' + str(outFilePath))
                raise
        elif fileType == 'fasta':
            self.fileType = 'fna'
        elif fileType == 'taxonomy':
            self.fileType = 'tax'
        else:
            raise Exception('The supported "fileType" is "stockholm", "fasta", or "taxonomy"\n')
        self.taxonomy = taxonomy


    def writeDnaRecord(self, dna, ncbid, intervalGi, intervalFrom, intervalTo, strand, gid):
        assert dna != None and ncbid != None and intervalGi != None and intervalFrom != None and intervalTo != None and strand != None
        assert self.seqType == 'dna', 'The object was created to store only the dna sequences'
        if int(strand) == 2:
            strandSign = '-'
        else:
            strandSign = '+'
        name = str('gi|' + str(intervalGi) + ':' + str(intervalFrom) + '-' + str(intervalTo) + strandSign
                   + '|gid:' + str(gid) + '|ncbid:' + str(ncbid))
        self._writeToFile(name, dna, ncbid)


    def writeProtRecord(self, prot, ncbid, seqLocWhole, gid):
        assert prot != None and ncbid != None and seqLocWhole != None
        assert self.seqType == 'prot', 'The object was created to store only the protein sequences'
        name = str('gi|' + str(seqLocWhole) + '|gid:' + str(gid) + '|ncbid:' + str(ncbid))
        self._writeToFile(name, prot, ncbid)


    def _writeToFile(self, name, seq, ncbid):
        if self.fileType == 'sto':
            entry = str(name + '\t' + seq + '\n')
        elif self.fileType == 'fna':
            entry = str('>' + name + '\n' + seq + '\n')
        elif self.fileType == 'tax':
            #pathToRoot = self.taxonomy.getPathToRootSemicolonSeparated(ncbid)
            #assert pathToRoot != None, str(name + ' ' + seq)
            #entry = str(name + '\t' + pathToRoot + '\n')
            entry = str(name + '\t' + ncbid + '\n')
        else:
            assert False, str('The file type is of the wrong type' + self.fileType)
        try:
            self.outFile.write(entry)
        except Exception:
            sys.stderr.write('Cannot write to file: ' + str(self.outFile.name) + '\n')
            raise


    def close(self):
        try:
            if self.fileType == 'sto':
                self.outFile.write('//')
            self.outFile.close()
        except Exception, e:
            sys.stderr.write(e + '\n' + 'Unable to close file: ' + self.outFile + '\n')



def main():
    parser = argparse.ArgumentParser(description='''Creates database for a marker gene''', epilog=''' ''')

    parser.add_argument('-m', '--marker-gene-name', action='store', nargs=1, required=True,
                        help='The name of a specific marker gene.',
                        dest='marker')

    parser.add_argument('-a', '--annotation-dir', action='store', nargs=1, required=True,
                        help='The name of the directory that contains annotation files.',
                        dest='annotationDir')

    parser.add_argument('-o', '--output-dir', action='store', nargs=1, required=True,
                        help='The name of the directory where the output files will be stored.',
                        dest='outDir')

    parser.add_argument('-t', '--taxonomyDb', action='store', nargs=1, required=True,
                        help='Taxonomy database file (SQLite).',
                        dest='taxonomyDb')

    parser.add_argument('-r', '--relax-gene-names', action='store_true',
                        help='If enabled, the script doesn`t control if the gene names are correct.',
                        dest='relaxGeneNames')

    parser.add_argument('-s', '--rec-skip', action='store', nargs=1,
                        help='The number of records that will be skipped at the beginning of the annotation file.',
                        dest='recSkip')

    parser.add_argument('-p', '--print-first-error', action='store_true',
                        help='The script stops after first error occurs',
                        dest='firstErrorStop')

    args = parser.parse_args()

    markerGeneName = str(args.marker[0])

    annotationDir = os.path.normpath(str(args.annotationDir[0]))

    outDir = os.path.normpath(str(args.outDir[0]))


    taxonomy = Taxonomy(os.path.normpath(str(args.taxonomyDb[0])), ['superkingdom','phylum','class','order','family','genus','species'])

    if args.recSkip:
        recSkip = int(args.recSkip[0])
    else:
        recSkip = 0

    if args.firstErrorStop:
        firstErrorStop = True
    else:
        firstErrorStop = False

    if args.relaxGeneNames:
        relaxGeneNames = True
    else:
        relaxGeneNames = False

    createGeneDb(markerGeneName, annotationDir, outDir, taxonomy, relaxGeneNames, recSkip, firstErrorStop)




if __name__ == "__main__":
    #if sys.platform == 'darwin' or sys.platform == 'win32':
    #    test()
    #else:
    main()


