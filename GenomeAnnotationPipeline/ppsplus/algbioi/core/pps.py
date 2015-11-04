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
import re
import subprocess

from algbioi.com import csv
from algbioi.com import common
from algbioi.com import fasta
from algbioi.com import taxonomy_ncbi
from algbioi.eval import accuracy
from algbioi.eval import confusion_matrix
from algbioi.core.taxonomy import Taxonomy
from algbioi.core import ref_seq


def computeTrainingAccuracy(workingDir, taWorkingDir, sampleSpecificDir, ppsTrainDataDir, outputDir, ppsInstallDir,
                            ppsScripts, ppsConfigFilePath, predictLogFileName, modelTaxonIdFilePath, databaseFile):
    """
        Computes the training accuracy for the PPS training data.
        This function doesn't consider training data used to train intermediate (misc?) nodes!
        The training data that correspond to the sample specific data is fragmented (via PPS) and
        contained in the training data of different lengths.

        @param workingDir: working directory of the PPS+ pipeline
        @param taWorkingDir: working directory for the accuracy computation
        @param sampleSpecificDir: directory containing the sample specific data
        @param ppsTrainDataDir: directory 'sampled_fasta' containing PPS training data
        @param outputDir: directory for output files
        @param ppsScripts: directory containing PPS scripts
        @param ppsConfigFilePath: the PPS configuration file
        @param ppsInstallDir: directory where PPS is installed
        @param predictLogFileName: logging file for PPS prediction
        @param modelTaxonIdFilePath: file containing all leaf ncbi taxon ids that are modelled
        @param databaseFile: ncbi taxonomy file in the sqlite3 format
    """
    for d in [workingDir, taWorkingDir, sampleSpecificDir,
              ppsTrainDataDir, outputDir, ppsInstallDir, ppsScripts, os.path.dirname(predictLogFileName)]:
        assert os.path.isdir(d), "Directory '%s' doesn't exist!" % d
    for f in [ppsConfigFilePath, databaseFile, modelTaxonIdFilePath]:
        assert os.path.isfile(f), "File '%s' doesn't exist!" % f

    # all directories that contain PPS training data
    trainDirList = [sampleSpecificDir]
    for d in os.listdir(ppsTrainDataDir):
        trainDirList.append(os.path.join(ppsTrainDataDir, d))

    # fasta file with all training sequences
    allTrainFastaFile = os.path.join(taWorkingDir, 'all_train_data.fna')
    out = csv.OutFileBuffer(allTrainFastaFile)
    seqIdToTruePred = {}

    # merge all training fasta files to one fasta file
    for d in trainDirList:
        dName = os.path.basename(d)
        for f in os.listdir(d):
            taxonId = int(os.path.basename(f).rsplit('.', 2)[0])
            for seqId, seq in fasta.fastaFileToDict(os.path.join(d, f)).iteritems():
                if d == sampleSpecificDir:
                    #label = int(str(str(seqId).rsplit('|', 1)[1]).split(':', 1)[1])
                    id = str(taxonId) + '|' + dName + '|' + seqId + '|label:' + str(taxonId)
                else:
                    id = str(taxonId) + '|' + dName + '|' + seqId
                out.writeText('>' + id + '\n' + seq + '\n')
                seqIdToTruePred[id] = taxonId
    out.close()

    # predict the merged file using the generated model
    if os.name == 'posix':
        predictCmd = str(os.path.join(ppsScripts, 'predict.rb') + ' ' + allTrainFastaFile + ' ' + ppsConfigFilePath)
        #print(predictCmd)
        logOut = open(predictLogFileName, 'w')
        predictProc = subprocess.Popen(predictCmd, shell=True, bufsize=-1, cwd=ppsInstallDir, stdout=logOut,
                                       stderr=subprocess.STDOUT)  # stdout=subprocess.STDOUT
        predictProc.wait()
        logOut.close()
        if predictProc.returncode != 0:
            raise Exception("PPS 'predict' training data returned with non-zero status: %s, cmd: %s" %
                            (predictProc.returncode, predictCmd))
    else:
        print("Can't run PPS on a non-posix system!")
        return

    # read in predicted train data
    seqIdToPred = csv.predToDict(allTrainFastaFile + '.nox.fna.out')

    # read fasta file
    seqIdToBp = fasta.getSequenceToBpDict(allTrainFastaFile)

    # leaf taxonIds that are modelled
    modelLeafTaxonIds = set(map(int, csv.getColumnAsList(modelTaxonIdFilePath)))

    taxonomyS = taxonomy_ncbi.TaxonomyNcbi(databaseFile, considerNoRank=True)
    notLeafTaxonIds = set()
    for id in modelLeafTaxonIds:
        notLeafTaxonIds.update(set(map(int, (taxonomyS.getParentsNcbidSet(id)))))
    taxonomyS.close()

    # get only sequences with true taxonId defined at leaf level that is modelled or lower
    seqIdToBp2 = {}
    seqIdToPred2 = {}
    seqIdToTruePred2 = {}
    seqIdToBpMisc = {}
    seqIdToPredMisc = {}
    seqIdToTruePredMisc = {}
    for seqId, bp in seqIdToBp.iteritems():
        label = int(str(str(seqId).rsplit('|', 1)[1]).split(':', 1)[1])
        if label not in notLeafTaxonIds:
            seqIdToBp2[seqId] = bp
            seqIdToPred2[seqId] = seqIdToPred[seqId]
            seqIdToTruePred2[seqId] = seqIdToTruePred[seqId]
        else:
            seqIdToBpMisc[seqId] = bp
            seqIdToPredMisc[seqId] = seqIdToPred[seqId]
            seqIdToTruePredMisc[seqId] = seqIdToTruePred[seqId]
    seqIdToBp = seqIdToBp2
    seqIdToPred = seqIdToPred2
    seqIdToTruePred = seqIdToTruePred2

    # accuracy for all, filter out sample specific data (whole length)
    seqIdToBpNoSampleSpec = {}
    for seqId, bp in seqIdToBp.iteritems():
        if str(seqId).split('|', 2)[1].strip() != os.path.basename(sampleSpecificDir).strip():
            seqIdToBpNoSampleSpec[seqId] = bp

    acc = accuracy.Accuracy(seqIdToBpNoSampleSpec, seqIdToPred, seqIdToTruePred, databaseFile)
    out = csv.OutFileBuffer(os.path.join(outputDir, 'train_accuracy_all.txt'))
    out.writeText(acc.getAccuracyPrint(taxonomy_ncbi.TAXONOMIC_RANKS[1:],
                                       minFracClade=None, minFracPred=None, overview=True))
    out.close()
    taxonomyA = acc.getTaxonomy()
    acc.close(closeTaxonomy=False)

    # accuracy for (misc) nodes
    acc = accuracy.Accuracy(seqIdToBpMisc, seqIdToPredMisc, seqIdToTruePredMisc, taxonomyA)
    out = csv.OutFileBuffer(os.path.join(outputDir, 'train_accuracy_misc.txt'))
    out.writeText(acc.getAccuracyPrint(taxonomy_ncbi.TAXONOMIC_RANKS[1:],
                                       minFracClade=None, minFracPred=None, overview=True))
    out.close()
    acc.close(closeTaxonomy=False)

    # generate the confusion matrices (for the "for all" scenario)
    cm = confusion_matrix.ConfusionMatrix(seqIdToBp, seqIdToPred, seqIdToTruePred, databaseFile,
                                          taxonomy_ncbi.TAXONOMIC_RANKS[1:])
    for rank in taxonomy_ncbi.TAXONOMIC_RANKS[1:]:
        cm.generateConfusionMatrix(rank, os.path.join(outputDir, 'train_accuracy_cmp_all'))
    taxonomyCM = cm.getTaxonomy()
    cm.close(closeTaxonomy=False)

    # accuracy for individual directories (seq lengths)
    # (the sample specific fragments are among PPS sampled fasta)
    for d in trainDirList:
        dName = os.path.basename(d)
        seqIdToBpSub = {}
        seqIdToPredSub = {}
        seqIdToTruePredSub = {}
        for seqId, bp in seqIdToBp.iteritems():
            if str(seqId).split('|', 2)[1].strip() == str(dName).strip():
                seqIdToBpSub[seqId] = seqIdToBp[seqId]
                seqIdToPredSub[seqId] = seqIdToPred[seqId]
                seqIdToTruePredSub[seqId] = seqIdToTruePred[seqId]

        # accuracy
        acc = accuracy.Accuracy(seqIdToBpSub, seqIdToPredSub, seqIdToTruePredSub, taxonomyA)
        out = csv.OutFileBuffer(os.path.join(outputDir, 'train_accuracy_' + dName + '.txt'))
        out.writeText(acc.getAccuracyPrint(taxonomy_ncbi.TAXONOMIC_RANKS[1:],
                                           minFracClade=None, minFracPred=None, overview=True))

        # confusion matrices
        cm = confusion_matrix.ConfusionMatrix(seqIdToBpSub, seqIdToPredSub, seqIdToTruePredSub, taxonomyCM,
                                              taxonomy_ncbi.TAXONOMIC_RANKS[1:])
        for rank in taxonomy_ncbi.TAXONOMIC_RANKS[1:]:
            cm.generateConfusionMatrix(rank, os.path.join(outputDir, 'train_accuracy_cmp_' + dName))
        cm.close(closeTaxonomy=False)

        out.close()
        acc.close(closeTaxonomy=False)
    taxonomyA.close()
    taxonomyCM.close()


def _trainAccuracyDataTest():
    workingDir = '/Users/ivan/Documents/nobackup/trainAccuracyTest/workingDir'
    taWorkingDir = '/Users/ivan/Documents/nobackup/trainAccuracyTest/workingDir/taWorkingDir'
    sampleSpecificDir = '/Users/ivan/Documents/nobackup/trainAccuracyTest/workingDir/sampleSpecificDir'
    ppsTrainDataDir = '/Users/ivan/Documents/nobackup/trainAccuracyTest/sampled_fasta'
    outputDir = '/Users/ivan/Documents/nobackup/trainAccuracyTest/outputDir'
    ppsInstallDir = '/Users/ivan/Documents/nobackup/trainAccuracyTest/ppsInstal'
    ppsScripts = '/Users/ivan/Documents/nobackup/trainAccuracyTest/ppsInstal/ppsScripts'
    ppsConfigFilePath = '/Users/ivan/Documents/nobackup/trainAccuracyTest/config_pps.txt'
    predictLogFileName = '/Users/ivan/Documents/nobackup/trainAccuracyTest/workingDir/predict_train_log.txt'
    databaseFile = '/Users/ivan/Documents/work/binning/taxonomy/20121122/ncbitax_sqlite.db'
    modelTaxonIdFilePath = '/Users/ivan/Documents/nobackup/trainAccuracyTest/workingDir/ncbids.txt'
    computeTrainingAccuracy(workingDir, taWorkingDir, sampleSpecificDir, ppsTrainDataDir, outputDir, ppsInstallDir,
                            ppsScripts, ppsConfigFilePath, predictLogFileName, modelTaxonIdFilePath, databaseFile)


def refToClades(refDir, taxonomyFile, rank='species', outFile=None):
    """
        Returns (stores) a list of all clades (at the given rank) sorted according to the abundance of
        the individual clades. Abundance in respect to the size of the reference data available.

        @param refDir: directory containing reference data (as needed for PPS)
        @param taxonomyFile: ncbi taxonomy in the sqlite3 format
        @param rank: consider clades at this rank
        @param outFile: tab sep file, first column taxon id, second column number of bp (can be None)
        @return: list of tuples (clade, bp)
    """
    taxonomy = taxonomy_ncbi.TaxonomyNcbi(taxonomyFile)
    cladeNcbiToBp = {}
    for fileName in os.listdir(refDir):
        size = os.path.getsize(os.path.join(refDir, fileName))
        ncbid = int(fileName.rsplit('.', 2)[0])
        current = ncbid
        while (current is not None) and (taxonomy.getRank(current) != rank):
            current = taxonomy.getParentNcbid(int(current))
        if current is not None:
            if current in cladeNcbiToBp:
                cladeNcbiToBp[current] += size
            else:
                cladeNcbiToBp[current] = size
        else:
            print('There is no ncbi taxon id defined at rank %s for ncbi taxon id %s' % (rank, ncbid))
    taxonomy.close()

    tuples = []
    for ncbid, size in cladeNcbiToBp.iteritems():
        tuples.append((ncbid, size))
    tuples.sort(key=lambda x: x[1], reverse=True)

    if outFile is not None:
        out = csv.OutFileBuffer(outFile)
        for t in tuples:
            out.writeText(str(t[0]) + '\t' + str(t[1]) + '\n')
        out.close()

    return tuples


def generateCladesForGeneralModel(refSeqDir, taxonomyDatabaseFile,
                                  rank, minTotalCount, minBpPerSpeciesCount, generalModelMaxClades, taxonIdListFile):
    """
        Generates the list of clades (file) to model for the general model

        @param refSeqDir: directory with reference data as needed for PPS
        @param taxonomyDatabaseFile: taxonomy file in the sqlite3 format
        @param rank: the clades will be considered at this rank
        @param minTotalCount: (see config)
        @param minBpPerSpeciesCount: (see config)
        @param generalModelMaxClades: maximum length of the list of the clades.
        @param taxonIdListFile: file to which the ncbi taxon ids will be stored
        @return: the number of the ncbi taxon ids stored in the file
    """
    cladeBpPairList = refToClades(refSeqDir, taxonomyDatabaseFile, rank)
    rs = ref_seq.RefSequences(refSeqDir, taxonomyDatabaseFile)

    cladeList = []
    count = 0
    for clade, bp in cladeBpPairList:
        if rs.isRefSufficient(int(clade), minTotalCount, minBpPerSpeciesCount):
            cladeList.append(int(clade))
            count += 1
        if count >= generalModelMaxClades:
            break

    out = csv.OutFileBuffer(taxonIdListFile)
    for clade in cladeList:
        out.writeText(str(clade) + '\n')
    out.close()
    rs.close()
    return len(cladeList)


def toRealNames(config, sequences):
    """
        Transforms a PPS file fileName.fas.PP.out that names sequences according to their ids to their real names.
    """
    outIdsPPSFile = str(config.get('inputIdsFastaFile') + '.PP.out')
    outNamesPPSFile = outIdsPPSFile + '.n'
    #os.path.normpath
    print outNamesPPSFile

    try:
        fr = open(os.path.normpath(outIdsPPSFile),'r')
        fw = open(os.path.normpath(outNamesPPSFile),'w')
    except Exception:
        print "Cannot open one of the files:", outIdsPPSFile, "or", outNamesPPSFile
        raise
    else:
        for line in fr:
            if re.match(r'^[0-9]+_[0-9]+[^0-9].*$', line):
                id = re.sub(r'^[0-9]+_([0-9]+)[^0-9].*$',r'\1' , line)
                rest = re.sub(r'^[0-9]+_[0-9]+([^0-9].*)$',r'\1' , line)
                seq = sequences.getSequence(int(id))
                fw.write(seq.name + rest) # seq.scaffold.name
            else:
                fw.write(line)
    finally:
        fr.close()
        fw.close()


def readPPSOutput(sequences, taxonomy, inputFastaIdsPPSFile, overwriteAllPlacements=False):
    """
        Reads the output file of PPS and for each sequence decides:
        if overwriteAllPlacements=True is, then the sequence is placed according to the PPS file regardless of its
        previous placement
        if overwriteAllPlacements=False then if a sequence is placed to a less specific rank, than PPS suggests then
        the sequence is placed according to the PPS file
    """

    infile = str(inputFastaIdsPPSFile + '.out')
    try:
        f = open(os.path.normpath(infile),'r')
    except Exception:
            print "Cannot open file:", infile
            raise
    else:
        #i = 0
        for line in f:
            line = common.noNewLine(line)
            if re.match(r'^[0-9]+_[0-9]+.*[^0-9]+[0-9]+[^0-9]*$', line):
                scaffoldId = int(re.sub(r'^([0-9]+)_[0-9]+.*[^0-9]+[0-9]+[^0-9]*$',r'\1' ,line))
                contigId = int(re.sub(r'^[0-9]+_([0-9]+).*[^0-9]+[0-9]+[^0-9]*$',r'\1' ,line))
                ncbid = int(re.sub(r'^[0-9]+_[0-9]+.*[^0-9]+([0-9]+)[^0-9]*$',r'\1' ,line))
                weight = None # the weight is not yet defined !!!
                if ncbid != 1:
                    #print line, ":", scaffoldId, contigId, ncbid
                    taxPathDictPPS = taxonomy.getPathToRoot(ncbid)
                    if taxPathDictPPS.keys() >= 1:
                        taxPathDictCurrent = sequences.getSequence(contigId).getTaxonomyPath()
                        if taxPathDictCurrent == None:
                            sequences.setTaxonomyPath(contigId, scaffoldId, taxPathDictPPS, weight)#weight = None !!!
                            #i += 1
                        else:
                            if ((overwriteAllPlacements) or (taxPathDictPPS.keys() > taxPathDictCurrent.keys())):
                                sequences.setTaxonomyPathOverride(contigId, scaffoldId, taxPathDictPPS, weight)#weight = None !!!
                                #i += 1
        #print "placed seq by PPS:", i

    finally:
        f.close()


def ppsOut2ppOut(inFile, outFile, taxonomicRanks, databaseFile):
    """
        Transforms a PPS output file into a file in the PP format.

        @param inFile: input file in the PPS format (first column: seq name, last column: ncbi taxon id)
        @param outFile: output file in the PP format
        @param taxonomicRanks: taxonomic ranks (starting from superkingdom)
        @param databaseFile: database file in the sqlite3 format
    """
    taxonomy = Taxonomy(databaseFile, taxonomicRanks)
    outBuff = csv.OutFileBuffer(outFile)
    namesList = csv.getColumnAsList(inFile, entryModifyFunction=None, colNum=0, sep='\t', comment='#')
    valCol = 1
    ncbidsList = csv.getColumnAsList(inFile, entryModifyFunction=None, colNum=valCol, sep='\t', comment='#')

    while True:  # this is not efficient!
        valCol += 1
        tmpList = csv.getColumnAsList(inFile, entryModifyFunction=None, colNum=valCol, sep='\t', comment='#')
        if len(tmpList) == len(namesList):
            ncbidsList = tmpList
        else:
            break

    header = str('#PPS file transformed to PP format, input file: ' + str(inFile) + '\n#ID' + '\t' + 'root')
    for rank in taxonomicRanks:
        header += str('\t' + rank)
    outBuff.writeText(str(header + '\n'))

    for i in range(len(namesList)):
        name = namesList[i]
        ncbid = ncbidsList[i]
        taxPathDict = taxonomy.getPathToRoot(int(ncbid))
        buff = str(name)
        if taxPathDict is None:
            buff += str('\t')
        else:
            buff += str('\t' + 'root')

        for rank in taxonomicRanks:
            if (taxPathDict is not None) and (rank in taxPathDict) and (not taxPathDict[rank].isCopy()):
                buff += str('\t' + taxPathDict[rank].name)
            else:
                buff += '\t'
        outBuff.writeText(str(buff + '\n'))
    outBuff.close()
    taxonomy.close()


class PP2PPSoutParser():

    def __init__(self,taxonomyNcbi, outBuffer):
        self.taxonomy = taxonomyNcbi
        self.out = outBuffer
        self.nameToNcbidDict = dict([])
        self.nameToNcbidDict['Spirochaetes (class)'] = 203691 # was changed to Spirochaetes
        self.nameToNcbidDict['Actinobacteria'] = 201174
        self.nameToNcbidDict['Actinobacteria (class)'] = 1760
        self.nameToNcbidDict['Fusobacteria (class)'] = 203490 #is considered to be: Fusobacteriia

    def parse(self, line):
        tokens = line.split('\t')
        id = str(tokens[0])
        length = len(tokens)
        label = None

        for i in range(len(tokens)-1):
            t = tokens[-(i+1)]
            if t != '':
                if t in self.nameToNcbidDict:
                    label = self.nameToNcbidDict[t]
                else:
                    label = self.taxonomy.getNcbid(t)
                    if label != None:
                        self.nameToNcbidDict[t] = label
                if label == None:
                    print 'CANNOT parse line:', line
                break

        if label != None:
            self.out.writeText(str(id + '\t' + str(label) + '\n'))

    def finalize(self):
        pass

def ppOut2PPSout():
    inFile = '/Users/ivan/Documents/work/binning/data/HumanGut/PP/TS29_scaff.file.0.5.txt'
    outFile = '/Users/ivan/Documents/work/binning/data/HumanGut/PP/TS29_scaff.file.0.5.PPS.txt'
    dbFile = '/Users/ivan/Documents/work/binning/taxonomy/20120828/ncbitax_sqlite.db' #DB
    taxonomy = taxonomy_ncbi.TaxonomyNcbi(dbFile)

    out = csv.OutFileBuffer(outFile)

    csv.forEachLine(inFile, PP2PPSoutParser(taxonomy, out))

    out.close()


def main01():
    #config = Config(open(os.path.normpath('/Users/ivan/Documents/work/binning/tests/CowRumen/01/config.cfg')), 'pPPS')
    #config = Config(open(os.path.normpath('/net/metagenomics/projects/PPSmg/tests/V35/config.cfg')), 'pPPS')
    #configMl = Config2(config, 'MLTreeMap')
    #configPPS = Config2(config, 'PPS')

    #read sequences
    #sequences = Sequences(config)

    #write ids file
    #sequences.writeSequences(config.get('inputIdsFastaFile'))

    #taxonomy = Taxonomy(config.get('databaseFile'), config.get('taxonomicRanks').split(','))

    taxonomicRanks = 'superkingdom,phylum,class,order,family,genus,species'.split(',')
    taxonomy = Taxonomy('/Users/ivan/Documents/work/binning/taxonomy/20120828/ncbitax_sqlite.db', taxonomicRanks)

    #ppsOut2ppOut('D:\\VM\\tmp\\simMC_AMD\\AMD.Arachne.genus', 'D:\\VM\\tmp\\simMC_AMD\\AMD.Arachne.genus.PP.out', taxonomy, config.get('taxonomicRanks').split(','))

    #ppsOut2ppOut('/Users/ivan/Documents/work/binning/data/CowRumen/cowRumenOrderNcbids.txt',
    #             '/Users/ivan/Documents/work/binning/data/CowRumen/cowRumenOrderNcbids.PP.txt', taxonomy, config.get('taxonomicRanks').split(','))

    #ppsOut2ppOut('/net/metagenomics/projects/PPSmg/data/V35/contigsMappedBlast1000LabelsSpecies.txt',
    #             '/net/metagenomics/projects/PPSmg/data/V35/contigsMappedBlast1000LabelsSpecies.PP.txt', taxonomy, config.get('taxonomicRanks').split(','))

    ppsOut2ppOut('/Users/ivan/Documents/work/binning/data/simMC/fromJohannes/contigs.genus.tax',
                 '/Users/ivan/Documents/work/binning/data/simMC/fromJohannes/contigs.genus.PP.tax', taxonomy, taxonomicRanks)


    #readPPSOutput(sequences, taxonomy, config.get('inputIdsFastaFile'))

    #sequences.writePlacements(str(config.get('inputIdsFastaFile') + '.pOUT'), config.get('taxonomicRanks').split(','))

    #toRealNames(config, sequences)
    taxonomy.close()


def collectChildren(taxonomy, ncbid):
    """
        Used in genomesToMask.
    """
    list = taxonomy.childrenNcbids(ncbid)
    if list == None:
        return [ncbid]
    else:
        resultList = []
        for i in list:
            #print 'i', i
            li = collectChildren(taxonomy, i)
            resultList.extend(li)
        return resultList


def genomesToMask():
    rank = 'genus' #which rank will be masked
    fileName = '/Users/ivan/Documents/work/binning/data/simMC/fromJohannes/contigs_genus_ncbids.txt'
    outFile = '/Users/ivan/Documents/work/binning/data/simMC/fromJohannes/genome_genus_masked.txt'
    outFile2 = '/Users/ivan/Documents/work/binning/data/simMC/fromJohannes/genome_ncbids_genus.txt'
    #outFile = '/Users/ivan/Documents/work/binning/data/V35/genome_species_masked.txt' #output file
    #outFile2 = '/Users/ivan/Documents/work/binning/data/V35/genome_ncbids_species.txt' #output file
    #fileName='/Users/ivan/Documents/work/binning/data/V35/genome_ncbids.txt' #list of all genome ncbids
    dbFile = '/Users/ivan/Documents/work/binning/taxonomy/20120828/ncbitax_sqlite.db' #DB
    out = csv.OutFileBuffer(outFile)
    out2 = csv.OutFileBuffer(outFile2)

    genomeNcbids = csv.getColumnAsList(fileName, entryModifyFunction=None, colNum=0, sep=None, comment='#')
    taxonomy = taxonomy_ncbi.TaxonomyNcbi(dbFile)

    maskNcbids = []
    #print len(genomeNcbids), genomeNcbids
    for ncbid in genomeNcbids:
        while taxonomy.getRank(ncbid) != rank:
            ncbid = taxonomy.getParentNcbid(ncbid)
            if int(ncbid) == 1:
                print 'root reached!'
                break
        maskNcbids.append(int(ncbid))

    #print len(Set(maskNcbids)), maskNcbids

    maskSet = set(maskNcbids)
    for i in maskSet:
        out2.writeText(str(str(i) + '\n'))

    resultList = []
    for ncbid in maskSet:
        list = collectChildren(taxonomy, ncbid)
        for i in list:
            out.writeText(str(str(i) + '\n'))
        print ncbid, list

    #print taxonomy.childrenNcbids(818) #997888,818


    out.close()
    out2.close()
    taxonomy.close()



if __name__ == "__main__":
  _trainAccuracyDataTest()
  #genomesToMask()
  #ppOut2PPSout()
  #main01()