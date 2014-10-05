#!/usr/bin/env python

import zlib
from algbioi.com.common import noNewLine
from algbioi.com.common import removeNonDna

#---------------------------------------------------
class Sequence():
    """Represents one aminoacid or protein sequence/contig

    Attributes:
        id -- unique id of the sequence/contig
        name
        getSeq() -- sequence
        seqBp -- number of bp of the sequence
        scaffold -- to which scaffold this sequence belongs
        getTaxPathDict() -- map: rank -> taxonomy node(ncbid, name, rank)
        ...
    """
    def __init__(self, id, name, seq):
        self.id = id
        self.name = noNewLine(name)
        seq = noNewLine(seq)
        self.seqBp = len(removeNonDna(seq))
        self._seqCompressed = zlib.compress(seq)
        self._taxPathDict = None
        self._placementWeight = None
        self._hash = hash(seq.upper())
        self._candidateTaxPathDictList = []
        self._candidateTaxPathDictWeightsList = []
        self._candidateTaxPathDictSourceList = [] #where does this prediction come from
        self._candidateTaxPathDictTagList = []
        self.scaffold = None
        self._removeNonDna = False

    def getSeq(self):
        seq = zlib.decompress(self._seqCompressed)
        if self._removeNonDna:
            return removeNonDna(seq)
        else:
            return seq

    def setRemoveNonDna(self, removeNonDnaChars):
        self._removeNonDna = removeNonDnaChars

    def setScaffold(self, scaffold):
        assert self.scaffold == None, 'Try to set the scaffold of a sequence twice'
        self.scaffold = scaffold

    #def setTaxonomyPath(self, taxPathDict):
    #    self._taxPathDict = taxPathDict
    #    self._placementWeight = None

    def setTaxonomyPath(self, taxPathDict, weight):
        self._taxPathDict = taxPathDict
        self._placementWeight = weight

    def delTaxonomyPath(self):
        assert self._taxPathDict != None
        self._taxPathDict = None
        self._placementWeight = None

    def setCandidateTaxonomyPath(self, taxPathDict, weight, source=None, tag=None):
        self._candidateTaxPathDictList.append(taxPathDict)
        self._candidateTaxPathDictWeightsList.append(weight)
        self._candidateTaxPathDictSourceList.append(source)
        self._candidateTaxPathDictTagList.append(tag)

    def getTaxonomyPath(self):
        return self._taxPathDict

    def getTaxonomyPathWeight(self):
        return self._placementWeight

    def getCandidateTaxPathDictList(self):
        return self._candidateTaxPathDictList

    def getCandidateTaxPathDictWeightsList(self):
        return self._candidateTaxPathDictWeightsList

    def getCandidateTaxPathSourceList(self):
        return self._candidateTaxPathDictSourceList

    def getCandidateTaxPathTagList(self):
        return self._candidateTaxPathDictTagList

    def getHash(self):
        return self._hash
