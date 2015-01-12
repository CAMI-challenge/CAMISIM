#!/usr/bin/env python

import zlib
from algbioi.com.common import noNewLine
from algbioi.com.common import removeNonDna


#---------------------------------------------------
class Scaffold():
    """Represents a set of contigs/sequences where the scaffold sequence can be known

    Attributes:
        id --
        name --
        contigs --
        ...
    """
    def __init__(self, id, name, contig, scaffoldSeq):
        self.id = id
        self.name = name
        self._taxPathDict = None
        self.contigs = []
        self._removeNonDna = False
        if (contig != None):
            self.contigs.append(contig)
        if (scaffoldSeq != None):
            seq = noNewLine(scaffoldSeq)
            self.seqBp = len(removeNonDna(seq))
            self._scaffCompressed = zlib.compress(seq)
            self._hash = hash(seq.upper())
            self._scaffDef = True
        else:
            self._scaffDef = False
            self._hash = None
            self.seqBp = 0


    def setRemoveNonDna(self, removeNonDnaChars):
        self._removeNonDna = removeNonDnaChars


    def getScaff(self):
        if self._scaffDef:
            seq = zlib.decompress(self._scaffCompressed)
            if self._removeNonDna:
                return removeNonDna(seq)
            else:
                return seq
        else:
            return None


    def removeScaffSeq(self):
        self._scaffDef = False
        self._hash = None
        self.seqBp = 0
        self._scaffCompressed = None


    def addNextContig(self, contig):
        self.contigs.append(contig)


    def getSeq(self):
        return self.getScaff()

    def getHash(self):
        assert self._hash != None, 'The hash code of the scaffold is not defined!'
        return self._hash

    def getScaffSeqDef(self):
        """
            @return: True if the sequence for this scaffold is defined.
        """
        return self._scaffDef
