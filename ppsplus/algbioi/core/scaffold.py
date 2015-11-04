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
