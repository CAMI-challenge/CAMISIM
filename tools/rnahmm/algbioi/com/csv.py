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
import re
import types

from algbioi.com.common import noNewLine


def getColumnAsList(fileName, entryModifyFunction=None, colNum=0, sep=None, comment='#'):
    """
    Returns a column of a file as a list.

    @param colNum: starts with 0, all entries will be taken from this column (default: 0)
    @param sep: the column separator (default: None ~ white space)
    @param comment: lines that starts with this substring are not considered (default: #)
    """
    lineParser = _ColumnEntryListBuffer(colNum, sep, comment, entryModifyFunction)
    return forEachLine(fileName, lineParser).retVal()


def filterOutLines(inFileName, outFileName, allowedEntriesSet, entryModifyFunction=None, colNum=0, sep=None, comment='#'):
    """
        From the input file filter out lines that does not contain entries from the allowedSet
        in a specific column and store the rest lines to the output file.

        @param allowedEntriesSet: the set of entries that are allowed in the respective column
        @param colNum: starts with 0, entries in this column will be considered (default: 0)
        @param sep: the column separator (default: None ~ white space)
        @param comment: lines that starts with this substring are just copied to the output file
    """
    outBuffer = OutFileBuffer(outFileName)
    lineCondition = _LineConditionFilterOutLines(allowedEntriesSet, colNum, sep, comment, entryModifyFunction)
    parser = _LineFilter(lineCondition, outBuffer)
    forEachLine(inFileName, parser)


def getMapping(inFileName, keyColNum, valColNum, sep=None, comment ='#'):
    """
        Transforms a tab separated file to a dictionary.

        @param inFileName: tab separated file
        @param keyColNum: the number of a column that represents keys (starting from 0)
        @param valColNum: the number of a column that represent values (starting from 0)
    """
    return forEachLine(inFileName, _MappingParser(keyColNum, valColNum, sep, comment)).getDict()


def getMappingTuple(inFileName, keyColNumTuple, valColNumTuple, sep=None, comment='#'):
    """
        Get a mapping where the key and the the value are tuples that can consist of several columns.
        E.g map: -> (col1, col2) -> (col5, col4, col3)

        @param inFileName: csv file
        @param keyColNumTuple: tuple of numbers that define key columns, e.g. (0,1)
        @param valColNumTuple: tuple of numbers that define value columns, e.g. (4,3,2)
        """
    return forEachLine(inFileName,
        _MappingTupleParser(keyColNumTuple, valColNumTuple, sep=sep, comment=comment)).getMapping()


class _MappingTupleParser():
    def __init__(self, keyColNumTuple, valColNumTuple, sep=None, comment='#'):
        self._map = {}
        self._keyColNumTuple = keyColNumTuple
        self._valColNumTuple = valColNumTuple
        self._sep = sep
        self._comment = comment

    def parse(self, line):
        if not isComment(line, self._comment):
            tokens = line.split(self._sep)
            k = []
            v = []
            for i in self._keyColNumTuple:
                k.append(tokens[i])
            for i in self._valColNumTuple:
                v.append(tokens[i])
            assert tuple(k) not in self._map, "Duplicate tuple keys."
            self._map[tuple(k)] = tuple(v)

    def finalize(self):
        pass

    def getMapping(self):
        return self._map


def forEachLine(filePath, parser):
    """
        For each line of the file call the parser, at the end call the finalize method of the parser if it`s defined.
    """
    try:
        f = open(os.path.normpath(filePath), 'r')
    except Exception:
        sys.stderr.write('Cannot open a file for reading: ' + filePath)
        raise
    else:
        try:
            for line in f:
                parser.parse(noNewLine(line))
        except Exception:
            sys.stderr.write('Cannot read from file: ' + filePath)
            raise
        finally:
            f.close()
    try:
        if isinstance(parser.finalize, types.MethodType):
            parser.finalize()
    except Exception:
        pass

    return parser


def isComment(line, comment):
    """
        Is the line a comment.
    """
    if re.match(str('[ \t]*' + comment), line):
        return True
    else:
        return False


class OutFileBuffer():
    """
        To append text to a file.
    """
    def __init__(self, outFilePath, bufferText = False, fileOpenMode='w'):
        self.outFilePath = outFilePath
        self.empty = True
        self.bufferText = bufferText
        self.textBuffer = ''
        self.opened = False
        try:
            self.outFile = open(os.path.normpath(self.outFilePath), fileOpenMode)
            self.opened = True
        except Exception:
            sys.stderr.write('Cannot open a file for writing: ' + outFilePath)
            raise

    def writeText(self, text):
        try:
            if not self.opened: #reopen to append
                self.outFile = open(os.path.normpath(self.outFilePath), 'a')
                self.opened = True
            self.outFile.write(text)
            if self.bufferText:
                self.textBuffer += text
            if self.empty:
                self.empty = False
        except Exception:
            sys.stderr.write('Cannot write to file: ' + self.outFilePath)
            self.close()
            raise

    def getTextBuffer(self):
        return self.textBuffer

    def isEmpty(self):
        return self.empty

    def close(self):
        self.outFile.close()
        self.opened = False


#HELPER FUNCTIONS AND CLASSES DEFINITION


class _MappingParser():
    def __init__(self, keyColNum, valColNum, sep, comment):
        self._dict = dict()
        self._keyColNum = keyColNum
        self._valColNum = valColNum
        self._sep = sep
        self._comment = comment
    def getDict(self):
        return self._dict
    def parse(self, line):
        if not isComment(line, self._comment):
            lineList = line.split(self._sep)
            if len(lineList) > self._keyColNum and len(lineList) > self._valColNum:
                key = lineList[self._keyColNum]
                val = lineList[self._valColNum]
                if key not in self._dict:
                    tmp = []
                    tmp.append(val)
                    self._dict[key] = tmp
                else:
                    self._dict[key].append(val)
            else:
                if len(lineList) > 0:
                    print str('TabSepFileFunctions:_MappingParser: line skipped: ' + line + ' doesn`t have enough entries\n')


class _LineConditionFilterOutLines():
    def __init__(self, allowedEntriesSet, colNum, sep, comment, entryModifyFunction=None):
        self.allowedEntriesSet = allowedEntriesSet
        self.colNum = colNum
        self.sep = sep
        self.comment = comment
        self.entryModifyFunction = entryModifyFunction

    def takeLine(self, line):
        if isComment(line, self.comment):
            return True
        else:
            lineList = line.split(self.sep)
            if len(lineList) > self.colNum:
                entry = lineList[self.colNum]
                if self.entryModifyFunction != None:
                    entry = self.entryModifyFunction(entry)
                if entry in self.allowedEntriesSet:
                    return True
        return False


class _LineFilter():
    def __init__(self, lineCondition, outFileBuffer):
        self.lineCondition = lineCondition
        self.outFileBuffer = outFileBuffer

    def parse(self, line):
        if self.lineCondition.takeLine(line):
            self.outFileBuffer.writeText(str(line + '\n'))

    def finalize(self):
        self.outFileBuffer.close()


def predToDict(predFilePath):
    """  Reads predictions. """
    return forEachLine(predFilePath, _PredParser()).getContigToPredDict()


class _PredParser():
    """
        To parse the prediction file.
    """
    def __init__(self):
        self._dict = dict()

    def parse(self, line):
        if not isComment(line, '#'):
            lineList = line.split()
            if len(lineList) >= 2:
                key = lineList[0]
                val = lineList[len(lineList) - 1]
                if key in self._dict:
                    sys.stderr.write('Consistency:PredParser: the contig "' + key + '" has already been assigned' )
                self._dict[key] = int(val)

    def getContigToPredDict(self):
        return self._dict


class _ColumnEntryListBuffer():
    """
        Parse a line and collects an entry in the specific column.
    """
    def __init__(self, colNum, sep, comment, entryModifyFunction=None):
        self.list = []
        self.colNum = colNum
        self.sep = sep
        self.comment = comment
        self.entryModifyFunction = entryModifyFunction

    def parse(self, line):
        if not isComment(line, self.comment): #the line is not a comment
            lineList = line.split(self.sep)
            if len(lineList) > self.colNum:
                entry = lineList[self.colNum]
                if self.entryModifyFunction is not None:
                    entry = self.entryModifyFunction(entry)
                self.list.append(entry)

    def retVal(self):
        return self.list


if __name__ == "__main__":
    pass