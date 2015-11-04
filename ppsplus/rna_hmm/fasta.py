
# Copyright (C) 2003, 2004, 2006 by  Thomas Mailund <mailund@birc.au.dk>
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307,
# USA.

"""
A parser for FASTA files.

Copyright (C) 2003, 2004, 2006 by  Thomas Mailund <mailund@birc.au.dk>
"""

class MalformedInput:
    "Exception raised when the input file does not look like a fasta file."
    pass

class FastaRecord:
    "Wrapper around a fasta record."
    def __init__(self, header, sequence):
        "Create a record with the given header and sequence."
        self.header = header
        self.sequence = sequence

    def __str__(self):
        result = ['>'+self.header]
        for i in xrange(0,len(self.sequence),60):
            result.append(self.sequence[i:i+60])
        return '\n'.join(result)
        

def _fasta_itr_from_file(file):
    "Provide an iteration through the fasta records in file."

    h = file.readline().strip()
    if h[0] != '>':
        raise MalformedInput()
    h = h[1:]

    seq = []
    for line in file:
        line = line.strip() # remove newline
        if not len(line):
            continue
        if line[0] == '>':
            yield FastaRecord(h,''.join(seq))

            h = line[1:]
            seq = []
            continue

        seq += [line]

    yield FastaRecord(h,''.join(seq))


def _fasta_itr_from_name(fname):
    "Provide an iteration through the fasta records in the file named fname. "
    f = open(fname)
    for rec in _fasta_itr_from_file(f):
        yield rec
    f.close()


def _fasta_itr(src):
    """Provide an iteration through the fasta records in file `src'.
    
    Here `src' can be either a file object or the name of a file.
    """
    if type(src) == str:
        return _fasta_itr_from_name(src)
    elif type(src) == file:
        return _fasta_itr_from_file(src)
    else:
        raise TypeError

def fasta_get_by_name(itr,name):
    "Return the record in itr with the given name."
    x = name.strip()
    for rec in itr:
        if rec.header.strip() == x:
            return rec
    return None

class fasta_itr:
    "An iterator through a sequence of fasta records."
    def __init__(self,src):
        "Create an iterator through the records in src."
        self.__itr = _fasta_itr(src)

    def __iter__(self):
        return self
    def next(self):
        return self.__itr.next()

    def __getitem__(self,name):
        return fasta_get_by_name(iter(self),name)

class fasta_slice:
    """Provide an iteration through the fasta records in file `src', from
    index `start' to index `stop'.

    Here `src' can be either a file object or the name of a file.
    """
    def __init__(self, src, start, stop):
        """Provide an iteration through the fasta records in file `src', from
        index `start' to index `stop'.

        Here `src' can be either a file object or the name of a file.
        """
        self.__itr = _fasta_itr(src)
        self.__current = 0
        self.__start = start
        self.__stop = stop

    def __iter__(self):
        return self

    def next(self):
        while self.__current < self.__start:
            # skip past first records until we get to `start'
            self.__itr.next()
            self.__current += 1

        if self.__current >= self.__stop:
            # stop after `stop'
            raise StopIteration

        self.__current += 1
        return self.__itr.next()

    def __getitem__(self,name):
        return fasta_get_by_name(iter(self),name)

def get_sequence(src,name):
    "Return the record in src with the given name."
    return fasta_itr(src)[name]


# TESTING...
if __name__ == '__main__':
    import sys
    if len(sys.argv) != 2:
        print "wrong programmer error"
        sys.exit(2)

    print 'iterating through all sequences in input file'
    for rec in fasta_itr(sys.argv[1]):
        print rec
    print

    #print 'input sequences (terminated with ^D)'
    #for rec in fasta_itr(sys.stdin):
    #    print rec
    #print

    print 'iterating through input, from the second sequence'
    for rec in fasta_slice(sys.argv[1], 1, 3):
        print rec
    print

    print 'the sequence for "bar"'
    print fasta_itr(sys.argv[1])["bar"]
    print fasta_slice(sys.argv[1],0,3)["bar"]
    print get_sequence(sys.argv[1],"bar")
    print


    
