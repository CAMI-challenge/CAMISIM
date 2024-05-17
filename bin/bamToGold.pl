#!/usr/bin/env perl

# abremges@cebitec.uni-bielefeld.de
# asczyrba@cebitec.uni-bielefeld.de

use strict;
use warnings;

use Getopt::Long;

my ($st, $r, $b, $n, $l, $c);

GetOptions("st=s" => \$st, "r=s" => \$r, "b=s" => \$b, "l=i" => \$l, "c=i" => \$c);
die "Usage: $0 -r <reference fasta (indexed with faidx)> -b <sorted bam file> -l <length cutoff> -c <min coverage>\n" unless ($r && $b && $l && $c);
if (!$st)
{
	$st = "samtools"
}
my $pos = -42;
my $seq = "";
my ($start, $stop);
my $seqname = "";
my $previous_sequence_name = "";
#print "samtools mpileup -B -Q 0 -f $r $b | cut -f1,2,3,4 | ";
open(PILEUP, "$st mpileup -B -Q 0 -f $r $b | cut -f1,2,3,4 | ") or die $!;

#      f1                  f2     f3      f4      f5      f6
# NZ_ARQX01000009.1       47061   T       1       ,       I
# [0]                      [1]   [2]     [3]

while(<PILEUP>)
{
	my @base = split;
	$seqname = $base[0];
	if ($seqname ne $previous_sequence_name or ($pos+1) != $base[1])  # contig done
	{
		if (length($seq) >= $l and $previous_sequence_name ne "")
		{
			$stop = $pos;
			print ">$previous_sequence_name\_from\_$start\_to\_$stop\_total\_".length($seq)."\n$seq\n";
		}
		$previous_sequence_name = $seqname;
		$start = $base[1];
		$seq = "";
	}
	$pos = $base[1];
	if ($base[3] >= $c)
	{
		$seq .= $base[2];
	}
}
# catch the last one
if (length($seq) >= $l)
{
	$stop = $pos;
	print ">$previous_sequence_name\_from\_$start\_to\_$stop\_total\_".length($seq)."\n$seq\n";
}
close(PILEUP);

#[fai_load] build FASTA index.
#[mpileup] 1 samples in 1 input files
#<mpileup> Set max per-file depth to 8000
#[fai_load] build FASTA index.
#[fai_build_core] different line length in sequence 'NZ_AGUD01000006.1'.
#[fai_load] fail to open FASTA index.

#My issue got resolved after I indexed the reference assembly by running the command samtools faidx ref.fasta.
