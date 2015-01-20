#! /usr/bin/perl

# Integrated test of esl-alimerge miniapp.
#
# Usage:     ./esl-alimerge.itest.pl <esl-alimerge binary> <tmpfile prefix>
# Example:   ./esl-alimerge.itest.pl ./esl-alimerge        foo
#
# EPN, Wed Nov 25 11:23:55 2009

$eslalimerge = shift;
$tmppfx      = shift;

if (! -x "$eslalimerge") { die "FAIL: didn't find esl-alimerge binary $eslalimerge"; }

open(ALIFILE, ">$tmppfx.1") || die "FAIL: couldn't open $tmppfx.1 for writing alifile";
print ALIFILE << "EOF";
# STOCKHOLM 1.0
seq1     aaAAAA.AAAA...Cc.cCCCCCC.C..GGGGGgggg
seq2     ..AAAAaAAAAaacCcccCCCCCCcCccGGGGG....
#=GC RF ..AAAA.AAAA...C...CCCCCC.c..GGGGG....
//
EOF
close ALIFILE;

open(ALIFILE, ">$tmppfx.2") || die "FAIL: couldn't open $tmppfx.2 for writing alifile";
print ALIFILE << "EOF";
# STOCKHOLM 1.0
#=GS seq3 sequence 3 is the best!
seq3     AAAAA..AAA.....CCC..CCCC....c..G..GG..GG.
seq4     AAAAAaaAAAaccccC-CccCCCCccccccgGggGGggGGg
#=GC RF AAAAA..AAA.....CCC..CCCC....c..G..GG..GG.
//
EOF
close ALIFILE;

open(LISTFILE, ">$tmppfx.list") || die "FAIL: couldn't open $tmppfx.list for writing list file";
print LISTFILE "$tmppfx.1\n";
print LISTFILE "$tmppfx.2\n";
close LISTFILE;

$output = `$eslalimerge -h`;
if ($? != 0)                                         { die "FAIL: esl-alimerge failed unexpectedly"; }
if ($output !~ /Usage: esl-alimerge/)                { die "FAIL: help output not right"; }

$output = `$eslalimerge --rna $tmppfx.1 $tmppfx.2 2>&1`;
if ($? != 0)                                                                              { die "FAIL: esl-alimerge failed unexpectedly"; }
if ($output !~ /sequence 3 is the best/)                                                  { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /aaAAAA\.A\.\.AAA\.\.\.\.\.Cc\.cCC\.\.CCCC\.\.\.\.C\.\.G\.\.GG\.\.GGgggg/) { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /\.\.AAAA\.AaaAAAaccccC\.\.\.\-CccCCCCccccccgGggGGggGGg\.\.\./)            { die "FAIL: alignments merged incorrectly"; }

$output = `$eslalimerge --dna $tmppfx.1 $tmppfx.2 2>&1`;
if ($? != 0)                                                                              { die "FAIL: esl-alimerge failed unexpectedly"; }
if ($output !~ /sequence 3 is the best/)                                                  { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /aaAAAA\.A\.\.AAA\.\.\.\.\.Cc\.cCC\.\.CCCC\.\.\.\.C\.\.G\.\.GG\.\.GGgggg/) { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /\.\.AAAA\.AaaAAAaccccC\.\.\.\-CccCCCCccccccgGggGGggGGg\.\.\./)            { die "FAIL: alignments merged incorrectly"; }

$output = `$eslalimerge --amino $tmppfx.1 $tmppfx.2 2>&1`;
if ($? != 0)                                                                              { die "FAIL: esl-alimerge failed unexpectedly"; }
if ($output !~ /sequence 3 is the best/)                                                  { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /aaAAAA\.A\.\.AAA\.\.\.\.\.Cc\.cCC\.\.CCCC\.\.\.\.C\.\.G\.\.GG\.\.GGgggg/) { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /\.\.AAAA\.AaaAAAaccccC\.\.\.\-CccCCCCccccccgGggGGggGGg\.\.\./)            { die "FAIL: alignments merged incorrectly"; }

$output = `$eslalimerge --rna --list $tmppfx.list 2>&1`;
if ($? != 0)                                                                              { die "FAIL: esl-alimerge failed unexpectedly"; }
if ($output !~ /sequence 3 is the best/)                                                  { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /aaAAAA\.A\.\.AAA\.\.\.\.\.Cc\.cCC\.\.CCCC\.\.\.\.C\.\.G\.\.GG\.\.GGgggg/) { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /\.\.AAAA\.AaaAAAaccccC\.\.\.\-CccCCCCccccccgGggGGggGGg\.\.\./)            { die "FAIL: alignments merged incorrectly"; }

system("$eslalimerge --rna -o $tmppfx.out $tmppfx.1 $tmppfx.2 > /dev/null");
if ($? != 0)                                                                              { die "FAIL: esl-alimerge failed unexpectedly"; }
$output = `cat $tmppfx.out`;
if ($output !~ /sequence 3 is the best/)                                                  { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /aaAAAA\.A\.\.AAA\.\.\.\.\.Cc\.cCC\.\.CCCC\.\.\.\.C\.\.G\.\.GG\.\.GGgggg/) { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /\.\.AAAA\.AaaAAAaccccC\.\.\.\-CccCCCCccccccgGggGGggGGg\.\.\./)            { die "FAIL: alignments merged incorrectly"; }

system("$eslalimerge --rna -v -o $tmppfx.out $tmppfx.1 $tmppfx.2 > /dev/null");
if ($? != 0)                                                                              { die "FAIL: esl-alimerge failed unexpectedly"; }
# don't worry about checking verbose output printed to stdout
$output = `cat $tmppfx.out`;
if ($output !~ /sequence 3 is the best/)                                                  { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /aaAAAA\.A\.\.AAA\.\.\.\.\.Cc\.cCC\.\.CCCC\.\.\.\.C\.\.G\.\.GG\.\.GGgggg/) { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /\.\.AAAA\.AaaAAAaccccC\.\.\.\-CccCCCCccccccgGggGGggGGg\.\.\./)            { die "FAIL: alignments merged incorrectly"; }
$output = `cat $tmppfx.out`;

$output = `$eslalimerge --rna --outformat pfam $tmppfx.1 $tmppfx.2 2>&1`;
if ($? != 0)                                                                              { die "FAIL: esl-alimerge failed unexpectedly"; }
if ($output !~ /sequence 3 is the best/)                                                  { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /aaAAAA\.A\.\.AAA\.\.\.\.\.Cc\.cCC\.\.CCCC\.\.\.\.C\.\.G\.\.GG\.\.GGgggg/) { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /\.\.AAAA\.AaaAAAaccccC\.\.\.\-CccCCCCccccccgGggGGggGGg\.\.\./)            { die "FAIL: alignments merged incorrectly"; }

$output = `$eslalimerge --rna --outformat a2m $tmppfx.1 $tmppfx.2 2>&1`;
if ($? != 0)                                                 { die "FAIL: esl-alimerge failed unexpectedly"; }
if ($output !~ /aaAAAAAAAACccCCCCCCCGGGGGgggg/)              { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /AAAAAaaAAAaccccC-CccCCCCccccCcgGggGGggGGg/)  { die "FAIL: alignments merged incorrectly"; }

$output = `$eslalimerge --rna --outformat psiblast $tmppfx.1 $tmppfx.2 2>&1`;
if ($? != 0)                                                          { die "FAIL: esl-alimerge failed unexpectedly"; }
if ($output !~ /aaAAAA-A--AAA-----Cc-cCC--CCCC----C--G--GG--GGgggg/)  { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /--AAAA-AaaAAAaccccC----CccCCCCccccCcgGggGGggGGg---/)  { die "FAIL: alignments merged incorrectly"; }

$output = `$eslalimerge --rna --outformat afa $tmppfx.1 $tmppfx.2 2>&1`;
if ($? != 0)                                                                              { die "FAIL: esl-alimerge failed unexpectedly"; }
if ($output !~ /aaAAAA\.A\.\.AAA\.\.\.\.\.Cc\.cCC\.\.CCCC\.\.\.\.C\.\.G\.\.GG\.\.GGgggg/) { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /\.\.AAAA\.AaaAAAaccccC\.\.\.\-CccCCCCccccccgGggGGggGGg\.\.\./)            { die "FAIL: alignments merged incorrectly"; }

$output = `$eslalimerge --rna --rfonly $tmppfx.1 $tmppfx.2 2>&1`;
if ($? != 0)                              { die "FAIL: esl-alimerge failed unexpectedly"; }
if ($output !~ /sequence 3 is the best/)  { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /AAAAAAAACCCCCCCCGGGGG/)   { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /AAAAAAAAC\-CCCCCcGGGGG/)  { die "FAIL: alignments merged incorrectly"; }

# repeat all the same tests (except a2m and psiblast output) but now in small memory mode
$output = `$eslalimerge --small --rna $tmppfx.1 $tmppfx.2 2>&1`;
if ($? != 0)                                                                              { die "FAIL: esl-alimerge failed unexpectedly"; }
if ($output !~ /sequence 3 is the best/)                                                  { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /aaAAAA\.A\.\.AAA\.\.\.\.\.Cc\.cCC\.\.CCCC\.\.\.\.C\.\.G\.\.GG\.\.GGgggg/) { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /\.\.AAAA\.AaaAAAaccccC\.\.\.\-CccCCCCccccccgGggGGggGGg\.\.\./)            { die "FAIL: alignments merged incorrectly"; }

$output = `$eslalimerge --small --dna $tmppfx.1 $tmppfx.2 2>&1`;
if ($? != 0)                                                                              { die "FAIL: esl-alimerge failed unexpectedly"; }
if ($output !~ /sequence 3 is the best/)                                                  { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /aaAAAA\.A\.\.AAA\.\.\.\.\.Cc\.cCC\.\.CCCC\.\.\.\.C\.\.G\.\.GG\.\.GGgggg/) { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /\.\.AAAA\.AaaAAAaccccC\.\.\.\-CccCCCCccccccgGggGGggGGg\.\.\./)            { die "FAIL: alignments merged incorrectly"; }

$output = `$eslalimerge --small --amino $tmppfx.1 $tmppfx.2 2>&1`;
if ($? != 0)                                                                              { die "FAIL: esl-alimerge failed unexpectedly"; }
if ($output !~ /sequence 3 is the best/)                                                  { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /aaAAAA\.A\.\.AAA\.\.\.\.\.Cc\.cCC\.\.CCCC\.\.\.\.C\.\.G\.\.GG\.\.GGgggg/) { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /\.\.AAAA\.AaaAAAaccccC\.\.\.\-CccCCCCccccccgGggGGggGGg\.\.\./)            { die "FAIL: alignments merged incorrectly"; }

$output = `$eslalimerge --small --rna --list $tmppfx.list 2>&1`;
if ($? != 0)                                                                              { die "FAIL: esl-alimerge failed unexpectedly"; }
if ($output !~ /sequence 3 is the best/)                                                  { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /aaAAAA\.A\.\.AAA\.\.\.\.\.Cc\.cCC\.\.CCCC\.\.\.\.C\.\.G\.\.GG\.\.GGgggg/) { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /\.\.AAAA\.AaaAAAaccccC\.\.\.\-CccCCCCccccccgGggGGggGGg\.\.\./)            { die "FAIL: alignments merged incorrectly"; }

system("$eslalimerge --small --rna -o $tmppfx.out $tmppfx.1 $tmppfx.2 > /dev/null");
if ($? != 0)                                                                              { die "FAIL: esl-alimerge failed unexpectedly"; }
$output = `cat $tmppfx.out`;
if ($output !~ /sequence 3 is the best/)                                                  { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /aaAAAA\.A\.\.AAA\.\.\.\.\.Cc\.cCC\.\.CCCC\.\.\.\.C\.\.G\.\.GG\.\.GGgggg/) { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /\.\.AAAA\.AaaAAAaccccC\.\.\.\-CccCCCCccccccgGggGGggGGg\.\.\./)            { die "FAIL: alignments merged incorrectly"; }

system("$eslalimerge --small --rna -v -o $tmppfx.out $tmppfx.1 $tmppfx.2 > /dev/null");
if ($? != 0)                                                                              { die "FAIL: esl-alimerge failed unexpectedly"; }
# don't worry about checking verbose output printed to stdout
$output = `cat $tmppfx.out`;
if ($output !~ /sequence 3 is the best/)                                                  { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /aaAAAA\.A\.\.AAA\.\.\.\.\.Cc\.cCC\.\.CCCC\.\.\.\.C\.\.G\.\.GG\.\.GGgggg/) { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /\.\.AAAA\.AaaAAAaccccC\.\.\.\-CccCCCCccccccgGggGGggGGg\.\.\./)            { die "FAIL: alignments merged incorrectly"; }
$output = `cat $tmppfx.out`;

$output = `$eslalimerge --small --rna --outformat pfam $tmppfx.1 $tmppfx.2 2>&1`;
if ($? != 0)                                                                              { die "FAIL: esl-alimerge failed unexpectedly"; }
if ($output !~ /sequence 3 is the best/)                                                  { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /aaAAAA\.A\.\.AAA\.\.\.\.\.Cc\.cCC\.\.CCCC\.\.\.\.C\.\.G\.\.GG\.\.GGgggg/) { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /\.\.AAAA\.AaaAAAaccccC\.\.\.\-CccCCCCccccccgGggGGggGGg\.\.\./)            { die "FAIL: alignments merged incorrectly"; }

$output = `$eslalimerge --small --rna --rfonly $tmppfx.1 $tmppfx.2 2>&1`;
if ($? != 0)                              { die "FAIL: esl-alimerge failed unexpectedly"; }
if ($output !~ /sequence 3 is the best/)  { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /AAAAAAAACCCCCCCCGGGGG/)   { die "FAIL: alignments merged incorrectly"; }
if ($output !~ /AAAAAAAAC\-CCCCCcGGGGG/)  { die "FAIL: alignments merged incorrectly"; }


print "ok\n"; 
unlink "$tmppfx.1";
unlink "$tmppfx.2";
unlink "$tmppfx.list";
unlink "$tmppfx.out";
exit 0;
