#! /usr/bin/perl

# Measures testsuite coverage (as percentage of source lines),
# using gcov.
#
# Usage: from testsuite directory:
#    ./coverage_report.pl
#
# This assumes you've already compiled the library. To recompile
# from scratch, do 
#    ./coverage_report.pl -c
#
# It assumes you have 'sloccount' installed, so it can count 
# ANSI C lines in files with no test driver. If you don't, use
#    ./coverage_report.pl -s
#
# SRE, Thu Mar  1 19:22:57 2007 (Janelia)
# SVN $Id: coverage_report.pl 231 2008-03-25 14:43:57Z eddys $
require  "getopts.pl";
$have_sloccount = 1;
&Getopts('cs');
if ($opt_c) { $do_recompile     = 1; }
if ($opt_s) { $have_sloccount   = 0; }

if ($ENV{'CC'}     ne "") { $CC     = $ENV{'CC'};  } else { $CC       = "gcc"; } 
$CFLAGS = "-g -Wall -fprofile-arcs -ftest-coverage";

printf("Code coverage test for Easel, using gcov:\n\n");

if ($do_recompile) { 
    print("Recompiling...      ");
    `(cd ..; make clean > /dev/null)`;                 if ($? != 0) { print "[make clean failed]\n"; exit; }
    `(cd ..; ./configure --enable-gcov > /dev/null)`;  if ($? != 0) { print "[configure failed]\n"; exit; }
    `(cd ..; make > /dev/null)`;                       if ($? != 0) { print "[make failed]\n"; exit; }
    print "ok.\n\n";
}


@modules = <../esl_*.c>;
unshift(@modules, "../easel.c");

$nmodules       = 0;
$npresent       = 0;
$ncompiled      = 0;
$nsuccess       = 0;
$nlines         = 0;
$nlines_covered = 0;
foreach $module (@modules) {
    $module =~ /^\.\.\/(\S+)/; 
    $basecfile = $1;
    $nmodules++;

    # create the eslDMATRIX_TESTDRIVE flag and dmatrix_utest program name from esl_dmatrix.c
    if ($module =~ /^\.\.\/(esl_)?(\S+).c/) { 
	$base     = $2;
	$progname = $base."_utest";
	$base     =~ tr/a-z/A-Z/;
	$flag     = "esl".$base."_TESTDRIVE";
    }

    printf("%-20s ", $basecfile);

    # one way to fail: there isn't a test driver at all
    `grep $flag $module`;
    if ($? != 0) { printf("%6.2f%% coverage  [driver ABSENT]\n", 0);  push @nodriverlist, $module; next; }
    $npresent++;

    `$CC $CFLAGS -I.. -L.. -o $progname -D$flag $module -leasel -lm  >& /dev/null`;
    if ($? != 0) { printf("%6.2f%% coverage   [compilation FAILED]\n", 0);       next; };
    $ncompiled++;
    
    `./$progname >& /dev/null`;
    if ($? != 0) { printf("%6.2f%% coverage   [test driver FAILED ]\n", 0);       next; };
    $nsuccess++;

    $output = `gcov $module`;
    if ($output =~ /File.*$module.*\nLines executed:\s*(\d+\.\d+)% of\s+(\d+)/) {
	$pct_cvg        = $1;
	$nlines         += $2;
	$nlines_covered += $1*$2/100;
	printf("%6.2f%% coverage\n", $pct_cvg);
    }
    else {die "failed to parse gcov output";}
}

if ($have_sloccount) {
    foreach $badmodule (@nodriverlist) {
	$output = `sloccount $badmodule`;
	if ($output =~ /ansic:\s+(\d+)/)  { $nlines_nodrivers += $1; }
	else { die("failed to parse sloccount output"); }
    }
}

printf("\nOf %d total modules in Easel:\n", $nmodules);
if ($npresent != $nmodules) {
    printf("   - %d have test drivers, %d do not\n", $npresent, $nmodules-$npresent);
} else {
    printf("   - All %d have test drivers\n", $npresent);
}
if ($ncompiled != $npresent) {
    printf("   - %d compiled, %d did not\n", $ncompiled, $npresent-$ncompiled);
} else {
    printf("   - All %d compiled\n", $ncompiled);
}
if ($nsuccess != $ncompiled) {
    printf("   - %d ran successfully, %d did not\n", $nsuccess, $ncompiled-$nsuccess);
} else {
    printf("   - All %d ran successfully\n", $nsuccess);
}

print "\n";
  printf("Total coverage (of .c's with test drivers): %.2f%%\n", 100.*$nlines_covered / $nlines);
if ($have_sloccount) {
    printf("Total coverage (including .c files without drivers yet): %.2f%%\n", 100.*$nlines_covered / ($nlines+$nlines_nodrivers));
}


