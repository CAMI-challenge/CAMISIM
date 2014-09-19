#! /usr/bin/perl

# Run the testsuite under Valgrind, to check for memory leakage.
#
# Usage: from testsuite directory:
#    testsuite/valgrind_report.pl
# This assumes you've already compiled the library. To recompile
# from scratch, do 
#    ./driver_report.pl -c
#
# SRE, Fri Mar  2 08:37:48 2007 [Janelia]
# SVN $Id: valgrind_report.pl 231 2008-03-25 14:43:57Z eddys $
require  "getopts.pl";
&Getopts('c');
if ($opt_c) { $do_recompile = 1; }

if ($ENV{'CC'}     ne "") { $CC     = $ENV{'CC'};     } else { $CC       = "gcc"; } 
if ($ENV{'CFLAGS'} ne "") { $CFLAGS = $ENV{'CFLAGS'}; } else { $CFLAGS   = "-g -Wall"; }

printf("Memory leak testing for Easel, using valgrind:\n\n");

if ($do_recompile) {
    print("Recompiling...       ");
    `(cd ..; make clean > /dev/null)`;                      if ($? != 0) { print "[make clean failed]\n"; exit; }
    `(cd ..; ./configure --enable-debugging > /dev/null)`;  if ($? != 0) { print "[configure failed]\n"; exit; }
    `(cd ..; make > /dev/null)`;                            if ($? != 0) { print "[make failed]\n"; exit; }
    print "ok.\n\n";
}

@modules = <../esl_*.c>;
unshift(@modules, "../easel.c");

$nmodules       = 0;
$npresent       = 0;
$ncompiled      = 0;
$nsuccess       = 0;
$nleaking       = 0;
foreach $module (@modules) {
    $module =~ /^\.\.\/(\S+)/; 
    $basecfile = $1;
    $nmodules++;

    # create the eslDMATRIX_TESTDRIVE flag and dmatrix_utest program name from esl_dmatrix.c
    if ($basecfile =~ /^(esl_)?(\S+).c/) { 
	$base     = $2;
	$progname = $base."_utest";
	$base     =~ tr/a-z/A-Z/;
	$flag     = "esl".$base."_TESTDRIVE";
    }

    printf("%-20s ", $basecfile);

    # one way to fail: there isn't a test driver at all
    `grep $flag $module`;
    if ($? != 0) { printf("                   [NO DRIVER]\n");      next; }
    $npresent++;

    `$CC $CFLAGS -I.. -L.. -o $progname -D$flag $module -leasel -lm  >& /dev/null`;
    if ($? != 0) { printf("                   [COMPILE FAILED]\n");       next; };
    $ncompiled++;
    push @proglist, $progname;
    
    $output = `valgrind ./$progname 2>&1`;
    if ($? != 0) { printf("                   [VALGRIND FAILED]\n");       next; };
    $nsuccess++;

    if ($output =~ /malloc\/free: in use at exit: (\S+) bytes in (\S+) blocks/)
    {
	if ($1 > 0) { 
	    $nleaking++;
	    print("[LEAK DETECTED ]\n");
	} else {
	    print("ok.\n");
	}
    } else { print "<< problem parsing valgrind output >>\n"; }                      
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
if ($nleaking == 0) {
    printf("None of %d .c's with running test drivers) show memory leaks\n", $nsuccess);
} else {
    printf("%d of %d .c's with running test drivers) are leaking.\n", $nleaking, $nsuccess);
}

unlink @proglist;


