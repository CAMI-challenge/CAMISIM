#! /usr/bin/perl

# Make sure that all drivers compile.
# (Eventually, we should also make sure they run! But that 
# means tracking their command-line arguments.)
#
# In all easel modules, we look for lines like
#   #ifdef eslFOO_(EXAMPLE|TESTDRIVE|REGRESSION|BENCHMARK|STATS)*
# that precede one or more example main()'s, then we compile the
# module with that #define.
#
# Usage: from testsuite directory:
#    ./driver_report.pl
# 
# This assumes you've already compiled the library. To recompile
# from scratch, do 
#    ./driver_report.pl -c
#
# Some regression tests run against the old squid library.
# squid must be installed and compiled; all drivers are linked
# against -lsquid.
# The path to squid's headers and libraries is assumed to be ~/src/squid.
# This location can be overridden by setting SQUIDSRC in the environment.
# 
# SRE, Fri Mar  2 10:01:44 2007 (Janelia)
# SVN $Id: driver_report.pl 326 2009-02-28 15:49:07Z eddys $

require  "getopts.pl";
&Getopts('c');
if ($opt_c) { $do_recompile = 1; }

if ($ENV{'CC'}       ne "") { $CC       = $ENV{'CC'};       } else { $CC        = "gcc"; } 
if ($ENV{'CFLAGS'}   ne "") { $CFLAGS   = $ENV{'CFLAGS'};   } else { $CFLAGS    = "-g -Wall"; }
if ($ENV{'SQUIDSRC'} ne "") { $SQUIDSRC = $ENV{'SQUIDSRC'}; } else { $SQUIDSRC  = "~/src/squid"; }
$progname = "drivertest";


print("Driver code compilation test for Easel:\n");
print("(Compiling with $CC $CFLAGS -L $SQUIDSRC -I $SQUIDSRC)\n\n");


if ($do_recompile) {
    print("Recompiling...      ");
    `(cd ..; make clean > /dev/null)`;                      if ($? != 0) { print "[make clean failed]\n"; exit; }
    `(cd ..; ./configure --enable-debugging > /dev/null)`;  if ($? != 0) { print "[configure failed]\n"; exit; }
    `(cd ..; make > /dev/null)`;                            if ($? != 0) { print "[make failed]\n"; exit; }
    print "ok.\n\n";
}

@modules = <../esl_*.c>;
unshift(@modules, "../easel.c");

$nmodules     = 0;
$ndrivers     = 0;
$nfailures    = 0;
$no_testdriver= 0;
$no_example   = 0;
foreach $module (@modules) {
    $module =~ /^\.\.\/(\S+)/; 
    $basecfile = $1;
    $nmodules++;
    printf("%-20s ", $basecfile);
    


    open(DRIVERFLAGS, qq/grep -E "#ifdef esl.*_(EXAMPLE|TESTDRIVE|REGRESSION|BENCHMARK|STATS)" $module | /) || die;
    $has_example = $has_testdriver = 0;
    $n = 0;
    while (<DRIVERFLAGS>) {
	/^#ifdef (esl\S+_(EXAMPLE|TESTDRIVE|REGRESSION|BENCHMARK|STATS)\d*)/;
	$flag = $1;
	$type = $2;
	if ($saw_flag{$flag}) { next; }

	if ($type eq "EXAMPLE")   { $has_example    = 1; }
	if ($type eq "TESTDRIVE") { $has_testdriver = 1; }

	if ($n == 0) { printf("%-30s ", $flag); }
	else         { printf("%20s %-30s ", "", $flag); }
	$n++;
	$ndrivers++;

        `$CC $CFLAGS -I.. -L.. -L $SQUIDSRC -I $SQUIDSRC -o drivertest -D$flag $module -leasel -lsquid -lm  >& /dev/null`;
 	if ($? != 0) { print("[FAILED]\n"); $nfailures++; }
	else         { print("ok.\n"); }
	$saw_flag{$flag} = 1;
    }
    if ($n == 0) { print "[NO DRIVERS PRESENT]\n"; }
    close DRIVERFLAGS;

    if (! $has_testdriver) { push @notestdrivelist, $module; $no_testdriver++;}
    if (! $has_example)    { push @noexamplelist,   $module; $no_example++;   } 
}

printf("\nOf %d total modules in Easel:\n", $nmodules);

if ($no_example == 0) { 
    printf("   - All %d have at least one example main()\n", $nmodules);
} else {
    printf("   - %d do not have an example main()\n", $no_example);
    foreach $module (@noexamplelist) {
	printf("        %s\n", $module);
    }
}
if ($no_testdriver == 0) { 
    printf("   - All %d have at least one test driver main()\n", $nmodules);
} else {
    printf("   - %d do not have a test driver main()\n", $no_testdriver);
    foreach $module (@notestdrivelist) {
	printf("        %s\n", $module);
    }
}

print "\n";
if ($nfailures == 0) {
    printf("All of the %d driver main()'s compile successfully\n", $ndrivers);
} else {
    printf("%d of the %d driver main()'s fail to compile\n", $nfailures, $ndrivers);
}

unlink $progname;
