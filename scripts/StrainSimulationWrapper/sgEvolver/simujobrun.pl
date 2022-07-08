#!/usr/bin/env perl

# This script automates the process of generating a set of evolved
# sequences using specific evolution parameters for seq-gen and
# sgEvolver, aligning those sequences, and scoring the alignments


use strict;
use POSIX;
use File::Basename;
use FindBin;
use lib "$FindBin::Bin/simulation_dir";
require simujobparams;

if( @ARGV < 3 ){
	die "Usage: simujobrun.pl <FastA input sequence> <GFF annotation file> <random seed> [debug]";
}

my $seq_filename = $ARGV[0];
my $gff = $ARGV[1];

$simujobparams::sgEvolver_random = $ARGV[2]+2;

# check whether we're debugging...
my $debug = 0;
$debug = 1 if @ARGV > 3 && $ARGV[3] eq "debug";

# an array of files that need to be deleted when cleaning up
my @delete_files = ();

# ROADMAP:
# 1) simulate evolution according to parameters
# 2) align genomes using the selected aligner
# 3) score alignments


		# write any command lines executed to a file
open( COMLINES, ">command_lines.txt"  );


my $tools_dir = dirname(__FILE__);

# run evolution with sgEvolver
my $sgevolver_cl = $tools_dir."/sgEvolver --stop-codon-bias=0.98 --ancestral-gff=$gff --accessory-gff=$gff --indel-size=1 --indel-freq=".$simujobparams::indel_rate." --small-ht-freq=$simujobparams::small_ht_rate --small-ht-size=$simujobparams::small_ht_size".
" --large-ht-freq=$simujobparams::large_ht_rate --inversion-freq=$simujobparams::inv_rate --large-ht-min=$simujobparams::large_ht_min --large-ht-max=$simujobparams::large_ht_max".
" --random-seed=$simujobparams::sgEvolver_random".
" --inversion-size=$simujobparams::inv_size $simujobparams::tree_filename $seq_filename $seq_filename $simujobparams::evolved_seqs_name $simujobparams::evolved_seqs_fname";
my $rval = executeCommand( $sgevolver_cl, "sgEvolver.out", "sgEvolver.err" );

		# die if sgEvolver failed
die "Failure in sgEvolver" if( $rval != 0 );

die "Failure in sgEvolver" unless -e "$simujobparams::evolved_seqs_name";
die "Failure in sgEvolver" unless -e "$simujobparams::evolved_seqs_fname";

my $break_cl = $tools_dir."/breakSimulatedGenomeOnAncestralContigs $ARGV[0] evolved_seqs.fas evolved.dat";
$rval = executeCommand( $break_cl, "break.out", "break.err" );
die "Failure in breakSimulatedGenomeOnAncestralContigs" if( $rval != 0 );


deleteFiles( @delete_files );



		# delete evolved sequence data
# push( @delete_files, $simujobparams::evolved_seqs_name );
# push( @delete_files, $simujobparams::evolved_seqs_fname );

deleteFiles( @delete_files );
exit(0);

sub executeCommand {
  my $command = shift;
  my $stdout_file = shift;
  my $stderr_file = shift;
  $command .= " >$stdout_file 2>$stderr_file";
  print "Executing $command\n";
  print COMLINES "$command\n";
  my $rval = system($command);
  `echo "Exited with code $rval" >> $stderr_file` if $rval != 0;
  return $rval;
}

sub deleteFiles {
  if( $debug != 1 ){
    foreach(@_){
      `rm -f $_`;
      print COMLINES "rm -f $_\n";
    }
  }
}


