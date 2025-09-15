#!/usr/bin/perl
use strict;

my $file = $ARGV[0];
my $out_file = $ARGV[1];
open(IN,$file) || die "Incorrect file $file. Exiting...\n";
open(my $OUT,'>',$out_file) || die "Incorrect file $out_file. Exiting...\n";

my ($seq, $name)=('','');
while(<IN>){
  chomp;
  my $line = $_;
  $seq.= uc($line) if(eof(IN));
  if (/\>(\S+)/ || eof(IN)){
    if($seq ne ''){
      my @seqgaps = split(/[N]{1,}/, $seq);
      if($#seqgaps > 0){
        my $ctgcount=0;
        foreach my $ctgseq (@seqgaps){
          $ctgcount++;
          my @fasta_names = split(/\s/,$name);
          my $fasta_name = "";
          my $i = 0;
          foreach my $spl (@fasta_names){
            $fasta_name .= $spl;
            if ($ctgcount > 1 and $i == 0){
                $fasta_name .= "_".$ctgcount;
            }
            $i++;
            $fasta_name .= " ";
          }
          print $OUT "$fasta_name (size=".length($ctgseq).")\n$ctgseq\n";
        }
      }else{
        print $OUT "$name\n$seq\n";
      }
    }
    $seq='';
    $name = $_;
  }else{
    $seq.=uc($line);
  }
}
