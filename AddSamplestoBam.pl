#!usr/bin/perl -w
# GJR 4/21/2012
#add read groups to sam file corresponding to barcode ids
# ALSO: changes overall quality scores of 255 to 99 to make compatible with GATK
#mod 2/8/2013 to fix problem grabbing id's that contain the "_" character

use strict;

my $idfile="barcodes.Rm.1.2.fixed.txt";
my $outfile="Rmendax_1_2_f.rg.sam";
my $idcol=1;
my $samfile="Rmendax_1_2_f.merge.sam";
#don't need this for output from bwa, but it also shouldn't affect it one way or the other
my $newqual=99;

my %barcodeIds;
my @barcodeArray;
open INFILE, "<$idfile";
while (<INFILE>) {
  chomp;
  my @vals = split "\t";
  my $id=$vals[$idcol];
  #optional, add regex to remove extraneous info
  #$id =~ s/_sequence//;
  $barcodeIds{$id} = 0;
  push(@barcodeArray,$id);
}
close INFILE;

open OUTFILE, ">$outfile";
open INFILE, "<$samfile";
my $switch=0;
while (<INFILE>) {
  if (/^\@/) {
    print OUTFILE "$_";
    next;
  }
  if ($_ !~ m/^\@/ and $switch == 0) {
    $switch = 1;
    for my $id (@barcodeArray) {
      print OUTFILE "\@RG\tID\:$id\tPL\:illumina\tSM\:$id\n";
    }
    print OUTFILE "\@RG\tID\:none\tPL\:illumina\tSM\:none\n";
  }
  chomp;
  #insert arbitrary quality score to make compatible with GATK
  s/\s255\s/\t$newqual\t/;
  print OUTFILE "$_";
  #modify regex to custom-grab sequence id
  /^\d+_([\w\d\-\.]+?)_/;
  my $id=$1;
  if (exists($barcodeIds{$id})) {
    print OUTFILE "\tRG\:Z\:$id\n";
  } else {
    print OUTFILE "\tRG\:Z\:none\n";
  }
}
close INFILE;
close OUTFILE;