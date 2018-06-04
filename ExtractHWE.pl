#!/usr/bin/perl -w
#GJR 11/18/2013



my $vcfFile1=$ARGV[0];
my $outfile=$ARGV[1];
my $thresh=0.01;
die "specify 1) first .vcf file 2) second .vcf file 3) outfile" unless @ARGV ==2;

open OUT, ">$outfile";
print OUT "Contig\tpVal\n";
open IN, "<$vcfFile1";
while (<IN>) {
  next if /^\#/;
  /HW\=(\d+)/;
  my $hw=$1;
  #convert from phred scale to probability
  $hw=10**(-$hw/10);
  my @vals= split "\t";
  my $loc=$vals[0]."_".$vals[1];
  print OUT "$loc\t$hw\n" if $hw > $thresh;
}
close IN;
close OUT;
