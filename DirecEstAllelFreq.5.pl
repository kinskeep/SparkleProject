#!/usr/bin/perl -w
#GJR 3/5/2014
#estimate allele frequencies of groups of individuals whose ids appear in a specified vcf file
#uses the direct method of freq = ( sum (i to n) [sum [(0,1,2) * Gi] ] ) / 2n
# where n = number of individuals in group and Gi is the ith vector of genotype probabilities
# gives the estimated frequency of the minor allele

use strict;

my $GenFile="Sparklef.keep.50.3.30.bi.75.filtered.HWE.DP14.recode.vcf";
my $outfile="DirectFreqEsts.Sparklef.5.txt";
#which column (0-indexed) of GroupFile to use
my $useCol=1;
my $GroupFile="everything.txt";
#minimum number of individuals in a category to estimate freq (below minN generated 'NA')
my $minN=5;

my %Idmap;
my @Cats;
my %CatHash;
open INFILE, "<$GroupFile";
while (<INFILE>) {
  chomp;
  my @vals=split "\t";
  my $pheno=$vals[$useCol];
  $Idmap{$vals[0]}=$pheno;
  push @Cats, $pheno unless exists $CatHash{$pheno};
  $CatHash{$pheno}=0;
}
close INFILE;




$"="\t";
open INFILE, "<$GenFile";
open OUT, ">$outfile";
print OUT "Locus\t@Cats";
my @Ids;
my @nAlleles=(0,1,2);
while (<INFILE>) {
  next if /^\#\#/;
  chomp;
  my @vals=split "\t";
  if (/^\#/) {
    @Ids=@vals[9..$#vals];
    next;
  }
  my $locus=$vals[0]."\_".$vals[1];
  my %CatSums=%CatHash;
  my %CatNs=%CatHash;
  my $IdInd=-1;
  for my $info (@vals[9..$#vals]) {
    $IdInd++;
    next if $info =~ m/\.\/\./;
    my @AllInfo=split ":", $info;
    my @PLs=split ",",$AllInfo[4];
    #convert from phred scale to probability
    @PLs = map {10**(-$_/10)} @PLs;
    #re-scale based on C*(p1+p2+p3)=1 (i.e., scaled probabilities sum to one)
    @PLs = map {1/(eval join '+', @PLs)*$_} @PLs;
    if (my $pheno = $Idmap{$Ids[$IdInd]}) {
      my @products=  map {$PLs[$_] * $nAlleles[$_]} 0..$#PLs;
      my $sum=eval join '+', @products;
      $CatSums{$pheno}+=$sum;
      $CatNs{$pheno}++;
    }
  }
  print OUT "\n$locus";
  for my $pheno (@Cats) {
    if ($CatNs{$pheno} < $minN) {
      print OUT "\tNA";
    } else {
      my $freq=$CatSums{$pheno}/(2*$CatNs{$pheno});
      print OUT "\t$freq";
    }
  }
} 
close OUT;
close INFILE;


#----------subroutines----------------


#not used in this script, but could be used to estimate freqs for multiple sets of groups


sub readFile {
  my $infile=shift;
  my %hashmap;
  my @cats;
  open INFILE, "<$infile";
  while (<INFILE>) {
    chomp;
    my @vals=split "\t";
    if ($.==1) {
      push @cats, {} for @vals[1..$#vals];
    }
    my @phens=@vals[1..$#vals];
    $hashmap{$vals[0]}=\@phens;
    my $count=0;
    for my $cat (@phens) {
      $cats[$count]->{$cat}=0;
      $count++;
    }
  }
  return (\%hashmap,\@cats);
}



