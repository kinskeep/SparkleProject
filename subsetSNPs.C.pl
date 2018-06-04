#!/usr/bin/perl -w
#GJR 11/18/2013
#MMD 12/4/2013

my $SNPfile=$ARGV[0];
my $FullVcfFile=$ARGV[1]; 
my $outfile=$ARGV[2];
my %goodSNPs;
die "specify 1) population .vcf file 2) full .vcf file 3) outfile" unless @ARGV ==3;

#------get list of SNPs to pull out------
open OUT, ">$outfile";
open IN, "<$SNPfile";
while (<IN>) {
  chomp;
  my $loc=$_;
  $goodSNPs{$loc}=1;
}
close IN;

#------filter SNPs in HWE out of full .vcf------
open IN, "<$FullVcfFile";
while (<IN>) {
  if (/^\#/) {
    print OUT $_;
    next;
    }
  my @vals= split "\t";
  my $loc=$vals[0]."_".$vals[1];
  print OUT $_ if $goodSNPs{$loc};
}
close IN;
close OUT;
