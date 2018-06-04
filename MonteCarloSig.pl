#!/usr/bin/perl -w
#GJR 9/5/2013

use strict;

#my $freqFile="Haw7VHaw30.AllDif.txt";
#my $outfile="Haw7VHaw30.sigLoci.txt";
#my $Popfile="/afs/crc.nd.edu/user/g/gragland/RAD_data/Scott_Selection/Bamova/Haw7vHaw30Ids.txt";
die 'need two arguments' unless $ARGV[1];
my $freqFile=shift @ARGV;
my $outfile=shift @ARGV;
my $Popfile=shift @ARGV;

my $GenFile="PJMclines.75.recode.vcf";
#my $hwFile="HweTest/HweCombinedPVal.0.05.txt";
my $nIter=1000;
my $minAsin=0.0;

#hash of ids and pop affiliations
my %Idmap = readFile($Popfile);

#hash of loci passing hw filter
#my %goodLoci;
#open IN, "<$hwFile";
#while (<IN>) {
#  my @vals = split "\t";
#  $goodLoci{$vals[0]}=0;
#}
#close IN;

my %Asins;
my %FreqDifs;
#store point estimates for p value estimation and triage
open IN, "<$freqFile";
while (<IN>) {
  next if $.==1;
  my @vals=split "\t";
  next if $vals[2]==0 or $vals[5]==0;
  #triage
  #next if abs($vals[5]) < $minAsin; 
  my $locus=$vals[0];
  $Asins{$locus}=$vals[5];
  $FreqDifs{$locus}=$vals[2];
}
close IN;

my @keys=keys %Asins;
my $as=@keys;

$"="\t";
open INFILE, "<$GenFile";
open OUTFILE, ">$outfile";
print OUTFILE "Locus\tmedFreq\tmedAsin\tpVal";
my @Ids;
my $n=0;
  my $st=time;
while (<INFILE>) {
  next if /^\#\#/;
  chomp;
  my @vals=split "\t";
  if (/^\#/) {
    @Ids=@vals[9..$#vals];
    next;
  }
  my $locus=$vals[0]."\_".$vals[1];
  #skip bad loci
#  next unless exists $goodLoci{$locus} and exists $FreqDifs{$locus};
  next unless exists $FreqDifs{$locus};
  my @genoProbs1;
  my @genoProbs2;
  my $IdInd=-1;
  for my $info (@vals[9..$#vals]) {
    $IdInd++;
    next if $info =~ m/\.\/\./;
    my @AllInfo=split ":", $info;
    my @PLs=split ",",$AllInfo[4];
#    my @order=order(@PLs);
    #convert from phred scale to probability
    @PLs = map {10**(-$_/10)} @PLs;
    #re-scale based on C*(p1+p2+p3)=1 (i.e., scaled probabilities sum to one)
    @PLs = map {1/(eval join '+', @PLs)*$_} @PLs;
    if (exists $Idmap{$Ids[$IdInd]}) {
      if ($Idmap{$Ids[$IdInd]} == 1) {push @genoProbs1,\@PLs} else {push @genoProbs2,\@PLs}
    }
  }
  #minimum number of individuals per pop
  next if @genoProbs1 < 10 or @genoProbs2 < 10;
  #randomly sample with replacement from pool
  my @RandFreqDifs;
  for (my $i = 0; $i < $nIter; $i++) {
    my $n1= @genoProbs1;
    my $n2= @genoProbs2;
    my @randPop1=randPop(\@genoProbs1,\@genoProbs2,$n1,$n2);
    my @randPop2=randPop(\@genoProbs2,\@genoProbs1,$n2,$n1);
    my @pop1=applyRand(@randPop1);
    my $q1=shift @pop1;
    my @pop2=applyRand(@randPop2);
    my $q2=shift @pop2;
    my $freq1=$q1/(2*@genoProbs1);
    my $freq2=$q2/(2*@genoProbs2);
    push @RandFreqDifs, abs($freq2-$freq1);
  }
  my $FreqDif=$FreqDifs{$locus};
  my $p=1 - getPercentile(abs($FreqDif),@RandFreqDifs);
  print OUTFILE "\n$locus\t$FreqDif\t".$Asins{$locus}."\t$p";
  $n++;
  if ($n==20) {
    my $etime=time-$st;
    my $dum;
  }
} 
close OUTFILE;
close INFILE;

#----------subroutines----------------

sub readFile {
  my $infile=shift;
  my %hashmap;
  open INFILE, "<$infile";
  while (<INFILE>) {
    my $len=length($_);
    next if $_ !~ m/\S/;
    /^(\S+)\s+?(\S+)/;
    my $typLine=$1.$2;
    my $lineLen=length($typLine);
    if ($len > ($lineLen*3)) {
      while (/(\S+)\s+?(\S+)/g) {
	$hashmap{$1}=$2;
      }
    } else {
      /(\S+)\s+?(\S+)/;
      $hashmap{$1}=$2;
    }
  }
  return %hashmap;
}

#reverse order of probsref1 and 2 input to create random pop1 and pop2
sub randPop {
  my $probsRef1=shift;
  my $probsRef2=shift;
  my $n1=shift;
  my $n2=shift;
  my @randPop;
  while (@randPop < @$probsRef1) {
    my $pop=int(rand(2));
    my $ind;
    if ($pop==0) {
      $ind=int(rand($n1));
      push @randPop, $probsRef1->[$ind];
    } else {
      $ind=int(rand($n2));
      push @randPop, $probsRef2->[$ind];
    }
  }
  return(@randPop);
}

#input is random array including point estimate as first value
sub getPercentile {
  my @order = rank(@_);
  my $percentile=$order[0]/@_;
  return($percentile);
}

#rank switches positions of @_ and @sorted in the loops compared to 'order'
sub rank {
  my @sorted = sort {$a <=> $b} @_;
  my @order;
  my %found;
  for my $val (@_) {
    my $count=-1;
    for my $ref (@sorted) {
      $count++;
      next if exists $found{$count};
      if ($val==$ref) {
	push @order, $count;
	$found{$count}=0;
	last;
      }
    }
  }
  return(@order);
}

sub order {
  my @sorted = sort {$a <=> $b} @_;
  my @order;
  my %found;
  for my $val (@sorted) {
    my $count=-1;
    for my $ref (@_) {
      $count++;
      next if exists $found{$count};
      if ($val==$ref) {
	push @order, $count;
	$found{$count}=0;
	last;
      }
    }
  }
  return(@order);
}

sub randGeno {
  my $sampleProb=rand();
  my $genref=shift;
  my @order=order(@$genref);
  my $genoCat=2;
  my $min=$genref->[$order[0]];
  my $mid=$genref->[$order[1]];
  my $max=$genref->[$order[2]];
  #$genoCat=1 if $sampleProb <= (1-$max) and $sampleProb >= $min;
  if ($sampleProb <= (1-$max) and $sampleProb >= $min) {
    $genoCat=1;
    my $dum;
  }
  $genoCat=0 if $sampleProb < $min;
  my $geno = $order[$genoCat];
  return($geno);
}

sub applyRand {
  my $q=0;
  my @pop;
  for my $genoref (@_) {
    my $geno=randGeno($genoref);
    $q=$q+$geno;
    push @pop,$geno;
  }
  return($q,@pop);
}
