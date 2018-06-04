#!/usr/bin/perl -w
#GJR 10/29/2012
# dia  =1, nondia =2
use strict;
#use lib "/afs/crc.nd.edu/user/g/gragland/PerlLibs/lib64/perl5/PDL";
#use lib "/afs/crc.nd.edu/user/g/gragland/PerlLibs/lib/perl5/x86_64-linux-thread-multi";
use PDL;
use PDL::Ufunc qw(pct);



my $GenFile="Sparklef.keep.50.3.30.bi.75.filtered.HWE.DP14.recode.vcf";
my $outfile="floridadif.txt";
#file mapping ids to population (either "1" or "2")
# frequency differences calculated as freq(population2) - freq(population1)
 #can be a subset of the total set of ids; all data not matching those ids will be skipped
my $Popfile="textfiles/florida.txt";
my $nIter=1000;




my %Idmap = readFile($Popfile);

$"="\t";
open INFILE, "<$GenFile";
open OUTFILE, ">$outfile";
print OUTFILE "Locus\tLowerFreq\tMedianFreq\tUpperFreq\tLowerAsin\tMedianAsin\tUpperAsin\tLowerTheta\tMedianTheta\tUpperTheta";
my @Ids;
while (<INFILE>) {
  next if /^\#\#/;
  chomp;
  my @vals=split "\t";
  if (/^\#/) {
    @Ids=@vals[9..$#vals];
    next;
  }
  my $locus=$vals[0]."\_".$vals[1];
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
  my @FreqDif;
  my @AsinDif;
  my @Thetas;
  for (my $i = 0; $i < $nIter; $i++) {
    my @pop1=applyRand(@genoProbs1);
    my $q1=shift @pop1;
    my @pop2=applyRand(@genoProbs2);
    my $q2=shift @pop2;
    my $freq1=$q1/(2*@genoProbs1);
    my $freq2=$q2/(2*@genoProbs2);
    my $theta=Theta(\@pop1,\@pop2);
    push @FreqDif, ($freq2-$freq1);
    if ($freq1 > 1 or $freq2 > 1) {
      my $dum;
    }
    push @AsinDif, (asin($freq2)-asin($freq1));
    push @Thetas, $theta;
  }
  my @Freq_CI=Percentiles(\@FreqDif);
  my @Asin_CI=Percentiles(\@AsinDif);
  my @Theta_CI=Percentiles(\@Thetas);
  print OUTFILE "\n$locus\t@Freq_CI\t@Asin_CI\t@Theta_CI";
} 
close OUTFILE;
close INFILE;


#----------subroutines----------------

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
      }
    }
  }
  return(@order);
}

sub log10 {
  my $n = shift;
  return log($n)/log(10);
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

sub Percentiles {
  my $datref=shift;
  my $datPdl = pdl $datref;
  my $lower = pct($datPdl, 0.025);
  my $med =  pct($datPdl, 0.5);
  my $upper =  pct($datPdl, 0.975);
  return($lower,$med,$upper);
}

sub asin {
  my $x=sqrt(shift);
  my $out=atan2($x,sqrt(1-$x*$x));
  return($out);
}

#Weir-Cockerham Fst
sub Theta {
  my $pop1r=shift;
  my $pop2r=shift;
  my $r=2;
  my $ns1=@$pop1r;
  my $ns2=@$pop2r;
  my $maxSize=(sort($ns1,$ns2))[1];
  my ($s11,$s12,$s21,$s22)=(0,0,0,0);
  for (my $i = 0; $i < $maxSize; $i++) {
    $s11++ if $i < $ns1 and $pop1r->[$i]==0;
    $s12++ if $i < $ns1 and $pop1r->[$i]==1;
    $s21++ if $i < $ns2 and $pop2r->[$i]==0;
    $s22++ if $i < $ns2 and $pop2r->[$i]==1;
  }  
  my $hTilda1=$s12/$ns1;
  my $hTilda2=$s22/$ns2;
  my $pTilda1=($s11*2 + $s12)/(2*$ns1);
  my $pTilda2=($s21*2 + $s22)/(2*$ns2);
  my $nBar=($ns1+$ns2)/2;
  my $C2=(sqrt(($ns1-$nBar)**2+($ns2-$nBar)**2)/$nBar)**2;
  my $nc=$nBar*(1-$C2/$r);
  my $pBar=(($ns1*$pTilda1)/($r*$nBar))+(($ns2*$pTilda2)/($r*$nBar));
  my $s2=(($ns1*($pTilda1-$pBar)**2)/(($r-1)*$nBar))+(($ns2*($pTilda2-$pBar)**2)/(($r-1)*$nBar));
  my $hBar=(($ns1*$hTilda1)/($r*$nBar))+(($ns2*$hTilda2)/($r*$nBar));
  my $a=($nBar/$nc)*($s2-(1/($nBar-1))*($pBar*(1-$pBar)-(($r-1)/$r)*$s2-(1/4)*$hBar));
  my $b=($nBar/($nBar-1))*($pBar*(1-$pBar)-(($r-1)/$r)*$s2-((2*$nBar-1)/(4*$nBar))*$hBar);
  my $c=(1/2)*$hBar;
  my $theta;
  if (($a+$b+$c) == 0) {$theta=-9999} else {$theta=$a/($a+$b+$c)}
  return($theta);
}



