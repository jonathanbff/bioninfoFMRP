#!/usr/bin/perl -w
use strict;

my %rHASH;
my %LOD_HASH;
my %DIST_HASH;

my $inFN = "LDPlot_Data.txt";
open(INFILE, $inFN) || die;
my $hdr = <INFILE>;
while(defined(my $line = <INFILE>)) { 
   chomp($line);
   $line =~ s/M\_//g;
   my @wds = split(/\t/, $line);
   my $posInfo = $wds[0] . "\/" . $wds[1];
   my $LOD = $wds[3];
   my $R2  = $wds[4];
   my $dist = $wds[7];
   $rHASH{$posInfo} = $R2;
   $LOD_HASH{$posInfo} = $LOD;
   $DIST_HASH{$posInfo} = $dist;
}
close(INFILE);
createScatterPlotFile(\%rHASH, \%DIST_HASH, "r2", "r2_scatter.txt");
createScatterPlotFile(\%LOD_HASH, \%DIST_HASH, "LOD", "LOD_scatter.txt");

sub createScatterPlotFile {
   my($hREF, $dREF, $ylab, $outFN) = @_;
   my %h = %$hREF;
   my %d = %$dREF;
   my @k = sort{$h{$b} <=> $h{$a}} keys(%h);
   my $numK = @k;
   open(OUTFILE, ">$outFN") || die("Cannot open $outFN for reading");
   print OUTFILE "DISTANCE\t$ylab\n";
   for(my $i = 0; $i < $numK; $i++) { 
      print OUTFILE "\"$k[$i]\"\t$d{$k[$i]}\t$h{$k[$i]}\n";
   }
   close(OUTFILE);
}

