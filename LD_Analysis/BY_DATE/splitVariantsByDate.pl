#!/usr/bin/perl -w
use strict;
my %monthHASH = getMonthHash("IDsByMonth.txt");
my %intToMonthHASH = getIntToMonthHash();
my @markerPos = getMarkerPos("markers_CUTOFF_0.01.info");
my $numMarkers = @markerPos;

my $inFN = "sequenceVariants_CUTOFF_0.01_ORIGGenotypes.ped";
open(INFILE, $inFN) || die;
my @arrREFARR;
for(my $i = 0; $i < 12; $i++) { 
   my @a;
   push(@arrREFARR, \@a);
}

while(defined(my $line = <INFILE>)) { 
   chomp($line);
   my @wds = split(/\t/, $line);
   my $p = $wds[0];
   $p =~ s/^P_//g;
   my @wds2 = split(/\-/, $p);
   my $id = $wds2[0];
   my $mo;
   if($id =~ /^DC_/) { 
      $mo = 4;
   } 
   else {
      $mo = $monthHASH{$id};
   }
   if(!defined($mo)) { 
      print "UNDEFINED FOR $id\n";
   }
   else {
      my $cREF = $arrREFARR[$mo - 1];
      my @currARR = @$cREF;
      push(@currARR, $line);
      $arrREFARR[$mo - 1] = \@currARR;
   }
}
close(INFILE);
for(my $i = 0; $i < 12; $i++) { 
   my $currMO = $intToMonthHASH{$i + 1};
   my $cREF = $arrREFARR[$i];
   my @cARR = @$cREF;
   printPEDFiles(\@cARR, \@markerPos, $currMO);
}

#-----------------------------------------------------------------------------
sub printPEDFiles {
   my($aREF, $mREF, $label) = @_;
   my @a = @$aREF;
   my @m = @$mREF;
   my $num = @a;
   my $numMarkers = @m;
   print "There are $num $label\n";
   if($num > 0) { 
      my $FREQ_OUT = $label . "_FREQS.txt";
      my $PED_OUT = "$label.ped";
      my $GENO_OUT = $label . "_GENOTYPES.txt";
   
      open(FREQ_FILE, ">$FREQ_OUT") || die("Cannot open $FREQ_OUT for writing");
      open(PED_FILE, ">$PED_OUT") || die("Cannot open $PED_OUT for writing");
      open(GENO_FILE, ">$GENO_OUT") || die("cannot open $GENO_OUT for writing");

      my %genotypeHASH;
      for(my $i = 0; $i < $num; $i++) { 
         print PED_FILE "$a[$i]\n";
         my $genotype="";
         my $currLine = $a[$i];
         $currLine =~ s/1 1/A/g;
         $currLine =~ s/2 2/C/g;
         $currLine =~ s/3 3/G/g;
         $currLine =~ s/4 4/T/g;
         $currLine =~ s/0 0/X/g;
         my @wds = split(/\t/, $currLine);
         for(my $i = 0; $i < $numMarkers; $i++) { 
            $genotype .= $wds[$i + 6];
         }
         if(!defined($genotypeHASH{$genotype})) { $genotypeHASH{$genotype} = 1; }
         else                                   { $genotypeHASH{$genotype}++;   }
      }
      my @gARR = sort{$genotypeHASH{$b} <=> $genotypeHASH{$a}} keys(%genotypeHASH);
      my $numG = @gARR;
      print GENO_FILE "Genotype\tCount\tFrequency\n";
      for(my $i = 0; $i < $numG; $i++) { 
         print GENO_FILE "$gARR[$i]\t$genotypeHASH{$gARR[$i]}\t" . $genotypeHASH{$gARR[$i]} / $num . "\n";
      }
      close(GENO_FILE);
      close(PED_FILE);

      print FREQ_FILE "MARKER\tA\tC\tG\tT\tUNDEF\n";
      for(my $i = 0; $i < $numMarkers; $i++) { 
         my %h;
         $h{"0 0"} = 0;
         $h{"1 1"} = 0;
         $h{"2 2"} = 0;
         $h{"3 3"} = 0;
         $h{"4 4"} = 0;
         for(my $j = 0; $j < $num; $j++) { 
            my $currS = $a[$j];
            my @wds = split(/\t/, $currS);
            my $currM = $wds[$i + 6];
            $h{$currM}++;
         }
         print FREQ_FILE "$m[$i]\t" . $h{"1 1"}/$num . "\t" . $h{"2 2"}/$num . "\t" . $h{"3 3"}/$num . "\t" . $h{"4 4"}/$num . "\t" . $h{"0 0"}/$num . "\n";
      }
      close(FREQ_FILE);
   }

}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
sub getMarkerPos { 
   my($inFN) = @_;
   my @m;
   open(INFILE, $inFN) || die("Cannot open $inFN for reading");
   while(defined(my $line = <INFILE>))  {
      chomp($line);
      my($id, $currMarker) = split(/\t/, $line);
      push(@m, $currMarker);
   }
   @m = sort {$a <=> $b} @m;
   return(@m);
}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
sub getMonthHash {
   my($inFN) = @_;
   my %h;
   open(INFILE, $inFN) || die("Cannot open $inFN for reading");
   while(defined(my $line = <INFILE>)) { 
      chomp($line);
      my($k, $v) = split(/\t/, $line);
      $h{$k} = $v;
   }
   close(INFILE);
   return(%h);
}
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
sub getIntToMonthHash { 
   my %h;
   $h{1} = "January";
   $h{2} = "February";
   $h{3} = "March";
   $h{4} = "April";
   $h{5} = "May";
   $h{6} = "June";
   $h{7} = "July";
   $h{8} = "August";
   $h{9} = "September";
   $h{10} = "October";
   $h{11} = "November";
   $h{12} = "December";
   return(%h);
} 
#-----------------------------------------------------------------------------
