#!/usr/bin/perl -w
use strict;

my @markerPos = getMarkerPos("markers_CUTOFF_0.01.info");
my $numMarkers = @markerPos;

my $inFN = "sequenceVariants_CUTOFF_0.01_ORIGGenotypes.ped";
open(INFILE, $inFN) || die;
my @cruise;
my @USA;
my @CHINA;
my @TAIWAN;
my @HKG;
my @IND;
my @EUROPE;
my @EAST_ASIA;
my @WEST_ASIA;
my @OTHER;

while(defined(my $line = <INFILE>)) { 
   chomp($line);
   if($line =~ /Cruise/) { 
      push(@cruise, $line);
   }
   else {
      if(($line =~ /USA/) || ($line =~ /KY/)) { 
         push(@USA, $line);
      }
      else {
         if($line =~ /CHN/) { 
            push(@CHINA, $line);
         }
         else {
            if($line =~ /TWN/) { 
               push(@TAIWAN, $line);
            }
            else {
               if($line =~ /HKG/) { 
                  push(@HKG, $line);
               }
               else {
                  if($line =~ /IND/) { 
                     push(@IND, $line);
                  }
                  else {
                     if(($line =~ /ESP/) || ($line =~ /GRC/) || ($line =~ /TUR/) || ($line =~ /SWE/) ||
                        ($line =~ /ITA/) || ($line =~ /FRA/) || ($line =~ /CZE/)) { 
                        push(@EUROPE, $line);
                     }
                     else {
                        if(($line =~ /KOR/) || ($line =~ /LKA/) || ($line =~ /NPL/) || ($line =~ /VNM/) || 
                           ($line =~ /MYS/)) { 
                            push(@EAST_ASIA, $line);
                        }
                        else {
                           if(($line =~ /KAZ/) || ($line =~ /PAK/) || ($line =~ /IRN/) || ($line =~ /ISR/)) { 
                              push(@WEST_ASIA, $line);
                           }
                           else {
                              push(@OTHER, $line);
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
}
close(INFILE);
printPEDFiles(\@cruise, \@markerPos, "CruiseA");
printPEDFiles(\@USA,\@markerPos,  "USA");
printPEDFiles(\@CHINA, \@markerPos, "CHINA");
printPEDFiles(\@TAIWAN, \@markerPos, "TAIWAN");
printPEDFiles(\@HKG, \@markerPos, "HONGKONG");
printPEDFiles(\@IND, \@markerPos, "INDIA");
printPEDFiles(\@EUROPE, \@markerPos, "EUROPE");
printPEDFiles(\@EAST_ASIA, \@markerPos, "EAST_ASIA");
printPEDFiles(\@WEST_ASIA, \@markerPos, "WEST_ASIA");
printPEDFiles(\@OTHER, \@markerPos, "OTHER");

#-----------------------------------------------------------------------------
sub printPEDFiles {
   my($aREF, $mREF, $label) = @_;
   my @a = @$aREF;
   my @m = @$mREF;
   my $num = @a;
   my $numMarkers = @m;
   print "There are $num $label\n";
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
