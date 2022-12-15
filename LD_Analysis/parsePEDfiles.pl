#!/usr/bin/perl -w 
use strict;
my $SEPARATOR  = "///";
my %multHASH = getMultipleHash("multipleHits.txt");
my %pedHASH  = getPedHash("sequenceVariants_CUTOFF_0.01.ped", \%multHASH);
createUpdatedPedFile(\%pedHASH, "sequenceVariants_CUTOFF_0.01_ORIGGenotypes.ped");
%pedHASH = updatePedHash(\%pedHASH);
createUpdatedPedFile(\%pedHASH, "sequenceVariants_CUTOFF_0.01_UPDATED.ped");

#==========================================================================
sub createUpdatedPedFile {
   my($pREF, $outFN) = @_;
   my %PEDHASH = %$pREF;
   my @pk = sort(keys(%PEDHASH));
   my $numPED = @pk;
   open(OUTFILE, ">$outFN") || die("Cannot open $outFN for writing");
   for(my $p = 0; $p < $numPED; $p++) { 
      my $currPed = $PEDHASH{$pk[$p]};
      print OUTFILE "$currPed\n";
   }
   close(OUTFILE);
}
#==========================================================================

#==========================================================================
sub updatePedHash {
   my($hREF) = @_;
   my %pHASH = %$hREF;
   my @pk = sort(keys(%pHASH));
   my $numPK = @pk;
   my $firstLine = $pHASH{$pk[0]};
   my @wds = split(/\t/, $firstLine);
   my $numWds = @wds;
   my $numMarkers = $numWds - 6; # 0 - Pedigree; 1 - Individual; 2 - Father; 3 - Mother; 4 - Sex; 5 - Affect Status; 6+ - Markers
   print "THERE ARE $numMarkers Markers\n";

   for(my $m = 0; $m < $numMarkers; $m++) { 
      my %cHash;
      $cHash{"1 1"} = 0;
      $cHash{"2 2"} = 0;
      $cHash{"3 3"} = 0;
      $cHash{"4 4"} = 0;
      for(my $p = 0; $p < $numPK; $p++) { 
         my $currPed = $pHASH{$pk[$p]}; 
         my @pWds = split(/\t/, $currPed);
         my $currGenotype = $pWds[$m + 6];
         if($currGenotype ne "0 0") { 
            $cHash{$currGenotype}++;
         }
      }
      my @countKeys = sort{$cHash{$b} <=> $cHash{$a}} keys(%cHash);
      my $majorAllele = $countKeys[0];
      my $minorAllele = $countKeys[1];

      for(my $p = 0; $p < $numPK; $p++) { 
         my $currPed = $pHASH{$pk[$p]}; 
         my @pWds = split(/\t/, $currPed);
         my $currGenotype = $pWds[$m + 6];
         if($currGenotype ne "0 0") { 
            if(!(($currGenotype eq $majorAllele) || ($currGenotype eq $minorAllele))) { 
               print "More than two genotypes for marker $m\n";
               $pWds[$m + 6] = "0 0";   ## Mark as unknown
               my $updatedPed = join("\t", @pWds);
               $pHASH{$pk[$p]} = $updatedPed;
            }
         }
      }
   }
   return(%pHASH);
}
#==========================================================================

#==========================================================================
sub getPedHash {
   my($inFN, $hREF) = @_;
   my %multiHASH = %$hREF;
   my %pedHASH;
   open(INFILE, $inFN) || die;
   while(defined(my $line = <INFILE>)) { 
      chomp($line);
      my @wds = split(/\t/, $line);
      my $id = $wds[1];
      if($id =~ /MULTIPLE/) {
         my $multi = $multiHASH{$id};
         if(!defined($multi)) { 
            die("Cannot find information for $multi\n");
         }
         my @wds2 = split(/\/\/\//, $multi);
         my $numMulti = @wds2;
         for(my $m = 0; $m < $numMulti; $m++) { 
            my $currID = $wds2[$m];
            my $currPED = "P_$currID";
            $wds[0] = $currPED;
            $wds[1] = $currID;
            my $currLine = join("\t", @wds);
            $pedHASH{$currID} = $currLine;
         } 
      }
      else {
         $pedHASH{$id} = $line;
      }
   }
   close(INFILE);
   return(%pedHASH);
}
#==========================================================================



#P_MT027064-USA-CA       MT027064-USA-CA 0       0       1       0       4 4     4 4     2 2     3 3     2 2     2 2     2 2     2 2
#2 2     3 3     2 2     2 2     2 2     4 4     2 2     1 1     2 2     4 4     2 2     2 2     1 1     1 1     2 2     3 3     2 2
#3 3     4 4     2 2     3 3     4 4     4 4     3 3     3 3     3 3     3 3     3 3     2 2     2 2     3 3     3 3     3 3     3 3
#4 4     3 3     2 2
#P_MT438742-USA-CA       MT438742-USA-CA 0       0       1       0       4 4     4 4     2 2     3 3     2 2     2 2     2 2     2 2
#2 2     3 3     2 2     2 2     2 2     4 4     2 2     1 1     2 2     4 4     2 2     2 2     1 1     1 1     2 2     3 3     4 4
#3 3     4 4     2 2     3 3     4 4     4 4     3 3     3 3     3 3     3 3     3 3     2 2     2 2     3 3     3 3     3 3     3 3
#4 4     3 3     2 2
#P_MT300186-USA-NC       MT300186-USA-NC 0       0       1       0       4 4     4 4     4 4     3 3     2 2     4 4     2 2     2 2
#2 2     3 3     2 2     4 4     2 2     4 4     2 2     1 1     2 2     4 4     2 2     2 2     1 1     3 3     2 2     4 4     2 2
#3 3     4 4     2 2     3 3     4 4     4 4     3 3     3 3     3 3     3 3     3 3     2 2     4 4     3 3     3 3     3 3     0 0
#0 0     0 0     0 0
#P_MT345887-USA-ID       MT345887-USA-ID 0       0       1       0       4 4     2 2     4 4     3 3     2 2     4 4     2 2     2 2
#2 2     3 3     2 2     4 4     2 2     4 4     2 2     1 1     2 2     4 4     2 2     2 2     1 1     3 3     2 2     4 4     2 2
#3 3     4 4     4 4     3 3     4 4     4 4     3 3     3 3     3 3     3 3     3 3     2 2     4 4     3 3     3 3     3 3     3 3
#4 4     3 3     2 2

#==========================================================================
sub getMultipleHash  {
   my($inFN) = @_;
   my %h;
   open(INFILE, $inFN) || die ("Cannot open $inFN for reading");
   while(defined(my $line = <INFILE>)) { 
      chomp($line);
      $line =~ s/\>//g;
      my($k,$v) = split(/\t/, $line);
      $h{$k} = $v;
   }
   close(INFILE);
   return(%h);   
}   
#==========================================================================

#multipleHits.txt";
# MULTIPLE_30     >MT434784-USA-NY///MT434814-USA-NY
# MULTIPLE_31     >MT350263-USA-WA///MT350264-USA-WA///MT350265-USA-WA///MT350266-USA-WA///MT350267-USA-WA///MT350272-USA-WA///MT350273-US
# A-WA///MT350276-USA-WA///MT350278-USA-WA///MT262896-USA-WA///MT262897-USA-WA///MT262898-USA-WA///MT262902-USA-WA///MT262903-USA-WA///MT2
# 62904-USA-WA///MT262906-USA-WA///MT262908-USA-WA///MT262909-USA-WA
# MULTIPLE_32     >MT276329-USA-FL///MT276330-USA-FL

