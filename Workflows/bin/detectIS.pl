#! /usr/bin/perl 

###### detectIS.pl 
##
##luigi.grassi@astrazeneca.com
##
##USAGE: perl detectIS.pl 
################################################################################

use warnings ;
use strict ;
use diagnostics ;
use FindBin; # locate this script to locate pm files
use lib "$FindBin::Bin/" ; 
use GeneralTools ;
use DetectSpltR ;
use DetectChmR ;
use Data::Dumper;
use List::Util qw( min max );
use Getopt::Long;

my $message_text  = "detectIS usage:\nperl detectIS.pl -h1 name_mate1_gnm.paf -h2 name_mate2_gnm.paf -v1 name_mate1_vir.paf -v2 name_mate2_vir.paf -o name\n-h1:  Aln results of R1 on host genome\n-h2:  Aln results of R2 on host genome\n-v1:  Aln results of R1 on virus/plasmid\n-v2:  Aln results of R2 on virus/plasmid\n-o:  Output prefix\nExtra options -mqual [1] -ovlwind [0.05]\n-mqual: minimum mapping quality [default 1, min:0 max:60]\n-ovlwind:  Ovl tolerability window (fraction of the read length) [default 0.05, min:0 max:1]\n-mspr [default 2]";

my $mqual=1 ;
my $mspr=2 ;
my $ovlwind=0.05 ;
my $GNM1 = '' ;
my $GNM2 = '' ;
my $PLSM1 = '' ;
my $PLSM2 = '' ;
my $OUTPREF = '' ;

GetOptions ('mspr=i' => \$mspr, 'mqual=i' => \$mqual, 'ovlwind=o' => \$ovlwind, 'h1=s' => \$GNM1, 'h2=s' => \$GNM2, 'v1=s'=> \$PLSM1, 'v2=s'=> \$PLSM2, 'o=s'=> \$OUTPREF) ;

( $GNM1 ne '' && $GNM2 ne '' && $PLSM1 ne '' && $PLSM2 ne '' && $OUTPREF ne '') || die $message_text ;

chomp $OUTPREF;
chomp $GNM1;
chomp $GNM2;
chomp $PLSM1;
chomp $PLSM2;

my $OUTPREFMD=$OUTPREF.'.md' ;
my $OUTPREFTXT=$OUTPREF.'.txt' ;
my $OUTPREFSPREADTXT=$OUTPREF.'_SRlist.txt' ;

open(OUT1, ">$OUTPREFMD") || die "ERROR: not possible to open output file $OUTPREFMD\n" ;
open(OUT2, ">$OUTPREFTXT") || die "ERROR: not possible to open output file $OUTPREFTXT\n" ;
open(OUT3, ">$OUTPREFSPREADTXT") || die "ERROR: not possible to open output file $OUTPREFSPREADTXT\n" ;

print OUT1 '\\titleAT'."\n\n" ;
print OUT1 '\\newpage'."\n\n" ;
print OUT1 "# detectIS Results\n" ;
print OUT1 "perl detectIS.pl -h1 $GNM1 -h2 $GNM2 -v1 $PLSM1 -v2 $PLSM2 -o $OUTPREF -mqual $mqual -ovlwind $ovlwind -mspr $mspr" ;
print OUT1 "\n\n----\n\n\n" ;

print OUT2 "IS\tTotSpReads\tR1R2SpReads\tR1SpReads\tR2SpReads\tChimReads\tSingleSplitRead\tInterval\n" ;

print OUT3 "IS\tReadID\tSRType\n" ;

###Detect potential Split Reads in genomic hits
my $splt1=DetectSpltR::detect_split_reads($GNM1, $mqual) ;
my $splt2=DetectSpltR::detect_split_reads($GNM2, $mqual) ;

###Remove fal Split Reads looking using plasmidic reads
my $ishash1=DetectSpltR::verify_split_reads($splt1, $PLSM1, $ovlwind) ;
my $ishash2=DetectSpltR::verify_split_reads($splt2, $PLSM2, $ovlwind) ;

###Merging the results from both reads and filtering by frequency
my ($ca1, $ca2, $ca3)=GeneralTools::merge_split_pairs($ishash1, $ishash2, $mspr) ;
my %ishash=%$ca1; # Key1:IS -> Element: number of fragment/read supporting it
my %isres=%$ca2; # Key1:IS -> Key2:Read_ID -> Element: Read_occurrence
my %iscoord=%$ca3; # Key IS -> Element: @ 0[Plm_chr] 1[Plm_pos] 2[Host_chr] 3[Host_pos]

if ( scalar (keys %ishash) > 0 ) { 
	my ($isfuss, $spreads) = DetectSpltR::verify_spreads_is($GNM1, $GNM2, $PLSM1, $PLSM2, \%isres, \%iscoord) ;
        my %hashres = %$isfuss;
        my %hspreads = %$spreads;
	foreach my $n(sort { $ishash{$b} <=> $ishash{$a} } keys %ishash) {
		if (exists $hashres{$n}) {
			print OUT1 "## $n\n" ;
			my $totsplit= $hashres{$n}{'R1'} + $hashres{$n}{'R2'} + $hashres{$n}{'R1R2'} ;
			print OUT1 "SPLIT READS: $totsplit (R1R2:$hashres{$n}{'R1R2'}; R1:$hashres{$n}{'R1'}; R2:$hashres{$n}{'R2'})\n";
			print OUT1 "CHIMERIC READS: $hashres{$n}{'CHIM'}\n";
			print OUT1 "SINGLE SPLIT READ: $hashres{$n}{'RU'}\n";
			print OUT1 "INTERVAL: $hashres{$n}{'INT'}\n";
			print OUT1 "\n\n----\n\n\n" ;
			print OUT2 "$n\t$totsplit\t$hashres{$n}{'R1R2'}\t$hashres{$n}{'R1'}\t$hashres{$n}{'R2'}\t$hashres{$n}{'CHIM'}\t$hashres{$n}{'RU'}\t$hashres{$n}{'INT'}\n" ;
			foreach my $spr (keys %{$hspreads{$n}}) {
				print OUT3 "$spr\t$n\t$hspreads{$n}{$spr}\n";
			}
		}
	}
}
else {
	print OUT1 "No split read identified! Looking only for chimeric reads\n" ;
	print OUT2 "No split read identified! Looking only for chimeric reads\n" ;

	my ($chm1, $chm2)=DetectChmR::detect_chimeric_reads($GNM1, $GNM2); #Detect genomic hits potentially in a chimeric pair
	DetectChmR::filter_chimeric_reads($PLSM1, $PLSM2, $chm1, $chm2); #Detect potential chimeric pairs using genomic and plasmidic hits
	
	my ($ChimReads2, $len2)=DetectChmR::count_and_collapse_chimeric_reads($PLSM1,$GNM2,$chm2) ; #Count hits dividing chromosomes in intervals 
	my ($ChimReads1, $len1)=DetectChmR::count_and_collapse_chimeric_reads($PLSM2,$GNM1,$chm1) ; #Count hits dividing chromosomes in intervals
	
	my %chmhash=();
		
	my %ch1 = %$ChimReads1;
	my %ch2 = %$ChimReads2;
	
        foreach my $k(keys %ch1) {
                $chmhash{$k}+=$ch1{$k} ;
        	
	}
	foreach my $k(keys %ch2) {
                $chmhash{$k}+=$ch2{$k} ;
        }
	
	foreach my $k(keys %chmhash) {
		if ($chmhash{$k}<2) {
               		delete $chmhash{$k} ;
        	}
	}
	if ( scalar (keys %chmhash) > 0 ) { #Only one IS supported by more than 1 chimeric pair
        	foreach my $n(sort { $chmhash{$b} <=> $chmhash{$a} } keys %chmhash) {
			my @tmp=split(/--/, $n) ;
                        my @tmp1=split(/:/, $tmp[0]) ;
                        my @tmp2=split(/:/, $tmp[1]) ;
                       	my $is1=$tmp1[1]*$len1*5 ;
                        my $is2=$is1+($len1*5);
                       	my $is3=$tmp2[1]*$len1*5 ;
                        my $is4=$is3+($len1*5);
                       	print OUT1 "## $tmp1[0]:$is1-$is2--$tmp2[0]:$is3-$is4\nCHIM:$chmhash{$n}\n" ;
                       	print OUT2 "## $tmp1[0]:$is1-$is2--$tmp2[0]:$is3-$is4\t0\t0\t0\t0\t$chmhash{$n}\t0\t0\n" ;
                        print OUT1 "\n\n----\n\n\n" ;
		}

        }
	else {
        	print OUT1 "No chimeric hit identified!\n ## No IS found!\n" ;
        	print OUT2 "No chimeric hit identified!\n ## No IS found!\n" ;
        }

}

close OUT1 ;
close OUT2 ;
close OUT3 ;
