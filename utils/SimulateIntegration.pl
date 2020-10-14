#!/usr/bin/perl
use strict;
use warnings;
use diagnostics ;

#perl SimulateIntegration.pl
#
################################################

(defined $ARGV[1]) || die "ERROR: please specify Chromosome and Plasmid fasta files\n" ;

my $genome = $ARGV[0];
my $plm = $ARGV[1];


my $chrseq='';
my $chrname='';
my $isgname='';
open(FILE, "<$genome") || die "ERROR: please chek $genome, I am not able to find it!\n" ;
while (my $line=<FILE> ) {
	chomp $line;
	if (length ($line)>0)  {
		if ($line=~/^>/) {
			$chrname=$line;
			my @tmp=split(/ /, $line) ;
                        $isgname=$tmp[0] ;
			$isgname=~s/>//;

		}
		else { 
			$chrseq.=$line ;
		}
	}	
}
close FILE ;


my $lastchr=(length($chrseq)) -1000 ;

my $plmseq='';
my $plname='';
open(FILEP, "<$plm") || die "ERROR: please chek $plm, I am not able to find it!\n" ;
while (my $line=<FILEP> ) {
        chomp $line;
        if (length ($line)>0)  {
                if ($line=~/^>/) {
			my @tmp=split(/ /, $line) ;
			$plname=$tmp[0] ;
			$plname=~s/>//;
		}
		else {
                        $plmseq.=$line ;
                }
        }       
}
close FILEP ;


my $V1=int(rand(length($plmseq)/2));
my $V2= (length($plmseq))-(int(rand(length($plmseq)/2)));
 
my $mcplm1=substr($plmseq, ($V1-1) ) ;
my $mcplm2=substr($plmseq, 0, $V2 ) ;
 
my $mcplm=$mcplm1.$plmseq.$plmseq.$mcplm2 ;

my $G1=1000+int(rand($lastchr)) ;

 

my $G2=$G1+length($mcplm);
print "$chrname $isgname:$G1--$plname:$V1 $plname:$V2--$isgname:$G2\n";

my $mcgnm1=substr($chrseq,0, $G1) ;
my $mcgnm2=substr($chrseq,$G1+length($mcplm)) ;

my $mcgnm=$mcgnm1.$mcplm.$mcgnm2;

print "$mcgnm\n";


