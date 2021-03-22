#!/usr/bin/perl
use strict;
use warnings;
use diagnostics ;

#perl SimulateIntegration.pl
#
#
#
################################################

(defined $ARGV[2]) || die "ERROR: please specify Chromosome fasta file, plasmid fasta file and number of inserted copies\n" ;

my $genome = $ARGV[0];
my $plm = $ARGV[1];
my $plmcp = $ARGV[2];

($plmcp=~ /^\d+\.?\d*$/ && $plmcp > 0 && $plmcp<=5) || die "ERROR: please number of inserted copies has to be a numeric value [0.1-5]\n" ; 

my $chrseq='';
my $chrname='';
my $isgname='';
open(FILE, "<$genome") || die "ERROR: please check the $genome file\n" ;
while (my $line=<FILE> ) {
	chomp $line;
	if (length ($line)>0)  {
		if ($line=~/^>/) {
			$chrname=$line;
			my @tmp=split(/\s+/, $line) ;
                        $isgname=$tmp[0] ;
			$isgname=~s/>//;

		}
		else { 
			$chrseq.=$line ;
		}
	}	
}
close FILE ;



my $plmseq='';
my $plname='';
open(FILEP, "<$plm") || die "ERROR: please check $plm file\n" ;
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


my $finalenght=int($plmcp*(length($plmseq))) ;
my $V1=int(rand(length($plmseq)));

my $mcplm1=substr($plmseq, $V1) ;
my $mcplm2=$mcplm1.$plmseq.$plmseq.$plmseq.$plmseq.$plmseq ;
my $mcplm=substr($mcplm2, 0, $finalenght ) ;
my $V2;
if ($finalenght <= length($mcplm1)) {
	$V2=$V1+($finalenght-1) ;
}
else {
	$V2=($finalenght-(length($mcplm1))-1) ;
	until ($V2<= length($plmseq)) {
		$V2=$V2-length($plmseq);
	} 
}

my $G1 = length($chrseq);
my $G2 = length($chrseq) ;
until ($G1<length($chrseq) && $G1>0) {
	$G1=int( (length($chrseq)/2) -rand(1000)+rand(1000) )  ;
}
until ($G2<length($chrseq) && $G2>$G1) {
	$G2=int($G1+(rand(1)*length($mcplm)));
}
print ">$isgname $isgname:$G1--$plname:$V1 $plname:$V2--$isgname:$G2 | " ;

my $mcgnm1=substr($chrseq,0, $G1) ;
my $mcgnm2=substr($chrseq,$G2) ;

my $mcgnm=$mcgnm1.$mcplm.$mcgnm2;
my $flen=length($mcgnm);
my $fleni=length($mcplm);
print "Total_len:$flen | Plasmid_len:$fleni\n";
print "$mcgnm\n";


