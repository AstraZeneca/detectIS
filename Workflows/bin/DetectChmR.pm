package DetectChmR;

use List::Util qw( min max );
use POSIX;
use FindBin; # locate this script
use lib "$FindBin::Bin/";
use GeneralTools ;


################################################################################
#########SUB
################################################################################

################################################################################
#####Passing a 2 genomic PAF file identifies potentially chimeric reads 
sub detect_chimeric_reads  
{
        my $paf1 = $_[0];
        my $paf2 = $_[1];
        my %chimreads1=() ;
        my %chimf1=() ;
        my %chimreads2=() ;
        my %chimf2=() ;
        open(PAF1 , "<", $paf1) || die "ERROR in opening $paf1: please check the file name\n" ;
        while (my $line=<PAF1>) {
                chomp $line;
                if (length ($line) >0) {
                        my @tmp =split("\t" , $line) ;
			if ( ($tmp[10]>=($tmp[1]*0.90)) && ($tmp[11]==60) && ($tmp[9]>($tmp[10]*0.99)) ) {
				$chimreads1{$tmp[0]}++ ;
			}
			else {
				$chimf1{$tmp[0]}++ ;
			}
                }
        }
        close PAF1;
	open(PAF2 , "<", $paf2) || die "ERROR in opening $paf2: please check the file name\n" ;
        while (my $line=<PAF2>) {
                chomp $line;
                if (length ($line) >0) {
                        my @tmp =split("\t" , $line) ;
                        if ( ($tmp[10]>=($tmp[1]*0.90)) && ($tmp[11]==60) && ($tmp[9]>($tmp[10]*0.99)) ) {
                                $chimreads2{$tmp[0]}++ ;
                        }
			else {
                                $chimf2{$tmp[0]}++ ;
                        }
                } 
        }
        close PAF2;
	my %chimrd1=();
	my %chimrd2=();
	foreach my $k(keys %chimreads1) {
		if ( (not exists $chimreads2{$k}) && (not exists $chimf2{$k}) && (not exists $chimf1{$k})) {
			$chimrd1{$k}=$chimreads1{$k};
		}
	}
	foreach my $k(keys %chimreads2) {
                if ( (not exists $chimreads1{$k}) && (not exists $chimf1{$k})  && (not exists $chimf2{$k}) ){
			$chimrd2{$k}=$chimreads2{$k};
		}
        }
	%chimreads1=();
	%chimreads2=();
	%chimf1=();
	%chimf2=();
	return (\%chimrd1, \%chimrd2);
}
################################################################################
	
sub filter_chimeric_reads  
{
        my $paf1 = $_[0];
        my $paf2 = $_[1]; 
        my $potchimr1= $_[2] ;
        my $potchimr2= $_[3] ;
	my %chimreads1=() ;
        my %chimreads2=() ;
	open(PAF1 , "<", $paf1) || die "ERROR in opening $paf1: please check the file name\n" ;
        while (my $line=<PAF1>) {
                chomp $line;
                if (length ($line) >0) {
                        my @tmp =split("\t" , $line) ;
                        if ( (exists $potchimr2->{ $tmp[0] })  && ($tmp[10]>=($tmp[1]*0.90)) && ($tmp[9]>($tmp[10]*0.90)) ) { 
                                $chimreads1{$tmp[0]}++ ;
                        }
			else {
				$chimf1{$tmp[0]}++ ;
			}
                }
        }
	close PAF1;
	open(PAF2 , "<", $paf2) || die "ERROR in opening $paf2: please check the file name\n" ;
        while (my $line=<PAF2>) {
                chomp $line;
                if (length ($line) >0) {
                        my @tmp =split("\t" , $line) ;
                	if ( (exists $potchimr1->{ $tmp[0] })  && ($tmp[10]>=($tmp[1]*0.90)) && ($tmp[9]>($tmp[10]*0.90)) ) { 
                                $chimreads2{$tmp[0]}++ ;

                        }
                        else {
                                $chimf2{$tmp[0]}++ ;
                        }
		}
        }
        close PAF2;
	foreach my $read (keys %$potchimr1) {
		if ( (exists $chimreads2{$read}) && (not exists $chimreads1{$read}) && (not exists $chimf1{$read}) &&  (not exists $chimf2{$read}))  {
			;
		}
		else {
			delete $potchimr1->{ $read } ;
		}
	}
	foreach my $read (keys %$potchimr2) {
                if ( (exists $chimreads1{$read}) && (not exists $chimreads2{$read}) && (not exists $chimf2{$read}) && (not exists $chimf1{$read}) ) {
                	;
		}
		else { 
		       delete $potchimr2->{ $read } ;  
                }
        }  
}


sub count_and_collapse_chimeric_reads  
{
        my $paf1 = $_[0];
        my $paf2 = $_[1];
        my $potchimr= $_[2] ;
	my %rpairs=();
	my @rlen=();
        open(PAF1 , "<", $paf1) || die "ERROR in opening $paf1: please check the file name\n" ;
        while (my $line=<PAF1>) {
                chomp $line;
                if (length ($line) >0) {
                        my @tmp =split("\t" , $line) ;
                        if (exists $potchimr->{ $tmp[0] })  {
                                push (@rlen, $tmp[1]); 
                                push (@{$rpairs{$tmp[0]}}, $tmp[5]); 
                                push (@{$rpairs{$tmp[0]}}, $tmp[6]); 
                                push (@{$rpairs{$tmp[0]}}, $tmp[7]); 
                                push (@{$rpairs{$tmp[0]}}, $tmp[8]); 
                        }
                }
        }
        close PAF1;
        open(PAF2 , "<", $paf2) || die "ERROR in opening $paf2: please check the file name\n" ;
        while (my $line=<PAF2>) {
                chomp $line;
                if (length ($line) >0) {
                        my @tmp =split("\t" , $line) ;
                	if (exists $potchimr->{ $tmp[0] })  {
                        	push (@rlen, $tmp[1]);              
				push (@{$rpairs{$tmp[0]}}, $tmp[5]); 
                        	push (@{$rpairs{$tmp[0]}}, $tmp[6]); 
                        	push (@{$rpairs{$tmp[0]}}, $tmp[7]); 
                                push (@{$rpairs{$tmp[0]}}, $tmp[8]);
			}

		}
        }
        close PAF2;
	my $len=0 ;
	if (@rlen >0 ) {
		$len=GeneralTools::average(\@rlen);
		@rlen=();
	}
	my %comb=();
        foreach my $read (keys %rpairs) {
			my $v1=floor(($rpairs{$read}[2]+(($rpairs{$read}[3]-$rpairs{$read}[2])/2))/($len*5)) ;
			my $v2=floor(($rpairs{$read}[6]+(($rpairs{$read}[7]-$rpairs{$read}[6])/2))/($len*5)) ;
			my $int1=$rpairs{$read}[0].":".$v1 ;
			my $int2=$rpairs{$read}[4].":".$v2;
			$comb{"$int1--$int2"}++ ;
	}
 	return(\%comb, $len);
}
1;
