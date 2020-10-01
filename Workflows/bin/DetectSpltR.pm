package DetectSpltR;

use List::Util qw( min max );
use POSIX;

################################################################################
#########SUB
################################################################################

################################################################################
#####Passing a genomic PAF file looks for potential split reads
sub detect_split_reads  
{ 
	my $paf = $_[0];
	my $mq = $_[1];
	my %spltreads=() ;
	my %blklist=() ;
	open(PAF , "<", $paf) || die "ERROR in opening $paf: please check the file name\n" ;
	while (my $line=<PAF>) {
        	chomp $line;
        	if (length ($line) >0) {
			my @tmp =split("\t" , $line) ;
			$blklist{$tmp[0]}++ ;
		}
	}
	close PAF;
	open(PAF , "<", $paf) || die "ERROR in opening $paf: please check the file name\n" ;
        while (my $line=<PAF>) {
                chomp $line;
                if (length ($line) >0) {
                        my @tmp =split("\t" , $line) ;
			if ( ($blklist{$tmp[0]}==1) && ($tmp[10]<($tmp[1])) && ($tmp[11]>=$mq) ) { #&& ($tmp[10]==$tmp[9]) ) { #Filter to consider potentialSplit Reads: Univocally mapped; mapping lenght shorter than read length and MQal==60
				$spltreads{$tmp[0]}[0]=$tmp[2]; #First mapped pos
				$spltreads{$tmp[0]}[1]=$tmp[3];	#Last mapped pos
				$spltreads{$tmp[0]}[2]=$tmp[4]; #STRAND
				$spltreads{$tmp[0]}[3]=$tmp[5];	#CHR
				$spltreads{$tmp[0]}[4]=$tmp[7];	#START
				$spltreads{$tmp[0]}[5]=$tmp[8];	#STOP
			}
		}
 	}
	%blklist=() ;
	return (\%spltreads); #1K: Read name; El: @   
} 
################################################################################

################################################################################
######Passing a Plasmidic PAF file verify the potentially split reads identified in the genome 
sub verify_split_reads
{
        my $splt1 = $_[0];
        my $paf = $_[1];
        my $ovlwin = $_[2];
        my %isarray = ();
        my %vircounts = ();

        open(PAF , "<", $paf) || die "ERROR in opening $paf: please check the file name\n" ;
        while (my $line=<PAF>) {
                chomp $line;
                if (length ($line) >0) {
                        my @tmp =split("\t" , $line) ;
                        if ( (exists $splt1->{ $tmp[0] }) ) {
				$vircounts{$tmp[0]}++ ;
			}
		}
	}
	close PAF ;
	open(PAF , "<", $paf) || die "ERROR in opening $paf: please check the file name\n" ;
        while (my $line=<PAF>) {
                chomp $line;
                if (length ($line) >0) {
                        my @tmp =split("\t" , $line) ;
                        if ( (exists $splt1->{ $tmp[0] })  && $vircounts{$tmp[0]}==1 ) {
                                my $rliminv=$tmp[1]*$ovlwin; #Overlap window shift allowed
                                if ( ($tmp[2]> ($splt1->{ $tmp[0] }[1] - $rliminv) && ($tmp[2]< ($splt1->{ $tmp[0] }[1] + $rliminv))) ) {
                                        if ( ($splt1->{ $tmp[0] }[2] eq "+") &&  ($tmp[4] eq "+") ){
                                                my $is=$tmp[5].":".$tmp[7]."--".$splt1->{ $tmp[0] }[3].":".$splt1->{ $tmp[0] }[5] ;
                                                $isarray{$is}{$tmp[0]}++;
                                        }
                                        elsif ( ($splt1->{ $tmp[0] }[2] eq "+") &&  ($tmp[4] eq "-") ){
                                                my $is=$tmp[5].":".$tmp[8]."--".$splt1->{ $tmp[0] }[3].":".$splt1->{ $tmp[0] }[5] ;
                                        	$isarray{$is}{$tmp[0]}++;
					}
                                        elsif ( ($splt1->{ $tmp[0] }[2] eq "-") &&  ($tmp[4] eq "+") ){
                                                my $is=$tmp[5].":".$tmp[7]."--".$splt1->{ $tmp[0] }[3].":".$splt1->{ $tmp[0] }[4] ;
                                                $isarray{$is}{$tmp[0]}++;
					}
                                        elsif ( ($splt1->{ $tmp[0] }[2] eq "-") &&  ($tmp[4] eq "-") ){
                                                my $is=$tmp[5].":".$tmp[8]."--".$splt1->{ $tmp[0] }[3].":".$splt1->{ $tmp[0] }[4] ;
                                                $isarray{$is}{$tmp[0]}++;
					}
                                }
                                elsif ( (($tmp[3] - $rliminv) < $splt1->{ $tmp[0] }[0] ) &&  (($tmp[3] + $rliminv) > $splt1->{ $tmp[0] }[0]) ) {
                                        if ( ($splt1->{ $tmp[0] }[2] eq "+") &&  ($tmp[4] eq "+") ){
                                        	my $is=$tmp[5].":".$tmp[8]."--".$splt1->{ $tmp[0] }[3].":".$splt1->{ $tmp[0] }[4] ;
                                               	$isarray{$is}{$tmp[0]}++;
					}
                                        elsif ( ($splt1->{ $tmp[0] }[2] eq "+") &&  ($tmp[4] eq "-") ){
                                                my $is=$tmp[5].":".$tmp[7]."--".$splt1->{ $tmp[0] }[3].":".$splt1->{ $tmp[0] }[4] ;
                                                $isarray{$is}{$tmp[0]}++;
					}
                                        elsif ( ($splt1->{ $tmp[0] }[2] eq "-") &&  ($tmp[4] eq "+") ){
                                                my $is=$tmp[5].":".$tmp[8]."--".$splt1->{ $tmp[0] }[3].":".$splt1->{ $tmp[0] }[5] ;
                                                $isarray{$is}{$tmp[0]}++;
					}
                                        elsif ( ($splt1->{ $tmp[0] }[2] eq "-") &&  ($tmp[4] eq "-") ){
                                                my $is=$tmp[5].":".$tmp[7]."--".$splt1->{ $tmp[0] }[3].":".$splt1->{ $tmp[0] }[5] ;
                                                $isarray{$is}{$tmp[0]}++;
					}

                                }
				else {
					delete $splt1->{ $tmp[0] } ;
				}
                        }
                }

        }
	close PAF ;
        $splt1=() ;
        return (\%isarray); #1K: IS; 2K: Supporting read
}


################################################################################
###### 2 genomic PAF file identifies potentially chimeric reads 
sub verify_spreads_is  
{
	my $pafg1 = $_[0]; #Genomic PAF Read1
        my $pafg2 = $_[1]; #Genomic PAF Read2
	my $pafv1 = $_[2]; #Viral PAF Read1
        my $pafv2 = $_[3]; #Viral PAF Read2
	my $ishash= $_[4]; #K1 IS; K2:Read_ID
	my $iscoord= $_[5]; #K: IS ->@ [VC] [VP] [GC] [GP]
	my %readstats = (); #K: Read_ID ->@ [G1 CHR] [G1 CO] [V1 CHR] [V1 CO] [G2 CHR] [G2 CO] [V2 CHR] [V2 CO] 
	my %readcoord = (); #K1:Genomic_Chr; K2:Genomic_Pos -> IS
	my %readtois1 = (); #K:Read_ID -> IS
	my %readtois = (); #K:Read_ID -> IS
	my %pmreads1 = (); #
	my %pmreads2 = ();
	my %chmreads = ();
	my @t1 = ();
	my @t2 = ();
	foreach my  $k (keys %{$ishash}) {
		foreach my $r (keys %{$ishash->{$k}}) {
			$readtois1{$r}{$k}++ ;
		}
	}
	
	foreach my  $k (keys %{$ishash}) {
		my $vl=0;
		$readcoord{$iscoord->{$k}[2]}{$iscoord->{$k}[3]}=$k;		
		foreach my $r (keys %{$ishash->{$k}}) {
			if (scalar (keys %{$readtois1{$r}}) ==1) {
				$vl++;
				$readtois{$r}=$k ;
				$readstats{$r}[0]=0;
				$readstats{$r}[1]=0;
				$readstats{$r}[2]=0;
				$readstats{$r}[3]=0;
				$readstats{$r}[4]=0;
                        	$readstats{$r}[5]=0;
                        	$readstats{$r}[6]=0;
                        	$readstats{$r}[7]=0;

			}
			else {
				delete $ishash->{$k}->{$r};
			}
		}
		if ($vl==0) {
			delete $ishash->{$k}; 
		}
	}
	open(PAF1 , "<", $pafg1) || die "ERROR in opening $pafg1: please check the file name\n" ;
        while (my $line=<PAF1>) {
                chomp $line;
                if (length ($line) >0) {
                        my @tmp =split("\t" , $line) ;
			if (exists $readstats{$tmp[0]}) {
				if ( ($tmp[5] eq $iscoord->{$readtois{$tmp[0]}}[2]) && (($tmp[7]>=($iscoord->{$readtois{$tmp[0]}}[3]-($tmp[1]*5))) && ($tmp[8]<=($iscoord->{$readtois{$tmp[0]}}[3]+($tmp[1]*5)))) ) { 
					$readstats{$tmp[0]}[0]=1;
				
					if ($tmp[7] == $iscoord->{$readtois{$tmp[0]}}[3]) {
                                        	$readstats{$tmp[0]}[1]=1;
						push (@{$t1{$readtois{$tmp[0]}}}, $tmp[8]) ;
                                	}
					elsif ($tmp[8] == $iscoord->{$readtois{$tmp[0]}}[3]) {
						$readstats{$tmp[0]}[1]=1;
						push (@{$t1{$readtois{$tmp[0]}}}, $tmp[7]) ;

					}
				}
			}
			elsif (exists $readcoord{$tmp[5]}) {
				foreach my $k(keys %{$readcoord{$tmp[5]}}) {
					if ( ($tmp[7]>=($k-($tmp[1]*5))) && ($tmp[8]<=($k+($tmp[1]*5))) ) { # && ($tmp[10]>=(0.9*$tmp[1]))) {
						$pmreads1{$tmp[0]}{$readcoord{$tmp[5]}{$k}}=$tmp[4];
					} 
				}
			}
		}
	}
	close PAF1 ;
	open(PAF2 , "<", $pafg2) || die "ERROR in opening $pafg2: please check the file name\n" ;
        while (my $line=<PAF2>) {
                chomp $line;
                if (length ($line) >0) {
                        my @tmp =split("\t" , $line) ;
                        if (exists $readstats{$tmp[0]}) {
                                if ( ($tmp[5] eq $iscoord->{$readtois{$tmp[0]}}[2]) && (($tmp[7]>=($iscoord->{$readtois{$tmp[0]}}[3]-($tmp[1]*5))) && ($tmp[8]<=($iscoord->{$readtois{$tmp[0]}}[3]+($tmp[1]*5)))) ) {        
					$readstats{$tmp[0]}[4]=1;
                                
                                	if ($tmp[7] == $iscoord->{$readtois{$tmp[0]}}[3]) {
                                        	$readstats{$tmp[0]}[5]=1;
						push (@{$t1{$readtois{$tmp[0]}}}, $tmp[8]) ;
                                	}
					elsif ($tmp[8] == $iscoord->{$readtois{$tmp[0]}}[3]) {
						$readstats{$tmp[0]}[5]=1;
						push (@{$t1{$readtois{$tmp[0]}}}, $tmp[7]) ;
					}		
                        	}
			}
                        elsif (exists $readcoord{$tmp[5]}) {
                                foreach my $k(keys %{$readcoord{$tmp[5]}}) {
                                        if ( ($tmp[7]>=($k-($tmp[1]*5))) && ($tmp[8]<=($k+($tmp[1]*5))) ) { #&& ($tmp[10]>=(0.9*$tmp[1]))) {
                                                $pmreads2{$tmp[0]}{$readcoord{$tmp[5]}{$k}}=$tmp[4];
                                        }
                                }
                        }
                }
        }
        close PAF2 ;
	open(PAF3 , "<", $pafv1) || die "ERROR in opening $pafv1: please check the file name\n" ;
        while (my $line=<PAF3>) {
                chomp $line;
                if (length ($line) >0) {
                        my @tmp =split("\t" , $line) ;
                        if (exists $readstats{$tmp[0]}) {
                        	if ( ($tmp[5] eq $iscoord->{$readtois{$tmp[0]}}[0]) && (($tmp[7]>=($iscoord->{$readtois{$tmp[0]}}[1]-($tmp[1]*5))) && ($tmp[8]<=($iscoord->{$readtois{$tmp[0]}}[1]+($tmp[1]*5)))) ) {
                                        $readstats{$tmp[0]}[2]=1;
                                	if ($tmp[7] == $iscoord->{$readtois{$tmp[0]}}[1]) {
                                        	$readstats{$tmp[0]}[3]=1;
                                		push (@{$t2{$readtois{$tmp[0]}}}, $tmp[8]) ;
					}
					elsif  ($tmp[8] == $iscoord->{$readtois{$tmp[0]}}[1]) {
						$readstats{$tmp[0]}[3]=1;
						push (@{$t2{$readtois{$tmp[0]}}}, $tmp[7]) ;
					}
				}

			}
                        elsif (exists $pmreads2{$tmp[0]}) {
				foreach my $is(keys %{$pmreads2{$tmp[0]}} ) {
					if ( ($tmp[5] eq $iscoord->{$is}[0]) && ($tmp[7]>=($iscoord->{$is}[1]-($tmp[1]*5))) && ($tmp[8]<=($iscoord->{$is}[1]+($tmp[1]*5))) && ($tmp[4] ne $pmreads2{$tmp[0]}{$is}) ) {
                                                $chmreads{$is}{$tmp[0]}++;
                                  
                        		}
                		}
        		}
		}
	}
        close PAF3 ;
	open(PAF4 , "<", $pafv2) || die "ERROR in opening $pafv2: please check the file name\n" ;
        while (my $line=<PAF4>) {
                chomp $line;
                if (length ($line) >0) {
                        my @tmp =split("\t" , $line) ;
                        if (exists $readstats{$tmp[0]}) {
                        	if ( ($tmp[5] eq $iscoord->{$readtois{$tmp[0]}}[0]) && (($tmp[7]>=($iscoord->{$readtois{$tmp[0]}}[1]-($tmp[1]*5))) && ($tmp[8]<=($iscoord->{$readtois{$tmp[0]}}[1]+($tmp[1]*5)))) ) {
                                        $readstats{$tmp[0]}[6]=1;
                                	if ($tmp[7] == $iscoord->{$readtois{$tmp[0]}}[1]) {
                                        	$readstats{$tmp[0]}[7]=1;
                                		push(@{$t2{$readtois{$tmp[0]}}}, $tmp[8]) ;
					}
					elsif ($tmp[8] == $iscoord->{$readtois{$tmp[0]}}[1])  {
						$readstats{$tmp[0]}[7]=1;
                                                push(@{$t2{$readtois{$tmp[0]}}}, $tmp[7]) ;
                                        }

				}
                        
			}
                        elsif (exists $pmreads1{$tmp[0]}) {
				foreach my $is(keys %{$pmreads1{$tmp[0]}} ) {
					if ( ($tmp[5] eq $iscoord->{$is}[0]) && ($tmp[7]>=($iscoord->{$is}[1]-($tmp[1]*5))) && ($tmp[8]<=($iscoord->{$is}[1]+($tmp[1]*5))) && ($tmp[4] ne $pmreads1{$tmp[0]}{$is}) ) {
                                                $chmreads{$is}{$tmp[0]}++;

                        		}
                		}
        		}
		}
	}
        close PAF4 ;
	my $verifis=();
	foreach my  $k (keys %{$ishash}) {
		my $chim= scalar keys(%{$chmreads{$k}}) ;
		my $R1R2=0;
		my $R1=0;
		my $R2=0;
		my $R0=0;
		foreach my $r (keys %{$ishash->{$k}}) {
			if ($readstats{$r}[0]==1 && $readstats{$r}[1]==1 && $readstats{$r}[2]==1 && $readstats{$r}[3]==1 && $readstats{$r}[4]==1 && $readstats{$r}[5]==1 && $readstats{$r}[6]==1 && $readstats{$r}[7]==1) {
				$R1R2++ ;
			}
			elsif ( ($readstats{$r}[0]==1 && $readstats{$r}[1]==1 && $readstats{$r}[2]==1 && $readstats{$r}[3]==1) && ( ($readstats{$r}[4]==1 && $readstats{$r}[6]==0) ||  ($readstats{$r}[4]==0 && $readstats{$r}[6]==1) ) ) {
                                $R1++ ;
                        }
			elsif ( ($readstats{$r}[4]==1 && $readstats{$r}[5]==1 && $readstats{$r}[6]==1 && $readstats{$r}[7]==1) && ( ($readstats{$r}[0]==1 && $readstats{$r}[2]==0) || ($readstats{$r}[0]==0 && $readstats{$r}[2]==1) ) ) {
				$R2++ ;
			}
			else {
				$R0++ ;
			}
		}
		if (exists $t2{$k} && $t1{$k} ) { 
			$verifis{$k}{"R1R2"}=$R1R2 ;
			$verifis{$k}{"R1"}=$R1 ;
			$verifis{$k}{"R2"}=$R2 ;
			$verifis{$k}{"CHIM"}=$chim ;
			$verifis{$k}{"RU"}=$R0 ;
			my @isc=split(/--/, $k) ;
                	my @iscp=split(/:/, $isc[0]) ;
                	my @iscg=split(/:/, $isc[1]) ;
			my $intp=max(@{$t2{$k}}) ;
			if ($intp<=$iscp[1]) {
				$intp=min(@{$t2{$k}}) ;
			}
			my $intg=max(@{$t1{$k}}) ;
                	if ($intg <= $iscg[1]) {
                        	$intg=min(@{$t1{$k}}) ;
                	}
			my $interval="$iscp[0]:$iscp[1]-$intp--$iscg[0]:$iscg[1]-$intg";
			$verifis{$k}{"INT"}=$interval;
		}
	
	}
	return(\%verifis) ;

}


################################################################################

1;
