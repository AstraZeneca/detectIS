package GeneralTools;

use List::Util qw( min max );
use POSIX;

################################################################################
#########SUB
################################################################################

sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty arrayn");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}

sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

sub merge_split_pairs  
{
        my $sp1 = $_[0];
        my $sp2 = $_[1];
       	my %ish1=(); # Key1:IS -> Element: number of fragment/read supporting it
	my %isr=(); # Key1:IS -> Key2:Read_ID -> Element: Read_occurrence
	my %isc=(); # Key IS -> Element: @ 0[Plm_chr] 1[Plm_pos] 2[Host_chr] 3[Host_pos]
 	foreach my $n(keys %{$sp1}) {
        	my @tmp=split(/--/, $n) ;
        	my @tmp1=split(/:/, $tmp[0]) ;
       		my @tmp2=split(/:/, $tmp[1]) ;
        	$isc{$n}[0]=$tmp1[0] ;
        	$isc{$n}[1]=$tmp1[1] ;
        	$isc{$n}[2]=$tmp2[0] ;
        	$isc{$n}[3]=$tmp2[1] ;
        	foreach my $j( keys(%{$sp1->{$n}}) ) {
                	$isr{$n}{$j}++ ;
        	}
	}
	foreach my $n(keys %{$sp2}) {
        	my @tmp=split(/--/, $n) ;
        	my @tmp1=split(/:/, $tmp[0]) ;
        	my @tmp2=split(/:/, $tmp[1]) ;
        	$isc{$n}[0]=$tmp1[0] ;
        	$isc{$n}[1]=$tmp1[1] ;
        	$isc{$n}[2]=$tmp2[0] ;
        	$isc{$n}[3]=$tmp2[1] ;
        	foreach my $j( keys(%{$sp2->{$n}}) ) {
                	$isr{$n}{$j}++;
        	}
	}
	foreach my $n(keys %isr) {
        	 $ish1{$n}=scalar (keys %{$isr{$n}}) ;
	}
	foreach my $n(keys %ish1) {  #Keeping IS with at least 2 supporting splt reads
        	if ($ish1{$n} <2) {
                	delete $ish1{$n} ;
                	delete $isr{$n} ;
                	delete $isc{$n} ;
        	}
	}
        return (\%ish1, \%isr, \%isc); #1K: Read name; El: @   
}

1;
