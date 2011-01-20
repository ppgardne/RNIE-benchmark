#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;

my ($help, $transtermfile, $fastafile, @warnings);
my $cmsearchThresh = 0;
&GetOptions( 
    "t|transterm=s"       => \$transtermfile,
    "f|fastafile=s"       => \$fastafile,
    "h|help"              => \$help 
    );

if( $help ) {
    &help();
    exit(1);
}

my %shortId2longId;
if(defined $fastafile && -s $fastafile){
    open(FA, "< $fastafile") or die  "FATAL: could not open for reading [$fastafile]\n[$!]";
    while(my $l = <FA>){
	next if $l !~ /^>/;
#Format is:
#>AAAB   501-542 D90227.1/10944-10985
	if($l=~/>(\S+)\s+\S+\s+(\S+)/){
	    $shortId2longId{$1}=$2;
	}
    }
    close(FA);
}
else {
    print "WARNING: fastafile not given or does not exist!\n";
}

# SEQUENCE D12501.1/1761-1799  (length 1038)

# g1                 10 - 100      + | 
# g2                200 - 250      + | 
# g3                300 - 350      + | 
# g4                400 - 490      + | 
#   TERM 1          468 - 482      + F    57  -4.2 -3.86346 | opp_overlap 475
#   TGGTCGCGCGTGCGC               GCTCAA AGC TCGAGC               TTTCTAGCACGATCC
#   TERM 2          505 - 530      + F   100 -16.1 -5.28176 | opp_overlap 505 503
#   CACGATCCCAAAAGA         GCCCGCACGGC GCAA GCCCTGCGGGC          TTTTTTGTGCCTCGA
# g5                600 - 690      + | 
#   TERM 3          693 - 705      + F    50  -4.7 -2.91264 | 
#   TAGGCCTCCGCGATC                GCACG CGA CGAGC                TGTGTCGATGGGAGG
# g6                750 - 800      + | 

if (defined $transtermfile && -s $transtermfile){
    my $name; 
    
    my %store;
    my $counter=10000;

    open(GFF, "> $transtermfile\.gff") or die  "FATAL: could not open for writing $transtermfile\.gff\n[$!]";
    open(TT, "< $transtermfile") or die  "FATAL: could not open for reading $transtermfile\n[$!]";
  TRANSTERM: while(<TT>){
	if(/^SEQUENCE\s+(\S+)/){
	    $name=$1;
	    $name=$shortId2longId{$name} if defined $shortId2longId{$name};
	}
	elsif(/TERM\s+\d+\s+(\d+)\s+\-\s+(\d+)\s+(\+|\-)\s+\w+\s+(\d+)/){
	    my ($start,$end,$strandString,$bits)=($1,$2,$3,$4);
	    #print "(name,start,end,strandString,bits)=($name,$start,$end,$strandString,$bits)\n";
	    #exit(0);
	    my ($maxOverlapPercent, $relativeStrand) = overlapCheck($name,$start,$end,\%store);
	    
	    if ($maxOverlapPercent>0.05){#more than 5% overlap, allows for some fuzzy boundaries...
		next TRANSTERM;
	    }
	    my $tag = 'terminator' . $counter;
	    push( @{ $store{$name} }, { 'start'  => $start,
					'end'    => $end,
					'bits'   => $bits,
					'tag'    => $tag
		  } );
	    
	    print GFF "$name\ttranstermHP\tterminator\t$start\t$end\t$bits\t$strandString\t\.\tID=$tag;Note=Predicted Rho independent terminator using transtermHP\n";
	    $counter++;
	}
    }
    close(TT);
    close(GFF);
}
else {
    die "FATAL: transterm file [$transtermfile] is either missing or empty!";
}

exit(0);

######################################################################
#fileRoot: strip the suffix off file names 
sub fileRoot {
    my $file = shift;
    $file =~ s/\.\S+$//;
    return $file;
}

######################################################################
#checkStrand: takes a star tand stop coord. Returns 1 if forward strand, -1 if reverse.
sub checkStrand{
    my ($a,$b)=@_;
    my $strand = 1;
    $strand = -1 if $b<$a;
    return $strand;
}

######################################################################
#overlapCheck: return true if the N/S-E overlaps with something in the %store hash:
sub overlapCheck{
    my ($name,$start,$end,$store)=@_;
    return (0,'') if not defined $store->{$name};
    my $strand = checkStrand($start,$end);
    ($start,$end) = reorder($start,$end);
    my ($maxOverlap,$relativeStrand)=(0,1);
    for (my $i=0; $i<scalar(@{$store->{$name}}); $i++){
	my $a = $store->{$name}[$i]{'start'};
	my $b = $store->{$name}[$i]{'end'};
	my $s = checkStrand($a,$b);
	($a,$b) = reorder($a,$b);
	if (overlap($a,$b,$start,$end)){
	    my $overlapExtent = overlapExtent($a,$b,$start,$end);
	    if ($maxOverlap < $overlapExtent){
		$maxOverlap = $overlapExtent;
		$relativeStrand = $strand * $s;
	    }
	}
    }
    return ($maxOverlap, $relativeStrand);
}

######################################################################
#Returns true if the coordinates for two regions ($x1, $y1) and ($x2, $y2) overlap:
# - assumes that $x1 < $y1 and $x2 < $y2.
sub overlap {
    my($x1, $y1, $x2, $y2) = @_;
    
    if ( ($x1<=$x2 && $x2<=$y1) || ($x1<=$y2 && $y2<=$y1) || ($x2<=$x1 && $x1<=$y2) || ($x2<=$y1 && $y1<=$y2)  ){
        return 1;
    }
    else {
        return 0;
    }
}

######################################################################
#Returns the extent of overlap between two regions A=($x1, $y1) and B=($x2, $y2):
# - assumes that $x1 < $y1 and $x2 < $y2.
#
# D = 2*|A n B|/(|A|+|B|)
#
sub overlapExtent {
    my($x1, $y1, $x2, $y2) = @_;
    
    ($x1, $y1) = reorder($x1, $y1);
    ($x2, $y2) = reorder($x2, $y2);
    # 1.
    # x1                   y1
    # |<---------A--------->|
    #    |<------B------>|
    #    x2             y2
    #    XXXXXXXXXXXXXXXXX
    # 2.  x1                     y1
    #     |<---------A----------->|
    # |<-------------B------>|
    # x2                    y2
    #     XXXXXXXXXXXXXXXXXXXX
    # 3. x1             y1
    #    |<------A------>|
    # |<---------B--------->|
    # x2                   y2
    #    XXXXXXXXXXXXXXXXX
    # 4. x1                    y1
    #    |<-------------A------>|
    #        |<---------B----------->|
    #        x2                     y2
    #        XXXXXXXXXXXXXXXXXXXX
    my $D=0;
    my $int=0;
    my $L1=$y1-$x1+1;
    my $L2=$y2-$x2+1;
    my $minL = min($L1,$L2);
    if ( ($x1<=$x2 && $x2<=$y1) && ($x1<=$y2 && $y2<=$y1) ){    #1.
	$D = $L2;
    }
    elsif ( ($x2<=$x1) && ($x1<=$y2 && $y2<=$y1) ){              #2.
	$D = $y2-$x1+1;
    }
    elsif ( ($x2<=$x1 && $x1<=$y2) && ($x2<=$y1 && $y1<=$y2) ){ #3.
	$D = $L1;
    }
    elsif ( ($x1<=$x2 && $x2<=$y1) && ($y1<=$y2) ){              #4.
	$D = $y1-$x2+1;
    }
    return $D/$minL;
}

######################################################################
#reorder: given 2 integers, return the smallest first & the largest last:
sub reorder {
    my ($x,$y)=@_;
    
    if ($y<$x){
	my $tmp = $x;
	$x = $y;
	$y = $tmp;
    }
    return ($x,$y);
}

######################################################################
#Max and Min
#max
sub max {
  return $_[0] if @_ == 1;
  $_[0] > $_[1] ? $_[0] : $_[1]
}

#min
sub min {
  return $_[0] if @_ == 1;
  $_[0] < $_[1] ? $_[0] : $_[1]
}

sub help {
    print STDERR <<EOF;

transtermHP2gff.pl: 

Usage:   transtermHP2gff.pl <options> 
Options:       -h|--help                     Show this help.

               -t|--transterm <file>         TranstermHP output file
               -f|--fastafile <file>         Fasta file used as input for TranstermHP
	                                          - header format should be:
	                                          >AAAB   501-542 D90227.1/10944-10985
EOF
}




