#!/software/bin/perl 

#Given multiple gff files, merge them together and combine entries with exactly the same coordinates. 

use warnings;
use strict;
use Getopt::Long;
use IO::File;

my (@gffs, $help);
my $overlapThreshold=1.0;
&GetOptions( 
    "h|help"              => \$help,
    "ot|overlapthresh=s"  => \$overlapThreshold,
    "g|gff=s@"            => \@gffs
    );

if( $help ) {
    &help();
    exit(1);
}
elsif (@gffs == 0){
    print "FATAL: no gff files given\n";
    &help();
    exit(1);
}

my $gffs = join(' ', @gffs);

my $fh = IO::File->new();
# sort features so they can be merged in coordinate
# order
$fh->open( "cat $gffs | sort -k4n |" );
my ($cur, $prev, $ext);
while($cur = <$fh>) {
    
    next if ($cur =~ /^#/);
    if(defined $prev){
	my @cur  = split(/\t/,  $cur);
	my @prev = split(/\t/, $prev);
	next if (not defined $cur[3] || not defined $cur[4] || not defined $cur[6]);
	$ext  = 0;
	my $new;
	$ext     = overlapExtent($cur[3], $cur[4], $prev[3], $prev[4]  ) if ($cur[6] eq $prev[6]);
	if ($ext>=$overlapThreshold){
	    $new  = merge(\@cur, \@prev) ;
	}
	else {
	    $new = $cur;
	}
	
	$cur  = $new;
	print $prev if $ext < $overlapThreshold;
    }
    $prev = $cur; 
}
$fh->close;

print $prev if $ext < $overlapThreshold;

exit(0);
######################################################################
#
sub merge {
    my ($gff1, $gff2)=@_;
    
    $gff1->[3]  = $gff2->[3] if ($gff1->[3] > $gff2->[3]);
    $gff1->[4]  = $gff2->[4] if ($gff1->[4] < $gff2->[4]);
    #Escape special chars
    if ((length($gff1->[8])>0) && (length($gff2->[8])>0) && ($gff1->[8] !~ /\Q$gff2->[8]\E/) && ($gff2->[8] !~ /\Q$gff1->[8]\E/)){
	chomp($gff1->[8]);
	$gff1->[8] .= ';' . $gff2->[8] ;
    }
    elsif (length($gff2->[8])>length($gff1->[8])) {
	$gff1->[8] = $gff2->[8];
    }
    
    my $new = join("\t", @{$gff1});
    return $new;
}

######################################################################
#Returns the extent of overlap between two regions A=($x1, $y1) and B=($x2, $y2)
#using the following metric:
#
# D1 = 2*|A n B|/(|A|+|B|)
#
#An alternative, more permissive, metric is:
#
# D2 = |A n B|/min(|A|,|B|)
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
    elsif ( ($x2<=$x1) && ($x1<=$y2 && $y2<=$y1) ){             #2.
	$D = $y2-$x1+1;
    }
    elsif ( ($x2<=$x1 && $x1<=$y2) && ($x2<=$y1 && $y1<=$y2) ){ #3.
	$D = $L1;
    }
    elsif ( ($x1<=$x2 && $x2<=$y1) && ($y1<=$y2) ){             #4.
	$D = $y1-$x2+1;
    }
#D1:
    return (2*$D)/($L1+$L2);
#D2:
#    return $D/$minL;
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

######################################################################
sub help {
    print STDERR <<EOF;

collapseGffs.pl: Given multiple gff files, merge them together and combine entries with exactly the same coordinates. 

Usage:   collapseGffs.pl -g <gff1> -g <gff2> -g <gff3> ... 
Options:       -h|--help                     Show this help.

               -ot|--overlapthresh <val>     Merge features overlapping by <val>\% or more. default=[$overlapThreshold]
	       -g|--gff <file>               Give GFF file names as input. For multiple gff files use additional -g\'s.

EOF
}

