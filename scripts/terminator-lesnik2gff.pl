#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;

my ($help, $lesnikfile, @warnings);
my $cmsearchThresh = 0;
&GetOptions( 
    "t=s"       => \$lesnikfile,
    "h|help"            => \$help 
    );

if( $help ) {
    &help();
    exit(1);
}

# 0                                                               1 2     3 4    5 6
# D12501.1/1761-1799                                              | -4.30 | 0.00 | -0.80 | -2.70 | -5.70 | -1.87 | 4.0 | 3.0 | 3.0 | 0.6 |      0  509 34 gcacggcgcaa gccc     tgc            gggc .  ttttt tgt gcct
# D12501.1/1761-1799                                              | -10.60 | 0.00 | -5.40 | -1.10 | -5.40 | -4.38 | 6.0 | 6.0 | 0.0 | 0.0 |     1  325 41 gcttctcggtg gcgccg   ctggga       cggcgc .  ttgct ctt gtgc
# D12501.1/1761-1799                                              | -8.90 | 0.00 | -5.40 | -1.10 | -5.40 | -2.68 | 5.0 | 8.0 | 0.0 | -2.0 |     1  325 41 gcttctcggtg gcgcc    gctgggac      ggcgc .  ttgct ctt gtgc
# D12501.1/1761-1799                                              | -7.50 | 0.00 | -2.20 | -4.20 | -4.30 | -2.98 | 7.0 | 7.0 | 5.0 | 0.0 |      1  534 44 aaaagcccgca gggcttg  cgccgtg     cgggctc .  ttttg gga tcgt
# D12501.1/1761-1799                                              | -9.90 | 0.00 | -2.20 | -4.20 | -4.30 | -5.38 | 8.0 | 5.0 | 5.0 | 1.0 |      1  534 44 aaaagcccgca gggcttgc gccgt      gcgggctc .  ttttg gga tcgt

# 1.  hairpin_score,    2
# 2.  spacer_score,     4
# 3.  proxy_T_score,    6
# 4.  dist_T_score,     8
# 5.  extra_T_score,    10
# 6.  dG_score, **      12
# 7.  stem_score,       14
# 8.  loop_score,       16
# 9.  A_score,          18
# 10. struct_score **   20

if (-s $lesnikfile){
    
    my %store;
    my $counter=10000;

    open(GFFss, "> $lesnikfile\.struct_score.gff") or die  "FATAL: could not open for writing $lesnikfile\.struct_score.gff\n[$!]";
    open(GFFdg, "> $lesnikfile\.dG_score.gff")     or die  "FATAL: could not open for writing $lesnikfile\.dG_score.gff\n[$!]";
    open(TT, "< $lesnikfile")                    or die  "FATAL: could not open for reading $lesnikfile\n[$!]";
  LESNIK: while(<TT>){
      next LESNIK if /^\#/;
      my @les = split(/\s+/, $_);
      #next LESNK if $les[22]; #Skip reverse strand hits
      my ($name, $start,$end,$dG_score,$struct_score)=($les[0],$les[23],$les[23]+$les[24],$les[12],$les[20]);
      my ($maxOverlapPercent, $relativeStrand) = overlapCheck($name,$start,$end,\%store);
      
      if ($maxOverlapPercent>0.05){#more than 5% overlap, allows for some fuzzy boundaries...
	  next LESNIK;
      }
      my $tag = 'terminator' . $counter;
      push( @{ $store{$name} }, { 'start'  => $start,
				  'end'    => $end,
				  'bits'   => $struct_score,
				  'tag'    => $tag
	    } );
      
      my $strandString = '+';
      $strandString = '-' if $les[22];
      my $transDGScore = $dG_score * (-1);
      print GFFss "$name\tRNAMotif\tterminator\t$start\t$end\t$struct_score\t$strandString\t\.\tID=$tag;Note=Predicted Rho independent terminator using Lesnik-RNAMotif;dG_score=$dG_score;struct_score=$struct_score\n";
      print GFFdg "$name\tRNAMotif\tterminator\t$start\t$end\t$transDGScore\t$strandString\t\.\tID=$tag;Note=Predicted Rho independent terminator using Lesnik-RNAMotif;dG_score=$dG_score;struct_score=$struct_score\n";
      $counter++;
      
  }
    close(TT);
    close(GFFss);
    close(GFFdg);
}
else {
    die "FATAL: lesnik file [$lesnikfile] is either missing or empty!";
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

terlminator-lesnik2gff.pl: 

Usage:   terlminator-lesnik2gff.pl <options> 
Options:       -h|--help                     Show this help.
               -t <file>                     RNAmotif with Lesnik desciptor output.

EOF
}




