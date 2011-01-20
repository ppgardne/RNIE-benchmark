#!/software/bin/perl 

#Given a terminator annotation and a genome annotation compute the
#distance of each terminator from the nearest 3' feature on the same
#strand.

use strict;
use warnings;

use Getopt::Long;
use IO::File;

my ($terminatorFeatures, $genomeFeatures, $help);
my $overlapThreshold=1.0;
&GetOptions( 
    "h|help"              => \$help,
    "t|terminators=s"     => \$terminatorFeatures, 
    "g|genome=s"          => \$genomeFeatures
    );

if( $help ) {
    &help();
    exit(1);
}
elsif ( (!-e $terminatorFeatures) || (!-e $genomeFeatures)){

    print "FATAL: either [$terminatorFeatures] or [$genomeFeatures] do not exist!\n\n";

    &help();
    exit(1);
}

my $fh = IO::File->new();
# slurp terminators ordered by coordinates
$fh->open( "cat $terminatorFeatures | sort -k4n |" );
my (%terminators, %dists);
while(my $ter = <$fh>) {
    next if ($ter =~ /^#/);
    my @ter = split(/\t/, $ter);
    next if (not defined $ter[3] || not defined $ter[4] || not defined $ter[6]);
    my $key = $ter[3] . '-' . $ter[4] . ':' . $ter[6];
    
    push(@{$terminators{$key}}, @ter);
    $dists{$key}       = 50000;
}
$fh->close;

$fh->open( "cat $genomeFeatures | sort -k4n |" );
my @features;
my ($minCoord, $maxCoord) = (0,0);
while(my $ft = <$fh>) {
    next if ($ft =~ /^#/);
    my @ft = split(/\t/, $ft);
    next if (not defined $ft[3] || not defined $ft[4] || not defined $ft[6]);
#    last if $ft[3] > 100000;
    chomp($ft);
    #print "ft:[$ft]\n";
    push(@features, \@ft);
    
    $maxCoord = $ft[4];
    #print "[$features[0][3]]\n";
    $minCoord = $features[0][3];
    
#    my ($i, $minTerm)=(0,0);
#    while($minTerm<$maxCoord){
    foreach my $key (sort {$terminators{$a}[3] <=> $terminators{$b}[3]} (keys %terminators)){

	#print "last if $minTerm<$maxCoord\n";
	
	###NEED SOMETHING LIKE THIS!!!!:::
	###last if $minTerm<$maxCoord;

	#print "$i:[$key:::$terminators{$key}]\n";
	#$minTerm=$terminators{$key}[3];
	
	next if not defined $terminators{$key}[3];
	next if $terminators{$key}[3] > $maxCoord + 1000;
	
	foreach my $ft2 (@features){
	    #print "ft2:[" . $ft2 . "] key:[$key] terminator{key}[6]:[$terminators{$key}[6]]\n";
	    #print "if::  [" . $terminators{$key}[6] . "] eq [" . $ft2->[6] . "]\n";
	    if ($terminators{$key}[6] eq $ft2->[6]){#Same strand
		my $dist = $terminators{$key}[3] - $ft2->[4] - 1;
		$dist    = $ft2->[3] - $terminators{$key}[4] - 1 if ($ft2->[6] eq '-');
		#print "key:[$key] ::: distHash:$dists{ $key }, distScalar:$dist\n" if(abs($dist)<1000);
		$dists{ $key } = minAbs($dists{ $key }, $dist);
		#print "key:[$key] ::: newDistHash:$dists{ $key }, distScalar:$dist\n" if(abs($dist)<1000);
	    }
	}
	
	#print "terminators{$key}[4]::[" . $terminators{$key}[4] . "]<minCoord::[$minCoord]\n";
	if ($terminators{$key}[4]<$minCoord){
	    #print "undef $key\n";
	    delete $terminators{$key};
	}
	
#	$i++;
    }
    
    my $shifted = shift(@features) if (scalar(@features)>50);
    #print "shifted::[$shifted->[3]-$shifted->[4]:$shifted->[6]]\n" if defined $shifted;
}
$fh->close;

#foreach my $key (sort {$a <=> $b} (keys %dists)){
foreach my $key ((keys %dists)){
    print "$key\t$dists{$key}\n";
}

exit(0);


######################################################################
#Max and Min
#return max in absolute value:
sub maxAbs {
  return $_[0] if @_ == 1;
  abs($_[0]) > abs($_[1]) ? $_[0] : $_[1]
}

#return min in absolute value:
sub minAbs {
  return $_[0] if @_ == 1;
  abs($_[0]) < abs($_[1]) ? $_[0] : $_[1]
}

######################################################################
sub help {
    print STDERR <<EOF;

terminatorDistances.pl: Given a terminator annotation and a genome
                        annotation compute the distance of each
                        terminator from the nearest 3\' feature on the
                        same strand.

Usage:   terminatorDistances.pl <gff1> <gff2> <gff3> ... 
Options:       -h|--help                     Show this help.

               -t|--terminators <file>       Annotation of predicted terminators 
	       -g|--genome      <file>       Annotation of all other features
	       

TODO:

EOF
}








