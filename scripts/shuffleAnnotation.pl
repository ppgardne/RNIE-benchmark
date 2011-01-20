#!/software/bin/perl 

#Take a gff file and either a genome length or genome sequence file and randomise the coordinates, keeping the same lengths and +/- strands. 

use warnings;
use strict;

use Getopt::Long;
use IO::File;

my ($gff, $genomeFile, $length, $help);
my $overlapThreshold=1.0;
&GetOptions( 
    "h|help"              => \$help,
    "l|length"            => \$length,
    "s|genome=s"          => \$genomeFile, 
    "g|gff=s"             => \$gff
    );

if( $help ) {
    &help();
    exit(1);
}
elsif ( (!-s $gff)){
    print "FATAL: [$gff] failed to stat!\n\n";
    &help();
    exit(1);
}
elsif ( defined $genomeFile && (!-s $genomeFile)){
    print "FATAL: [$genomeFile] failed to stat!\n\n";
    &help();
    exit(1);
}
elsif ( not defined $genomeFile && not defined $length){
    print "FATAL: you need to provide either a genome sequence file (-s <file>) of a genome length (-l <num>)!\n\n";
    &help();
    exit(1);
}

my $fh = IO::File->new();
if (not defined $length){
    
    my $seqstat = 'seqstat';
    $fh->open( "$seqstat $genomeFile |" );
    while(my $ss = <$fh>) {
	if($ss=~/Smallest:\s+(\d+)/){
	    $length=$1;
	}
	elsif($ss=~/Number of sequences:\s+(\d+)/){
	    die "FATAL: too many sequences in [$genomeFile]! Should be just 1." if $1>1;
	}
    }
    $fh->close();
}

$fh->open( $gff );
my $outFile = "$$.gff"; 
open(UT, "> $outFile");
while(my $g = <$fh>) {
    
    my @g = split(/\t/, $g);
    if ($g=~/^#/ || not defined $g[3] || not defined $g[4] || not defined $g[6] || @g != 9){
	print UT $g;
	next;
    }
    my $new3 = int(rand($length))+1;
    my $new4 = $new3 + ($g[4] - $g[3]);
    ($g[3], $g[4]) = ($new3, $new4);
    
    my $newG = join("\t", @g);
    print UT $newG;
}
close(UT);
$fh->close();

system("cat $$.gff | sort -k4n");

unlink("$$.gff");
######################################################################
sub help {
    print STDERR <<EOF;

$0: Given a gff file and either a genome length or
                      genome sequence file and permute the
                      coordinates, keeping the same length and strand
                      distributions.

Usage:   $0 -g <gff1> -s <genome seq. file>
         $0 -g <gff1> -l <genome length>

Options:       -h|--help                     Show this help.
	       -g|--gff <file>               A GFF file
	       -s|--genome <file>            A genome sequence file
               -l|--length <num>             The genome length

Dependencies:
seqstat from the SQUID components of older HMMER/INFERNAL package.

EOF
}
