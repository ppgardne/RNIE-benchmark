#!/usr/bin/perl 

use warnings;
use strict;
use Getopt::Long;

my ($help, $rniegfffile);
my $cmsearchThresh = 0;
&GetOptions( 
    "r|rnie=s"       => \$rniegfffile,
    "h|help"         => \$help 
    );

my @models = qw(
seed-a-14
seed-b-14
);
# super43-45-seed-1
# terminator32-40-1
# terminator40-43-1
# terminator43-45-1
# );


if (defined $rniegfffile && -s $rniegfffile){
    
    my $outfilePrefix = $rniegfffile;
    $outfilePrefix =~ s/\.gff//i;
    
    open(IN,  "< $rniegfffile") or die  "FATAL: could not open for reading [$rniegfffile]\n[$!]";
    open(UTb, "> $outfilePrefix\.bits.gff") or die  "FATAL: could not open for writing [$outfilePrefix\.bits.gff]\n[$!]";
    open(UTe, "> $outfilePrefix\.E_value.gff") or die  "FATAL: could not open for writing [$outfilePrefix\.E_value.gff]\n[$!]";
    open(UTs, "> $outfilePrefix\.sumBits.gff") or die  "FATAL: could not open for writing [$outfilePrefix\.sumBits.gff]\n[$!]";
    
    my %modFH;
    foreach my $mod (@models){
	open($modFH{$mod}, "> $outfilePrefix\.$mod\.bits.gff") or die  "FATAL: could not open for writing [$outfilePrefix\.$mod\.bits.gff]\n[$!]";
    }
    
    while(my $l = <IN>){
	next if $l =~ /^#/;
	my @l = split(/\t/, $l);
	next if $l[6] =~ /\-/; #Only counting the positive strand
	my $val;
	$val = fetchTag($l[8], 'Bits');
	printFeature(\*UTb, \@l, $val*1) if $val ne 'NA';
	
	$val = fetchTag($l[8], 'E_value');
	printFeature(\*UTe, \@l, 1000/$val) if $val ne 'NA';

	$val = fetchTag($l[8], 'sumBits');
	printFeature(\*UTs, \@l, $val*1) if $val ne 'NA';
	
	foreach my $mod (@models){
	    $val = fetchTag($l[8], $mod);
	    if ($val ne 'NA'){
		$l[5] = $val*1;
		my $printMe = join("\t", @l);
		print {$modFH{$mod}} "$printMe";
	    }
	    #printFeature({$modFH{$mod}}, \@l, $val*1) if $val ne 'NA';
	}
    }
    
    close(IN);
    close(UTb);
    close(UTe);
    close(UTs);
    foreach my $mod (@models){
	close($modFH{$mod});
    }
    
}

exit(0);

sub fetchTag {
    my ($comment, $tag) = @_;
    my $sep = "\,\:\;\\s";
    if($comment=~/[$sep]$tag\s*=\s*(\S+?)[$sep\n]/){
	return $1;
    }
    else {
	return 'NA';
    }
}

sub printFeature {
    my ($fh, $gff, $val)=@_;
    
    $gff->[5] = $val;
    my $printMe = join("\t", @{$gff});
    print $fh "$printMe";
}
