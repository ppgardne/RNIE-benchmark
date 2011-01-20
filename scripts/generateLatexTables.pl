#!/usr/bin/perl 

#Read stats on the control genome seqs. 

use warnings;
use strict;
use Env;

my $genomesDir = $ENV{RNIE} . "/benchmark/genomes";
my @genomeFiles = glob("$genomesDir/*embl");

my %table;
foreach my $f (@genomeFiles){
#genomes/AE000511.1.embl
    if($f=~/$genomesDir\/(\S+)\.(\d+)\.embl/){    
	my $acc = $1; 
	my $version = $2;
	my ($species, $phyla, $length, $cds, $gc)=('','', 0, 0, 0);
	open(GE, "< $f");
	while(<GE>){
	    if(/^OS\s+(.*)\s*$/){
		$species = $1;
		$species = species2shortspecies($species);
		print "$acc.$version\t$species";
	    }
	    elsif(/^OC\s+Bacteria;\s+(\S+);/){
		$phyla = $1;
	    }
	    elsif(/FT\s+CDS\s+/){
		$cds++;
	    }
	}
	close(GE);
	
	open(GE, "esl-seqstat -c $f |");
	while(<GE>){
	    if(/Smallest\:\s+(\d+)/){
		$length=$1;
	    }
	    elsif(/^residue:\s+[C|G]\s+\d+\s+(\d\.\d+)/){
		$gc+=$1;
	    }
	}
	close(GE);

	print "\tlength=$length\n";
	
	my $numNative   = `wc -l $genomesDir/$acc\.$version-genomeMode-rnie.gff`;
	my $numShuffled = `wc -l $genomesDir/$acc\.$version-shuffled-genomeMode-rnie.gff`;
	$numNative    =~ s/\s+\S+$//;
	$numShuffled =~ s/\s+\S+$//;

	my $numNativeGeneMode   = `wc -l $genomesDir/$acc\.$version-geneMode-rnie.gff`;
	#                                                      -shuffled-geneMode-rnie.gff
	my $numShuffledGeneMode = `wc -l $genomesDir/$acc\.$version-shuffled-geneMode-rnie.gff`;
	$numNativeGeneMode    =~ s/\s+\S+$//;
	$numShuffledGeneMode =~ s/\s+\S+$//;

	$table{$acc} =           { 'species'   => $species,
				   'phyla'     => $phyla,
				   'length'    => $length,
				   'gc'        => $gc,
				   'cds'       => $cds,
				   'native'    => $numNative,
				   'shuffled'  => $numShuffled,
				   'nativeS'   => $numNativeGeneMode,
				   'shuffledS' => $numShuffledGeneMode
	};
    }
    
}

#my @train = qw(AL009126 U00096);
my @test  = qw(
AE000516
AP009493
AE015928
AE001363
AE017126
AE000513
AL009126
AM180355
AE009951
CP001147
U00096
AE000511
AE014613
AE016823
AF222894
CP000771
CP000975
);

open(TB, "> table.tex");

print TB "

\\hspace{-10mm}
\\begin{table}
\\caption[]{Control genomes. }\\label{table:1}
\\begin{tabular}{lllcrcrrrr}
\\hline
\\hline
Species & EMBL      & Phylum & Genome & Number & G\+C    & \\multicolumn{4}{c}{Number predictions}               \\\\
        & accession &        & size   & CDSs   & content & \\multicolumn{2}{c}{Genome} & \\multicolumn{2}{c}{Gene} \\\\
        &           &        & (MB)   &        &         & Nat. & Shuff.             & Nat. & Shuff.             \\\\
\\hline
";

#foreach my $embl (@train, @test){
foreach my $embl (@test){
    
    printf TB     "\\emph\{%s\} & $embl & %s & %0.2f & %d & %0.2f & %d & %d & %d & %d \\\\ \n", $table{$embl}{'species'}, $table{$embl}{'phyla'}, ($table{$embl}{'length'})/1000000, $table{$embl}{'cds'}, $table{$embl}{'gc'}, $table{$embl}{'native'}, $table{$embl}{'shuffled'}, $table{$embl}{'nativeS'}, $table{$embl}{'shuffledS'};
    printf STDOUT "\\emph\{%s\} & $embl & %s & %0.2f & %d & %0.2f & %d & %d & %d & %d \\\\ \n", $table{$embl}{'species'}, $table{$embl}{'phyla'}, ($table{$embl}{'length'})/1000000, $table{$embl}{'cds'}, $table{$embl}{'gc'}, $table{$embl}{'native'}, $table{$embl}{'shuffled'}, $table{$embl}{'nativeS'}, $table{$embl}{'shuffledS'};
    
}

print TB "
\\hline
\\hline
  \\end{tabular}
\\end{table}
";

close(TB);

exit(0);


######################################################################
#species2shortspecies: Given a species string eg. "Homo sapiens
#                      (human)" generate a nicely formated short name
#                      with no whitespace eg. "H.sapiens".
sub species2shortspecies {
    my $species = shift;
    my $shortSpecies;
    
    if ($species=~/(.*)\s+sp\./){
	$shortSpecies = $1;
    }
    elsif ($species=~/metagenome/i or $species=~/uncultured/i){
	$species=~s/metagenome/metag\./gi;
	$species=~s/uncultured/uncult\./gi;
	my @w = split(/\s+/,$species);
	if(scalar(@w)>2){
	    foreach my $w (@w){
		$shortSpecies .= substr($w, 0, 5) . '.';
	    }
	}
	else {
	    $shortSpecies = $species;
	    $shortSpecies =~ s/\s+/_/g;
	}
    }#lots of conditions here. Need else you get some ridiculous species names.
    elsif($species=~/^(\S+)\s+(\S{4,})/ ){
	$shortSpecies = substr($1,0,1) . "." . $2; 
    }
    else {
	$shortSpecies = $species;
    }
    
    $shortSpecies =~ s/\s+/_/g;
    $shortSpecies =~ s/[\'\(\)\:\/]//g;
#    $shortSpecies = substr($shortSpecies,0,20) if (length($shortSpecies) > 20);
    
#   H.P 
    return $shortSpecies;
}

