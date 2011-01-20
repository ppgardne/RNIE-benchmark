#!/software/bin/perl 

use warnings;
use strict;

my $emblFile=shift;

die "FATAL: [$emblFile] does not exist!" if (not -s $emblFile);

open(EM, "< $emblFile") or die "FATAL: failed to open [$emblFile]!\n[$!]" ;
my(%store);#$type,$s,$e,$strand,$gene,$product,$notes);
my $seqId = 'undefined';
while(my $l=<EM>){
    
    if($l=~/^ID\s+(\S+);/){
	$seqId=$1;
    }
    elsif($l=~/^FT\s{3}(\S+)\s+(\d+)\.\.(\d+)/){
	next if $1 eq 'source';
	printGff($seqId,\%store) if (defined $seqId && defined $store{'s'} && defined $store{'e'} && defined $store{'strand'}); 
	undef %store;
	($store{'type'},$store{'s'},$store{'e'},$store{'strand'})=($1,$2,$3,'+');
    }
    elsif($l=~/^FT\s{3}(\S+)\s+complement\((\d+)\.\.(\d+)\)/){
	printGff($seqId,\%store) if (defined $seqId && defined $store{'s'} && defined $store{'e'} && defined $store{'strand'}); 
	undef %store;
	($store{'type'},$store{'s'},$store{'e'},$store{'strand'})=($1,$2,$3,'-');
    }
    elsif($l=~/^FT\s+\/gene="(.*)/){
	$store{'gene'}=$1;
	$store{'gene'}=~s/\"$//;
    }
    elsif($l=~/^FT\s+\/product="(.*)/){
	$store{'product'}=$1;
	$store{'product'}=~s/\"$//;
    }
    elsif($l=~/^FT\s+\/note="(.*)/){
	$store{'note'}=$1;
	$store{'note'}=~s/\"$//;
    }
    
    
}
close(EM);

exit(0);

sub printGff {
    my ($seqId,$store)=@_;
    $store->{'type'}='gene' if (not defined $store->{'type'});
    my $phase='.';
    $phase=0 if $store->{'type'}=~/CDS/;
    my ($attributes,$spacer)=('','');
    $attributes .= 'Name='    . $store->{'gene'}                     if (defined $store->{'gene'});
    $spacer = ';' if length($attributes)>0;
    $attributes .= $spacer . 'Product="' . $store->{'product'} . '"' if (defined $store->{'product'});
    $spacer = ';' if length($attributes)>0;
    $attributes .= $spacer . 'Note="'    . $store->{'notes'}   . '"' if (defined $store->{'notes'});
    $attributes .= ';';

    print "$seqId\tEMBL\t".$store->{'type'}."\t".$store->{'s'}."\t".$store->{'e'}."\t.\t".$store->{'strand'}."\t.\t$attributes\n";
    
}
