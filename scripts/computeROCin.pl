#!/usr/bin/perl 

#read in gff files, return dat files, sorted on score and annotated as
#either FP or TP, tag FN's to the end of the file.
#

use warnings;
use strict;



my @files = qw(
Embedded-rnie.bits.gff
Embedded-rnie.E_value.gff
Embedded-rnie.sumBits.gff
Embedded-2features.transterm.gff
Embedded-4features.transterm.gff
Embedded-9features.transterm.gff
Embedded-10features.transterm.gff
Embedded.rnall-barrick.dG.gff
Embedded.rnall-barrick.hbG.gff
Embedded.rnall-wan.dG.gff
Embedded.rnall-wan.hbG.gff
Embedded.terminator-lesnik.out.dG_score.gff
Embedded.terminator-lesnik.out.struct_score.gff
);

my @rnieModels = qw(
terminator32-40
terminator40-43
terminator43-45
);

#Coords of terminators within trueEmbedded.fa 
my $trueEmbedded = 'trueEmbedded.fa';
open(TE, "< $trueEmbedded");
my %trueEmbeddedCoords;
while(<TE>){
    if(/^>(\S+)\s+(\d+)\-(\d+)$/){
	#print "trueEmbeddedCoords{$1}=($2, $3)\n";
	push(@{ $trueEmbeddedCoords{$1} }, ($2, $3));

    }
}
close(TE);

foreach my $file (@files){
    my $tfile = 'true' . $file;
    my %fileHandles;
    my %seenTp;
    my %seen = %trueEmbeddedCoords;
    open(TF, "< $tfile");
    open(UT, "> $tfile\.dat1");
    open(UT2, "> $tfile\.dat2") if $tfile !~ /trueEmbedded.transterm/;
    open(UT3, "> $tfile\.dat3") if $tfile =~ /trueEmbedded.rnall/;
    if( $tfile =~ /all-trueEmbedded.gff/){
	open(UT3, "> $tfile\.dat3");
	###
	foreach my $m (@rnieModels){
	    open(my $fh, "> $tfile\.$m\.dat1");
	    $fileHandles{$m}=$fh;
	}
    }
    
    my $count=0;
    #print "tfile:[$tfile]\n";
    while(<TF>){
	my @gff = split(/\t/, $_);
	next if $gff[6] !~ /\+/;
	next if $tfile =~ /all-trueEmbedded.gff/ && $gff[5]<5;
	my $ol = 0;
	#print "overlap($gff[3],$gff[4],trueEmbeddedCoords{$gff[0]}=|$trueEmbeddedCoords{$gff[0]}|)\n";
	$ol = overlap($gff[3],$gff[4],$trueEmbeddedCoords{$gff[0]}[0], $trueEmbeddedCoords{$gff[0]}[1]) if $trueEmbeddedCoords{$gff[0]}[0] and $trueEmbeddedCoords{$gff[0]}[1];
      	my $score2='NA';
	my $score3='NA';                                    #E_value=4.34e-03;sumBits=73.92;Note
	if( $tfile !~ /trueEmbedded.transterm/ && ($gff[8]=~/E_value\=(\S+)\;sumBits\=(\S+)\;Note/ or $gff[8]=~/dG_score\=(\S+)/ or $gff[8]=~/Twt=(\S+);hbG=(\s*\S+\s*);Thresh_dG/)){
	    $score2=$1;
	    $score3=$2           if $2;
	    $score2=$score2*(-1) if $tfile =~ /trueEmbedded.terminator-lesnik.out.gff/;
	    $score2=$score2*(1)  if $tfile =~ /all-trueEmbedded.gff/;
	    $score3 = 'NA' if $score3 !~ /NA/ && $score3 == 9999;
	}
	delete($seen{$gff[0]});
	
	my $score1 = $gff[5];
	$score1=$score1*(-1) if $tfile =~ /trueEmbedded.rnall/;
	$score3=$score3*(-1) if $score3 !~ /NA/ && $tfile =~ /trueEmbedded.rnall/;
	#print $_ if $score3 !~ /NA/ && $score3<0 && $tfile =~ /trueEmbedded.rnall/;
	#print "score1:[$score1] score2:[$score2] score3:[$score3]\n" if $tfile =~ /trueEmbedded.rnaII.gff/;
	
	my $tag = 'FP';
	if($ol&&!$seenTp{$gff[0]}){
	    $tag = 'TP';
	    $seenTp{$gff[0]}=1;
	}
	#print "ol:[$ol] tag:[$tag] [$gff[3]],[$gff[4]],[$trueEmbeddedCoords{$gff[0]}[0]], [$trueEmbeddedCoords{$gff[0]}[1]]\n" if $count<5;

	printf UT  "%0.2f\t$tag\t%d\n",   $score1, $gff[4]-$gff[3];
	printf UT2 "%0.8f\t$tag\t%d\n",   $score2, $gff[4]-$gff[3] if defined($score2) && $score2 !~ /NA/;
	printf UT3 "%0.3f\t$tag\t%d\n",   $score3, $gff[4]-$gff[3] if defined($score3) && $score3 !~ /NA/;
	
	if( $tfile =~ /all-trueEmbedded.gff/ && $gff[8]=~/maxModel=(\S+),allModels/){
		my $model = $1;
		printf { $fileHandles{$model} } "%0.2f\t$tag\t%d\n", $gff[5], $gff[4]-$gff[3] if defined $fileHandles{$model};
	}
	$count++;
    }
    close(TF);
    
    foreach my $k (keys %seen){
	printf UT "NA\tFN\t%d\n",  $seen{$k}[1]-$seen{$k}[0];
	printf UT2 "NA\tFN\t%d\n", $seen{$k}[1]-$seen{$k}[0] if $tfile !~ /trueEmbedded.transterm/;
	printf UT3 "NA\tFN\t%d\n", $seen{$k}[1]-$seen{$k}[0] if $tfile =~ /all-trueEmbedded.gff/ or $tfile =~ /trueEmbedded.rnall/;
    }
    
    close(UT);
    close(UT2) if $tfile !~ /trueEmbedded.transterm/;
    close(UT3) if $tfile =~ /trueEmbedded.rnall/;
    if( $tfile =~ /all-trueEmbedded.gff/){
	close(UT3);
	foreach my $m (@rnieModels){
	    close($fileHandles{$m});
	}
    }
}

foreach my $file (@files){
    my $ffile = 'false' . $file;

    my %fileHandles;
    open(FF, "< $ffile");
    open(UT, "> $ffile\.dat1");
    open(UT2, "> $ffile\.dat2") if $ffile !~ /falseEmbedded.transterm/;
    open(UT3, "> $ffile\.dat3") if $ffile =~ /falseEmbedded.rnall/;
    
    if( $ffile =~ /all-falseEmbedded.gff/){
	open(UT3, "> $ffile\.dat3");
	foreach my $m (@rnieModels){
	    open(my $fh, "> $ffile\.$m\.dat1");
	    $fileHandles{$m}=$fh;
	}

    }

    #print "ffile:[$ffile]\n";
    
    while(<FF>){
	
	my @gff = split(/\t/, $_);
	next if $gff[6] !~ /\+/;
	next if $ffile =~ /all-falseEmbedded.gff/ && $gff[5]<5;
	
	my $score2='NA';
	my $score3='NA';
#	if( $ffile !~ /falseEmbedded.transterm/ && ($gff[8]=~/E_value\=(\S+)\;Note/                 or $gff[8]=~/dG_score\=(\S+)/)){
	if( $ffile !~ /falseEmbedded.transterm/ && ($gff[8]=~/E_value\=(\S+)\;sumBits\=(\S+)\;Note/ or $gff[8]=~/dG_score\=(\S+)/ or $gff[8]=~/Twt=(\S+);hbG=(\s*\S+\s*);Thresh_dG/)){
	    $score2=$1;
	    $score3=$2      if $2;
	    $score2=$score2*(-1) if $ffile =~ /falseEmbedded.terminator-lesnik.out.gff/;
	    $score2=$score2*(1)  if $ffile =~ /all-falseEmbedded.gff/;
	    $score3 = 'NA' if $score3 !~ /NA/ && $score3 == 9999;
	}
	my $score1 = $gff[5];
	$score1=$score1*(-1) if $ffile =~ /falseEmbedded.rnall/;
	$score3=$score3*(-1) if $score3 !~ /NA/ && $ffile =~ /falseEmbedded.rnall/;
	
	#print "score1:[$score1] score2:[$score2] score3:[$score3]\n" if $ffile =~ /falseEmbedded.rnaII.gff/;

	printf UT  "%0.2f\tFP\t%d\n", $score1, $gff[4]-$gff[3];
	printf UT2 "%0.8f\tFP\t%d\n", $score2, $gff[4]-$gff[3] if defined($score2) && $score2 !~ /NA/;
	printf UT3 "%0.3f\tFP\t%d\n", $score3, $gff[4]-$gff[3] if defined($score3) && $score3 !~ /NA/;
	
	if( $ffile =~ /all-falseEmbedded.gff/ && $gff[8]=~/maxModel=(\S+),allModels/){
	    printf { $fileHandles{$1} } "%0.2f\tFP\t%d\n", $gff[5], $gff[4]-$gff[3];
	}
	
    }
    close(TF);
    
    close(UT);
    close(UT2) if $ffile !~ /falseEmbedded.transterm/;
    close(UT3) if $ffile =~ /falseEmbedded.rnall/;
    
    if( $ffile =~ /all-falseEmbedded.gff/){
	close(UT3);

	foreach my $m (@rnieModels){
	    close($fileHandles{$m});
	}
    }
}

my @datFiles;
for (my $i=0; $i<scalar(@falseFiles); $i++){
    my $rev1 = '';
    $rev1 = 'r' if $trueFiles[$i] =~ /all-trueEmbedded.gff/ || $trueFiles[$i] =~ /trueEmbedded.transterm/;
    my $rev2 = '';
    $rev2 = 'r' if $trueFiles[$i] =~ /trueEmbedded.terminator-lesnik.out.gff/;
    my $rev3 = '';
    $rev3 = 'r' if $trueFiles[$i] =~ /all-trueEmbedded.gff/;
    
    system("cat $trueFiles[$i]\.dat1 $falseFiles[$i]\.dat1 | sort -n$rev1 > $trueFiles[$i]\.1.dat");
    system("cat $trueFiles[$i]\.dat2 $falseFiles[$i]\.dat2 | sort -n$rev2 > $trueFiles[$i]\.2.dat") if $trueFiles[$i] !~ /trueEmbedded.transterm/;
    system("cat $trueFiles[$i]\.dat3 $falseFiles[$i]\.dat3 | sort -n$rev3 > $trueFiles[$i]\.3.dat") if $trueFiles[$i] =~ /trueEmbedded.rnall/;
    push(@datFiles, "$trueFiles[$i]\.1.dat");
    push(@datFiles, "$trueFiles[$i]\.2.dat") if $trueFiles[$i] !~ /trueEmbedded.transterm/;    
    push(@datFiles, "$trueFiles[$i]\.3.dat") if $trueFiles[$i] =~ /trueEmbedded.rnall/;    

    if( $trueFiles[$i]  =~ /all-trueEmbedded.gff/){
	system("cat $trueFiles[$i]\.dat3 $falseFiles[$i]\.dat3 | sort -n$rev3 > $trueFiles[$i]\.3.dat");
	push(@datFiles, "$trueFiles[$i]\.3.dat") if $trueFiles[$i]  =~ /all-trueEmbedded.gff/;
	
	foreach my $m (@rnieModels){
	    system("cat $trueFiles[$i]\.$m\.dat1 $falseFiles[$i]\.$m\.dat1 | sort -n$rev1 > $trueFiles[$i]\.$m\.1\.dat");
	    push(@datFiles, "$trueFiles[$i]\.$m\.1.dat");
	}
    }
}

my $totalTrue     = 981;
my $totalTrueNucs = 42419;
my $totalNucs     = 1022438 + 102243350;
$totalNucs /= 1000000; #in MB...
#Compute the stats:
foreach my $dfile (@datFiles){
    my $maxMcc=0;
    open(DF, "< $dfile");
    open(DO,  "> $dfile\.rin");
    open(DO2, "> $dfile\.rin2") if$dfile=~/all-trueEmbedded.gff\S+1.dat/;
    my $cnt=0;
    open(DF, "< $dfile");
    my ($tp,$fp,$fn)            =(0,0,0);
    my ($tpNucs,$tnNucs,$fpNucs,$fnNucs)=(0,0,0,0);
    my ($lastSens,$lastPpv,$lastScore) = (0,0,0);
    my ($lastMccNucs,$lastScoreNucs) = (0,0);
    while(<DF>){
	my @df = split(/\t/, $_);
	next if $df[1] =~ /FN$/;
	######TP,TN,FP,FN
	$tp++           if $df[1] =~ /TP$/;
	$tpNucs+=$df[2] if $df[1] =~ /TP$/;
	
	$fp++           if $df[1] =~ /FP$/;
	$fpNucs+=$df[2] if $df[1] =~ /FP$/;
	
	$fn     = $totalTrue     - $tp;
	$fnNucs = $totalTrueNucs - $tpNucs; 
	
	$tnNucs = $totalNucs*1000000 - $fnNucs - $fpNucs - $tpNucs;
	######SENS,PPV,MCC,FPR
	my $sens = calcSENS($tp, $fn);
	my $ppv  = calcPPV( $tp, $fp);
	my $mcc  = calcMCC($tp, $tnNucs, $fp, $fn);
	my $fpr  = $fp/$totalNucs;
	$maxMcc=$mcc if $maxMcc<$mcc;
	($lastScore,$lastScoreNucs)    =($df[0],$df[0]) if $cnt == 0;
	#Only print if significantly different from previous computed vals.
	if( abs($lastSens - $sens)>0.0001 && abs($lastPpv - $ppv)>0.0001 && $lastScore!=$df[0]){
	    printf DO "%0.6f\t%0.6f\t%0.6f\t$df[0]\t%0.6f\n", $sens, $ppv, $fpr, $mcc;
	    $lastSens  = $sens; 
	    $lastPpv   = $ppv;
	    $lastScore = $df[0];
	}
	
	if($dfile=~/all-trueEmbedded.gff\S+1.dat/){
	    #print "tpNucs,tnNucs,fpNucs,fnNucs = $tpNucs, $tnNucs, $fpNucs, $fnNucs\n" if $cnt < 5;
	    my $sensNucs = calcSENS($tpNucs, $fnNucs);
	    my $ppvNucs  = calcPPV( $tpNucs, $fpNucs);
	    my $mccNucs  = calcMCC($tpNucs, $tnNucs, $fpNucs, $fnNucs);
	    my $fprNucs  = $fpNucs/$totalNucs;
	    #printf "%0.6f\t%0.6f\t%0.6f\t$df[0]\t%0.6f\n", $sensNucs, $ppvNucs, $fprNucs, $mccNucs if $cnt < 5;
	    
	    if(abs($lastMccNucs-$mccNucs)>0.0001 && $lastScoreNucs!=$df[0]){
		printf DO2 "%0.6f\t%0.6f\t%0.6f\t$df[0]\t%0.6f\n", $sensNucs, $ppvNucs, $fprNucs, $mccNucs;
		$lastMccNucs   = $mccNucs; 
		$lastScoreNucs = $df[0];
	    }
	}

	$cnt++;
    }

    printf "maxMcc:[%0.8f]\t$dfile\.rin\n", $maxMcc;
    close(DF);
    close(DO);
    close(DO2) if $dfile=~/all-trueEmbedded.gff\S+1.dat/;
}
 
system("R CMD BATCH --no-save plotROC.R");
system("egrep \47Print\|min\\(\|max\\(\|\^\\[1\\]\47 plotROC.Rout"); 

exit(0);
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

#Compute the Mathew's correlation coefficient from TP, TN, FP & FN counts:
sub calcMCC {
    my ($tp, $tn, $fp, $fn) = @_;
    #Array contains TP, TN, FP and FN:
    
    my $denom = sqrt( ($tp + $fp)*($tp + $fn)*($tn + $fp)*($tn + $fn) );
    my $mcc = 0.00;
    $mcc = ($tp * $tn - $fp * $fn)/$denom if $denom > 0;
    return $mcc;
}

#Compute the Sensitivity from TP & FN counts:
sub calcSENS {
    my ($tp, $fn) = @_;
    #TPR = TP / P = TP / (TP + FN)
    
    my $denom = ($tp + $fn);
    my $sens = 0.00;
    $sens = $tp/$denom if $denom > 0;
    return $sens;
}

#Compute the PPV from TP & FP counts:
sub calcPPV {
    my ($tp, $fp) = @_;
    #PPV = TP / (TP + FP)
    
    my $denom = ($tp + $fp);
    my $ppv = 1.00;
    $ppv = $tp/$denom if $denom > 0;
    return $ppv;
}
