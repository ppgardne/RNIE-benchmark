#!/usr/bin/perl 

#read in gff files, return dat files, sorted on score and annotated as
#either FP or TP, tag FN's to the end of the file.
#

use warnings;
use strict;

#USE THE RNIE PARAMETER!
my $rnieHome = '';
$rnieHome = $ENV{'RNIE'} . '/benchmark/training-test/' if defined $ENV{'RNIE'};

# Embedded-fast-rnie.bits
# Embedded-fast-rnie.E_value
# Embedded-slow-rnie.bits
# Embedded-slow-rnie.E_value
# Embedded-slow-rnie.super43-45-seed-1.bits
# Embedded-slow-rnie.terminator32-40-1.bits
# Embedded-slow-rnie.terminator40-43-1.bits
# Embedded-slow-rnie.terminator43-45-1.bits

my @files = qw(
Embedded-genomeMode-rnie.bits
Embedded-geneMode-rnie.bits
Embedded-fast-rnie.bits
Embedded-slow-rnie.bits
Embedded-2features.transterm
Embedded-4features.transterm
Embedded-9features.transterm
Embedded-10features.transterm
Embedded.rnall-barrick.dG
Embedded.rnall-barrick.hbG
Embedded.rnall-wan.dG
Embedded.rnall-wan.hbG
Embedded.terminator-lesnik.out.dG_score
Embedded.terminator-lesnik.out.struct_score
);

my $totalNucs     = 431453 + 43145150;

foreach my $fileInfix (@files){ 
    my ($roc, $rocNucs) = gff2accuracy($totalNucs, $rnieHome. 'trueEmbedded.gff', $rnieHome.'true' . $fileInfix . '.gff', $rnieHome.'false' . $fileInfix . '.gff');
    printDat($fileInfix . '-accuracy.dat', $roc);
    printDat($fileInfix . '-accuracyNucs.dat', $rocNucs);
}
system("R CMD BATCH --no-save $rnieHome/../scripts/plotROC.R");

exit(0);

######################################################################
#Print an accuracy datastructure out to file:
sub printDat {
    my ($fileName, $roc)=@_;
    
    open(DAT, "> $rnieHome$fileName") or die "FATAL: failed to open [$rnieHome$fileName] for writing!\n[$!]";
    printf DAT "%0.6s\t%0.6s\t%0.6s\t%0.6s\t%0.6s\n", 'sens', 'ppv', 'mcc', 'fpr/KB', 'score';
    foreach my $acc (@{$roc}){
	printf DAT "%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\n", $acc->{'sens'}, $acc->{'ppv'}, $acc->{'mcc'}, 1000*$acc->{'fpr'}, $acc->{'score'};
    }
    close(DAT);
}

######################################################################
#Given 2 (or more) gff files compute accuracy statistics for the method:
# gff2accuracy($trueAnnotationGffFile, $predictionGffFile1, $predictionGffFile2, ...)
# return ROC arrays (scores, sensitivity, ppv, fpr, mcc, )
sub gff2accuracy {
    my $totalNucs = shift; 
    my @files = @_;
    
    my $trueFile = shift(@files);
    
    
    
    open(TR, "< $trueFile") or die "FATAL: failed to open file [$trueFile] for reading!\n[$!]";
    my %trueCoords;
    my ($totalTrue, $totalTrueNucs) = (0,0); 
    while(<TR>){
	my @gff = split(/\t/, $_);
	if($gff[0]=~/^\S+$/ && $gff[3]=~/^\d+$/ && $gff[4]=~/^\d+$/){
	    push(@{ $trueCoords{$gff[0]} }, {'start' => $gff[3], 'end' => $gff[4], 'strand' => $gff[6]});
	    $totalTrue++;
	    $totalTrueNucs += (abs($gff[3] - $gff[4])+1);
	    #print "TRUE:trueCoords{$gff[0]} }, {\'start\' => $gff[3], \'end\' => $gff[4], \'strand\' => $gff[6]\n";
	}
    }
    close(TR);
    
    my ($tp,$fp,$fn)                      = (0,0,0);
    my ($tpNucs,$tnNucs,$fpNucs,$fnNucs)  = (0,0,0,0);
    my ($lastSens,$lastPpv,$lastScore)    = (0,0,0);   #Only record values that differ significantly from previous records
    my ($lastMccNucs,$maxMcc,$maxMccNucs) = (0,0,0);    
    my (@roc, @rocNucs);
    my $cnt = 0;
    my $predFiles = join(' ', @files);
    #MYAHAHAHA:
    open( PF, "cat $predFiles | /usr/bin/sort -k6nr |" ) or die "FATAL: failed to open pipe [cat $rnieHome$predFiles | /usr/bin/sort -k7n |]!\n[$!]";
    while(my $in = <PF>){
	my @gff = split(/\t/, $in);
	($lastScore)    =($gff[5]) if $cnt == 0;
	die "FATAL: check the file sort on [@files], lastScore[$lastScore] < nowScore[$gff[5]]!" if $lastScore < $gff[5];
	my $olNucs = 0;
	foreach my $trueElement ( @{ $trueCoords{$gff[0]} } ){
	    next if $trueElement->{'strand'} ne $gff[6];
	    #print "overlapExtentNucs($gff[3],$gff[4]," . $trueElement->{'start'} .",". $trueElement->{'end'}.")";
	    $olNucs += overlapExtentNucs($gff[3],$gff[4],$trueElement->{'start'}, $trueElement->{'end'}) if ($trueElement->{'start'} and $trueElement->{'end'});
	}
	
	if($olNucs){
	    $tp++;
	    $tpNucs+=$olNucs;
	}
	else {
	    $fp++;
	    $fpNucs+=(abs($gff[3] - $gff[4])+1);
	}
	
	$fn     = $totalTrue     - $tp;
	$fnNucs = $totalTrueNucs - $tpNucs; 
	$tnNucs = $totalNucs - $fnNucs - $fpNucs - $tpNucs;
	
	######PER FEATURE SENS,PPV,MCC,FPR
	my $sens = calcSENS($tp, $fn);
	my $ppv  = calcPPV( $tp, $fp);
	my $mcc  = calcMCC($tp, $tnNucs, $fp, $fn);
	my $fpr  = $fp/$totalNucs;
	$maxMcc=$mcc if $maxMcc<$mcc;
	
	######PER NUCLEOTIDE SENS,PPV,MCC,FPR
	my $sensNucs = calcSENS($tpNucs, $fnNucs);
	my $ppvNucs  = calcPPV( $tpNucs, $fpNucs);
	my $mccNucs  = calcMCC($tpNucs, $tnNucs, $fpNucs, $fnNucs);
	my $fprNucs  = $fpNucs/$totalNucs;
	$maxMccNucs=$mcc if $maxMccNucs<$mccNucs;
	
	#print "PRED:$gff[0] }, {\'start\' => $gff[3], \'end\' => $gff[4], \'strand\' => $gff[6]\n";
	#print "sens[$sens] ppv[$ppv] mcc[$mcc] fpr[$fpr]\n";
	#Only print if significantly different from previous computed vals.
	if( (abs($lastSens - $sens)>0.001 || abs($lastPpv - $ppv)>0.001 || abs($lastMccNucs-$mccNucs)>0.001) && ($lastScore!=$gff[5])){
	    
	    #SAVE STATE!
	    push(@roc,     {'sens' => $sens,     'ppv' => $ppv,     'mcc' => $mcc,     'fpr' => $fpr,     'score' => $gff[5]} );
	    push(@rocNucs, {'sens' => $sensNucs, 'ppv' => $ppvNucs, 'mcc' => $mccNucs, 'fpr' => $fprNucs, 'score' => $gff[5]} );
	    #print "   sens[$sens]        ppv[$ppv]        mcc[$mcc]        fpr[$fpr]\n";
	    #print "sensNuc[$sensNucs] ppvNuc[$ppvNucs] mccNuc[$mccNucs] fprNuc[$fprNucs]\n";
	    $lastSens  = $sens; 
	    $lastPpv   = $ppv;
	    $lastMccNucs   = $mccNucs; 
	    $lastScore = $gff[5];
	    
	}
	
	
	$cnt++;
    }
    close(PF);
    
    for (my $i=0; $i<(@files); $i++){
	$files[$i] = fullPath2FileName($files[$i]);
    }
    
    printf "max(MCC) = [%0.6f] max(MCCNucs) = [%0.6f] [@files]\n", $maxMcc, $maxMccNucs;
    
    return (\@roc, \@rocNucs);
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
# D = |A n B|
#
sub overlapExtentNucs {
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
    return $D;
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
#Compute the Mathew's correlation coefficient from TP, TN, FP & FN counts:
sub calcMCC {
    my ($tp, $tn, $fp, $fn) = @_;
    #Array contains TP, TN, FP and FN:
    
    my $denom = sqrt( ($tp + $fp)*($tp + $fn)*($tn + $fp)*($tn + $fn) );
    my $mcc = 0.00;
    $mcc = ($tp * $tn - $fp * $fn)/$denom if $denom > 0;
    $mcc =  1.0 if $mcc >  1.0;
    $mcc = -1.0 if $mcc < -1.0;
    return $mcc;
}

######################################################################
#Compute the Sensitivity from TP & FN counts:
sub calcSENS {
    my ($tp, $fn) = @_;
    #TPR = TP / P = TP / (TP + FN)
    
    my $denom = ($tp + $fn);
    my $sens = 0.00;
    $sens = $tp/$denom if $denom > 0;
    $sens = 1.0 if $sens>1.0;
    return $sens;
}

######################################################################
#Compute the PPV from TP & FP counts:
sub calcPPV {
    my ($tp, $fp) = @_;
    #PPV = TP / (TP + FP)
    
    my $denom = ($tp + $fp);
    my $ppv = 1.00;
    $ppv = $tp/$denom if $denom > 0;
    $ppv = 1.0 if $ppv>1.0;
    return $ppv;
}

######################################################################
#fullPath2FileName
sub fullPath2FileName {
    my $name = shift;
    my @dirs = split(/\//, $name);
    $name = pop(@dirs);
    return $name;
}
