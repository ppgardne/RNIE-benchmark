#cd $RNIE/benchmark/genomes
#R CMD BATCH --no-save ../scripts/plotDists.R

# # AE000511.1   Hpylori    Helicobacter pylori 26695
# # AE000513.1   Dradiodur  Deinococcus radiodurans R1
# # AE000516.2   Mtubercul  Mycobacterium tuberculosis CDC1551
# # AE001363.1   Cpneum	  Chlamydophila pneumoniae CWL029
# # AE009951.2   Fnucleat	  Fusobacterium nucleatum subsp. nucleatum ATCC 25586
# # AE014613.1   Styphi	  Salmonella enterica subsp. enterica serovar Typhi str. Ty2
# # AE015928.1   Bthetaio	  Bacteroides thetaiotaomicron VPI-5482
# # AE016823.1   Linterr	  Leptospira interrogans serovar Copenhageni str. Fiocruz L1-130
# # AE017126.1   Pmarinus	  Prochlorococcus marinus subsp. marinus str. CCMP1375
# # AF222894.1   Uparvum	  Ureaplasma parvum serovar 3 str. ATCC 700970
# # AL009126.3   Bsubtilis  Bacillus subtilis subsp. subtilis str. 168
# # AM180355.1   Cdiff	  Clostridium difficile 630
# # AP009493.1   Sgris	  Streptomyces griseus subsp. griseus NBRC 13350
# # CP000771.1   Fnodo	  Fervidobacterium nodosum Rt17-B1
# # CP000975.1   Minfern	  Methylacidiphilum infernorum V4
# # CP001147.1   Tyellow	  Thermodesulfovibrio yellowstonii DSM 11347
# # U00096.2     Ecoli	  Escherichia coli str. K-12 substr. MG1655

#grep ^"# # " plotDists.R | perl -lane 'print "rnie$F[3]\t<-read.table(\"$F[2].fasta-rnie.dists\",header = F, sep = \"\\t\")"'
#grep ^"# # " plotDists.R | perl -lane 'print "rnie$F[3]\[,2\]"' | tr "\n" ","
#grep ^"# # " plotDists.R | perl -lane 'print "rnie$F[3]H<-\thist(rnie$F[3]\[,2\],breaks=breaks,plot=F)"' 
#grep ^"# # " plotDists.R | perl -lane 'print "points(rnie$F[3]H\$mids, rnie$F[3]H\$counts, pch=\47o\47, col=\47red\47); \tlines(rnie$F[3]H\$mids, rnie$F[3]H\$counts,col=\47red\47)"' 

rnieHpylori     <-read.table("AE000511.1.fasta-rnie.dists",header = F, sep = "\t")
rnieDradiodur   <-read.table("AE000513.1.fasta-rnie.dists",header = F, sep = "\t")
rnieMtubercul   <-read.table("AE000516.2.fasta-rnie.dists",header = F, sep = "\t")
rnieCpneum      <-read.table("AE001363.1.fasta-rnie.dists",header = F, sep = "\t")
rnieFnucleat    <-read.table("AE009951.2.fasta-rnie.dists",header = F, sep = "\t")
rnieStyphi      <-read.table("AE014613.1.fasta-rnie.dists",header = F, sep = "\t")
rnieBthetaio    <-read.table("AE015928.1.fasta-rnie.dists",header = F, sep = "\t")
rnieLinterr     <-read.table("AE016823.1.fasta-rnie.dists",header = F, sep = "\t")

rniePmarinus    <-read.table("AE017126.1.fasta-rnie.dists",header = F, sep = "\t")
rnieUparvum     <-read.table("AF222894.1.fasta-rnie.dists",header = F, sep = "\t")
rnieBsubtilis   <-read.table("AL009126.3.fasta-rnie.dists",header = F, sep = "\t")
rnieCdiff       <-read.table("AM180355.1.fasta-rnie.dists",header = F, sep = "\t")
rnieSgris       <-read.table("AP009493.1.fasta-rnie.dists",header = F, sep = "\t")
rnieFnodo       <-read.table("CP000771.1.fasta-rnie.dists",header = F, sep = "\t")
rnieMinfern     <-read.table("CP000975.1.fasta-rnie.dists",header = F, sep = "\t")
rnieTyellow     <-read.table("CP001147.1.fasta-rnie.dists",header = F, sep = "\t")
rnieEcoli       <-read.table("U00096.2.fasta-rnie.dists",header = F, sep = "\t")



mx<-ceiling(max(c(rnieHpylori[,2],rnieDradiodur[,2],rnieMtubercul[,2],rnieCpneum[,2],rnieFnucleat[,2],rnieStyphi[,2],rnieBthetaio[,2],rnieLinterr[,2],rniePmarinus[,2],rnieUparvum[,2],rnieBsubtilis[,2],rnieCdiff[,2],rnieSgris[,2],rnieFnodo[,2],rnieMinfern[,2],rnieTyellow[,2],rnieEcoli[,2]))+1)
mn<-  floor(min(c(rnieHpylori[,2],rnieDradiodur[,2],rnieMtubercul[,2],rnieCpneum[,2],rnieFnucleat[,2],rnieStyphi[,2],rnieBthetaio[,2],rnieLinterr[,2],rniePmarinus[,2],rnieUparvum[,2],rnieBsubtilis[,2],rnieCdiff[,2],rnieSgris[,2],rnieFnodo[,2],rnieMinfern[,2],rnieTyellow[,2],rnieEcoli[,2]))-1)
breaks <- seq(-500,500,by=15)
breaks <- c(mn, breaks, mx)
rnieHpyloriH<-  hist(rnieHpylori[,2],breaks=breaks,plot=F)
rnieDradiodurH<-hist(rnieDradiodur[,2],breaks=breaks,plot=F)
rnieMtuberculH<-hist(rnieMtubercul[,2],breaks=breaks,plot=F)
rnieCpneumH<-   hist(rnieCpneum[,2],breaks=breaks,plot=F)
rnieFnucleatH<- hist(rnieFnucleat[,2],breaks=breaks,plot=F)
rnieStyphiH<-   hist(rnieStyphi[,2],breaks=breaks,plot=F)
rnieBthetaioH<- hist(rnieBthetaio[,2],breaks=breaks,plot=F)

rnieLinterrH<-  hist(rnieLinterr[,2],breaks=breaks,plot=F)
rniePmarinusH<- hist(rniePmarinus[,2],breaks=breaks,plot=F)
rnieUparvumH<-  hist(rnieUparvum[,2],breaks=breaks,plot=F)
rnieBsubtilisH<-hist(rnieBsubtilis[,2],breaks=breaks,plot=F)
rnieCdiffH<-    hist(rnieCdiff[,2],breaks=breaks,plot=F)
rnieSgrisH<-    hist(rnieSgris[,2],breaks=breaks,plot=F)
rnieFnodoH<-    hist(rnieFnodo[,2],breaks=breaks,plot=F)
rnieMinfernH<-  hist(rnieMinfern[,2],breaks=breaks,plot=F)
rnieTyellowH<-  hist(rnieTyellow[,2],breaks=breaks,plot=F)
rnieEcoliH<-    hist(rnieEcoli[,2],breaks=breaks,plot=F)

######################################################################
#Negative controls:
#1. permuted terminator coordinates:
#grep ^"# # " plotDists.R | perl -lane 'print "rnie$F[3]Perm\t<-read.table(\"$F[2].fasta-rnie-permuted.dists\",header = F, sep = \"\\t\")"'
#grep ^"# # " plotDists.R | perl -lane 'print "rnie$F[3]Perm\[,2\]"' | tr "\n" ","
#2. shuffled genome annotations:
#grep ^"# # " plotDists.R | perl -lane 'print "rnie$F[3]Shuff\t<-read.table(\"$F[2]-shuffled.fasta-rnie.dists\",header = F, sep = \"\\t\")"' | egrep -v 'AE000513|AE000516|AF222894|AP009493'
#grep ^"# # " plotDists.R | perl -lane 'print "rnie$F[3]Shuff\[,2\]\t$F[2]"' | egrep -v 'AE000513|AE000516|AF222894|AP009493' | cut -f 1 | tr "\n" ","  

rnieHpyloriShuff        <-read.table("AE000511.1-shuffled.fasta-rnie.dists",header = F, sep = "\t")
rnieCpneumShuff <-read.table("AE001363.1-shuffled.fasta-rnie.dists",header = F, sep = "\t")
rnieFnucleatShuff       <-read.table("AE009951.2-shuffled.fasta-rnie.dists",header = F, sep = "\t")
rnieStyphiShuff <-read.table("AE014613.1-shuffled.fasta-rnie.dists",header = F, sep = "\t")
rnieBthetaioShuff       <-read.table("AE015928.1-shuffled.fasta-rnie.dists",header = F, sep = "\t")
rnieLinterrShuff        <-read.table("AE016823.1-shuffled.fasta-rnie.dists",header = F, sep = "\t")

rniePmarinusShuff       <-read.table("AE017126.1-shuffled.fasta-rnie.dists",header = F, sep = "\t")
rnieBsubtilisShuff      <-read.table("AL009126.3-shuffled.fasta-rnie.dists",header = F, sep = "\t")
rnieCdiffShuff  <-read.table("AM180355.1-shuffled.fasta-rnie.dists",header = F, sep = "\t")
rnieFnodoShuff  <-read.table("CP000771.1-shuffled.fasta-rnie.dists",header = F, sep = "\t")
rnieMinfernShuff        <-read.table("CP000975.1-shuffled.fasta-rnie.dists",header = F, sep = "\t")
rnieTyellowShuff        <-read.table("CP001147.1-shuffled.fasta-rnie.dists",header = F, sep = "\t")
rnieEcoliShuff  <-read.table("U00096.2-shuffled.fasta-rnie.dists",header = F, sep = "\t")

#######

#grep ^"# # " plotDists.R |  perl -lane 'print "../genomes/$F[2].transterm.dists"' | tr "\n" " "
#cat ../genomes/AE000511.1.transterm.dists ../genomes/AE000513.1.transterm.dists ../genomes/AE000516.2.transterm.dists ../genomes/AE001363.1.transterm.dists ../genomes/AE009951.2.transterm.dists ../genomes/AE014613.1.transterm.dists ../genomes/AE015928.1.transterm.dists ../genomes/AE016823.1.transterm.dists ../genomes/AE017126.1.transterm.dists ../genomes/AF222894.1.transterm.dists ../genomes/AL009126.3.transterm.dists ../genomes/AM180355.1.transterm.dists ../genomes/AP009493.1.transterm.dists ../genomes/CP000771.1.transterm.dists ../genomes/CP000975.1.transterm.dists ../genomes/CP001147.1.transterm.dists ../genomes/U00096.2.transterm.dists > ../genomes/ALL.transterm.dists
#cat ../genomes/AE000511.1-shuffled.transterm.dists ../genomes/AE000513.1-shuffled.transterm.dists ../genomes/AE000516.2-shuffled.transterm.dists ../genomes/AE001363.1-shuffled.transterm.dists ../genomes/AE009951.2-shuffled.transterm.dists ../genomes/AE014613.1-shuffled.transterm.dists ../genomes/AE015928.1-shuffled.transterm.dists ../genomes/AE016823.1-shuffled.transterm.dists ../genomes/AE017126.1-shuffled.transterm.dists ../genomes/AF222894.1-shuffled.transterm.dists ../genomes/AL009126.3-shuffled.transterm.dists ../genomes/AM180355.1-shuffled.transterm.dists ../genomes/AP009493.1-shuffled.transterm.dists ../genomes/CP000771.1-shuffled.transterm.dists ../genomes/CP000975.1-shuffled.transterm.dists ../genomes/CP001147.1-shuffled.transterm.dists ../genomes/U00096.2-shuffled.transterm.dists > ../genomes/ALL-shuffled.transterm.dists
#grep ^"# # " plotDists.R |  egrep -iv 'coli|subtilis' | perl -lane 'print "../genomes/$F[2].transterm.dists"' | tr "\n" " "
#cat ../genomes/AE000511.1.transterm.dists ../genomes/AE000513.1.transterm.dists ../genomes/AE000516.2.transterm.dists ../genomes/AE001363.1.transterm.dists ../genomes/AE009951.2.transterm.dists ../genomes/AE014613.1.transterm.dists ../genomes/AE015928.1.transterm.dists ../genomes/AE016823.1.transterm.dists ../genomes/AE017126.1.transterm.dists ../genomes/AF222894.1.transterm.dists ../genomes/AM180355.1.transterm.dists ../genomes/AP009493.1.transterm.dists ../genomes/CP000771.1.transterm.dists ../genomes/CP000975.1.transterm.dists ../genomes/CP001147.1.transterm.dists                                                                           > ../genomes/SOME.transterm.dists
#cat ../genomes/AE000511.1-shuffled.transterm.dists ../genomes/AE000513.1-shuffled.transterm.dists ../genomes/AE000516.2-shuffled.transterm.dists ../genomes/AE001363.1-shuffled.transterm.dists ../genomes/AE009951.2-shuffled.transterm.dists ../genomes/AE014613.1-shuffled.transterm.dists ../genomes/AE015928.1-shuffled.transterm.dists ../genomes/AE016823.1-shuffled.transterm.dists ../genomes/AE017126.1-shuffled.transterm.dists ../genomes/AF222894.1-shuffled.transterm.dists ../genomes/AM180355.1-shuffled.transterm.dists ../genomes/AP009493.1-shuffled.transterm.dists ../genomes/CP000771.1-shuffled.transterm.dists ../genomes/CP000975.1-shuffled.transterm.dists ../genomes/CP001147.1-shuffled.transterm.dists                                                                           > ../genomes/SOME-shuffled.transterm.dists

transtermALL            <-read.table("ALL.transterm.dists",header = F, sep = "\t")
transtermSOME           <-read.table("SOME.transterm.dists",header = F, sep = "\t")
transtermALLShuff       <-read.table("ALL-shuffled.transterm.dists",header = F, sep = "\t")
transtermSOMEShuff      <-read.table("SOME-shuffled.transterm.dists",header = F, sep = "\t")


breaks <- seq(-800,800,by=2)
breaks <- c(mn, breaks, mx)

rnieAllH<-    hist(c(rnieHpylori[,2],rnieDradiodur[,2],rnieMtubercul[,2],rnieCpneum[,2],rnieFnucleat[,2],rnieStyphi[,2],rnieBthetaio[,2],rnieLinterr[,2],rniePmarinus[,2],rnieUparvum[,2],rnieBsubtilis[,2],rnieCdiff[,2],rnieSgris[,2],rnieFnodo[,2],rnieMinfern[,2],rnieTyellow[,2],rnieEcoli[,2]),breaks=breaks,plot=F)
rnieSomeH<-    hist(c(rnieHpylori[,2],rnieDradiodur[,2],rnieMtubercul[,2],rnieCpneum[,2],rnieFnucleat[,2],rnieStyphi[,2],rnieBthetaio[,2],rnieLinterr[,2],rniePmarinus[,2],rnieUparvum[,2],rnieCdiff[,2],rnieSgris[,2],rnieFnodo[,2],rnieMinfern[,2],rnieTyellow[,2]),breaks=breaks,plot=F)

rnieAllPermH<-    hist(c(rnieHpyloriPerm[,2],rnieDradiodurPerm[,2],rnieMtuberculPerm[,2],rnieCpneumPerm[,2],rnieFnucleatPerm[,2],rnieStyphiPerm[,2],rnieBthetaioPerm[,2],rnieLinterrPerm[,2],rniePmarinusPerm[,2],rnieUparvumPerm[,2],rnieBsubtilisPerm[,2],rnieCdiffPerm[,2],rnieSgrisPerm[,2],rnieFnodoPerm[,2],rnieMinfernPerm[,2],rnieTyellowPerm[,2],rnieEcoliPerm[,2]),breaks=breaks,plot=F)
rnieSomePermH<-    hist(c(rnieHpyloriPerm[,2],rnieDradiodurPerm[,2],rnieMtuberculPerm[,2],rnieCpneumPerm[,2],rnieFnucleatPerm[,2],rnieStyphiPerm[,2],rnieBthetaioPerm[,2],rnieLinterrPerm[,2],rniePmarinusPerm[,2],rnieUparvumPerm[,2],rnieCdiffPerm[,2],rnieSgrisPerm[,2],rnieFnodoPerm[,2],rnieMinfernPerm[,2],rnieTyellowPerm[,2]),breaks=breaks,plot=F)

rnieAllShuffH<-    hist(c(rnieHpyloriShuff[,2],rnieCpneumShuff[,2],rnieFnucleatShuff[,2],rnieStyphiShuff[,2],rnieBthetaioShuff[,2],rnieLinterrShuff[,2],rniePmarinusShuff[,2],rnieBsubtilisShuff[,2],rnieCdiffShuff[,2],rnieFnodoShuff[,2],rnieMinfernShuff[,2],rnieTyellowShuff[,2],rnieEcoliShuff[,2]),breaks=breaks,plot=F)
rnieSomeShuffH<-   hist(c(rnieHpyloriShuff[,2],rnieCpneumShuff[,2],rnieFnucleatShuff[,2],rnieStyphiShuff[,2],rnieBthetaioShuff[,2],rnieLinterrShuff[,2],rniePmarinusShuff[,2],rnieCdiffShuff[,2],rnieFnodoShuff[,2],rnieMinfernShuff[,2],rnieTyellowShuff[,2]),breaks=breaks,plot=F)

transtermALLH<-    hist(transtermALL[,2],breaks=breaks,plot=F)
transtermSOMEH<-   hist(transtermSOME[,2],breaks=breaks,plot=F)
transtermALLShuffH<-    hist(transtermALLShuff[,2],breaks=breaks,plot=F)
transtermSOMEShuffH<-   hist(transtermSOMEShuff[,2],breaks=breaks,plot=F)

totals<-c(
sum(transtermSOMEShuffH$counts), 
sum(transtermALLShuffH$counts), 
sum(transtermSOMEH$counts), 
sum(transtermALLH$counts), 
sum(rnieSomeShuffH$counts), 
sum(rnieAllShuffH$counts), 
sum(rnieSomeH$counts), 
sum(rnieAllH$counts)
)

totalsM<-mat.or.vec(length(totals), length(totals))
for (i in seq(1,length(totals))){
    totalsM[i,i]<-log10(totals[i])
}

pdf(file="geneProximity.pdf")

op<-par(mfrow=c(1,1),cex=0.75,las=1,ylog='TRUE',mai=c(1.02, 0.82, 0.82, 0.42))
#plot(rnieAllH$mids, log10(rnieAllH$counts+1e-10), xlim=c(-100,500), ylim=c(0,3), pch='x', col='black', main='Terminator proximity to annotated genic features', xlab='Minimum distance from a 3\47 end of all genic features (nucs)', ylab='Frequency', yaxt="n")
plot(lowess(rnieAllH$mids[2:(length( rnieAllH$mids)-1)],log10(rnieAllH$counts[2:(length( rnieAllH$mids)-1)]+1), f = 1/30), xlim=c(-100,500), ylim=c(0,3), type='l', col='black', main='Terminator proximity to annotated genic features', xlab='Minimum distance from a 3\47 end of all genic features (nucs)', ylab='Frequency', yaxt="n", lwd=3)
axis(2,at=c(0,log10(5),1,log10(50),2,log10(500),3),labels=c(0,5,10,50,100,500,1000), las=1)
#lines(lowess(rnieAllH$mids[2:(length( rnieAllH$mids)-1)],log10(rnieAllH$counts[2:(length( rnieAllH$mids)-1)]+1), f = 1/30),col='black',lwd=3)

#points(rnieSomeH$mids, log10(rnieSomeH$counts+1e-10), pch='+', col='grey'); 
lines(lowess(rnieSomeH$mids[2:(length( rnieAllH$mids)-1)], log10(rnieSomeH$counts[2:(length( rnieAllH$mids)-1)]+1), f = 1/30),col='grey',lwd=3)
lines(lowess(rnieAllShuffH$mids[2:(length( rnieAllH$mids)-1)],  log10(rnieAllShuffH$counts[2:(length( rnieAllH$mids)-1)]+1), f = 1/30),col='black', lty=2,lwd=3)
lines(lowess(rnieSomeShuffH$mids[2:(length( rnieAllH$mids)-1)], log10(rnieSomeShuffH$counts[2:(length( rnieAllH$mids)-1)]+1), f = 1/30),col='grey', lty=2,lwd=3)

#points(transtermALLH$mids,  log10(transtermALLH$counts+1e-10),      pch='x', col='green');          
lines(lowess(transtermALLH$mids[2:(length( rnieAllH$mids)-1)],  log10(transtermALLH$counts[2:(length( rnieAllH$mids)-1)]+1), f = 1/30),col='green',lwd=3)
#points(transtermSOMEH$mids, log10(transtermSOMEH$counts+1e-10),     pch='+', col='lightgreen');     
lines(lowess(transtermSOMEH$mids[2:(length( rnieAllH$mids)-1)], log10(transtermSOMEH$counts[2:(length( rnieAllH$mids)-1)]+1), f = 1/30),col='lightgreen',lwd=3)
lines(lowess(transtermALLShuffH$mids[2:(length( rnieAllH$mids)-1)],  log10(transtermALLShuffH$counts[2:(length( rnieAllH$mids)-1)]+1), f = 1/30),col='green', lty=2,lwd=3)
lines(lowess(transtermSOMEShuffH$mids[2:(length( rnieAllH$mids)-1)], log10(transtermSOMEShuffH$counts[2:(length( rnieAllH$mids)-1)]+1), f = 1/30),col='lightgreen', lty=2,lwd=3)

leg<-legend(123, 3 ,c(
"A. RNIE (all genomes)",
"B. RNIE (no E. coli & B. subtilis)",
"C. RNIE on permuted genomes (all genomes)",
"D. RNIE on permuted genomes (no E. coli & B. subtilis)",
"E. TransTermHP (all genomes)",
"F. TransTermHP (no E. coli & B. subtilis)",
"G. TransTermHP on permuted genomes (all genomes)",
"H. TransTermHP on permuted genomes (no E. coli & B. subtilis)"),
col=c("black","grey","black","grey","green","lightgreen","green","lightgreen"),pch=c('x','+','','','x','+','',''),lty=c(1,1,2,2,1,1,2,2),lwd=c(3,3,3,3,3,3,3,3),ncol=1,cex=0.75)
lines(c(0,0), c(-500,500),col='black', lty=3)

names<-c(
"TransTerm.noEC/BS.perm",
"TransTerm.perm",
"TransTerm.noEC/BS",
"TransTerm",
"RNIE.noEC/BS.perm",
"RNIE.perm",
"RNIE.noEC/BS",
"RNIE"
)
cols<-c(
"lightgreen",
"green",
"lightgreen",
"green",
"grey",
"black",
"grey",
"black"
)
#densities<-c(1000,1000,50,50,1000,1000,50,50)
densities<-c(50,50,1000,1000,50,50,1000,1000)
par(fig=c(0.47,0.89,0.48,0.70),mai=c(0.2,0.2,0.09,0),las=1, new=TRUE,cex=0.6)
#c(x1, x2, y1, y2)
#mai: bottom, left, top, right)
#pars$mai 1.02 0.82 0.82 0.42
barplot(totalsM, xlab="",names.arg=c("H","G","F","E","D","C","B","A"), col=cols, density=densities, horiz='TRUE', main="Total number of predictions", xaxt="n")
#,names.arg=names
axis(1,at=c(0,log10(5),1,log10(50),2,log10(500),3,log10(5000)),labels=c(0,"",10,"",100,"",1000,""), las=1)

dev.off()

######################################################################














