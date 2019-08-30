	Starting 1/12/2018.
	Experiment 3 is to visualize the genomic distribution of diversity enhancement in core
	using the genome itself as the source of information, not geographical or environmental variables.
	Goal is to produce "heatmaps" showing the enhancement across each chromosome for all blocklengths,
	effectively mapping the genomewide correlation between diversity at a position and diversity
	in the genome at large.
	
	As with analysis 2, note 
	that the Arabidopsis data has been modified to 40 populations from the 41 before
	(two had identical coordinates, which I missed previously).  Also that the outlying
	Sorghum population in India has been removed, reducing the number of populations from 23
	to 22.
	
	No binning procedure is needed because genomic data are categorical. Most data is co-located with
	/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper
	
	

***POPULUS***

	Run complete analysis on Ceres and blip, make sure results are congruent. <--They were.
	Calculate mean M+ enrichment across core sizes, for each locus, for each blocklength. This uses the
	SUMfiles as input and produces the +STATS files.


**********BASHPARALLEL-FASTER-RUNINSCREEN**********
mysd() {
       #calculate sum, mean, sd, n, mean,using Welford's algorithm and AWK
       #here, sum is the cumulative enrichment across all core sizes, mean is the average
       #supply space delimited string of values as argument like: mysd $foo, where foo="1 2 3 4.5"
       awk '{
         sum = 0;                   # Initialize running sum (for mean calculation)
         M = 0;
         S = 0;
         for (k=1; k <= NF; k++) { 
              sum += $k;                # Update running sum
              x = $k;
              oldM = M;
              M = M + ((x - M)/k);
              S = S + (x - M)*(x - oldM);
         }
         var = S/(NF - 1);
         print sum "\t" sum/NF "\t" sqrt(var) "\t" NF;
       }' <<< "$1"; # the first argument when calling mysd 
}
export -f mysd;

main() {
        l=$1;
        ofile="$p""+STATS.R"$e"Tgen.loc"$l".b"$b".txt"; #temporary outfile for each parallel process
        sumfileIN="$p""SUM.b"$b".R"$e"Tgen.txt";

        boo=$(grep -F $l "$sumfileIN" | awk '{ print $13 - $10 }' | tr "\n" " "); #extract all core sizes for this locus, calculate target enrichment
        st=$(mysd "$boo"); #calculate sum, mean, sd, and n using function mysd()
        echo "$b"$'\t'"$l"$'\t'"$st" > "$ofile"; #write stats to output file
}
export -f main;

#parameters
t="Populus";
elist="gen";
bmin=2; #min blocklength
bmax=50; #max blocklength
bstep=2; #step rate
nnode=7; #number of cluster nodes
blist=$(seq $bmin $bstep $bmax); #blocklengthrange
p=$(pwd)"/"; 
export p;

#parallelize on locus, within loop on blocklength
for e in $elist;
  do echo $e;
    export e;
    for b in $blist;
      do echo "b=$b";
        export b;
        sumfileIN="$p""SUM.b"$b".R"$e"Tgen.txt";
        uloc=$(cut -d$'\t' -f2 "$sumfileIN" | tail -n +2 | uniq); #find unique loci from all loci in column 2
        numloc=$(echo "$uloc" | wc -l);
        
        #parallel step
        #sends $gnuN records to each node via parallel step #1, this is piped to parallel step #2, which starts 96 jobs per node from the instruction set it has received
        echo "processing $sumfileIN...";
        gnuN=$(echo "$numloc / $nnode" | bc); #number of lines of input file to send to gnu parallel per node
        echo "$uloc" | parallel --sshloginfile ~/machines --jobs 1 --env main --env mysd --env e --env b --env p --pipe -N"$gnuN" /home/reevesp/bin/parallel -j96 --env main --env mysd --env e --env b --env p main;
        #run on Mac
        #echo "$uloc" | parallel --env main --env mysd --env e --env b --env p main;

        #consolidate results for each locus into one file
        ofile="+STATS.R"$e"Tgen.b"$b".txt"; #temporary outfile
        > "$p""$ofile";
        for l in $uloc;
          do echo "concatenating loc file $l to +STATS.R"$e"Tgen.b"$b".txt ...";
            cat "$p""+STATS.R"$e"Tgen.loc"$l".b"$b".txt" >> "$p""$ofile";
          done;
          #test whether all were concatenated
          s1=$(ls "$p" | grep -c "loc"); #count the number of files with "loc" in their name using grep. can't wc -l because too many arguments
          s2=$(wc -l "$p""$ofile" |  awk '{print $1}');
         if [ $s1 == $s2 ]; 
           then find "$p" -maxdepth 1 -name "*loc*" -delete; #remove files when there are too many for rm
           else echo "s1="$s1", s2="$s2". Consolidated file "$ofile" no good. Quitting..." > err.txt;
           exit 1;
         fi;
      done; #$b
      
    #concatenate files for this Rgeo/Rpco
    statsfile="$p""+STATS."$t".R"$e"Tgen.txt";
    echo blocklength$'\t'locus$'\t'sum$'\t'mean$'\t'sd$'\t'n > "$statsfile";
    >"$statsfile"TMP;
    #quick cat and sort (some out files from ceres are not sorted correctly by locus)
    cat "$p""+STATS.R"$e"Tgen.b"*".txt" >> "$statsfile"TMP; #cat together the whole mess
    sortkey1=$(cut -d$'\t' -f2 "$statsfile"TMP | cut -d'.' -f3); #get a column containing only the locus index number, field 3 of column 2
    paste -d$'\t' "$statsfile"TMP <(echo "$sortkey1") | sort -t$'\t' -n -k1,1 -n -k7,7 | cut -d$'\t' -f1-6 > "$statsfile"TMP2; #add sort key to last column, sort by blocklength then locus index, remove last column
    cat "$statsfile"TMP2 >> "$statsfile";
    rm "$statsfile"TMP;
    rm "$statsfile"TMP2;

    #test whether all individual +STATS files have been added to the final concatenated file using number of lines
    s1=$(wc -l +STATS.R"$e"Tgen.b*.txt | grep total | awk '{print $1}');
    s2=$(wc -l "$statsfile" | awk '{print $1}');
    s2=$(( $s2 - 1 )); #decount header
    if [ $s1 == $s2 ]; 
      then rm +STATS.R"$e"Tgen.b*.txt; #remove files if number of lines in concat file is the same as the sum of all input files
      else echo "the number of lines ain't the same. something is 'crewed.";
    fi;

  done; #$e
**********BASHEND**********
	Takes about 6 hours on Mac. ~9 minutes on blip.

	add the proportion of NG,SS,NS sites from PopLocusWeights.txt to the +STATS summary file,
	generating the +v2STATS file.  Requires the PopLocusWeights.txt file as input.  

**********BASH**********
#use md5 on mac md5sum on linux
#modify weightsfile to include only the blocklengths analyzed
weightsfile="PopLocusWeights.txt";
bmin=2; #min blocklength
bmax=50; #max blocklength
bstep=2; #step rate
blist=$(seq $bmin $bstep $bmax); #blocklengthrange

wtfhead=$(head -1 "$weightsfile");
wtfbod=$(for i in $blist;
  do awk -F"\t" -v i=$i '$2 == i {print}' "$weightsfile"; #accumulate each blocklength needed
  done;)
wtf=$(echo "$wtfhead"; echo "$wtfbod");

#echo "$wtf" > wtf.txt; 

#wtf=$(awk -F"\t" '$2 != "1" {print}' "$weightsfile"); #exclude blocklength=1 from the weightsfile for use with blocklength 2-200
elist="gen";
for e in $elist;
  do statsfile="+STATS.Populus.R"$e"Tgen.txt";
    s=$(cut -d$'\t' -f2 $statsfile | tail -n +2 | md5sum); #get md5 of locus ID column from the +STATS file
    w1=$(echo "$wtf" | cut -d$'\t' -f4 | sed 's/_/./g'); #get the locus ID less blocklength index column, replace _ with . for md5 check
    w2=$(echo "$wtf" | cut -d$'\t' -f3); #get the blocklength index column
    w=$(paste -d'.' <(echo "$w1") <(echo "$w2") | tail -n +2 | md5sum); #get md5 of locus ID column from the Weights file
    if [ "$s" == "$w" ]; #verify that files are in same order using a checksum of the unique locus ID column
      then boo=$(echo "$wtf" | sed 's/blocklength/blocklength2/g' | sed 's/locus/locus2/g'); #modify the header in the Weights file so that there are no redundant column names after combining
        paste -d$'\t' "$statsfile" <(echo "$boo") > "+v2STATS.Populus.R"$e"Tgen.txt"; #paste the statsfile and the weights file together
      else echo "md5s do not match, $e, aborting.";
        kill -INT $$; #terminate the script, return to the shell
    fi;
  done;
**********BASHEND**********
	Takes a couple seconds

	Plot summed enrichment by genomic position across blocklengths as a heatmap. Use Rscript
	PopGenomicGeography.r. Uses the +v2STATS files as input.


	Using the file +*SNPtoBlockMap.txt, calculate the sum of M+ enrichment value across blocklengths,
	for each SNP, i.e. for each position in the genome.  These data will be extracted from the 
	+v2STATS* files in several steps.  The first makes a file containing a table with the max 
	enrichment value at each site for each block, *EnrichAcrossBlocks*. Step 1:

**********BASHPARALLEL**********
myp() {
  e=$1;
  b=$2;
  v2name="$p""+v2STATS."$v2p".R"$e"Tgen.txt";
        
        echo "e=$e  b=$b";
        bout="$p""rr.R"$e"Tgen.b""$b"".tmp"; #temporary outfile

        #cut out the locus id and sum columns for current blocklenth
        gg=$(grep ^$b$'\t' "$v2name" | cut -d$'\t' -f1-3);
        
        #cut column with locus names for current blocklength from SNPtoBlockMap
        if [ $b = 1 ]; then
          cc=$(cut -d$'\t' -f4 "$tname" | tail -n +2 | sed 's/_/./g' | sed 's/$/./g');
        else col=$(( $b + 5 )); #calculate the column to cut out for each blocklength
          cc=$(cut -d$'\t' -f$col "$tname" | tail -n +2 | sed 's/_/./g' | sed 's/$/./g');
        fi;

        #search for locus names in $gg, retrieve enrichment value
        rr="b$b";

#put another parallel here if necessary
        for l in $cc;
          do if [[ $l = "--." ]]; then
                r="--"; #test for empty enrichment value
              else
                r=$(grep -m 1 ^"$b"$'\t'"$l" <<<"$gg" | awk '{print $3}');
              fi;
            rr+=$'\n'$r; #add current locus enrichment value to list
          done;

         #write out the result
         echo "$rr" > "$bout";
         #truncate -s -1 "$bout"; #remove the hanging newline, use 'gtruncate' only on osx
}
export -f myp;

#parameters
v2p="Populus"; export v2p;
t="Pop";
minbl=2; #1,2
maxbl=50; #50,200
blstep=2; #step
blist=$(seq $minbl $blstep $maxbl); #for even numbered runs 2-50
elist="gen";

p=$(pwd)"/"; export p;
tname="$p""+""$t""SNPtoBlockMap.txt"; export tname;

#begin parallel
parallel --sshloginfile ~/machines --jobs 24 --env myp --env p --env v2p --env tname myp ::: $elist ::: $blist;
#parallel --jobs 24 --env myp --env p --env v2p --env tname myp ::: $elist ::: $blist; #parallel command for head node only
#parallel --env myp --env p --env v2p --env tname myp ::: $elist ::: $blist; #parallel command for Mac

  
#reconstitute
echo "reconstituting..."
for e in $elist;
  do echo "e=$e";
    comm=""; #create a list of tmp files to paste together
    for b in $blist;
      do comm+="rr.R"$e"Tgen.b$b.tmp ";
    done;
    pos=$(cut -d$'\t' -f4 "$tname" | sed 's/'^.*_'//g' | sed 's/locus/pos/g');
    h=$(cut -d$'\t' -f1-6 "$tname" | tr $'\t' ' ');
    h2=$(paste -d' ' <(echo "$h") <(echo "$pos"));
    paste -d' ' <(echo "$h2") $comm > "+"$t"EnrichAcrossBlocks.R"$e"Tgen.txt";
    rm $comm; #clean up
  done;#e
  
**********BASHEND**********
	Takes ~1.5hrs for maxbl=50 on Mac. Like 5 minutes on blip.

	Process each line (SNP position) of *EnrichAcrossBlocks*.txt to get sum across all blocklengths
	of cumulative enrichment across all core sizes of max enrichment value from 10 reps of M+. Step 2:

**********BASH**********
mysd() {
       #calculate sum, mean, sd, n, mean,using Welford's algorithm and AWK
       #here, sum is the cumulative enrichment across all blocklengths, mean is the average
       #supply space delimited string of values as argument like: mysd $foo, where foo="1 2 3 4.5"
       #returns space delimited string like "sum mean sd n"
       awk '{
         sum = 0;                   # Initialize running sum (for mean calculation)
         M = 0;
         S = 0;
         for (k=1; k <= NF; k++) { 
              sum += $k;                # Update running sum
              x = $k;
              oldM = M;
              M = M + ((x - M)/k);
              S = S + (x - M)*(x - oldM);
         }
         var = S/(NF - 1);
         print sum " " sum/NF " " sqrt(var) " " NF;
       }' <<< "$1"; # the first argument when calling mysd 
}
export -f mysd;

myp() {
    l=$1

    p=$(echo "$l" | cut -d' ' -f8-); #cut out the portion of the line with max enrichment values

    #skip loci that have missing blocks
    if [[ $p != *"--"* ]]; then
      m=$(mysd "$p"); #calculate sum, mean, sd, n
      chr=$(echo "$l" | cut -d' ' -f5);
      pos=$(echo "$l" | cut -d' ' -f7);
      i=$(echo "$l" | cut -d' ' -f1); #get genome index
      #loc=$chr"."$pos"."$(echo "$l" | cut -d' ' -f6);
      loc=$chr"."$pos"."$i;
      echo "snp=$i";

      echo $i" "$chr" "$pos" "$loc" "$m > "$path"$i".tmp"; #add new values to growing string, write to tmp file
    fi;
}
export -f myp;

#parameters
elist="gen";
t="Pop";
path=$(pwd)"/"; export path;

#parallel
for e in $elist;
do echo $e;
  #nnode=8;
  #numloc=$(wc -l +PopEnrichAcrossBlocks.RgenTgen.txt | awk '{print $1}');
  #numloc=$(( $numloc - 1 ));
  #gnuN=$(echo "$numloc / $nnode" | bc); #number of lines of input file to send to gnu parallel per node

  #parallel command
  #for cluster
  #tail -n +2 "+"$t"EnrichAcrossBlocks.R"$e"Tgen.txt" | parallel  --sshloginfile ~/machines --jobs 1 --env myp --env mysd --env path --pipe --round-robin /home/reevesp/bin/parallel -j96  --env myp --env mysd --env path myp;
  #for head node only
  tail -n +2 "+"$t"EnrichAcrossBlocks.R"$e"Tgen.txt" | parallel --jobs 1 --env myp --env mysd --env path --pipe --round-robin /home/reevesp/bin/parallel -j96  --env myp --env mysd --env path myp;
  #for mac
  #tail -n +2 "+"$t"EnrichAcrossBlocks.R"$e"Tgen.txt" | parallel  --env myp --env mysd --env path myp;
  
  #assemble results
  echo "concatenating...";
  ofile="+"$t"SUMEnrichAcrossBlocks.R"$e"Tgen.txt";
  o="genomeindex chr pos loc sum mean sd n";
  echo "$o" > "$ofile";
  
  find "$path" -name "*.tmp" -print0 | xargs -0 cat | sort -n -t' ' -k1 >> "$ofile";
  
  #clean up
  echo "clean up...";
  find "$path" -name "*.tmp" -print0 | xargs -0 rm;
done;
**********ENDBASH**********
	Takes ~3 minutes on Mac, seconds on blip.

	Use R to make some pretty plots of significant enrichment across the genome. Switch setwd between
	hclust and kmeans.

**********BEGINR**********
#install.packages("fitdistrplus")
#install.packages("plotrix")
options(error = recover)
rm(list=ls()) 
setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Populus/Experiment2/+revision\ with\ hclust/+results/hclust/+PopRgenRuns/rarefaction/+analysis")
library(fitdistrplus)
library(plotrix)

t<-"Pop"
usezero="no" #this sets whether to include significant sites only when they are also > or < 0
             #to use the zero criterion, say "yes", to just use the ecdf around the mean, say "no"
pvals<-c(0.001,0.05)
devvals<-c("pdf","svg","png")
for (pp in pvals)
{
	for (zg in devvals)
	{
		zgg<-get(zg) #e.g convert string "pdf" to function pdf now named 'zgg'
		if ( zg == "pdf" | zg == "svg" ) { zgg(file=paste(t,"M+GenomeScan",pp, ".", zg, sep="")) }
		else { 
		par(mar=c(1,1,1,1))
		zgg(file=paste(t,"M+GenomeScan",pp, ".", zg, sep=""),width=2000,height=2000,res=300) } #call zgg function with relevant parameters
		
		par(mfrow=c(2,1))

			i="gen"
			g <- read.table(paste("+",t,"SUMEnrichAcrossBlocks.R",i,"Tgen.txt", sep=""), header=TRUE, sep=" ")
			mycdf <- ecdf(g$sum) #calculate the cdf of the empirical distribution
			p <- mycdf(g$sum) #get the probability of each observation (or lower)
			g$pecdf <- p #add new column to g, containing the probability from the empirical cdf
			g$index <- 1:nrow(g) #add a column with the sequential index
			h <- g[which(g$pecdf > 1-pp),] #get snps that are highly enriched
			if ( usezero == "yes" ) {
				hgeo <- h[which(h$sum>0),] #remove snps where sum<0, i.e. less than random
			}
			else {
				hgeo <- h # do not remove snps based on positive or negative
			}
			l <- g[which(g$pecdf < pp),] #get snps that are lowly enriched
			if ( usezero == "yes" ) {
				lgeo <- l[which(l$sum<0),] #remove snps where sum>0, i.e. more than random
			}
			else {
				lgeo <- l # do not remove snps based on positive or negative
			}
			highposgeo<-hgeo$index #get the sequential position of the highly enriched
			lowposgeo<-lgeo$index #get the sequential position of the deficient

			plot(g$index, g$sum, type="l", main=paste(t,i,pp,sep=""))
			color.scale.lines(g$index, g$sum, col=factor(g$chr))
			points(g$index[highposgeo], g$sum[highposgeo], col = "red", cex=0.7) #mark significant values
			points(g$index[lowposgeo], g$sum[lowposgeo], col = "blue", cex=0.7) #mark significant values

	
			#create BED formatted output for hgeo, hpco. verify first whether there exist any snps that are
			#both significantly enriched and better than random (or significantly suppressed and worse than random.)
			#if not, do not print an output table.
			if ( length(hgeo$pos) != 0 ) {
			hgeobed <- data.frame(chr=paste("Chr", hgeo$chr, sep=""), chromStart=hgeo$pos, chromEnd=hgeo$pos+1)	
			write.table(hgeobed, file=paste(t,"hGen",pp,".bed", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
			}
		
			#create BED formatted output for lgeo, lpco
			if ( length(lgeo$pos) != 0 ) {
			lgeobed <- data.frame(chr=paste("Chr", lgeo$chr, sep=""), chromStart=lgeo$pos, chromEnd=lgeo$pos+1)	
			write.table(lgeobed, file=paste(t,"lGen",pp,".bed", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
			}
		
		dev.off()
	}
}
**********ENDR**********

	Compare results from blip and ceres for populus.  They are highly correlated, thus combining resources is OK.
ceres=$(cut -d' ' -f5 +PopSUMEnrichAcrossBlocks.RgenTgen.txt); #do this, then cd to blip results
blip=$(cut -d' ' -f5 +PopSUMEnrichAcrossBlocks.RgenTgen.txt);
paste -d' ' <(echo "$blip") <(echo "$ceres") > blipcerescompare.txt;

	Compare results from blip 2-50all and 2-50even only:
a=$(cut -d' ' -f4 +PopSUMEnrichAcrossBlocks.RgenTgen.txt);
b=$(cut -d' ' -f4 ../Popb2-50all/+PopSUMEnrichAcrossBlocks.RgenTgen.txt);
echo "$a" | wc -l; # number of positions included in 2-50even = 32113
echo "$b" | wc -l; # number of positions included in 2-50all = 32060
c=$(comm -3 <(echo "$a") <(echo "$b")); # get positions excluded from 2-50all but present in 2-50even
grep -v -f <(echo "$c") +PopSUMEnrichAcrossBlocks.RgenTgen.txt > b2-50evenComparer.txt# remove lines from even not found in all

all=$(cut -d' ' -f5 ../Popb2-50all/+PopSUMEnrichAcrossBlocks.RgenTgen.txt);
even=$(cut -d' ' -f5 b2-50evenComparer.txt);
paste -d' ' <(echo "$even") <(echo "$all") > evenallcompare.txt;




	Determine the identity of each important SNP

	Download Populus trichocarpa annotation v2.2
wget http://www.plantgdb.org/download/Download/xGDB/PtGDB/Ptrichocarpa_156_gene.gff3.bz2
bzip -d Ptrichocarpa_156_gene.gff3.bz2

	Convert to BED format using bedops, extract only genes from GFF file
cat Ptrichocarpa_156_gene.gff3 | gff2bed | grep -w gene > Pt22genes.bed

	Create background data file, a bed file containing only those genes for which there are SNPs.
	Create files listing genes enriched by M+ analyses for downstream GO over-representation analysis.
	Uses bedops. For whatever reason, this script has to be run in two parts or it fails.
**********BASH**********
***Part 1***
i="Pop"; #taxon identifier
ag="Pt22genes.bed"; #bed file containing all annotated genes for given genome release

#Create a bed file for each SNP locus where we have a measurement (only needs to be done for geo or pco, since SNPs are the same):
chr=$(cut -d' ' -f2 "+"$i"SUMEnrichAcrossBlocks.RgenTgen.txt" | tail -n +2 | sed 's/^/scaffold_/g');
pos=$(cut -d' ' -f3 "+"$i"SUMEnrichAcrossBlocks.RgenTgen.txt" | tail -n +2);
pos2=$(awk '{$1 = $1 + 1; print}' <<<"$pos"); 
paste -d$'\t' <(echo "$chr") <(echo "$pos") <(echo "$pos2") > $i"AllSNPs.bed";

#Sort to make sure it will work for bedops --element-of
/usr/local/bin/sort-bed $i"AllSNPs.bed" > tmp.txt;
mv tmp.txt $i"AllSNPs.bed";

#Calculate the intersection of PopAllSNPs.bed and PT22genes.bed, the list of genes for which there are SNPs in the study (the background):
/usr/local/bin/bedops --element-of 1 "$ag" $i"AllSNPs.bed" > $i"AllGenesWithSNPs.bed";
  
#Make background bed files into a simple gene list.
cut -d$'\t' -f4 $i"AllGenesWithSNPs.bed" | sort -u > $i"BackgroundGenesList.txt";
#clean up
rm $i"AllSNPs.bed";
rm $i"AllGenesWithSNPs.bed";

***Part 2***
#Compute lists of genes enriched by M+
for k in "Gen";
  do echo $k;
    #Verify that bed files are sorted correctly using bedops:sort-bed, e.g.
    for j in "0.001" "0.05";
      do echo "  $j";
        for h in "h" "l";
          do echo "    $h";
            /usr/local/bin/sort-bed $i$h$k$j.bed | sed 's/^Chr/scaffold_/g'> tmp.txt;
            mv tmp.txt $i$h$k$j.bed;
    
            #Determine nearest gene to highly enriched snps using closest-features
            /usr/local/bin/closest-features --closest --dist $i$h$k$j.bed "$ag" > $i$h"$k"EnrichedGenes$j.bed;
    
            #Extract hits that are within genes, not just close
            grep "|0$" $i$h"$k"EnrichedGenes$j.bed > $i$h"$k"EnrichedWithinGenes$j.bed;
            rm $i$h"$k"EnrichedGenes$j.bed;
    
            #Get gene identifiers for downstream GO analysis
            cut -d$'\t' -f6 $i$h"$k"EnrichedWithinGenes$j.bed | sort -u > $i$h"$k"EnrichedGenesList$j.txt;
            rm $i$h"$k"EnrichedWithinGenes$j.bed;
          done;
      done;
  done;
**********ENDBASH**********

	Use Python GOATOOLS to analyze over-representation:
easy_install goatools;
easy_install fisher;
easy_install statsmodels;
wget http://geneontology.org/ontology/go-basic.obo;
wget http://www.geneontology.org/ontology/subsets/goslim_generic.obo;

	Download a raw gene association file (GAF) from http://amigo.geneontology.org/amigo/search/annotation
	Choose species in 'Species' menu, remove 'not', 'contributes_to', and "colocalizes_with'
	in 'Annotation Qualifier' menu.  Call it PopRawGAF.txt. This GAF will include GO terms for 
	Biological Process, Molecular Function, and Cellular Component ontologies.
	Process the file into an "association" input file, called PopGAF.txt, for GOATOOLS script find_enrichment.py.
	
	In principle this does not need to be repeated for Final Analysis 3.  The recreated PopGAF.txt file
	contains the same relationships between gene name and GO term as the old PopGAF.txt file (Final Analysis 2), although
	it is not in the same sort order, which is curious.  The old PopGAF.txt file from Final Analysis 2 is
	sorted alphabetically, although I don't recall doing that and it is not in the code. SEE ARABIDOPSIS FOR A
	PARALLEL VERSION.
**********BASH**********
t=Pop;
cut -d$'\t' -f3 "$t"RawGAF.txt | sort -u > GAFgenestmp.txt;
raw=$(cut -d$'\t' -f3,5 "$t"RawGAF.txt);
>tmp.txt;
lf=$(wc -l GAFgenestmp.txt | awk '{print $1}'); #get total line number for progress bar
z=0; #counter
while read -r line;
do
    z=$(($z+1));
    echo -en "\r"$(($lf-$z)); #display countdown progress 
    f="";
    f=$(grep ^"$line"$'\t' <<<"$raw" | cut -d$'\t' -f2 | tr '\n' ';');
    if [[ "$f" != "" ]]; then
      echo "$line"$'\t'"$f" >> tmp.txt;
    fi;
done < GAFgenestmp.txt;
sed 's/;$//g' tmp.txt > "$t"GAF.txt;
rm GAFgenestmp.txt;
rm tmp.txt;
**********ENDBASH**********
	Takes 1.5 hours on Mac.

	Run GOATOOLS find_enrichment.py to get the enrichment status, over (e), or under (p),
	representation (column 3 in output).
	The --no_propagate_counts selects only the least inclusive GO term, i.e no parent terms.
	This eliminates multiple significant results along a parent-child path, but has the
	undesirable consequence of making significant under-representation (p) a questionable
	result.

	For Populus only, need to add a 'g' to the end of each line in source and population input
	files.  These are the '...EnrichedGenesList..." files.
for p in 0.05 0.001;
  do echo $p;
  for e in Gen;
    do echo "  $e";
      for h in "h" "l";
        do echo "    $h";
          sed 's/$/g/g' Pop"$h""$e"EnrichedGenesList"$p".txt > tmp.txt
          mv tmp.txt Pop"$h""$e"EnrichedGenesList"$p".txt;
        done;
    done;
  done;
sed 's/$/g/g' PopBackgroundGenesList.txt > tmp.txt
mv tmp.txt PopBackgroundGenesList.txt;






	Move the relevant items (inc. *GAF.txt, go-basic.obo) to a nested folder, "statistical tests", then run
	the GOATOOLS script find_enrichment.py:
**********BASH**********
mv *EnrichedGenesList* "statistical tests";
mv PopBackgroundGenesList.txt "statistical tests";
cd "statistical tests";

t=Pop;
comm="--alpha 0.05 --pval 0.05 --obo go-basic.obo --no_propagate_counts --method holm";

for p in 0.05 0.001;
  do echo $p;
  for e in Gen;
    do echo "  $e";
      for h in "h" "l";
        do echo "    $h";
          find_enrichment.py $(echo "$comm" "$t""$h""$e"EnrichedGenesList"$p".txt "$t"BackgroundGenesList.txt "$t"GAF.txt) > +"$t""$h""$e""$p"GOAout.txt;
        done;
    done;
  done;
**********ENDBASH**********

	Examine first part of output files to determine if any GO terms are significantly over-
	or under-represented.
head -20 +*GOA*;

	FINAL RESULT:
There are no Holm-Bonferroni corrected, significantly over- or under-represented GO terms for Populus among genes that
are well-collected or poorly-collected by genetic information.



	Make a pie chart showing GOslim category representation for genes enriched by M+ at the
	0.05 level. Do this in a folder called "pie charts".

	Extract the three GO categories (biological_process, cellular_component, molecular_function)
	from goslim_plant.obo.  Use PERL modules go-perl > go-filter-subset.pl.  THIS DOES NOT NEED TO BE REDONE.
	USE EXTRACTED SUBSETS FROM FINAL ANALYSIS 2.

for i in biological_process cellular_component molecular_function;
  do
    /Users/shrub/perl5/bin/go-filter-subset.pl -namespace "$i" goslim_plant.obo > goslim_plant_"$i".obo;
  done;

	For whatever reason (actually, it is because go-filter-subset.pl includes 'part-of's, not
	just 'is_a's as members of a category), two "cellular_component"s remain in goslim_plant_biological_process.obo.
	Ten "biological_process"s remain in goslim_plant_molecular_function.obo.
	Remove them manually.
	



	Make a gene association file (GAF) for M+ well-collected and poorly-collected genes.
	In files: Pop[h,l][Gen]EnrichedGenesList0.05.txt

**********BASHPARALLEL**********
myp() {
    line=$1;
    o=$(mktemp /tmp/tmp.XXXXXXXX);
    f="";
    f=$(grep ^"$line"$'\t' in.txt | cut -d$'\t' -f2 | tr '\n' ';');
    if [[ "$f" != "" ]]; then
      echo "$line"$'\t'"$f" > $o;
    fi;
}
export -f myp;

t=Pop;
rm /tmp/tmp.*; #clear tmp directory of tmp. files generated by this script
cut -d$'\t' -f3,5 "$t"RawGAF.txt > in.txt; #make a simplified GAF file of all genes to query
  for e in Gen;
    do echo "  $e";
      for h in "h" "l";
        do echo "    $h";
          cat $t$h$e"EnrichedGenesList0.05.txt" | parallel --env myp myp;

          cat /tmp/tmp.* | sort | sed 's/;$//g' > ./$t$h$e"0.05GAF.txt";
          rm /tmp/tmp.*;
        done;
    done;
rm in.txt;
**********ENDBASH**********

	Count the number of occurrences of each GO slim term in the gene association files
	derived from the EnrichedGenesList(s). These are the files Pop[h,l][Gen]0.05GAF.txt.
	Do this for each major GO category (bp, cc, mf). Output files are like Pop[h,l][Gen]piechart_[bp,cc,mf].txt.
	Goal here is to calculate the frequency of the function in the well/poorly collected genes and
	compare that to the frequency of the same functions when all annotated genes are considered.
**********BASH**********
  
t="Pop";
b=$(wc -l "$t"GAF.txt | awk '{print $1}'); #the total number of annotated genes
for j in biological_process cellular_component molecular_function;
  do echo $j;
    map_to_slim.py --slim_out=direct --association_file=$t"GAF.txt" go-basic.obo goslim_plant_"$j".obo > "$t"slim_"$j".txt; #get the goslim_plant terms associated with all annotated genes using goatools map_to_slim.py
    a=$(awk -F$'\t' '{print $2}' "$t"slim_"$j".txt | tail -n +5 | sed '/^$/d' | tr ";" "\n" | wc -l  | awk '{print $1}'); #number of function hits, i.e. total number of GO terms found in slim for all collected genes.

    for e in Gen;
      do echo "  $e";
        for h in "h" "l";
          do echo "    $h";
            map_to_slim.py --slim_out=direct --association_file=$t$h$e"0.05GAF.txt" go-basic.obo goslim_plant_"$j".obo > "$t$h$e"slim_"$j".txt; #get the goslim_plant terms associated with the genes using goatools map_to_slim.py

            go=$(awk '{print $2}' "$t$h$e"slim_"$j".txt | tail -n +5 | tr ";" "\n" | sort -u); #list of unique GO terms for well/poorly collected genes
            o=$(awk '{print $2}' "$t$h$e"slim_"$j".txt | tail -n +5 | sed '/^$/d' | tr ";" "\n" | wc -l  | awk '{print $1}'); #number of function hits, i.e. total number of GO terms found in slim for well/poorly collected genes.
            n=$(wc -l "$t$h$e"slim_"$j".txt | awk '{print $1}'); #number of lines = number of genes for well/poorly collected genes
            n=$(( $n - 4 )); #subtract off header lines
            echo "GOterm ObsTotNumGenes ObsTotNumFunctionHits ObsGOcount obsfreq TotNumAnnotatedGenes ExpTotNumFuncHits ExpGOcount expfreq ratio absratio plratio diff absdiff" > Piechart"$t$h$e"_"$j".txt; #for each unique GO term, count how many times it occurs in slim.txt, divide by number of genes to get the proportion of genes representing that function in the M+ enriched/purified gene set.

            for i in $go;
              do f=$(grep $i "$t$h$e"slim_"$j".txt | wc -l | awk '{print $1}'); #get the number of lines the GO term is found in in the slim for well/poorly collected genes
                x=$(grep $i "$t"slim_"$j".txt | wc -l | awk '{print $1}'); #get the number of lines the GO term is found in in the slim for all genes
                g=$(echo "scale=4;$f/$o" | bc); #freq of function among all enriched functions
                y=$(echo "scale=4;$x/$a" | bc); #freq of function among all functions
                z=$(echo "$g - $y" | bc); #diff btw freq of function among all functions, and among enriched functions (negative means less common in enriched subset)
                zz=${z#-}; #absolute value of z
                rrx=$(echo "scale=4;$g/$y" | bc); #for log2 plotting
                #use below 4 lines for non-log2 axes
                if (( $(bc <<< "$g >= $y") )); 
                  then rr=$(echo "scale=4;$g/$y" | bc); #over represented results in positive fold enrichment
                  else rr=$(echo "scale=4;-$y/$g" | bc); #under represented results in negative fold enrichment
                fi; #ratio of frequency in enriched regions / frequency in all regions
                rrr=${rr#-}; #absolute value of rr
                 echo "$i $n $o $f $g $b $a $x $y $rr $rrr $rrx $z $zz" >> Piechart"$t$h$e"_"$j".txt;
              done;
              
            rm "$t$h$e"slim_"$j".txt;
          done;
      done;
  done;
**********ENDBASH**********

	Get verbal description of slimmed GO terms for eventual plotting. Include the cellular_component
	even though it is just the location the protein is found. Filter out singletons, where there is a unique
	observation of a GO term in the set of GO terms from enriched regions:

**********BASH**********
t="Pop";
for e in Gen;
  do echo "  $e";
    for h in "h" "l";
      do echo "    $h";
        > o.txt; #freq output
        > p.txt; #ratio output
        for j in biological_process molecular_function cellular_component;
          do echo $j;
            #filter out singletons, GO terms with only 1 observation in the enriched regions
            sin=$(grep -v "GO:"[0-9]*" "[0-9]*" "[0-9]*" "1" " Piechart$t$h$e"_"$j.txt | tail -n +2);

            #extract columns holding freq and ratio metrics for reformatting
            f=$(echo "$sin" | cut -d' ' -f1); #get GO slim terms
            p=$(echo "$sin" | cut -d' ' -f13); #get diff in freq of function btw enriched functions and all functions
            pp=$(echo "$sin" | cut -d' ' -f14); #get abs of diff in freq of function btw enriched functions and all functions
            q=$(echo "$sin" | cut -d' ' -f10); #get ratio between enriched/all
            qq=$(echo "$sin" | cut -d' ' -f11); #get abs of ratio between enriched/all
            pl=$(echo "$sin" | cut -d' ' -f12); #get ratio between enriched/all to be plotted on log2 axis
            c=$(for i in $f;
              do sed -n -e '/id: '$i'/,$p' goslim_plant_"$j".obo | head -2 | tail -1 | cut -d' ' -f2-;
              done;);
            d=$(for i in $f; do echo "$j"; done;); #create a header column with bp, cc, mf category
            paste -d$'\t' <(echo "$d") <(echo "$f") <(echo "$c") <(echo "$p")  <(echo "$pp") >> o.txt; #write to temporary freq output file
            paste -d$'\t' <(echo "$d") <(echo "$f") <(echo "$c") <(echo "$q")  <(echo "$qq") <(echo "$pl")>> p.txt; #write to temporary ratio output file
         done;
          echo Gocat$'\t'Goterm$'\t'name$'\t'diff$'\t'absdiff > PlotRpieFREQ$t$h$e.txt; #make header for freq output file
          echo Gocat$'\t'Goterm$'\t'name$'\t'ratio$'\t'absratio$'\t'plratio > PlotRpieRATIO$t$h$e.txt; #make header for ratio output file
          sort -t$'\t' -nr -k5,5 o.txt >> PlotRpieFREQ$t$h$e.txt;
          sort -t$'\t' -nr -k5,5 p.txt >> PlotRpieRATIO$t$h$e.txt;
         rm o.txt p.txt;
      done;
  done;

**********ENDBASH**********
	
	Use R to construct charts, change setwd as necessary:
**********RSCRIPT**********
#install.packages("ggplot2")
library(ggplot2)
options(error = recover)
rm(list=ls()) 
setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Populus/Experiment2/+revision\ with\ hclust/+results/hclust/+PopRgenRuns/rarefaction/+analysis/+PopRgenGOWork/pie\ charts")

t<-"Pop"
rvals<-c("RATIO", "FREQ")
cvals<-c("h", "l")
evals<-c("Gen")
for (rr in rvals)
{
  for (cc in cvals)
  {
    for (ee in evals)
    {
      pdf(file=paste("GOchart",rr,t,cc,ee,".pdf", sep=""))
      #par(mfrow=c(1,2))
  
      g <- read.table(paste("PlotRpie",rr,t,cc,ee,".txt", sep=""), header=TRUE, sep="\t")
      h <- g[g$Goterm!="GO:0008150",] #Remove overarching GO categories biological_process, molecular_function, cellular_component
      h <- h[h$Goterm!="GO:0003674",]
      h <- h[h$Goterm!="GO:0005575",]
  
      par(mar = c(4,20,4,2) + 0.1)
      if(rr == "FREQ") {
        barplot(h$diff, main=paste("GOchart",rr,t,cc,ee,".pdf", sep=""), names.arg=paste(h$Goterm, h$name, sep=", "), las=2, horiz=TRUE, cex.names=0.6, xlim=c(min(h$diff),max(h$diff)))
        } 
      else {
        #barplot(h$plratio, log="x", main=paste("GOchart",rr,t,cc,ee,".pdf", sep=""), names.arg=paste(h$Goterm, h$name, sep=", "), las=2, horiz=TRUE, cex.names=0.6, xlim=c(min(h$plratio),max(h$plratio)))
        hh=h #transfer h table to a new variable, which will be modified
        hh$id=paste(hh$Goterm,hh$name,sep=" ") #assemble label term
        hh$id <- factor(hh$id, levels=hh$id) #fix order by making id a factor with levels

        ggout<-ggplot(hh,aes(id,plratio,width=0.8)) + 
        geom_bar(stat="identity",color="black",fill="dark grey",size=0.25) + 
        scale_y_continuous(trans='log2', breaks=c(0.125,0.25,0.5,1,2,4,8)) + 
        coord_flip() + 
        theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), text=element_text(size=12, family="ArialMT"))
        
        print(ggout)
        }      
      dev.off()
    }
  }
}
**********ENDR**********




@@@@@@@@@@@@@@@@@@@@ END of tested processing for PopulusRgen @@@@@@@@@@@@@@@@@@@@@@@@@@



























***SORGHUM***
	Calculate mean M+ enrichment across core sizes, for each locus, for each blocklength. This uses the
	SUMfiles as input and produces the +STATS files.  Do this on the cluster.



**********BASHPARALLEL-FASTER-RUNINSCREEN**********
mysd() {
       #calculate sum, mean, sd, n, mean,using Welford's algorithm and AWK
       #here, sum is the cumulative enrichment across all core sizes, mean is the average
       #supply space delimited string of values as argument like: mysd $foo, where foo="1 2 3 4.5"
       awk '{
         sum = 0;                   # Initialize running sum (for mean calculation)
         M = 0;
         S = 0;
         for (k=1; k <= NF; k++) { 
              sum += $k;                # Update running sum
              x = $k;
              oldM = M;
              M = M + ((x - M)/k);
              S = S + (x - M)*(x - oldM);
         }
         var = S/(NF - 1);
         print sum "\t" sum/NF "\t" sqrt(var) "\t" NF;
       }' <<< "$1"; # the first argument when calling mysd 
}
export -f mysd;

main() {
        l=$1;
        ofile="$p""+STATS.R"$e"Tgen.loc"$l".b"$b".txt"; #temporary outfile for each parallel process
        sumfileIN="$p""SUM.b"$b".R"$e"Tgen.txt";

        boo=$(grep -F $l "$sumfileIN" | awk '{ print $13 - $10 }' | tr "\n" " "); #extract all core sizes for this locus, calculate target enrichment
        st=$(mysd "$boo"); #calculate sum, mean, sd, and n using function mysd()
        echo "$b"$'\t'"$l"$'\t'"$st" > "$ofile"; #write stats to output file
}
export -f main;

#parameters
t="Sorghum";
elist="gen";
bmin=2; #min blocklength
bmax=50; #max blocklength
bstep=2; #step rate
nnode=9; #number of cluster nodes
blist=$(seq $bmin $bstep $bmax); #blocklengthrange
p=$(pwd)"/"; 
export p;

#parallelize on locus, within loop on blocklength
for e in $elist;
  do echo $e;
    export e;
    for b in $blist;
      do echo "b=$b";
        export b;
        sumfileIN="$p""SUM.b"$b".R"$e"Tgen.txt";
        uloc=$(cut -d$'\t' -f2 "$sumfileIN" | tail -n +2 | uniq); #find unique loci from all loci in column 2
        numloc=$(echo "$uloc" | wc -l);
        
        #parallel step
        #sends $gnuN records to each node via parallel step #1, this is piped to parallel step #2, which starts 96 jobs per node from the instruction set it has received
        echo "processing $sumfileIN...";
        gnuN=$(echo "$numloc / $nnode" | bc); #number of lines of input file to send to gnu parallel per node
        echo "$uloc" | parallel --sshloginfile ~/machines --jobs 1 --env main --env mysd --env e --env b --env p --pipe -N"$gnuN" /home/reevesp/bin/parallel -j96 --env main --env mysd --env e --env b --env p main;
        #run on head node only
        #echo "$uloc" | parallel --jobs 1 --env main --env mysd --env e --env b --env p --pipe -N"$gnuN" /home/reevesp/bin/parallel -j96 --env main --env mysd --env e --env b --env p main;
        #run on ceres, one node
        #echo "$uloc" | parallel --jobs 1 --env main --env mysd --env e --env b --env p --pipe -N"$gnuN" parallel -j96 --env main --env mysd --env e --env b --env p main;
        #run on Mac
        #echo "$uloc" | parallel --env main --env mysd --env e --env b --env p main;

        #consolidate results for each locus into one file
        ofile="+STATS.R"$e"Tgen.b"$b".txt"; #temporary outfile
        > "$p""$ofile";
        #find "$p" -name "+STATS.R"$e"Tgen.loc*.b"$b".txt" | while read line; do cat "$line" >> "$p""$ofile"; done; #fast cat, unsorted by $uloc
        for l in $uloc;
          do echo "concatenating loc file $l to +STATS.R"$e"Tgen.b"$b".txt ...";
            cat "$p""+STATS.R"$e"Tgen.loc"$l".b"$b".txt" >> "$p""$ofile";
          done;
          #test whether all were concatenated
          s1=$(ls "$p" | grep -c "loc"); #count the number of files with "loc" in their name using grep. can't wc -l because too many arguments
          s2=$(wc -l "$p""$ofile" |  awk '{print $1}');
         if [ $s1 == $s2 ]; 
           then find "$p" -maxdepth 1 -name "*loc*" -delete; #remove files when there are too many for rm
           else echo "s1="$s1", s2="$s2". Consolidated file "$ofile" no good. Quitting..." > err.txt;
           exit 1;
         fi;
      done; #$b
      
    #concatenate files for this Rgeo/Rpco
    statsfile="$p""+STATS."$t".R"$e"Tgen.txt";
    echo blocklength$'\t'locus$'\t'sum$'\t'mean$'\t'sd$'\t'n > "$statsfile";
    >"$statsfile"TMP;
    #quick cat and sort (some out files from ceres are not sorted correctly by locus)
    cat "$p""+STATS.R"$e"Tgen.b"*".txt" >> "$statsfile"TMP; #cat together the whole mess
    sortkey1=$(cut -d$'\t' -f2 "$statsfile"TMP | cut -d'.' -f3); #get a column containing only the locus index number, field 3 of column 2
    paste -d$'\t' "$statsfile"TMP <(echo "$sortkey1") | sort -t$'\t' -n -k1,1 -n -k7,7 | cut -d$'\t' -f1-6 > "$statsfile"TMP2; #add sort key to last column, sort by blocklength then locus index, remove last column
    cat "$statsfile"TMP2 >> "$statsfile";
    rm "$statsfile"TMP;
    rm "$statsfile"TMP2;

    #test whether all individual +STATS files have been added to the final concatenated file using number of lines
    s1=$(wc -l +STATS.R"$e"Tgen.b*.txt | grep total | awk '{print $1}');
    s2=$(wc -l "$statsfile" | awk '{print $1}');
    s2=$(( $s2 - 1 )); #decount header
    if [ $s1 == $s2 ]; 
      then rm +STATS.R"$e"Tgen.b*.txt; #remove files if number of lines in concat file is the same as the sum of all input files
      else echo "the number of lines ain't the same. something is 'crewed.";
    fi;

  done; #$e
**********BASHEND**********
	Takes ~45 minutes on blip. 


	add the proportion of NG,SS,NS sites from SorLocusWeights.txt to the +STATS summary file,
	generating the +v2STATS file.  Requires the SorLocusWeights.txt file as input.  

**********BASH**********
#use md5 on mac md5sum on linux
#modify weightsfile to include only the blocklengths analyzed
weightsfile="SorLocusWeights.txt";
bmin=2; #min blocklength
bmax=50; #max blocklength
bstep=2; #step rate
blist=$(seq $bmin $bstep $bmax); #blocklengthrange

wtfhead=$(head -1 "$weightsfile");
wtfbod=$(awk -v blist="$blist" 'BEGIN{split(blist,t); for (i in t) vals[t[i]]} ($2 in vals)' "$weightsfile");
#below is slower, but more human-readable alternative
#wtfbod=$(for i in $blist;
#  do awk -F"\t" -v i=$i '$2 == i {print}' "$weightsfile"; #accumulate each blocklength needed
#  done;)
wtf=$(echo "$wtfhead"; echo "$wtfbod");

echo "$wtf" > wtf.txt; 

#wtf=$(awk -F"\t" '$2 != "1" {print}' "$weightsfile"); #exclude blocklength=1 from the weightsfile for use with blocklength 2-200
elist="gen";
for e in $elist;
  do statsfile="+STATS.Sorghum.R"$e"Tgen.txt";
    s=$(cut -d$'\t' -f2 $statsfile | tail -n +2 | md5sum); #get md5 of locus ID column from the +STATS file
    w1=$(echo "$wtf" | cut -d$'\t' -f4 | sed 's/_/./g'); #get the locus ID less blocklength index column, replace _ with . for md5 check
    w2=$(echo "$wtf" | cut -d$'\t' -f3); #get the blocklength index column
    w=$(paste -d'.' <(echo "$w1") <(echo "$w2") | tail -n +2 | md5sum); #get md5 of locus ID column from the Weights file
    if [ "$s" == "$w" ]; #verify that files are in same order using a checksum of the unique locus ID column
      then boo=$(echo "$wtf" | sed 's/blocklength/blocklength2/g' | sed 's/locus/locus2/g'); #modify the header in the Weights file so that there are no redundant column names after combining
        paste -d$'\t' "$statsfile" <(echo "$boo") > "+v2STATS.Sorghum.R"$e"Tgen.txt"; #paste the statsfile and the weights file together
      else echo "md5s do not match, $e, aborting.";
        kill -INT $$; #terminate the script, return to the shell
    fi;
  done;
**********BASHEND**********




#Plot summed enrichment by genomic position across blocklengths as a heatmap. Use Rscript
#SorGenomicGeography.r. Uses the +v2STATS files as input.



	Using the file +*SNPtoBlockMap.txt, calculate the sum of M+ enrichment value across blocklengths,
	for each SNP, i.e. for each position in the genome.  These data will be extracted from the 
	+v2STATS* files in several steps.  The first makes a file containing a table with the max 
	enrichment value at each site for each block, *EnrichAcrossBlocks*. Step 1:

**********BASHPARALLEL**********
myp() {
  e=$1;
  b=$2;
  v2name="$p""+v2STATS."$v2p".R"$e"Tgen.txt";
        
        echo "e=$e  b=$b";
        bout="$p""rr.R"$e"Tgen.b""$b"".tmp"; #temporary outfile

        #cut out the locus id and sum columns for current blocklenth
        gg=$(grep ^$b$'\t' "$v2name" | cut -d$'\t' -f1-3);
        
        #cut column with locus names for current blocklength from SNPtoBlockMap
        if [ $b = 1 ]; then
          cc=$(cut -d$'\t' -f4 "$tname" | tail -n +2 | sed 's/_/./g' | sed 's/$/./g');
        else col=$(( $b + 5 )); #calculate the column to cut out for each blocklength
          cc=$(cut -d$'\t' -f$col "$tname" | tail -n +2 | sed 's/_/./g' | sed 's/$/./g');
        fi;

        #search for locus names in $gg, retrieve enrichment value
        rr="b$b";

#put another parallel here if necessary
        for l in $cc;
          do if [[ $l = "--." ]]; then
                r="--"; #test for empty enrichment value
              else
                r=$(grep -m 1 ^"$b"$'\t'"$l" <<<"$gg" | awk '{print $3}');
              fi;
            rr+=$'\n'$r; #add current locus enrichment value to list
          done;

         #write out the result
         echo "$rr" > "$bout";
         #truncate -s -1 "$bout"; #remove the hanging newline, use 'gtruncate' only on osx
}
export -f myp;

#parameters
v2p="Sorghum"; export v2p;
t="Sor";
minbl=2; #1,2
maxbl=50; #50,200
blstep=2; #step
blist=$(seq $minbl $blstep $maxbl); #for even numbered runs 2-50
elist="gen";

p=$(pwd)"/"; export p;
tname="$p""+""$t""SNPtoBlockMap.txt"; export tname;

#begin parallel
parallel --sshloginfile ~/machines --jobs 24 --env myp --env p --env v2p --env tname myp ::: $elist ::: $blist;
#parallel --jobs 24 --env myp --env p --env v2p --env tname myp ::: $elist ::: $blist; #parallel command for head node only
#parallel --jobs 25 --env myp --env p --env v2p --env tname myp ::: $elist ::: $blist; #parallel command for ceres
#parallel --env myp --env p --env v2p --env tname myp ::: $elist ::: $blist; #parallel command for Mac

  
#reconstitute
echo "reconstituting..."
for e in $elist;
  do echo "e=$e";
    comm=""; #create a list of tmp files to paste together
    for b in $blist;
      do comm+="rr.R"$e"Tgen.b$b.tmp ";
    done;
    pos=$(cut -d$'\t' -f4 "$tname" | sed 's/'^.*_'//g' | sed 's/locus/pos/g');
    h=$(cut -d$'\t' -f1-6 "$tname" | tr $'\t' ' ');
    h2=$(paste -d' ' <(echo "$h") <(echo "$pos"));
    paste -d' ' <(echo "$h2") $comm > "+"$t"EnrichAcrossBlocks.R"$e"Tgen.txt";
    rm $comm; #clean up
  done;#e
  
**********BASHEND**********
	Takes ~19 hours.




	Process each line (SNP position) of *EnrichAcrossBlocks*.txt to get sum across all blocklengths
	of cumulative enrichment across all core sizes of max enrichment value from 10 reps of M+.
	Produces a file called SUM*EnrichAcrossBlocks... Step 2:

**********BASH**********
mysd() {
       #calculate sum, mean, sd, n, mean,using Welford's algorithm and AWK
       #here, sum is the cumulative enrichment across all blocklengths, mean is the average
       #supply space delimited string of values as argument like: mysd $foo, where foo="1 2 3 4.5"
       #returns space delimited string like "sum mean sd n"
       awk '{
         sum = 0;                   # Initialize running sum (for mean calculation)
         M = 0;
         S = 0;
         for (k=1; k <= NF; k++) { 
              sum += $k;                # Update running sum
              x = $k;
              oldM = M;
              M = M + ((x - M)/k);
              S = S + (x - M)*(x - oldM);
         }
         var = S/(NF - 1);
         print sum " " sum/NF " " sqrt(var) " " NF;
       }' <<< "$1"; # the first argument when calling mysd 
}
export -f mysd;

myp() {
    l=$1

    p=$(echo "$l" | cut -d' ' -f8-); #cut out the portion of the line with max enrichment values

    #skip loci that have missing blocks
    if [[ $p != *"--"* ]]; then
      m=$(mysd "$p"); #calculate sum, mean, sd, n
      chr=$(echo "$l" | cut -d' ' -f5);
      pos=$(echo "$l" | cut -d' ' -f7);
      i=$(echo "$l" | cut -d' ' -f1); #get genome index
      #loc=$chr"."$pos"."$(echo "$l" | cut -d' ' -f6);
      loc=$chr"."$pos"."$i;
      echo "snp=$i";

      echo $i" "$chr" "$pos" "$loc" "$m > "$path"$i".tmp"; #add new values to growing string, write to tmp file
    fi;
}
export -f myp;

#parameters
elist="gen";
t="Sor";
path=$(pwd)"/"; export path;

#parallel
for e in $elist;
do echo $e;
  #nnode=8;
  #numloc=$(wc -l +AtEnrichAcrossBlocks.RgenTgen.txt | awk '{print $1}');
  #numloc=$(( $numloc - 1 ));
  #gnuN=$(echo "$numloc / $nnode" | bc); #number of lines of input file to send to gnu parallel per node

  #parallel command
  #for cluster
  tail -n +2 "+"$t"EnrichAcrossBlocks.R"$e"Tgen.txt" | parallel  --sshloginfile ~/machines --jobs 1 --env myp --env mysd --env path --pipe --round-robin /home/reevesp/bin/parallel -j96  --env myp --env mysd --env path myp;
  #for head node only
  #tail -n +2 "+"$t"EnrichAcrossBlocks.R"$e"Tgen.txt" | parallel --jobs 1 --env myp --env mysd --env path --pipe --round-robin /home/reevesp/bin/parallel -j96  --env myp --env mysd --env path myp;
  #for ceres N=26
  #tail -n +2 "+"$t"EnrichAcrossBlocks.R"$e"Tgen.txt" | parallel  --jobs 26 --env myp --env mysd --env path myp;
  #for mac
  #tail -n +2 "+"$t"EnrichAcrossBlocks.R"$e"Tgen.txt" | parallel  --env myp --env mysd --env path myp;
  
  #assemble results
  echo "concatenating...";
  ofile="+"$t"SUMEnrichAcrossBlocks.R"$e"Tgen.txt";
  o="genomeindex chr pos loc sum mean sd n";
  echo "$o" > "$ofile";
  
  find "$path" -name "*.tmp" -print0 | xargs -0 cat | sort -n -t' ' -k1 >> "$ofile";
  
  #clean up
  echo "clean up...";
  find "$path" -name "*.tmp" -print0 | xargs -0 rm;
done;
**********ENDBASH**********
	Takes ~10 minutes, about 7 hours on ceres N=26

	Use R to make some pretty plots of significant enrichment across the genome.

#For Rgen use this
**********BEGINR**********
#install.packages("fitdistrplus")
#install.packages("plotrix")
options(error = recover)
rm(list=ls()) 
setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/ComparisonFigs/RecalculateM+2.50/SorghumRecalc2/+analysis")
library(fitdistrplus)
library(plotrix)

t<-"Sor"
usezero="no" #this sets whether to include significant sites only when they are also > or < 0
             #to use the zero criterion, say "yes", to just use the ecdf around the mean, say "no"
pvals<-c(0.001,0.05)
devvals<-c("pdf","svg","png")
for (pp in pvals)
{
	for (zg in devvals)
	{
		zgg<-get(zg) #e.g convert string "pdf" to function pdf now named 'zgg'
		if ( zg == "pdf" | zg == "svg" ) { zgg(file=paste(t,"M+GenomeScanRgen",pp, ".", zg, sep="")) }
		else { 
		    par(mar=c(1,1,1,1))
		    zgg(file=paste(t,"M+GenomeScanRgen",pp, ".", zg, sep=""),width=2000,height=2000,res=300) } #call zgg function with relevant parameters
		
		par(mfrow=c(2,1))

			i="gen"
			g <- read.table(paste("+",t,"SUMEnrichAcrossBlocks.R",i,"Tgen.2.50.txt", sep=""), header=TRUE, sep=" ")
			mycdf <- ecdf(g$sum) #calculate the cdf of the empirical distribution
			p <- mycdf(g$sum) #get the probability of each observation (or lower)
			g$pecdf <- p #add new column to g, containing the probability from the empirical cdf
			g$index <- 1:nrow(g) #add a column with the sequential index
			h <- g[which(g$pecdf > 1-pp),] #get snps that are highly enriched
			if ( usezero == "yes" ) {
				hgeo <- h[which(h$sum>0),] #remove snps where sum<0, i.e. less than random
			}
			else {
				hgeo <- h # do not remove snps based on positive or negative
			}
			l <- g[which(g$pecdf < pp),] #get snps that are lowly enriched
			if ( usezero == "yes" ) {
				lgeo <- l[which(l$sum<0),] #remove snps where sum>0, i.e. more than random
			}
			else {
				lgeo <- l # do not remove snps based on positive or negative
			}
			highposgeo<-hgeo$index #get the sequential position of the highly enriched
			lowposgeo<-lgeo$index #get the sequential position of the deficient

			plot(g$index, g$sum, type="l", main=paste(t,i,pp,sep=""))
			color.scale.lines(g$index, g$sum, col=factor(g$chr))
			points(g$index[highposgeo], g$sum[highposgeo], col = "red", cex=0.7) #mark significant values
			points(g$index[lowposgeo], g$sum[lowposgeo], col = "blue", cex=0.7) #mark significant values

	
			#create BED formatted output for hgeo, hpco. verify first whether there exist any snps that are
			#both significantly enriched and better than random (or significantly suppressed and worse than random.)
			#if not, do not print an output table.
			if ( length(hgeo$pos) != 0 ) {
			hgeobed <- data.frame(chr=paste("Chr", hgeo$chr, sep=""), chromStart=hgeo$pos, chromEnd=hgeo$pos+1)	
			write.table(hgeobed, file=paste(t,"hGen",pp,".bed", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
			}
		
			#create BED formatted output for lgeo, lpco
			if ( length(lgeo$pos) != 0 ) {
			lgeobed <- data.frame(chr=paste("Chr", lgeo$chr, sep=""), chromStart=lgeo$pos, chromEnd=lgeo$pos+1)	
			write.table(lgeobed, file=paste(t,"lGen",pp,".bed", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
			}
		
		dev.off()
	}
}
**********ENDR**********


For Rgeo and Rpco use this
**********BEGINR**********
#install.packages("fitdistrplus")
#install.packages("plotrix")
options(error = recover)
rm(list=ls()) 
setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/ComparisonFigs/RecalculateM+2.50/SorghumRecalc2/+analysis")
library(fitdistrplus)
library(plotrix)

t<-"Sor"
usezero="no" #this sets whether to include significant sites only when they are also > or < 0
             #to use the zero criterion, say "yes", to just use the ecdf around the mean, say "no"
pvals<-c(0.001,0.05)
for (pp in pvals)
{
	pdf(file=paste(t,"M+GenomeScan",pp,".pdf", sep=""))
	par(mfrow=c(2,1))

		i="geo"
		g <- read.table(paste("+",t,"SUMEnrichAcrossBlocks.R",i,"Tgen.2.50.txt", sep=""), header=TRUE, sep=" ")
		mycdf <- ecdf(g$sum) #calculate the cdf of the empirical distribution
		p <- mycdf(g$sum) #get the probability of each observation (or lower)
		g$pecdf <- p #add new column to g, containing the probability from the empirical cdf
		g$index <- 1:nrow(g) #add a column with the sequential index
		h <- g[which(g$pecdf > 1-pp),] #get snps that are highly enriched
		if ( usezero == "yes" ) {
			hgeo <- h[which(h$sum>0),] #remove snps where sum<0, i.e. less than random
		}
		else {
			hgeo <- h # do not remove snps based on positive or negative
		}
		l <- g[which(g$pecdf < pp),] #get snps that are lowly enriched
		if ( usezero == "yes" ) {
			lgeo <- l[which(l$sum<0),] #remove snps where sum>0, i.e. more than random
		}
		else {
			lgeo <- l # do not remove snps based on positive or negative
		}
		highposgeo<-hgeo$index #get the sequential position of the highly enriched
		lowposgeo<-lgeo$index #get the sequential position of the deficient

		plot(g$index, g$sum, type="l", main=paste(t,i,pp,sep=""))
		color.scale.lines(g$index, g$sum, col=factor(g$chr))
		points(g$index[highposgeo], g$sum[highposgeo], col = "red", cex=0.7) #mark significant values
		points(g$index[lowposgeo], g$sum[lowposgeo], col = "blue", cex=0.7) #mark significant values

		i="pco"
		g <- read.table(paste("+",t,"SUMEnrichAcrossBlocks.R",i,"Tgen.2.50.txt", sep=""), header=TRUE, sep=" ")
		mycdf <- ecdf(g$sum) #calculate the cdf of the empirical distribution
		p <- mycdf(g$sum) #get the probability of each observation (or lower)
		g$pecdf <- p #add new column to g, containing the probability from the empirical cdf
		g$index <- 1:nrow(g) #add a column with the sequential index
		h <- g[which(g$pecdf > 1-pp),] #get snps that are highly enriched
		if ( usezero == "yes" ) {
		hpco <- h[which(h$sum>0),] #remove snps where sum<0, i.e. less than random
		}
		else {
			hpco <- h
		}
		l <- g[which(g$pecdf < pp),] #get snps that are lowly enriched
		if ( usezero == "yes" ) {
		lpco <- l[which(l$sum<0),] #remove snps where sum>0, i.e. more than random
		}
		else {
			lpco <- l
		}
		highpospco<-hpco$index #get the sequential position of the highly enriched
		lowpospco<-lpco$index #get the sequential position of the deficient

		plot(g$index, g$sum, type="l", main=paste(t,i,pp,sep=""))
		color.scale.lines(g$index, g$sum, col=factor(g$chr))
		points(g$index[highpospco], g$sum[highpospco], col = "red", cex=0.7) #mark significant values
		points(g$index[lowpospco], g$sum[lowpospco], col = "blue", cex=0.7) #mark significant values

		#calculate intersection of well-collected significant snps
		isect <- hgeo[is.element(hgeo$pos, intersect(hgeo$pos,hpco$pos)),]
		write.table(isect, file=paste(t,"HighCommonToGeoPco",pp,".txt", sep=""), row.names=FALSE, quote=FALSE)

		#calculate intersection of poorly-collected significant snps
		isect <- lgeo[is.element(lgeo$pos, intersect(lgeo$pos,lpco$pos)),]
		write.table(isect, file=paste(t,"LowCommonToGeoPco",pp,".txt", sep=""), row.names=FALSE, quote=FALSE)
	
		#create BED formatted output for hgeo, hpco. verify first whether there exist any snps that are
		#both significantly enriched and better than random (or significantly suppressed and worse than random.)
		#if not, do not print an output table.
		if ( length(hgeo$pos) != 0 ) {
		hgeobed <- data.frame(chr=paste("Chr", hgeo$chr, sep=""), chromStart=hgeo$pos, chromEnd=hgeo$pos+1)	
		write.table(hgeobed, file=paste(t,"hGeo",pp,".bed", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
		}
		if ( length(hpco$pos) != 0 ) {
		hpcobed <- data.frame(chr=paste("Chr", hpco$chr, sep=""), chromStart=hpco$pos, chromEnd=hpco$pos+1)	
		write.table(hpcobed, file=paste(t,"hPco",pp,".bed", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
		}
		
		#create BED formatted output for lgeo, lpco
		if ( length(lgeo$pos) != 0 ) {
		lgeobed <- data.frame(chr=paste("Chr", lgeo$chr, sep=""), chromStart=lgeo$pos, chromEnd=lgeo$pos+1)	
		write.table(lgeobed, file=paste(t,"lGeo",pp,".bed", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
		}
		if ( length(lpco$pos) != 0 ) {
		lpcobed <- data.frame(chr=paste("Chr", lpco$chr, sep=""), chromStart=lpco$pos, chromEnd=lpco$pos+1)	
		write.table(lpcobed, file=paste(t,"lPco",pp,".bed", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
		}

	dev.off()
}
**********ENDR**********






	Determine the identity of each important Sorghum SNP

	Download Sorghum bicolor annotation v1.4
wget ftp://ftp.jgi-psf.org/pub/JGI_data/Sorghum_bicolor/v1.0/Sbi/annotation/Sbi1.4/Sbi1.4.gff3.gz
gunzip Sbi1.4.gff3.gz

	Convert to BED format using bedops, extract only genes from GFF file
cat Sbi1.4.gff3 | gff2bed | grep -w gene > Sb14genes.bed

	Create background data file, a bed file containing only those genes for which there are SNPs.
	Create files listing genes enriched by M+ analyses for downstream GO over-representation analysis.
	Uses bedops. For whatever reason, this script has to be run in two parts or it fails.
**********BASH**********
***Part 1***
i="Sor"; #taxon identifier
ag="Sb14genes.bed"; #bed file containing all annotated genes for given genome release

#Create a bed file for each SNP locus where we have a measurement (only needs to be done for geo or pco, since SNPs are the same):
chr=$(cut -d' ' -f2 "+"$i"SUMEnrichAcrossBlocks.RgenTgen.2.50.txt" | tail -n +2 | sed 's/^/chromosome_/g');
pos=$(cut -d' ' -f3 "+"$i"SUMEnrichAcrossBlocks.RgenTgen.2.50.txt" | tail -n +2);
pos2=$(awk '{$1 = $1 + 1; print}' <<<"$pos"); 
paste -d$'\t' <(echo "$chr") <(echo "$pos") <(echo "$pos2") > $i"AllSNPs.bed";

#Sort to make sure it will work for bedops --element-of
/usr/local/bin/sort-bed $i"AllSNPs.bed" > tmp.txt;
mv tmp.txt $i"AllSNPs.bed";

#Calculate the intersection of SorAllSNPs.bed and Sb14genes.bed, the list of genes for which there are SNPs in the study (the background):
/usr/local/bin/bedops --element-of 1 "$ag" $i"AllSNPs.bed" > $i"AllGenesWithSNPs.bed";
  
#Make background bed files into a simple gene list.
cut -d$'\t' -f4 $i"AllGenesWithSNPs.bed" | sort -u > $i"BackgroundGenesList.txt";
#clean up
rm $i"AllSNPs.bed";
rm $i"AllGenesWithSNPs.bed";

***Part 2***
#Compute lists of genes enriched by M+
for k in "Gen";
  do echo $k;
    #Verify that bed files are sorted correctly using bedops:sort-bed, e.g.
    for j in "0.001" "0.05";
      do echo "  $j";
        for h in "h" "l";
          do echo "    $h";
            /usr/local/bin/sort-bed $i$h$k$j.bed | sed 's/^Chr/chromosome_/g' > tmp.txt;
            mv tmp.txt $i$h$k$j.bed;
    
            #Determine nearest gene to highly enriched snps using closest-features
            /usr/local/bin/closest-features --closest --dist $i$h$k$j.bed "$ag" > $i$h"$k"EnrichedGenes$j.bed;
    
            #Extract hits that are within genes, not just close
            grep "|0$" $i$h"$k"EnrichedGenes$j.bed > $i$h"$k"EnrichedWithinGenes$j.bed;
            rm $i$h"$k"EnrichedGenes$j.bed;
    
            #Get gene identifiers for downstream GO analysis
            cut -d$'\t' -f6 $i$h"$k"EnrichedWithinGenes$j.bed | sort -u > $i$h"$k"EnrichedGenesList$j.txt;
            rm $i$h"$k"EnrichedWithinGenes$j.bed;
          done;
      done;
  done;
**********ENDBASH**********

	Use Python GOATOOLS to analyze over-representation:
easy_install goatools;
easy_install fisher;
easy_install statsmodels;
wget http://geneontology.org/ontology/go-basic.obo;
wget http://www.geneontology.org/ontology/subsets/goslim_generic.obo;

	Download a raw gene association file (GAF) from http://amigo.geneontology.org/amigo/search/annotation
	Choose species in 'Species' menu, remove 'not', 'contributes_to', and "colocalizes_with'
	in 'Annotation Qualifier' menu.  Call it SorRawGAF.txt. This GAF will include GO terms for 
	Biological Process, Molecular Function, and Cellular Component ontologies.
	Process the file into an "association" input file, called SorGAF.txt, for GOATOOLS script find_enrichment.py
	
	THIS DOESNT NEED TO BE REPEATED. USE FROM FINAL ANALYSIS 2. SEE ARABIDOPSIS FOR A
	PARALLEL VERSION.
**********BASH**********
t=Sor;
cut -d$'\t' -f3 "$t"RawGAF.txt | sort -u > GAFgenestmp.txt;
raw=$(cut -d$'\t' -f3,5 "$t"RawGAF.txt);
>tmp.txt;
lf=$(wc -l GAFgenestmp.txt | awk '{print $1}'); #get total line number for progress bar
z=0; #counter
while read -r line;
do
    z=$(($z+1));
    echo -en "\r"$(($lf-$z)); #display countdown progress 
    f="";
    f=$(grep ^"$line"$'\t' <<<"$raw" | cut -d$'\t' -f2 | tr '\n' ';');
    if [[ "$f" != "" ]]; then
      echo "$line"$'\t'"$f" >> tmp.txt;
    fi;
done < GAFgenestmp.txt;
sed 's/;$//g' tmp.txt > "$t"GAF.txt;
rm GAFgenestmp.txt;
rm tmp.txt;
**********ENDBASH**********
	Takes a few hours.

	Run GOATOOLS find_enrichment.py to get the enrichment status, over (e), or under (p),
	representation (column 3 in output).
	The --no_propagate_counts selects only the least inclusive GO term, i.e no parent terms.
	This eliminates multiple significant results along a parent-child path, but has the
	undesirable consequence of making significant under-representation (p) a questionable
	result.

	Move the relevant items (inc. *GAF.txt, go-basic.obo) to a nested folder, "statistical tests", then run
	the GOATOOLS script find_enrichment.py:
**********BASH**********
mv *EnrichedGenesList* "statistical tests";
mv SorBackgroundGenesList.txt "statistical tests";
cd "statistical tests";

t=Sor;
comm="--alpha 0.05 --pval 0.05 --obo go-basic.obo --no_propagate_counts --method holm";
for p in 0.05 0.001;
  do echo $p;
  for e in Gen;
    do echo "  $e";
      for h in "h" "l";
        do echo "    $h";
          find_enrichment.py $(echo "$comm" "$t""$h""$e"EnrichedGenesList"$p".txt "$t"BackgroundGenesList.txt "$t"GAF.txt) > +"$t""$h""$e""$p"GOAout.txt;
        done;
    done;
  done;
**********ENDBASH**********

	Examine first part of output files to determine if any GO terms are significantly over-
	or under-represented.
head -20 +*GOA*;

	FINAL RESULT:
	There are no Holm-Bonferroni corrected, significantly over- or under-represented GO terms among Sorghum genes that are well-
	or poorly- collected using genetic data.




	Make a pie chart showing GOslim category representation for genes enriched by M+ at the
	0.05 level.

	Extract the three GO categories (biological_process, cellular_component, molecular_function)
	from goslim_plant.obo.  Use PERL modules go-perl > go-filter-subset.pl. NO NEED TO REPEAT.

for i in biological_process cellular_component molecular_function;
  do
    /Users/shrub/perl5/bin/go-filter-subset.pl -namespace "$i" goslim_plant.obo > goslim_plant_"$i".obo;
  done;

	For whatever reason (actually, it is because go-filter-subset.pl includes 'part-of's, not
	just 'is_a's as members of a category), two "cellular_component"s remain in goslim_plant_biological_process.obo.
	Ten "biological_process"s remain in goslim_plant_molecular_function.obo.
	Remove them manually.
	
	Make a gene association file (GAF) for M+ well-collected and poorly-collected genes.
	Input files: Sor[h,l][Geon]EnrichedGenesList0.05.txt

**********BASHPARALLEL**********
myp() {
    line=$1;
    o=$(mktemp /tmp/tmp.XXXXXXXX);
    f="";
    f=$(grep ^"$line"$'\t' in.txt | cut -d$'\t' -f2 | tr '\n' ';');
    if [[ "$f" != "" ]]; then
      echo "$line"$'\t'"$f" > $o;
    fi;
}
export -f myp;

t=Sor;
rm /tmp/tmp.*; #clear tmp directory of tmp. files generated by this script
cut -d$'\t' -f3,5 "$t"RawGAF.txt > in.txt; #make a simplified GAF file of all genes to query
  for e in Gen;
    do echo "  $e";
      for h in "h" "l";
        do echo "    $h";
          cat $t$h$e"EnrichedGenesList0.05.txt" | parallel --env myp myp;

          cat /tmp/tmp.* | sort | sed 's/;$//g' > ./$t$h$e"0.05GAF.txt";
          rm /tmp/tmp.*;
        done;
    done;
rm in.txt;
**********ENDBASH**********

	Count the number of occurrences of each GO slim term in the gene association files
	derived from the EnrichedGenesList(s). These are the files Sor[h,l][Pco,Geo]0.05GAF.txt.
	Do this for each major GO category (bp, cc, mf). Output files are like PiechartSor[h,l][Pco,Geo]_[bp,cc,mf].txt.
	Goal here is to calculate the frequency of the function in the well/poorly collected genes and
	compare that to the frequency of the same functions when all annotated genes are considered.
	Requires go-basic.obo and SorGAF.txt as input.
**********BASH**********
  
t="Sor";
b=$(wc -l "$t"GAF.txt | awk '{print $1}'); #the total number of annotated genes
for j in biological_process cellular_component molecular_function;
  do echo $j;
    map_to_slim.py --slim_out=direct --association_file=$t"GAF.txt" go-basic.obo goslim_plant_"$j".obo > "$t"slim_"$j".txt; #get the goslim_plant terms associated with all annotated genes using goatools map_to_slim.py
    a=$(awk -F$'\t' '{print $2}' "$t"slim_"$j".txt | tail -n +5 | sed '/^$/d' | tr ";" "\n" | wc -l  | awk '{print $1}'); #number of function hits, i.e. total number of GO terms found in slim for all collected genes.

    for e in Gen;
      do echo "  $e";
        for h in "h" "l";
          do echo "    $h";
            map_to_slim.py --slim_out=direct --association_file=$t$h$e"0.05GAF.txt" go-basic.obo goslim_plant_"$j".obo > "$t$h$e"slim_"$j".txt; #get the goslim_plant terms associated with the genes using goatools map_to_slim.py

            go=$(awk '{print $2}' "$t$h$e"slim_"$j".txt | tail -n +5 | tr ";" "\n" | sort -u); #list of unique GO terms for well/poorly collected genes
            o=$(awk '{print $2}' "$t$h$e"slim_"$j".txt | tail -n +5 | sed '/^$/d' | tr ";" "\n" | wc -l  | awk '{print $1}'); #number of function hits, i.e. total number of GO terms found in slim for well/poorly collected genes.
            n=$(wc -l "$t$h$e"slim_"$j".txt | awk '{print $1}'); #number of lines = number of genes for well/poorly collected genes
            n=$(( $n - 4 )); #subtract off header lines
            echo "GOterm ObsTotNumGenes ObsTotNumFunctionHits ObsGOcount obsfreq TotNumAnnotatedGenes ExpTotNumFuncHits ExpGOcount expfreq ratio absratio plratio diff absdiff" > Piechart"$t$h$e"_"$j".txt; #for each unique GO term, count how many times it occurs in slim.txt, divide by number of genes to get the proportion of genes representing that function in the M+ enriched/purified gene set.

            for i in $go;
              do f=$(grep $i "$t$h$e"slim_"$j".txt | wc -l | awk '{print $1}'); #get the number of lines the GO term is found in in the slim for well/poorly collected genes
                x=$(grep $i "$t"slim_"$j".txt | wc -l | awk '{print $1}'); #get the number of lines the GO term is found in in the slim for all genes
                g=$(echo "scale=4;$f/$o" | bc); #freq of function among all enriched functions
                y=$(echo "scale=4;$x/$a" | bc); #freq of function among all functions
                z=$(echo "$g - $y" | bc); #diff btw freq of function among all functions, and among enriched functions (negative means less common in enriched subset)
                zz=${z#-}; #absolute value of z
                rrx=$(echo "scale=4;$g/$y" | bc); #for log2 plotting
                #use below 4 lines for non-log2 axes
                if (( $(bc <<< "$g >= $y") )); 
                  then rr=$(echo "scale=4;$g/$y" | bc); #over represented results in positive fold enrichment
                  else rr=$(echo "scale=4;-$y/$g" | bc); #under represented results in negative fold enrichment
                fi; #ratio of frequency in enriched regions / frequency in all regions
                rrr=${rr#-}; #absolute value of rr
                 echo "$i $n $o $f $g $b $a $x $y $rr $rrr $rrx $z $zz" >> Piechart"$t$h$e"_"$j".txt;
              done;
              
            rm "$t$h$e"slim_"$j".txt;
          done;
      done;
  done;
**********ENDBASH**********

	Get verbal description of slimmed GO terms for eventual plotting. Include the cellular_component
	even though it is just the location the protein is found:
t="Sor";
for e in Gen;
  do echo "  $e";
    for h in "h" "l";
      do echo "    $h";
        > o.txt; #freq output
        > p.txt; #ratio output
        for j in biological_process molecular_function cellular_component;
          do echo $j;
            #filter out singletons, GO terms with only 1 observation in the enriched regions
            sin=$(grep -v "GO:"[0-9]*" "[0-9]*" "[0-9]*" "1" " Piechart$t$h$e"_"$j.txt | tail -n +2);

            #extract columns holding freq and ratio metrics for reformatting
            f=$(echo "$sin" | cut -d' ' -f1); #get GO slim terms
            p=$(echo "$sin" | cut -d' ' -f13); #get diff in freq of function btw enriched functions and all functions
            pp=$(echo "$sin" | cut -d' ' -f14); #get abs of diff in freq of function btw enriched functions and all functions
            q=$(echo "$sin" | cut -d' ' -f10); #get ratio between enriched/all
            qq=$(echo "$sin" | cut -d' ' -f11); #get abs of ratio between enriched/all
            pl=$(echo "$sin" | cut -d' ' -f12); #get ratio between enriched/all to be plotted on log2 axis
            c=$(for i in $f;
              do sed -n -e '/id: '$i'/,$p' goslim_plant_"$j".obo | head -2 | tail -1 | cut -d' ' -f2-;
              done;);
            d=$(for i in $f; do echo "$j"; done;); #create a header column with bp, cc, mf category
            paste -d$'\t' <(echo "$d") <(echo "$f") <(echo "$c") <(echo "$p")  <(echo "$pp") >> o.txt; #write to temporary freq output file
            paste -d$'\t' <(echo "$d") <(echo "$f") <(echo "$c") <(echo "$q")  <(echo "$qq") <(echo "$pl")>> p.txt; #write to temporary ratio output file
         done;
          echo Gocat$'\t'Goterm$'\t'name$'\t'diff$'\t'absdiff > PlotRpieFREQ$t$h$e.txt; #make header for freq output file
          echo Gocat$'\t'Goterm$'\t'name$'\t'ratio$'\t'absratio$'\t'plratio > PlotRpieRATIO$t$h$e.txt; #make header for ratio output file
          sort -t$'\t' -nr -k5,5 o.txt >> PlotRpieFREQ$t$h$e.txt;
          sort -t$'\t' -nr -k5,5 p.txt >> PlotRpieRATIO$t$h$e.txt;
         rm o.txt p.txt;
      done;
  done;
	
	Use R to construct charts:
**********RSCRIPT**********
#install.packages("ggplot2")
library(ggplot2)
options(error = recover)
rm(list=ls()) 
setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/ComparisonFigs/RecalculateM+2.50/SorghumRecalc2/+analysis/+SorGOWork/pie\ charts")

t<-"Sor"
rvals<-c("RATIO", "FREQ")
cvals<-c("h", "l")
evals<-c("Gen")
for (rr in rvals)
{
  for (cc in cvals)
  {
    for (ee in evals)
    {
      pdf(file=paste("GOchart",rr,t,cc,ee,".pdf", sep=""))
      #par(mfrow=c(1,2))
  
      g <- read.table(paste("PlotRpie",rr,t,cc,ee,".txt", sep=""), header=TRUE, sep="\t")
      h <- g[g$Goterm!="GO:0008150",] #Remove overarching GO categories biological_process, molecular_function
      h <- h[h$Goterm!="GO:0003674",]
      h <- h[h$Goterm!="GO:0005575",]
  
      par(mar = c(4,20,4,2) + 0.1)
      if(rr == "FREQ") {
        barplot(h$diff, main=paste("GOchart",rr,t,cc,ee,".pdf", sep=""), names.arg=paste(h$Goterm, h$name, sep=", "), las=2, horiz=TRUE, cex.names=0.6, xlim=c(min(h$diff),max(h$diff)))
        } 
      else {
        #barplot(h$plratio, log="x", main=paste("GOchart",rr,t,cc,ee,".pdf", sep=""), names.arg=paste(h$Goterm, h$name, sep=", "), las=2, horiz=TRUE, cex.names=0.6, xlim=c(min(h$plratio),max(h$plratio)))
        hh=h #transfer h table to a new variable, which will be modified
        hh$id=paste(hh$Goterm,hh$name,sep=" ") #assemble label term
        hh$id <- factor(hh$id, levels=hh$id) #fix order by making id a factor with levels

        ggout<-ggplot(hh,aes(id,plratio,width=0.8)) + 
        geom_bar(stat="identity",color="black",fill="dark grey",size=0.25) + 
        scale_y_continuous(trans='log2', breaks=c(0.125,0.25,0.5,1,2,4,8)) + 
        coord_flip() + 
        theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), text=element_text(size=12, family="ArialMT"))
        
        print(ggout)
        }      
      dev.off()
    }
  }
}
**********ENDR**********


@@@@@@@@@@@@@@@@@@@@ END of tested processing for SorRgen @@@@@@@@@@@@@@@@@@@@@@@@@@



























***ARABIDOPSIS***
	
	Calculate mean M+ enrichment across core sizes, for each locus, for each blocklength. This uses the
	SUMfiles as input and produces the +STATS files.  Do this on the cluster.



**********BASHPARALLEL-FASTER-RUNINSCREEN**********
mysd() {
       #calculate sum, mean, sd, n, mean,using Welford's algorithm and AWK
       #here, sum is the cumulative enrichment across all core sizes, mean is the average
       #supply space delimited string of values as argument like: mysd $foo, where foo="1 2 3 4.5"
       awk '{
         sum = 0;                   # Initialize running sum (for mean calculation)
         M = 0;
         S = 0;
         for (k=1; k <= NF; k++) { 
              sum += $k;                # Update running sum
              x = $k;
              oldM = M;
              M = M + ((x - M)/k);
              S = S + (x - M)*(x - oldM);
         }
         var = S/(NF - 1);
         print sum "\t" sum/NF "\t" sqrt(var) "\t" NF;
       }' <<< "$1"; # the first argument when calling mysd 
}
export -f mysd;

main() {
        l=$1;
        ofile="$p""+STATS.R"$e"Tgen.loc"$l".b"$b".txt"; #temporary outfile for each parallel process
        sumfileIN="$p""SUM.b"$b".R"$e"Tgen.txt";

        boo=$(grep -F $l "$sumfileIN" | awk '{ print $13 - $10 }' | tr "\n" " "); #extract all core sizes for this locus, calculate target enrichment
        st=$(mysd "$boo"); #calculate sum, mean, sd, and n using function mysd()
        echo "$b"$'\t'"$l"$'\t'"$st" > "$ofile"; #write stats to output file
}
export -f main;

#parameters
t="Arabidopsis";
elist="gen";
bmin=2; #min blocklength
bmax=50; #max blocklength
bstep=2; #step rate
nnode=8; #number of cluster nodes
blist=$(seq $bmin $bstep $bmax); #blocklengthrange
p=$(pwd)"/"; 
export p;

#parallelize on locus, within loop on blocklength
for e in $elist;
  do echo $e;
    export e;
    for b in $blist;
      do echo "b=$b";
        export b;
        sumfileIN="$p""SUM.b"$b".R"$e"Tgen.txt";
        uloc=$(cut -d$'\t' -f2 "$sumfileIN" | tail -n +2 | uniq); #find unique loci from all loci in column 2
        numloc=$(echo "$uloc" | wc -l);
        
        #parallel step
        #sends $gnuN records to each node via parallel step #1, this is piped to parallel step #2, which starts 96 jobs per node from the instruction set it has received
        echo "processing $sumfileIN...";
        gnuN=$(echo "$numloc / $nnode" | bc); #number of lines of input file to send to gnu parallel per node
        echo "$uloc" | parallel --sshloginfile ~/machines --jobs 1 --env main --env mysd --env e --env b --env p --pipe -N"$gnuN" /home/reevesp/bin/parallel -j96 --env main --env mysd --env e --env b --env p main;
        #run on head node only
        #echo "$uloc" | parallel --jobs 1 --env main --env mysd --env e --env b --env p --pipe -N"$gnuN" /home/reevesp/bin/parallel -j96 --env main --env mysd --env e --env b --env p main;
        #run on ceres, one node
        #echo "$uloc" | parallel --jobs 1 --env main --env mysd --env e --env b --env p --pipe -N"$gnuN" parallel -j96 --env main --env mysd --env e --env b --env p main;
        #run on Mac
        #echo "$uloc" | parallel --env main --env mysd --env e --env b --env p main;

        #consolidate results for each locus into one file
        ofile="+STATS.R"$e"Tgen.b"$b".txt"; #temporary outfile
        > "$p""$ofile";
        #find "$p" -name "+STATS.R"$e"Tgen.loc*.b"$b".txt" | while read line; do cat "$line" >> "$p""$ofile"; done; #fast cat, unsorted by $uloc
        for l in $uloc;
          do echo "concatenating loc file $l to +STATS.R"$e"Tgen.b"$b".txt ...";
            cat "$p""+STATS.R"$e"Tgen.loc"$l".b"$b".txt" >> "$p""$ofile";
          done;
          #test whether all were concatenated
          s1=$(ls "$p" | grep -c "loc"); #count the number of files with "loc" in their name using grep. can't wc -l because too many arguments
          s2=$(wc -l "$p""$ofile" |  awk '{print $1}');
         if [ $s1 == $s2 ]; 
           then find "$p" -maxdepth 1 -name "*loc*" -delete; #remove files when there are too many for rm
           else echo "s1="$s1", s2="$s2". Consolidated file "$ofile" no good. Quitting..." > err.txt;
           exit 1;
         fi;
      done; #$b
      
    #concatenate files for this Rgeo/Rpco
    statsfile="$p""+STATS."$t".R"$e"Tgen.txt";
    echo blocklength$'\t'locus$'\t'sum$'\t'mean$'\t'sd$'\t'n > "$statsfile";
    >"$statsfile"TMP;
    #quick cat and sort (some out files from ceres are not sorted correctly by locus)
    cat "$p""+STATS.R"$e"Tgen.b"*".txt" >> "$statsfile"TMP; #cat together the whole mess
    sortkey1=$(cut -d$'\t' -f2 "$statsfile"TMP | cut -d'.' -f3); #get a column containing only the locus index number, field 3 of column 2
    paste -d$'\t' "$statsfile"TMP <(echo "$sortkey1") | sort -t$'\t' -n -k1,1 -n -k7,7 | cut -d$'\t' -f1-6 > "$statsfile"TMP2; #add sort key to last column, sort by blocklength then locus index, remove last column
    cat "$statsfile"TMP2 >> "$statsfile";
    rm "$statsfile"TMP;
    rm "$statsfile"TMP2;

    #test whether all individual +STATS files have been added to the final concatenated file using number of lines
    s1=$(wc -l +STATS.R"$e"Tgen.b*.txt | grep total | awk '{print $1}');
    s2=$(wc -l "$statsfile" | awk '{print $1}');
    s2=$(( $s2 - 1 )); #decount header
    if [ $s1 == $s2 ]; 
      then rm +STATS.R"$e"Tgen.b*.txt; #remove files if number of lines in concat file is the same as the sum of all input files
      else echo "the number of lines ain't the same. something is 'crewed.";
    fi;

  done; #$e
**********BASHEND**********
	Takes ~15 minutes on blip. 


	add the proportion of NG,SS,NS sites from AtLocusWeights.txt to the +STATS summary file,
	generating the +v2STATS file.  Requires the AtLocusWeights.txt file as input.  

**********BASH**********
#use md5 on mac md5sum on linux
#modify weightsfile to include only the blocklengths analyzed
weightsfile="AtLocusWeights.txt";
bmin=2; #min blocklength
bmax=50; #max blocklength
bstep=2; #step rate
blist=$(seq $bmin $bstep $bmax); #blocklengthrange

wtfhead=$(head -1 "$weightsfile");
wtfbod=$(for i in $blist;
  do awk -F"\t" -v i=$i '$2 == i {print}' "$weightsfile"; #accumulate each blocklength needed
  done;)
wtf=$(echo "$wtfhead"; echo "$wtfbod");

#echo "$wtf" > wtf.txt; 

#wtf=$(awk -F"\t" '$2 != "1" {print}' "$weightsfile"); #exclude blocklength=1 from the weightsfile for use with blocklength 2-200
elist="gen";
for e in $elist;
  do statsfile="+STATS.Arabidopsis.R"$e"Tgen.txt";
    s=$(cut -d$'\t' -f2 $statsfile | tail -n +2 | md5sum); #get md5 of locus ID column from the +STATS file
    w1=$(echo "$wtf" | cut -d$'\t' -f4 | sed 's/_/./g'); #get the locus ID less blocklength index column, replace _ with . for md5 check
    w2=$(echo "$wtf" | cut -d$'\t' -f3); #get the blocklength index column
    w=$(paste -d'.' <(echo "$w1") <(echo "$w2") | tail -n +2 | md5sum); #get md5 of locus ID column from the Weights file
    if [ "$s" == "$w" ]; #verify that files are in same order using a checksum of the unique locus ID column
      then boo=$(echo "$wtf" | sed 's/blocklength/blocklength2/g' | sed 's/locus/locus2/g'); #modify the header in the Weights file so that there are no redundant column names after combining
        paste -d$'\t' "$statsfile" <(echo "$boo") > "+v2STATS.Arabidopsis.R"$e"Tgen.txt"; #paste the statsfile and the weights file together
      else echo "md5s do not match, $e, aborting.";
        kill -INT $$; #terminate the script, return to the shell
    fi;
  done;
**********BASHEND**********




#Plot summed enrichment by genomic position across blocklengths as a heatmap. Use Rscript
#AtGenomicGeography.r. Uses the +v2STATS files as input.



	Using the file +*SNPtoBlockMap.txt, calculate the sum of M+ enrichment value across blocklengths,
	for each SNP, i.e. for each position in the genome.  These data will be extracted from the 
	+v2STATS* files in several steps.  The first makes a file containing a table with the max 
	enrichment value at each site for each block, *EnrichAcrossBlocks*. Step 1:

**********BASHPARALLEL**********
myp() {
  e=$1;
  b=$2;
  v2name="$p""+v2STATS."$v2p".R"$e"Tgen.txt";
        
        echo "e=$e  b=$b";
        bout="$p""rr.R"$e"Tgen.b""$b"".tmp"; #temporary outfile

        #cut out the locus id and sum columns for current blocklenth
        gg=$(grep ^$b$'\t' "$v2name" | cut -d$'\t' -f1-3);
        
        #cut column with locus names for current blocklength from SNPtoBlockMap
        if [ $b = 1 ]; then
          cc=$(cut -d$'\t' -f4 "$tname" | tail -n +2 | sed 's/_/./g' | sed 's/$/./g');
        else col=$(( $b + 5 )); #calculate the column to cut out for each blocklength
          cc=$(cut -d$'\t' -f$col "$tname" | tail -n +2 | sed 's/_/./g' | sed 's/$/./g');
        fi;

        #search for locus names in $gg, retrieve enrichment value
        rr="b$b";

#put another parallel here if necessary
        for l in $cc;
          do if [[ $l = "--." ]]; then
                r="--"; #test for empty enrichment value
              else
                r=$(grep -m 1 ^"$b"$'\t'"$l" <<<"$gg" | awk '{print $3}');
              fi;
            rr+=$'\n'$r; #add current locus enrichment value to list
          done;

         #write out the result
         echo "$rr" > "$bout";
         #truncate -s -1 "$bout"; #remove the hanging newline, use 'gtruncate' only on osx
}
export -f myp;

#parameters
v2p="Arabidopsis"; export v2p;
t="At";
minbl=2; #1,2
maxbl=50; #50,200
blstep=2; #step
blist=$(seq $minbl $blstep $maxbl); #for even numbered runs 2-50
elist="gen";

p=$(pwd)"/"; export p;
tname="$p""+""$t""SNPtoBlockMap.txt"; export tname;

#begin parallel
parallel --sshloginfile ~/machines --jobs 24 --env myp --env p --env v2p --env tname myp ::: $elist ::: $blist;
#parallel --jobs 24 --env myp --env p --env v2p --env tname myp ::: $elist ::: $blist; #parallel command for head node only
#parallel --jobs 25 --env myp --env p --env v2p --env tname myp ::: $elist ::: $blist; #parallel command for ceres
#parallel --env myp --env p --env v2p --env tname myp ::: $elist ::: $blist; #parallel command for Mac

  
#reconstitute
echo "reconstituting..."
for e in $elist;
  do echo "e=$e";
    comm=""; #create a list of tmp files to paste together
    for b in $blist;
      do comm+="rr.R"$e"Tgen.b$b.tmp ";
    done;
    pos=$(cut -d$'\t' -f4 "$tname" | sed 's/'^.*_'//g' | sed 's/locus/pos/g');
    h=$(cut -d$'\t' -f1-6 "$tname" | tr $'\t' ' ');
    h2=$(paste -d' ' <(echo "$h") <(echo "$pos"));
    paste -d' ' <(echo "$h2") $comm > "+"$t"EnrichAcrossBlocks.R"$e"Tgen.txt";
    rm $comm; #clean up
  done;#e
  
**********BASHEND**********
	Takes ~8 hours.




	Process each line (SNP position) of *EnrichAcrossBlocks*.txt to get sum across all blocklengths
	of cumulative enrichment across all core sizes of max enrichment value from 10 reps of M+.
	Produces a file called SUM*EnrichAcrossBlocks... Step 2:

**********BASH**********
mysd() {
       #calculate sum, mean, sd, n, mean,using Welford's algorithm and AWK
       #here, sum is the cumulative enrichment across all blocklengths, mean is the average
       #supply space delimited string of values as argument like: mysd $foo, where foo="1 2 3 4.5"
       #returns space delimited string like "sum mean sd n"
       awk '{
         sum = 0;                   # Initialize running sum (for mean calculation)
         M = 0;
         S = 0;
         for (k=1; k <= NF; k++) { 
              sum += $k;                # Update running sum
              x = $k;
              oldM = M;
              M = M + ((x - M)/k);
              S = S + (x - M)*(x - oldM);
         }
         var = S/(NF - 1);
         print sum " " sum/NF " " sqrt(var) " " NF;
       }' <<< "$1"; # the first argument when calling mysd 
}
export -f mysd;

myp() {
    l=$1

    p=$(echo "$l" | cut -d' ' -f8-); #cut out the portion of the line with max enrichment values

    #skip loci that have missing blocks
    if [[ $p != *"--"* ]]; then
      m=$(mysd "$p"); #calculate sum, mean, sd, n
      chr=$(echo "$l" | cut -d' ' -f5);
      pos=$(echo "$l" | cut -d' ' -f7);
      i=$(echo "$l" | cut -d' ' -f1); #get genome index
      #loc=$chr"."$pos"."$(echo "$l" | cut -d' ' -f6);
      loc=$chr"."$pos"."$i;
      echo "snp=$i";

      echo $i" "$chr" "$pos" "$loc" "$m > "$path"$i".tmp"; #add new values to growing string, write to tmp file
    fi;
}
export -f myp;

#parameters
elist="gen";
t="At";
path=$(pwd)"/"; export path;

#parallel
for e in $elist;
do echo $e;
  #nnode=8;
  #numloc=$(wc -l +AtEnrichAcrossBlocks.RgenTgen.txt | awk '{print $1}');
  #numloc=$(( $numloc - 1 ));
  #gnuN=$(echo "$numloc / $nnode" | bc); #number of lines of input file to send to gnu parallel per node

  #parallel command
  #for cluster
  tail -n +2 "+"$t"EnrichAcrossBlocks.R"$e"Tgen.txt" | parallel  --sshloginfile ~/machines --jobs 1 --env myp --env mysd --env path --pipe --round-robin /home/reevesp/bin/parallel -j96  --env myp --env mysd --env path myp;
  #for head node only
  #tail -n +2 "+"$t"EnrichAcrossBlocks.R"$e"Tgen.txt" | parallel --jobs 1 --env myp --env mysd --env path --pipe --round-robin /home/reevesp/bin/parallel -j96  --env myp --env mysd --env path myp;
  #for ceres N=26
  #tail -n +2 "+"$t"EnrichAcrossBlocks.R"$e"Tgen.txt" | parallel  --jobs 26 --env myp --env mysd --env path myp;
  #for mac
  #tail -n +2 "+"$t"EnrichAcrossBlocks.R"$e"Tgen.txt" | parallel  --env myp --env mysd --env path myp;
  
  #assemble results
  echo "concatenating...";
  ofile="+"$t"SUMEnrichAcrossBlocks.R"$e"Tgen.txt";
  o="genomeindex chr pos loc sum mean sd n";
  echo "$o" > "$ofile";
  
  find "$path" -name "*.tmp" -print0 | xargs -0 cat | sort -n -t' ' -k1 >> "$ofile";
  
  #clean up
  echo "clean up...";
  find "$path" -name "*.tmp" -print0 | xargs -0 rm;
done;
**********ENDBASH**********
	Takes ~10 minutes, about 7 hours on ceres N=26

	Use R to make some pretty plots of significant enrichment across the genome.

**********BEGINR**********
#install.packages("fitdistrplus")
#install.packages("plotrix")
options(error = recover)
rm(list=ls()) 
setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Arabidopsis/Experiment2/+revision\ with\ hclust/+results/hclust/+AtRgenRuns/rarefaction/+analysis")
library(fitdistrplus)
library(plotrix)

t<-"At"
usezero="no" #this sets whether to include significant sites only when they are also > or < 0
             #to use the zero criterion, say "yes", to just use the ecdf around the mean, say "no"
pvals<-c(0.001,0.05)
devvals<-c("pdf","svg","png")
for (pp in pvals)
{
	for (zg in devvals)
	{
		zgg<-get(zg) #e.g convert string "pdf" to function pdf now named 'zgg'
		if ( zg == "pdf" | zg == "svg" ) { zgg(file=paste(t,"M+GenomeScan",pp, ".", zg, sep="")) }
		else { 
		par(mar=c(1,1,1,1))
		zgg(file=paste(t,"M+GenomeScan",pp, ".", zg, sep=""),width=2000,height=2000,res=300) } #call zgg function with relevant parameters
		
		par(mfrow=c(2,1))

			i="gen"
			g <- read.table(paste("+",t,"SUMEnrichAcrossBlocks.R",i,"Tgen.txt", sep=""), header=TRUE, sep=" ")
			mycdf <- ecdf(g$sum) #calculate the cdf of the empirical distribution
			p <- mycdf(g$sum) #get the probability of each observation (or lower)
			g$pecdf <- p #add new column to g, containing the probability from the empirical cdf
			g$index <- 1:nrow(g) #add a column with the sequential index
			h <- g[which(g$pecdf > 1-pp),] #get snps that are highly enriched
			if ( usezero == "yes" ) {
				hgeo <- h[which(h$sum>0),] #remove snps where sum<0, i.e. less than random
			}
			else {
				hgeo <- h # do not remove snps based on positive or negative
			}
			l <- g[which(g$pecdf < pp),] #get snps that are lowly enriched
			if ( usezero == "yes" ) {
				lgeo <- l[which(l$sum<0),] #remove snps where sum>0, i.e. more than random
			}
			else {
				lgeo <- l # do not remove snps based on positive or negative
			}
			highposgeo<-hgeo$index #get the sequential position of the highly enriched
			lowposgeo<-lgeo$index #get the sequential position of the deficient

			plot(g$index, g$sum, type="l", main=paste(t,i,pp,sep=""))
			color.scale.lines(g$index, g$sum, col=factor(g$chr))
			points(g$index[highposgeo], g$sum[highposgeo], col = "red", cex=0.7) #mark significant values
			points(g$index[lowposgeo], g$sum[lowposgeo], col = "blue", cex=0.7) #mark significant values

	
			#create BED formatted output for hgeo, hpco. verify first whether there exist any snps that are
			#both significantly enriched and better than random (or significantly suppressed and worse than random.)
			#if not, do not print an output table.
			if ( length(hgeo$pos) != 0 ) {
			hgeobed <- data.frame(chr=paste("Chr", hgeo$chr, sep=""), chromStart=hgeo$pos, chromEnd=hgeo$pos+1)	
			write.table(hgeobed, file=paste(t,"hGen",pp,".bed", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
			}
		
			#create BED formatted output for lgeo, lpco
			if ( length(lgeo$pos) != 0 ) {
			lgeobed <- data.frame(chr=paste("Chr", lgeo$chr, sep=""), chromStart=lgeo$pos, chromEnd=lgeo$pos+1)	
			write.table(lgeobed, file=paste(t,"lGen",pp,".bed", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)
			}
		
		dev.off()
	}
}
**********ENDR**********





	Determine the identity of each important SNP

	Download TAIR 10 annotation
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff 
	Convert to BED format using bedops, extract only genes from GFF file
cat TAIR10_GFF3_genes.gff | gff2bed | grep -w gene > TAIR10genes.bed

	Create background data file, a bed file containing only those genes for which there are SNPs.
	Create files listing genes enriched by M+ analyses for downstream GO over-representation analysis.
	Uses bedops. For whatever reason, this script has to be run in two parts or it fails.
**********BASH**********
***Part 1***
i="At"; #taxon identifier
ag="TAIR10genes.bed"; #bed file containing all annotated genes for given genome release

#Create a bed file for each SNP locus where we have a measurement (only needs to be done for geo or pco, since SNPs are the same):
chr=$(cut -d' ' -f2 "+"$i"SUMEnrichAcrossBlocks.RgenTgen.txt" | tail -n +2 | sed 's/^/Chr/g');
pos=$(cut -d' ' -f3 "+"$i"SUMEnrichAcrossBlocks.RgenTgen.txt" | tail -n +2);
pos2=$(awk '{$1 = $1 + 1; print}' <<<"$pos"); 
paste -d$'\t' <(echo "$chr") <(echo "$pos") <(echo "$pos2") > $i"AllSNPs.bed";

#Sort to make sure it will work for bedops --element-of
/usr/local/bin/sort-bed $i"AllSNPs.bed" > tmp.txt;
mv tmp.txt $i"AllSNPs.bed";

#Calculate the intersection of AtAllSNPs.bed and TAIR10genes.bed, the list of genes for which there are SNPs in the study (the background):
/usr/local/bin/bedops --element-of 1 "$ag" $i"AllSNPs.bed" > $i"AllGenesWithSNPs.bed";
  
#Make background bed files into a simple gene list.
cut -d$'\t' -f4 $i"AllGenesWithSNPs.bed" | sort -u > $i"BackgroundGenesList.txt";
#clean up
rm $i"AllSNPs.bed";
rm $i"AllGenesWithSNPs.bed";

***Part 2***
#Compute lists of genes enriched by M+
for k in "Gen";
  do echo $k;
    #Verify that bed files are sorted correctly using bedops:sort-bed, e.g.
    for j in "0.001" "0.05";
      do echo "  $j";
        for h in "h" "l";
          do echo "    $h";
            /usr/local/bin/sort-bed $i$h$k$j.bed > tmp.txt;
            mv tmp.txt $i$h$k$j.bed;
    
            #Determine nearest gene to highly enriched snps using closest-features
            /usr/local/bin/closest-features --closest --dist $i$h$k$j.bed "$ag" > $i$h"$k"EnrichedGenes$j.bed;
    
            #Extract hits that are within genes, not just close
            grep "|0$" $i$h"$k"EnrichedGenes$j.bed > $i$h"$k"EnrichedWithinGenes$j.bed;
            rm $i$h"$k"EnrichedGenes$j.bed;
    
            #Get gene identifiers for downstream GO analysis
            cut -d$'\t' -f6 $i$h"$k"EnrichedWithinGenes$j.bed | sort -u > $i$h"$k"EnrichedGenesList$j.txt;
            rm $i$h"$k"EnrichedWithinGenes$j.bed;
          done;
      done;
  done;
**********ENDBASH**********

	Use Python GOATOOLS to analyze over-representation:
easy_install goatools;
easy_install fisher;
easy_install statsmodels;
wget http://geneontology.org/ontology/go-basic.obo;
wget http://www.geneontology.org/ontology/subsets/goslim_generic.obo;

	Download a raw gene association file (GAF) from http://amigo.geneontology.org/amigo/search/annotation
	Choose species in 'Species' menu, remove 'not', 'contributes_to', and "colocalizes_with'
	in 'Annotation Qualifier' menu.  For Arabidopsis, this was done in three downloads, one for
	Cellular Component, one for Molecular Function, and one for Biological Process due to limits
	on the number of lines that could be downloaded at once. Use cat * and call it ATRawGAF.txt. This GAF will include GO terms for 
	Biological Process, Molecular Function, and Cellular Component ontologies.
	Process the file into an "association" input file, called AtGAF.txt, for GOATOOLS script find_enrichment.py


	THIS DOESNT NEED TO BE REPEATED. USE FROM FINAL ANALYSIS 2. SEE BELOW FOR A
	PARALLEL VERSION
**********BASH**********
t=At;
cut -d$'\t' -f3 "$t"RawGAF.txt | sort -u > GAFgenestmp.txt;
raw=$(cut -d$'\t' -f3,5 "$t"RawGAF.txt);
>tmp.txt;
lf=$(wc -l GAFgenestmp.txt | awk '{print $1}'); #get total line number for progress bar
z=0; #counter
while read -r line;
do
    z=$(($z+1));
    echo -en "\r"$(($lf-$z)); #display countdown progress 
    f="";
    f=$(grep ^"$line"$'\t' <<<"$raw" | cut -d$'\t' -f2 | tr '\n' ';');
    if [[ "$f" != "" ]]; then
      echo "$line"$'\t'"$f" >> tmp.txt;
    fi;
done < GAFgenestmp.txt;
sed 's/;$//g' tmp.txt > "$t"GAF.txt;
rm GAFgenestmp.txt;
rm tmp.txt;
**********ENDBASH**********

**********BASHPARALLEL**********
myp() {
    line=$1;
    o=$(mktemp /tmp/tmp.XXXXXXXX);
    f="";
    f=$(grep ^"$line"$'\t' in.txt | cut -d$'\t' -f2 | tr '\n' ';');
    if [[ "$f" != "" ]]; then
      echo "$line"$'\t'"$f" > $o;
    fi;
}
export -f myp;

t=At;
cut -d$'\t' -f3 "$t"RawGAF.txt | sort -u > GAFgenestmp.txt;
cut -d$'\t' -f3,5 "$t"RawGAF.txt > in.txt;
find /tmp/ -name "tmp.*" -print0 | xargs -0 rm;
#rm /tmp/tmp.*;
cat GAFgenestmp.txt | parallel --env myp myp;


#cat /tmp/tmp.* | sort | sed 's/;$//g' > ./"$t"GAF.txt;
find /tmp/ -name "tmp.*" -print0 | xargs -0 cat > ./"$t"GAFtmp.txt;
sort ./"$t"GAFtmp.txt | sed 's/;$//g' > ./"$t"GAF.txt; #clean up GAFtmp
rm GAFgenestmp.txt;
rm ./"$t"GAFtmp.txt;
find /tmp/ -name "tmp.*" -print0 | xargs -0 rm;
#rm /tmp/tmp.*;
**********ENDBASH**********
	Takes a minute or so on glitch. 10 minutes on Mac.

	One gene, FIP1[V], is missing because grep can't handle the brackets.  So, get it manually:
f=$(grep ^'FIP1\[V\]'$'\t' in.txt | cut -d$'\t' -f2 | tr '\n' ';' | sed 's/;$//g');
cat AtGAF.txt <(echo 'FIP1[V]'$'\t'"$f") | sort > tmp.txt;
mv tmp.txt AtGAF.txt;
rm in.txt;





@@@@@@@@@@@@@@@@@@@@ END of tested processing for AtRgen @@@@@@@@@@@@@@@@@@@@@@@@@@








	Run GOATOOLS find_enrichment.py to get the enrichment status, over (e), or under (p),
	representation (column 3 in output).
	The --no_propagate_counts selects only the least inclusive GO term, i.e no parent terms.
	This eliminates multiple significant results along a parent-child path, but has the
	undesirable consequence of making significant under-representation (p) a questionable
	result because higher level terms are "shorted" their descendant terms.

	Move the relevant items (inc. *GAF.txt, go-basic.obo) to a nested folder, "statistical tests", then run
	the GOATOOLS script find_enrichment.py:
**********BASH**********
mv *EnrichedGenesList* "statistical tests";
mv AtBackgroundGenesList.txt "statistical tests";
cd "statistical tests";

t=At;
comm="--alpha 0.05 --pval 0.05 --obo go-basic.obo --no_propagate_counts --method holm";

for p in 0.05 0.001;
  do echo $p;
  for e in Geo Pco;
    do echo "  $e";
      for h in "h" "l";
        do echo "    $h";
          find_enrichment.py $(echo "$comm" "$t""$h""$e"EnrichedGenesList"$p".txt "$t"BackgroundGenesList.txt "$t"GAF.txt) > +"$t""$h""$e""$p"GOAout.txt;
        done;
    done;
  done;
**********ENDBASH**********

	Examine first part of output files to determine if any GO terms are significantly over-
	or under-represented.
head -20 +*GOA*;

	FINAL RESULT:

For MaxBL=50:
	There are four significantly over-represented GO terms for Arabidopsis among genes that are
	well-collected by environmental data and Hclust.  Three of these GO terms are close but no cigar over-represented
	among genes that are well-collected by geographic data. The three significant results for +AthPco0.05GOAout.txt
	are attributable to the same cluster of genes: AT4G03440 AT4G03450 AT4G03460 AT4G03470 AT4G03480 AT4G03490 AT4G03500.
	These can be found in /Volumes/J22/M+ Study/Analysis/Final analysis for 2nd paper/Arabidopsis/Experiment2/+revision with hclust/+results/hclust/+AtHclustGO work/ToMaxBL50/gene identification/README.txt.
	Similar findings occurred for MaxBL=200 Hclust and Rcut.
	Using Hclust:
	==> +AthPco0.001GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:0009117	BP	e	nucleotide metabolic process	2/33	4/24259	1.07e-05	n.a.	2	0.0183
	==> +AthPco0.05GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:2000031	BP	e	regulation of salicylic acid mediated signaling pathway	7/1424	10/24259	2.43e-07	n.a.	7	0.000413
	GO:0071446	BP	e	cellular response to salicylic acid stimulus	7/1424	11/24259	6.34e-07	n.a.	7	0.00108
	GO:0031347	BP	e	regulation of defense response	7/1424	17/24259	2.73e-05	n.a.	7	0.0464
	
For MaxBL=200:
	Using Rcut:
	There are three significantly over-represented GO terms among Arabidopsis genes that are well-collected
	using geographic data. There are no significantly over- or under-represented GO terms among 
	Arabidopsis genes that are well- or poorly- collected using environmental data and Rcut.
	==> +AthGeo0.05GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:2000031	BP	e	regulation of salicylic acid mediated signaling pathway	7/896	10/24203	1.01e-08	n.a.	7	1.73e-05
	GO:0071446	BP	e	cellular response to salicylic acid stimulus	7/896	11/24203	2.7e-08	n.a.	7	4.59e-05
	GO:0031347	BP	e	regulation of defense response	7/896	17/24203	1.31e-06	n.a.	7	0.00222

	Using Rcut2:
	There is one significantly under-represented GO term among Arabidopsis genes that are well-collected
	using environmental data.  There are two significantly over-represented GO terms among Arabidopsis genes
	that are poorly collected using environmental data and Rcut2.
	==> +AthPco0.05GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:0008150	BP	p	biological_process	419/2232	5740/24203	4.67e-09	n.a.	419	7.95e-06
	==> +AtlPco0.05GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:0006952	BP	e	defense response	27/848	235/24203	6.76e-08	n.a.	27	0.000115
	GO:0008150	BP	e	biological_process	254/848	5740/24203	2.23e-05	n.a.	254	0.0379

	Using Hclust:
	There are three significantly over-represented GO terms for Arabidopsis among genes that are
	well-collected by geographic data and Hclust.  The same three GO terms are significantly over-represented
	among genes that are well-collected by environmental data. These same findings occurred for Rcut.
	There is one significantly over-represented GO term among genes that are poorly collected
	by geographic data.
	==> +AthGeo0.05GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:2000031	BP	e	regulation of salicylic acid mediated signaling pathway	7/1236	10/24203	9.33e-08	n.a.	7	0.000159
	GO:0071446	BP	e	cellular response to salicylic acid stimulus	7/1236	11/24203	2.45e-07	n.a.	7	0.000417
	GO:0031347	BP	e	regulation of defense response	7/1236	17/24203	1.1e-05	n.a.	7	0.0187
	==> +AthPco0.05GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:2000031	BP	e	regulation of salicylic acid mediated signaling pathway	7/1304	10/24203	1.35e-07	n.a.	7	0.00023
	GO:0071446	BP	e	cellular response to salicylic acid stimulus	7/1304	11/24203	3.53e-07	n.a.	7	0.000601
	GO:0031347	BP	e	regulation of defense response	7/1304	17/24203	1.56e-05	n.a.	7	0.0266
	==> +AtlGeo0.05GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:0000272	BP	e	polysaccharide catabolic process	5/1576	6/24203	6.6e-06	n.a.	5	0.0112

	Using Kmeans:
	There are two significantly over-represented GO terms for Arabidopsis among genes in regions that are
	well-collected by environmental information:
	==> +AthPco0.05GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:2000031	BP	e	regulation of salicylic acid mediated signaling pathway	6/1652	10/24203	1.66e-05	n.a.	6	0.0282
	GO:0009117	BP	e	nucleotide metabolic process	4/1652	4/24203	2.16e-05	n.a.	4	0.0368
	There are three significantly over-represented GO terms for Arabidopsis among genes in regions that are
	poorly-collected using geographic information:
	==> +AtlGeo0.001GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:0000272	BP	e	polysaccharide catabolic process	4/44	6/24203	1.42e-10	n.a.	4	2.42e-07
	GO:0006032	BP	e	chitin catabolic process	4/44	9/24203	1.19e-09	n.a.	4	2.02e-06
	GO:0016998	BP	e	cell wall macromolecule catabolic process	4/44	13/24203	6.71e-09	n.a.	4	1.14e-05

	Above result with p<0.001 consistent with below result with p<0.05:
	==> +AtlGeo0.05GOAout.txt <==
	GO	NS	enrichment	name	ratio_in_study	ratio_in_pop	p_uncorrected	depth	study_count	p_holm
	GO:0000272	BP	e	polysaccharide catabolic process	5/1859	6/24203	1.49e-05	n.a.	5	0.0254




	Make a pie chart showing GOslim category representation for genes enriched by M+ at the
	0.05 level.

	Extract the three GO categories (biological_process, cellular_component, molecular_function)
	from goslim_plant.obo.  Use PERL modules go-perl > go-filter-subset.pl.

for i in biological_process cellular_component molecular_function;
  do
    /Users/shrub/perl5/bin/go-filter-subset.pl -namespace "$i" goslim_plant.obo > goslim_plant_"$i".obo;
  done;

	For whatever reason (actually, it is because go-filter-subset.pl includes 'part-of's, not
	just 'is_a's as members of a category), two "cellular_component"s remain in goslim_plant_biological_process.obo.
	Ten "biological_process"s remain in goslim_plant_molecular_function.obo.
	Remove them manually.
	
	Make a gene association file (GAF) for M+ well-collected and poorly-collected genes.
	Input files: At[h,l][Geo,Pco]EnrichedGenesList0.05.txt, AtRawGAF.txt

**********BASHPARALLEL**********
myp() {
    line=$1;
    o=$(mktemp /tmp/tmp.XXXXXXXX);
    f="";
    f=$(grep ^"$line"$'\t' in.txt | cut -d$'\t' -f2 | tr '\n' ';');
    if [[ "$f" != "" ]]; then
      echo "$line"$'\t'"$f" > $o;
    fi;
}
export -f myp;

t=At;
rm /tmp/tmp.*; #clear tmp directory of tmp. files generated by this script
cut -d$'\t' -f3,5 "$t"RawGAF.txt > in.txt; #make a simplified GAF file of all genes to query
  for e in Geo Pco;
    do echo "  $e";
      for h in "h" "l";
        do echo "    $h";
          cat $t$h$e"EnrichedGenesList0.05.txt" | parallel --env myp myp;

          cat /tmp/tmp.* | sort | sed 's/;$//g' > ./$t$h$e"0.05GAF.txt";
          rm /tmp/tmp.*;
        done;
    done;
rm in.txt;
**********ENDBASH**********

	Count the number of occurrences of each GO slim term in the gene association files
	derived from the EnrichedGenesList(s). These are the files At[h,l][Pco,Geo]0.05GAF.txt.
	Do this for each major GO category (bp, cc, mf). Output files are like At[h,l][Pco,Geo]piechart_[bp,cc,mf].txt.
	Goal here is to calculate the frequency of the function in the well/poorly collected genes and
	compare that to the frequency of the same functions when all annotated genes are considered.
**********BASH**********
  
t="At";
b=$(wc -l "$t"GAF.txt | awk '{print $1}'); #the total number of annotated genes
for j in biological_process cellular_component molecular_function;
  do echo $j;
    map_to_slim.py --slim_out=direct --association_file=$t"GAF.txt" go-basic.obo goslim_plant_"$j".obo > "$t"slim_"$j".txt; #get the goslim_plant terms associated with all annotated genes using goatools map_to_slim.py
    a=$(awk -F$'\t' '{print $2}' "$t"slim_"$j".txt | tail -n +5 | sed '/^$/d' | tr ";" "\n" | wc -l  | awk '{print $1}'); #number of function hits, i.e. total number of GO terms found in slim for all collected genes.

    for e in Geo Pco;
      do echo "  $e";
        for h in "h" "l";
          do echo "    $h";
            map_to_slim.py --slim_out=direct --association_file=$t$h$e"0.05GAF.txt" go-basic.obo goslim_plant_"$j".obo > "$t$h$e"slim_"$j".txt; #get the goslim_plant terms associated with the genes using goatools map_to_slim.py

            go=$(awk '{print $2}' "$t$h$e"slim_"$j".txt | tail -n +5 | tr ";" "\n" | sort -u); #list of unique GO terms for well/poorly collected genes
            o=$(awk '{print $2}' "$t$h$e"slim_"$j".txt | tail -n +5 | sed '/^$/d' | tr ";" "\n" | wc -l  | awk '{print $1}'); #number of function hits, i.e. total number of GO terms found in slim for well/poorly collected genes.
            n=$(wc -l "$t$h$e"slim_"$j".txt | awk '{print $1}'); #number of lines = number of genes for well/poorly collected genes
            n=$(( $n - 4 )); #subtract off header lines
            echo "GOterm ObsTotNumGenes ObsTotNumFunctionHits ObsGOcount obsfreq TotNumAnnotatedGenes ExpTotNumFuncHits ExpGOcount expfreq ratio absratio plratio diff absdiff" > Piechart"$t$h$e"_"$j".txt; #for each unique GO term, count how many times it occurs in slim.txt, divide by number of genes to get the proportion of genes representing that function in the M+ enriched/purified gene set.

            for i in $go;
              do f=$(grep $i "$t$h$e"slim_"$j".txt | wc -l | awk '{print $1}'); #get the number of lines the GO term is found in in the slim for well/poorly collected genes
                x=$(grep $i "$t"slim_"$j".txt | wc -l | awk '{print $1}'); #get the number of lines the GO term is found in in the slim for all genes
                g=$(echo "scale=4;$f/$o" | bc); #freq of function among all enriched functions
                y=$(echo "scale=4;$x/$a" | bc); #freq of function among all functions
                z=$(echo "$g - $y" | bc); #diff btw freq of function among all functions, and among enriched functions (negative means less common in enriched subset)
                zz=${z#-}; #absolute value of z
                rrx=$(echo "scale=4;$g/$y" | bc); #for log2 plotting
                #use below 4 lines for non-log2 axes
                if (( $(bc <<< "$g >= $y") )); 
                  then rr=$(echo "scale=4;$g/$y" | bc); #over represented results in positive fold enrichment
                  else rr=$(echo "scale=4;-$y/$g" | bc); #under represented results in negative fold enrichment
                fi; #ratio of frequency in enriched regions / frequency in all regions
                rrr=${rr#-}; #absolute value of rr
                 echo "$i $n $o $f $g $b $a $x $y $rr $rrr $rrx $z $zz" >> Piechart"$t$h$e"_"$j".txt;
              done;
              
            rm "$t$h$e"slim_"$j".txt;
          done;
      done;
  done;
**********ENDBASH**********


	Get verbal description of slimmed GO terms for eventual plotting. Exclude the cellular_component
	as it is just the location the protein is found:
t="At";
for e in Geo Pco;
  do echo "  $e";
    for h in "h" "l";
      do echo "    $h";
        > o.txt; #freq output
        > p.txt; #ratio output
        for j in biological_process molecular_function;
          do echo $j;
            #filter out singletons, GO terms with only 1 observation in the enriched regions
            sin=$(grep -v "GO:"[0-9]*" "[0-9]*" "[0-9]*" "1" " Piechart$t$h$e"_"$j.txt | tail -n +2);

            #extract columns holding freq and ratio metrics for reformatting
            f=$(echo "$sin" | cut -d' ' -f1); #get GO slim terms
            p=$(echo "$sin" | cut -d' ' -f13); #get diff in freq of function btw enriched functions and all functions
            pp=$(echo "$sin" | cut -d' ' -f14); #get abs of diff in freq of function btw enriched functions and all functions
            q=$(echo "$sin" | cut -d' ' -f10); #get ratio between enriched/all
            qq=$(echo "$sin" | cut -d' ' -f11); #get abs of ratio between enriched/all
            pl=$(echo "$sin" | cut -d' ' -f12); #get ratio between enriched/all to be plotted on log2 axis
            c=$(for i in $f;
              do sed -n -e '/id: '$i'/,$p' goslim_plant_"$j".obo | head -2 | tail -1 | cut -d' ' -f2-;
              done;);
            d=$(for i in $f; do echo "$j"; done;); #create a header column with bp, cc, mf category
            paste -d$'\t' <(echo "$d") <(echo "$f") <(echo "$c") <(echo "$p")  <(echo "$pp") >> o.txt; #write to temporary freq output file
            paste -d$'\t' <(echo "$d") <(echo "$f") <(echo "$c") <(echo "$q")  <(echo "$qq") <(echo "$pl") >> p.txt; #write to temporary ratio output file
         done;
          echo Gocat$'\t'Goterm$'\t'name$'\t'diff$'\t'absdiff > PlotRpieFREQ$t$h$e.txt; #make header for freq output file
          echo Gocat$'\t'Goterm$'\t'name$'\t'ratio$'\t'absratio$'\t'plratio > PlotRpieRATIO$t$h$e.txt; #make header for ratio output file
          sort -t$'\t' -nr -k5,5 o.txt >> PlotRpieFREQ$t$h$e.txt;
          sort -t$'\t' -nr -k5,5 p.txt >> PlotRpieRATIO$t$h$e.txt;
         rm o.txt p.txt;
      done;
  done;
	
	Use R to construct charts, change hclust and kmeans as necessary in setwd:
**********RSCRIPT**********
#install.packages("ggplot2")
#library(ggplot2)
options(error = recover)
rm(list=ls()) 
setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Arabidopsis/Experiment2/+revision\ with\ hclust/+results/hclust/+AtHclustGO\ work/ToMaxBL50/pie\ charts")

t<-"At"
rvals<-c("RATIO", "FREQ")
cvals<-c("h", "l")
evals<-c("Gen")
for (rr in rvals)
{
  for (cc in cvals)
  {
    for (ee in evals)
    {
      pdf(file=paste("GOchart",rr,t,cc,ee,".pdf", sep=""))
      #par(mfrow=c(1,2))
  
      g <- read.table(paste("PlotRpie",rr,t,cc,ee,".txt", sep=""), header=TRUE, sep="\t")
      h <- g[g$Goterm!="GO:0008150",] #Remove overarching GO categories biological_process, molecular_function
      h <- h[h$Goterm!="GO:0003674",]
      h <- h[h$Goterm!="GO:0005575",]
  
      par(mar = c(4,20,4,2) + 0.1)
      if(rr == "FREQ") {
        barplot(h$diff, main=paste("GOchart",rr,t,cc,ee,".pdf", sep=""), names.arg=paste(h$Goterm, h$name, sep=", "), las=2, horiz=TRUE, cex.names=0.6, xlim=c(min(h$diff),max(h$diff)))
        } 
      else {
        #barplot(h$plratio, log="x", main=paste("GOchart",rr,t,cc,ee,".pdf", sep=""), names.arg=paste(h$Goterm, h$name, sep=", "), las=2, horiz=TRUE, cex.names=0.6, xlim=c(min(h$plratio),max(h$plratio)))
        hh=h #transfer h table to a new variable, which will be modified
        hh$id=paste(hh$Goterm,hh$name,sep=" ") #assemble label term
        hh$id <- factor(hh$id, levels=hh$id) #fix order by making id a factor with levels

        ggout<-ggplot(hh,aes(id,plratio,width=0.8)) + 
        geom_bar(stat="identity",color="black",fill="dark grey",size=0.25) + 
        scale_y_continuous(trans='log2', breaks=c(0.125,0.25,0.5,1,2,4,8)) + 
        coord_flip() + 
        theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), text=element_text(size=12, family="ArialMT"))
        
        print(ggout)
        }      
      dev.off()
    }
  }
}
**********ENDR**********

	For Kmeans, because sometimes well-collected regions do not exist for pco because they are all worse than random.
**********RSCRIPT**********
options(error = recover)
rm(list=ls()) 
setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Arabidopsis/Experiment2/+revision\ with\ Rcut/+AtRcutGO\ work/pie\ charts")

t<-"At"
cvals<-c("h", "l")
evals<-c("Pco","Geo")
for (cc in cvals)
{
  for (ee in evals)
  {
    pdf(file=paste("GOchart",t,cc,ee,".pdf", sep=""))
    #par(mfrow=c(1,2))

    if ( ( cc=="h" ) & ( ee=="Pco" ) ) {
      #do nothing
    }
    else {
      g <- read.table(paste("PlotRpie",t,cc,ee,".txt", sep=""), header=TRUE, sep="\t")
      h <- g[g$Goterm!="GO:0008150",] #Remove overarching GO categories biological_process, molecular_function
      h <- h[h$Goterm!="GO:0003674",]

      par(mar = c(4,20,4,2) + 0.1)
      barplot(h$diff, main=paste("GOchart",t,cc,ee,".pdf", sep=""), names.arg=paste(h$Goterm, h$name, sep=", "), las=2, horiz=TRUE, cex.names=0.6, xlim=c(min(h$diff),max(h$diff)))
    }
    dev.off()
  }
}
**********ENDR**********






*******POSTANALYSIS TRICKS********
	Calculate the number of MF and BP GO terms that are enriched in well- and poorly-conserved genomic regions.
	Start in the folder 'pie charts':

t="At";
for e in Geo Pco;
do for h in h l;
  do echo $h$e;
  awk -F$'\t' '$6 >= 1.5' PlotRpieRATIO"$t$h$e".txt | tail -n +2 | wc -l;
  done;
done;









**********TOOLS**********

	Check whether +STATS* files have the same number of lines.
	They should:
wc -l +STATS*;

	Check whether +STATS* files have the same loci in column 2:
g=$(cut -d$'\t' -f2 +STATS.*.RgeoTgen.txt);
p=$(cut -d$'\t' -f2 +STATS.*.RpcoTgen.txt);
diff <(echo "$g") <(echo "$p");


	Check whether SUM.RgeoTgen and SUM.RpcoTgen have the same number of lines.
	They should:

v="1";
#(for i in $v;
(for i in {1..1000}; 
do f=$(wc -l SUM.b$i.R*Tgen.txt| awk '{print $1}' | head -2 | sort -u);
  if (( $(echo "$f" | wc -l) == 2 )); then
    ls -l SUM.b$i.R*Tgen.txt;
    echo;
  fi;
done;) > xlinenum.txt;

#If same (or different), figure out if Rgeo and Rpco list different loci:
#(for i in $v;
(for i in {1..1000}; 
do echo $i;
  g=$(cut -d$'\t' -f2 SUM.b$i.RgeoTgen.txt | sort -u | md5sum);
  p=$(cut -d$'\t' -f2 SUM.b$i.RpcoTgen.txt | sort -u | md5sum);
  if [[ $g != $p ]]; then
    ls -l SUM.b$i.R*Tgen.txt;
    echo;
  fi;
done;) > xdiffloci.txt;

If there are some files that have the same number of lines, but differ in the loci listed
in those lines, figure out where the differences are:

for i in "6";  #cycle through list of offending files
do g=$(cut -d$'\t' -f2 SUM.b$i.RgeoTgen.txt);
  p=$(cut -d$'\t' -f2 SUM.b$i.RpcoTgen.txt);
  diff <(echo "$g") <(echo "$p") > diff$i.txt;
done;




If the same, compute +STATS* files.

If different, determine where the locus names differ between +STATS* files
sg=$(cut -d$'\t' -f2 +STATS.Sorghum.RgeoTgen.txt);
sp=$(cut -d$'\t' -f2 +STATS.Sorghum.RpcoTgen.txt);
diff <(echo "$sg") <(echo "$sp") > f.txt

Determine which blocklengths have a scrambled order relative to what's expected:

#(for i in $v;
(for i in {1..1000}; 
do g=$(cut -d$'\t' -f2 SUM.b$i.RgeoTgen.txt | md5sum);
  gsort=$(cut -d$'\t' -f2 SUM.b$i.RgeoTgen.txt | sort -t'.' -n -k1,1 -k2,2 | md5sum);
  p=$(cut -d$'\t' -f2 SUM.b$i.RpcoTgen.txt | md5sum);
  psort=$(cut -d$'\t' -f2 SUM.b$i.RpcoTgen.txt | sort -t'.' -n -k1,1 -k2,2 | md5sum);
  if [[ $g != $gsort ]]; then
    echo "SUM.b$i.RgeoTgen.txt";
  fi;
  if [[ $p != $psort ]]; then
    echo "SUM.b$i.RpcoTgen.txt";
  fi;
done; ) > xscrambled.txt

#Verify that all loci are included, done on Mac:
(nq=63; #(AtHclust=117,AtKmeans=89 ,AtRcut=27) number of env+geo+pco characters, to be cut off
        #(SorHclust=63,SorKmeans=53)
t="../../script7outputHclust/SorM";
#for i in $v;
for i in {1..1000};
do totloc=$(awk '{print $1}' "$t".b"$i".RgenTenv.var | tail -n +4 | ghead -n -$nq | sort -u);
  g=$(cut -d$'\t' -f2 SUM.b$i.RgeoTgen.txt | tail -n +2 | sort -u);
  p=$(cut -d$'\t' -f2 SUM.b$i.RpcoTgen.txt | tail -n +2 | sort -u);
  if [[ $g != $totloc ]]; then
    echo "problem: SUM.b$i.RgeoTgen.txt";
  fi;
  if [[ $p != $totloc ]]; then
    echo "problem: SUM.b$i.RpcoTgen.txt";
  fi;
done; ) > xallloci.txt

If all loci are included, but order is scrambled, the file can be sorted:

****
myp() {
      i=$1;
      echo $i.$e;
      sortkey1=$(cut -d$'\t' -f2 SUM.b"$i".R"$e"Tgen.txt | cut -d'.' -f3);
      paste -d$'\t' SUM.b"$i".R"$e"Tgen.txt <(echo "$sortkey1") | sort -t$'\t' -n -k28 -n -k3,3 > tmp"$i"."$e".txt 
      #Determine whether blocklengths are scrambled again
      g=$(cut -d$'\t' -f2 tmp"$i"."$e".txt | md5sum);
      gsort=$(cut -d$'\t' -f2 tmp"$i"."$e".txt | sort -t'.' -n -k1,1 -k2,2 | md5sum);
      if [[ $g != $gsort ]]; then
        echo "tmp"$i"."$e".txt is still screwed";
      fi;
}
export -f myp;

slist=$(grep SUM xscrambled.txt | cut -d'.' -f2 | sed 's/b//g' | sort -nu); #get list of scrambled output files
elist="geo pco";
for e in $elist;
  do export e; 
    echo "$slist" | parallel --env myp --env e myp;
  done;
****




#If no longer scrambled, mv back to original file name, eliminating last column (added for sorting):
****
mya() {
      i=$1;
      echo $i.$e; 
      rev tmp"$i"."$e".txt | cut -d$'\t' -f2- | rev > SUM.b"$i".R"$e"Tgen.txt;
      rm tmp"$i"."$e".txt;
}
export -f mya;

slist=$(grep SUM xscrambled.txt | cut -d'.' -f2 | sed 's/b//g' | sort -nu); #get list of scrambled output files
elist="geo pco";
for e in $elist;
  do export e; 
    echo "$slist" | parallel --env mya --env e mya;
  done;
****



#Determine whether all loci contain the same number of reported lines in the SUM file.
#Here, the grep -v " 22 " refers to the number of core sizes searched for Sorghum, 2-23.
#This shows loci where there are not 22 summary lines:

cut -d$'\t' -f2 SUM.b1.RgeoTgen.txt | sort -n -t'.' -k3,3 | uniq -c | grep -v " 21 "; # returns '1 locus' if all loci have the correct number of lines
cut -d$'\t' -f2 SUM.b1.RpcoTgen.txt | sort -n -t'.' -k3,3 | uniq -c | grep -v " 21 ";


#repeat needed analyses:
a=$(cat rgeobad.txt);
for i in $a; do z=$(echo "$i" | sed 's/-/ /g'); ./8dSorghumSlice1RgeoCeresLociSubsetBeginEnd.sh $z; done;
for i in $a; do z=$(echo "$i" | sed 's/-/ /g'); ./8dAtSlice1RpcoLociSubsetBeginEnd.sh $z; done;


#take old sum file, remove bad lines from above:
v="10.948500.375681 10.952995.375683 10.953073.375684 10.953103.375685 1.949028.1064 5.953135.229205"
grep -e 10.948500.375681 -e 10.952995.375683 -e 10.953073.375684 -e 10.953103.375685 -e 1.949028.1064 -e 5.953135.229205 SUM.b1.RgeoTgen.txtORIG > SUM.b1.RgeoTgen.txtORIG2; 


#concatenate
a=$(tail -n +2 newSUM.txt);
cat SUM.b1.RgeoTgen.txtORIG2 <(echo "$a") > SUM.b1.RgeoTgen.txt;
mv SUM.b1.RgeoTgen.txt ../compare;


#find rgeo missing loci by finding non-consecutive locus id number, will not find missing loci at end of file: 

l=$(cut -d$'\t' -f2 SUM.b1.RgeoTgen.txt | cut -d'.' -f3 | sort -nu); conseq "$l" > torepeatrgeo.txt;
l=$(cut -d$'\t' -f2 SUM.b1.RpcoTgen.txt | cut -d'.' -f3 | sort -nu); conseq "$l" > torepeatrpco.txt;
	
	ll=$(cut -d$'\t' -f2 tmp1.geo.txt | cut -d'.' -f3 | uniq)
	
	none









#compute +STATS* files

#additional tools
#determine loci with fewer than expected lines in the SUM file, there should be (#pops - 1)
#lines per locus:

cut -d$'\t' -f2 SUM.b1.RgeoTgen.txt | sort | uniq -c | grep -v "22 " #assumes Sorghum with #pops=23


