	Starting 2/7/2019, in folder /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/monolignolAnalysis
	Use list of monolignol loci provided by Hannah Tetreault in LigninRegionsToPipe.txt
	as a single set of reference loci, calculate diversity retained at all genomic loci (target),
	including the monolignol references.
	
	This analysis will require appending the output of the haplotypista run that takes LigninRegionsToPipe.txt
	as input and prints out M+ files for blocklengths 2 to 50 by 2.
	
	
{  This first approach doesn't work, skip it...

	#in the folder monolignolDeployment (the destination folder, $depath), have a set of folders labeled
	#2,4,6,8,...50, each with the unmodified .dat and .var files output by various applescripts
	#(see README Final Analysis2.rev expt 2 with hclust.txt), and a copy of m+1 for ceres.
	#For current purposes you can just copy the .dat and .var and m+1 from the Final Analysis 3 deployment:
	mkdir monolignolDeployment;
	cd monolignolDeployment;
	seq 2 2 50 | parallel mkdir {};
	seq 2 2 50 | parallel cp "/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Sorghum/Experiment2/+revision\ with\ hclust/m+deployment/hclust/+forRgenruns/rarefaction/RgenDeployment/"{}"/*[vdm][a+][rt1]" {};

	#Use haplotypista to generate M+ data sets for blocklength 2 to 50 by 2 for the set of monolignol
	#genes in file LigninRegionsToPipe.txt

	#make a pop id file to use for M+ output switch -v of haplotypista
	cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/monolignolAnalysis/monolignolHaplotypista;
	cat SorghumHaplotypistaInputMOD3.txt | tail -n +4 | cut -d. -f1 > Sorpopid.txt;

	#Pipe the lignin gene regions from /Volumes/Pbot/NP02-PREEVES/Edrive/Pat\'s\ Files/Pat\'s/Research\ 2017/Hannah\ Tetreault/MonolignolGenes_locations_PR.xlsx 
	#to haplotypista with the following command line:
	cp /Volumes/attica/Desktop/Attica-Pats\ Folder/Scripts/Haplotype\ analysis/haplotypista/haplotypista .;
	cat LigninRegionsToPipe.txt | ./haplotypista -i SorghumHaplotypistaInputMOD3.txt -o SorLig -l SorLiglog.txt -b 1 50 -m ? -p 2 -v Sorpopid.txt;

	#calculate the number of monolignol loci remaining for each blocklength
	for i in $(seq 2 2 50);
	  do echo -n "$i ";
		grep -a " 2 1 0 1 5" "SorLig.b"$i".var" | sort | uniq | wc -l;
	  done;
			Result:
			2      420
			4      176
			6      104
			8       63
			10       43
			12       31
			14       22
			16       17
			18        9
			20        5
			22        4
			24        4
			26        2
			28        2
			30        2
			32        2
			34        1
			36        0
			38        0
			40        0
			42        0
			44        0
			46        0
			48        0
			50        0

	#This is weird.  Because haplotypista deletes any data that doesn't fall within the regions of interest at the outset,
	you end up with many fewer remaining loci than if you look to see if a haplotype block created from the full data set
	contains the region of interest.  So, maybe the old way of extracting the reference loci from var files created by
	haplotypista for the whole data set is the best way to get the reference set.  But, they should be copied out
	as a separate data set and appended to the end of the var file. That way, they can be used for maximization collectively
	but used as target loci separately.  Start over...
}

#setup
cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/monolignolAnalysis/monolignolDeployment
seq 2 2 50 | parallel mkdir {};
seq 2 2 50 | parallel cp "/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Sorghum/Experiment2/+revision\ with\ hclust/m+deployment/hclust/+forRgenruns/rarefaction/RgenDeploymentPaper3/"{}"/*[vdm][a+][rt1]" {};

t="Sor";
export t;

#make working copies of the original var files
for b in $(seq 2 2 50);
  do cp ./$b/"$t"SNP.b"$b".RgenTenv.var ./$b/"$t"SNP.b"$b".RgenTenv.varORIG;
  done;

#turn all mstrat variable codings to ignore, 10115
for b in $(seq 2 2 50);
  do echo "$b";
    inf=./$b/"$t"SNP.b"$b".RgenTenv.var; #define input file
    cat "$inf" | tr "\t" " " | sed 's/2 1 0 1 5/1 0 1 1 5/g' | sed 's/2 0 1 1 5/1 0 1 1 5/g' > tmp.txt;
    mv tmp.txt "$inf";
  done;


#identify snps that belong to monolignol regions, change those to 21015
#There is a much faster version of this in READMEGOAnalysis.txt
####BEGIN BASH####
mypp() {
     b=$1; #blocklength
     tmpfile="$wd/monolignolDeployment/$b/$t"SNP.b"$b".RgenTenv.var;
     rangefile="$wd/LigninRegionsToPipe.txt";  #define file containing genomic ranges to treat as reference loci
     
     #define coding
     referencea='2 1 0 1 5';
     target='2 0 1 1 5';
     ignore='1 0 1 1 5';

     #remove regions of interest from the list of potential target loci
#     tttf=$(cat "$tmpfile"); #load var file into memory
#     vin=$(cat "$rangefile"); #genomic ranges to set as reference loci into memory
     locusstr=$(grep '\.' "$tmpfile" | cut -d$'\t' -f1 | uniq | cut -d' ' -f1); #read all genetic locus names in the var file into a space delimited string, works because only genetic loci have . character in them
     locusstr2="";
     cchrold=0;

     for i in $locusstr;
       do lp=$(echo "$i" | cut -d. -f2); #the locus position of potential target
         cchr=$(echo "$i" | cut -d. -f1); #determine chromosome for this locus
         if (( $cchr != $cchrold )); #only waste time with the grep enclosed if the chromosome has changed
         then
           vintouse=$(grep "^$cchr"'\.' "$rangefile"); #only consider ranges that come from the same chromosome
#           vintouse=$(echo "$vin" | grep "^$cchr"'\.'); #only consider ranges that come from the same chromosome
         fi;
         zz="targ"; #assume it is a target locus
         
         #echo "$i";
         
         for j in $vintouse;
           do rlow=$(echo "$j" | cut -d: -f1 | cut -d. -f2); #the left end of the genomic range to exclude
             rhigh=$(echo "$j" | cut -d: -f2 | cut -d. -f2); #the right end of the genomic range to exclude
             
             if (( ($rlow <= $lp) && ($lp <= $rhigh) )); #test whether current locus position is a reference locus
             then
               sed -i "s/$i $ignore/$i $referencea/" "$tmpfile"; #it is a reference locus
#               tttf=$(echo "$tttf" | sed "s/$i $ignore/$i $referencea/"); #it is a reference locus
               zz="ref"; #reset to reference
               break;
             fi;
           done;
#         if [[ $zz != "ref" ]]; then locusstr2+=$(echo "$i "); fi; #add locus to list of target loci if it is not a reference locus
         cchrold=$cchr;
       done;
       
#       echo "$tttf" > "$tmpfile"1;
       cp "$tmpfile" "$tmpfile"1;
       echo "b$b";
}
export -f mypp;

t="Sor";
export t;
wd=$(pwd);
export wd

seq 2 2 50 | parallel --env t --env wd --env mypp mypp;

####END BASH####
	Takes ~12 hours on compute-0-9. There is a much faster version of this in READMEGOAnalysis.txt (like 5 minutes)






#calculate the number of monolignol loci remaining for each blocklength
for i in $(seq 2 2 50);
  do echo -n "$i ";
    grep " 2 1 0 1 5" "$i/SorSNP.b"$i".RgenTenv.var1" | sort | uniq | wc -l;
  done;
		Result:
		2 429
		4 187
		6 109
		8 74
		10 57
		12 37
		14 36
		16 28
		18 34
		20 26
		22 17
		24 17
		26 13
		28 12
		30 12
		32 15
		34 13
		36 10
		38 10
		40 9
		42 10
		44 10
		46 10
		48 2
		50 7

#convert var file names and clean up temporary var files
for i in $(seq 2 2 50);
  do echo "$i";
    cp "$i/SorSNP.b"$i".RgenTenv.var1" "$i/SorSNP.b"$i".Rlig.var";
    rm "$i/SorSNP.b"$i".RgenTenv.var" "$i/SorSNP.b"$i".RgenTenv.var1" "$i/SorSNP.b"$i".RgenTenv.varORIG";
  done;

#make single space delimiter in dat file, set up as new file name, envgeopcolig, which will
#get the reference loci duplicated onto its end
for i in $(seq 2 2 50);
  do echo "$i";
    sed 's/  */ /g' "$i/SorSNP.b"$i".envgeopco.dat" > "$i/SorSNP.b"$i".envgeopcolig.dat";
  done;

#copy out all "2 1 0 1 5" (reference) loci and paste them at the end of the var file,
#then find them in dat file and replicate them at the end as well
for i in $(seq 2 2 50);
  do echo "$i";
    a=$(grep " 2 1 0 1 5" "$i/SorSNP.b"$i".Rlig.var" | sed 's/^/lig/g'); #get all var file lines that are 21015 reference
    sed 's/ 2 1 0 1 5/ 1 0 1 1 5/g' "$i/SorSNP.b"$i".Rlig.var" > "$i/SorSNP.b"$i".Rlig2.var"; #make a copy of original var file, setting existing 21015 reference to 10115 ignore
    echo "$a" >> "$i/SorSNP.b"$i".Rlig2.var"; # add set of 21015 reference loci to end of file

    blist=$(grep -n ' 2 1 0 1 5' "$i/SorSNP.b"$i".Rlig.var" | cut -d: -f1); #get the line numbers of reference loci from var file, these are column numbers in dat file
    allc="";
    for b in $blist;
      do echo col$b;
        c=$(cut -d' ' -f$b "$i/SorSNP.b"$i".envgeopcolig.dat");
        allc=$(paste -d' ' <(echo "$allc") <(echo "$c")); #get all columns containing reference data pasted together
      done;
      allc=$(echo "$allc" | sed 's/^ //g'); #remove line leading space
    paste -d' ' "$i/SorSNP.b"$i".envgeopcolig.dat" <(echo "$allc") > "$i/SorSNP.b"$i".envgeopcolig2.dat"
  done;
#Takes about 10 minutes on compute-0-9

#clean up by making a new folder, monolignolDeployment above the work directory that contains just the
#files needed for execution, relabeled following conventions (e.g. *Rlig.var not *Rlig2.var
cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/monolignolAnalysis;
cp -r /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/monolignolAnalysis/work/monolignolDeployment .;
cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/monolignolAnalysis/monolignolDeployment;
for i in $(seq 2 2 50);
  do echo "$i";
    mv "$i"/SorSNP.b"$i".Rlig2.var "$i"/SorSNP.b"$i".Rlig.var;
    mv "$i"/SorSNP.b"$i".envgeopcolig2.dat "$i"/SorSNP.b"$i".envgeopco.dat;
    rm "$i"/SorSNP.b"$i".envgeopcolig.dat;
    rm "$i"/*.sh;
  done;






#Run 9dSlicedBatchGeneratorRaref.sh to create slurm submit files. It will probably not be necessary to slice
#the data to fit in a 48 hour run time, since there are so many fewer loci to optimize,
#but let's just stick with that framework for consistency.
#For instructions on how to run 9dSlicedBatchGeneratorRaref.sh, see its header.
#In folder "/Volumes/J22/M+ Study/Analysis/Final analysis for 3rd paper/ForSlicedSearch" issue:

cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/ForSlicedSearch;
./9dSlicedBatchGeneratorRaref.sh;

#remove null bytes for some reason polluting all the var files
cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/monolignolAnalysis/monolignolDeployment
for i in $(seq 2 2 50);
  do echo $i;
    gsed -i 's/\x0//g' $i/*.var; #use gsed on Mac only, sed otherwise
  done;

#transfer to ceres and make some mods:
rsync -aP /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/monolignolAnalysis/monolignolDeployment/* pat.reeves@login.scinet.science:"/home/pat.reeves/mplusrunsSor/monolignol";
cd /home/pat.reeves/mplusrunsSor/monolignol/+SorSlicedTools;
sed -i 's/short/short,medium,long,mem,mem768/' slurmArraySorb*Sliced.sh; #add more node possibilities to all slurm files
sed -i 's/--mem=0/--mem=124G/' slurmArraySorb*Sliced.sh;
sed -i 's:gcc/64/::g' slurmArraySorb*Sliced.sh;
sed -i 's:_libs::g' slurmArraySorb*Sliced.sh;
cd /home/pat.reeves/mplusrunsSor/monolignol;
for i in $(seq 2 2 50);
  do echo $i;
    cp /home/pat.reeves/m+1/m+1 $i;
  done;
cd /home/pat.reeves/mplusrunsSor/monolignol/+SorSlicedTools;
chgrp -Rv proj-patellifolia ~/mplusrunsSor/monolignol; #chgrp for all folders and files so they can be written to by my or ann's account
chmod -Rv g+w ~/mplusrunsSor/monolignol; #allow write on all objects for all members of group
cd /home/pat.reeves/mplusrunsSor/monolignol/+SorSlicedTools;
grep ": sbatch" *Sliced.sh | cut -d: -f3 | sort -t'b' -k3,3n; #print all the sbatch commands


#or move to blip, from blip issue:
mkdir /home/reevesp/mplusrunsSor/monolignol;
rsync -aP shrub@10.177.9.250:"/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/monolignolAnalysis/monolignolDeployment/*" /home/reevesp/mplusrunsSor/monolignol;

#convert ceres SLURM files to work on blip, add a set that maximize on raw allele account ("slurmArray*.L.sh")
cd /home/reevesp/mplusrunsSor/monolignol/+SorSlicedTools;
for i in $(seq 2 2 50);
  do echo $i;
    sed -i 's/-p short/-p CLUSTER/' slurmArraySorb"$i"Sliced.sh;
    sed -i 's/-n 40/-n 24/' slurmArraySorb"$i"Sliced.sh;
    sed -i 's/--mem=0/--mem=44G/' slurmArraySorb"$i"Sliced.sh;
    sed -i 's/-t 48:00:00/-t 0/' slurmArraySorb"$i"Sliced.sh;
    sed -i 's/^module/#module/g' slurmArraySorb"$i"Sliced.sh;
    #sed -i 's:#rm -fr $TMPDIR:rm -fr $TMPDIR/*:' slurmArraySorb"$i"Sliced.sh;
    cat slurmArraySorb"$i"Sliced.sh | sed 's/\.sh/\.L\.sh/g' | sed 's/m+1/m+1L/g' > slurmArraySorb"$i"Sliced.L.sh; #make a copy for executing with non-standardized optimality criterion (raw allele count)
  done;
  
#modify b2 runs on blip to not exceed memory or only use compute-0-9
sed -i 's/-n 24/-n 12/' slurmArraySorb2Sliced.*sh; #run on all nodes with 12 cores
#sed -i 's/-n 24/-n 40/' slurmArraySorb2Sliced.*sh; #run on compute-0-9
#sed -i 's/--mem=0/--mem=160G/' slurmArraySorb2Sliced.*sh; #run on compute-0-9

#get the proper m+1, and m+1L which optimizes on allele count
cd /home/reevesp/mplusrunsSor/monolignol;
for i in $(seq 2 2 50);
  do echo $i;
    cp /home/reevesp/Mplus/m+1 $i;
    cp /home/reevesp/Mplus/m+1L $i;
  done;

#make a copy of the 9d* execution scripts and modify for use with m+1L
cd /home/reevesp/mplusrunsSor/monolignol;
for i in $(seq 2 2 50);
  do echo $i;
    ls $i/9d* | parallel "sed 's/m+1/m+1L/g' {} > {}.tmp";
    ls $i/9d*.tmp | cut -d'.' -f1-2 | parallel "mv {}.sh.tmp {}.L.sh";
    chmod a+x $i/*.sh;
  done;
 
#run on blip, takes ~17 hrs 
cd /home/reevesp/mplusrunsSor/monolignol/+SorSlicedTools;
grep ": sbatch" *Sliced.sh | cut -d: -f3 | sort -t'b' -k3,3n; #print all the sbatch commands
grep ": sbatch" *Sliced.L.sh | cut -d: -f3 | sort -t'b' -k3,3n; #print all the sbatch commands

#in a screen, run the 9dSlicedSumMonitor.sh by cut and paste

#you can verify that the new SUM files differ from the old ones by the addition of the reference loci, which
#are now also treated as target.  This will allow the +v2STATS file to be calculated correctly later.

for i in $(seq 2 2 50);
  do echo $i;
    diff <(cut -d$'\t' -f2 ./48/SUM.b"$i".RligTgen.txt) <(cut -d$'\t' -f2 ../monolignolOLD/work/SUM.b"$i".RligTgen.txt)  | grep '<' | sort -u;
  done;




###ANALYSIS OF MONOLIGNOL LOCI###


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
elist="lig";
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
	Takes ~50 minutes on blip. 


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
elist="lig";
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
	Takes a minute or so.



#Plot summed enrichment by genomic position across blocklengths as a heatmap. Use Rscript
#SorGenomicGeography.r. Uses the +v2STATS files as input.



	Using the file +*SNPtoBlockMap.txt, calculate the sum of M+ enrichment value across blocklengths,
	for each SNP, i.e. for each position in the genome.  These data will be extracted from the 
	+v2STATS* files in several steps.  The first makes a file containing a table with the max 
	enrichment value at each site for each block, *EnrichAcrossBlocks*. Step 1:

**********BASHPARALLEL**********
myqq() {
  l=$(echo $1 | cut -d' ' -f2); #extract locus pos
  indx=$(echo $1 | cut -d' ' -f1); #extract the line index of the locus
  if [[ $l = "--." ]]; then
        r="--"; #test for empty enrichment value
      else
        r=$(grep -m 1 ^"$b"$'\t'"$l" "$ggfs" | cut -d$'\t' -f3);
      fi;
    echo "$indx $r";
}
export -f myqq;

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
        cc=$(echo "$cc" | nl -n ln -s' ' | sed 's/  */ /g'); #number lines of cc



        #search for locus names in $gg, retrieve enrichment value
        echo " b$b" > "$bout";

        #second parallel, across lines of SNPtoBlockMap
        ggfs=$(mktemp);
        echo "$gg" > "$ggfs";
        export ggfs;
        export b;
        (echo "$cc" | /home/reevesp/bin/parallel -j24 --pipe -N96 --round-robin --env b --env ggfs --env myqq /home/reevesp/bin/parallel -j2 --env b --env ggfs --env myqq myqq) >> "$bout";
        cp "$bout" "$ggfs"XXX; #copy outfile to tmp location
        sort -t' ' -k1,1n "$ggfs"XXX | cut -d' ' -f2 > "$bout"; #sort by line index, then remove line index
        rm "$ggfs";
        rm "$ggfs"XXX;
}
export -f myp;

#parameters
v2p="Sorghum"; export v2p;
t="Sor";
minbl=2; #1,2
maxbl=50; #50,200
blstep=2; #step
blist=$(seq $minbl $blstep $maxbl); #for even numbered runs 2-50
elist="lig";

p=$(pwd)"/"; export p;
tname="$p""+""$t""SNPtoBlockMap.txt"; export tname;

#first parallel, distribute across envgeopco models and blocklengths
parallel --sshloginfile ~/machines --jobs 1 --env myp --env myqq --env p --env v2p --env tname myp ::: $elist ::: $blist;
#parallel --jobs 40 --env myp --env p --env v2p --env tname myp ::: $elist ::: $blist; #parallel command for head node only
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
    #rm $comm; #clean up
  done;#e
  
**********BASHEND**********
	Since redesigning parallel steps, takes 25 minutes on 9 nodes.
	
	***WATCH FOR THIS ERROR: Signal SIGCHLD received, but no signal handler set.***
	If you get it, repeat the affected block length and reconstitute.
	




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
elist="lig";
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

**********BEGINR**********
#install.packages("fitdistrplus")
#install.packages("plotrix")
options(error = recover)
rm(list=ls()) 
setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/monolignolAnalysis/monolignolAnalysis/+results/rarefaction/+analysis")
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
		if ( zg == "pdf" | zg == "svg" ) { zgg(file=paste(t,"M+GenomeScan",pp, ".", zg, sep="")) }
		else { 
		par(mar=c(1,1,1,1))
		zgg(file=paste(t,"M+GenomeScan",pp, ".", zg, sep=""),width=2000,height=2000,res=300) } #call zgg function with relevant parameters
		
		par(mfrow=c(2,1))

			i="lig"
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
chr=$(cut -d' ' -f2 "+"$i"SUMEnrichAcrossBlocks.RgenTgen.txt" | tail -n +2 | sed 's/^/chromosome_/g');
pos=$(cut -d' ' -f3 "+"$i"SUMEnrichAcrossBlocks.RgenTgen.txt" | tail -n +2);
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
setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Sorghum/Experiment2/+revision\ with\ hclust/+results/hclust/+SorRgenRuns/rarefaction/+analysis/+SorRgenGOWork/pie\ charts")

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




 