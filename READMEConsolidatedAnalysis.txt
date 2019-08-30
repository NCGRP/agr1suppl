#!/bin/bash

#	This script contains the consolidated processing of M+ results. The taxon, reference
#	model (Rgen,Rgeo,Rpco,Rlig,Rgo,Rlsa), and range of block lengths can be set.
#	cd to directory containing SUM* files for Rgen/Rgeo/Rpco/Rlig/Rgo/Rlsa.
#	You must have the files *LocusWeights.txt and +*SNPtoBlockMap.txt in the same folder.
#	Execute as a bash script like ./READMEConsolidatedAnalysis.txt 2 50.

#parameters
t="Sorghum"
v2p="Sorghum"; export v2p;
ts="Sor";
elist="gen geo pco lig lsa regulat"; #model for reference loci: "gen geo pco lig lsa regulat"

de=""; #estimator of allelic diversity at target loci, use nothing ("") for m+ estimator, use "ALLELECNT" for allele counts
export de;
mx=""; #optimality criterion used during M+ search. "L" for allele count, nothing ("") for m+ criterion
export mx;

bmin=$1; export bmin; #minimum blocklength
bmax=$2; export bmax; #maximum blocklength
if [[ $bmin == $bmax ]]; then dlc="off"; else dlc="on"; fi; #turn off "do line count" check if only one block length is considered

bstep=2; #step rate
nnode=9; #number of cluster nodes
blist=$(seq $bmin $bstep $bmax); #blocklengthrange

p=$(pwd)"/"; export p;

weightsfile="$p""SorRedoLocusWeights.txt";
tname="$p""+""$ts""RedoSNPtoBlockMap.txt"; export tname;



#	Calculate mean M+ enrichment across core sizes, for each locus, for each blocklength. This uses the
#	SUMfiles as input and produces the +STATS files.  Do this on the cluster.
#	
#**********BASHPARALLEL-FASTER-RUNINSCREEN**********
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
        sumfileIN="$p""SUM.b"$b".R"$e"Tgen.txt"$mx;
        ofile="$p""+STATS.R"$e"Tgen.loc"$l".b"$b".txt"$mx$de; #temporary outfile for each parallel process

        if [[ $de = "ALLELECNT" ]];
          then boo=$(grep -F $l "$sumfileIN" | awk '{ print $25 - $22 }' | tr "\n" " "); #extract all core sizes for this locus, calculate target enrichment using allele count
          else boo=$(grep -F $l "$sumfileIN" | awk '{ print $13 - $10 }' | tr "\n" " "); #extract all core sizes for this locus, calculate target enrichment using m+ diversity criterion
        fi;

        st=$(mysd "$boo"); #calculate sum, mean, sd, and n using function mysd()
        echo "$b"$'\t'"$l"$'\t'"$st" > "$ofile"; #write stats to output file
}
export -f main;

#parallelize on locus, within loop on blocklength
for e in $elist;
  do echo $e;
    export e;
    for b in $blist;
      do echo "b=$b";
        export b;
        sumfileIN="$p""SUM.b"$b".R"$e"Tgen.txt"$mx;
        uloc=$(cut -d$'\t' -f2 "$sumfileIN" | tail -n +2 | uniq); #find unique loci from all loci in column 2
        numloc=$(echo "$uloc" | wc -l);
        
        #parallel step
        #sends $gnuN records to each node via parallel step #1, this is piped to parallel step #2, which starts 96 jobs per node from the instruction set it has received
        echo "processing $sumfileIN...";
        gnuN=$(echo "$numloc / $nnode" | bc); #number of lines of input file to send to gnu parallel per node
        echo "$uloc" | parallel --sshloginfile ~/machinesALLCPUSPEC --jobs 1 --env main --env mysd --env mx --env de --env e --env b --env p --pipe -N"$gnuN" /home/reevesp/bin/parallel -j96 --env main --env mysd --env mx --env de --env e --env b --env p main;
        #run on head node only
        #echo "$uloc" | parallel --jobs 1 --env main --env mysd --env mx --env de --env e --env b --env p --pipe -N"$gnuN" /home/reevesp/bin/parallel -j96 --env main --env mysd --env mx --env de --env e --env b --env p main;
        #run on ceres, one node
        #echo "$uloc" | parallel --jobs 1 --env main --env mysd --env mx --env de --env e --env b --env p --pipe -N"$gnuN" parallel -j96 --env main --env mysd --env mx --env de --env e --env b --env p main;
        #run on Mac
        #echo "$uloc" | parallel --env main --env mysd --env mx --env de --env e --env b --env p main;

        #consolidate results for each locus into one file
        ofile="+STATS.R"$e"Tgen.b"$b".txt"$mx$de; #temporary outfile
        > "$p""$ofile";
        > "$p""$ofile"TMP;
       
        #fast cat
        find "$p" -maxdepth 1 -name "+STATS.R"$e"Tgen.loc*.b"$b".txt"$mx$de -print0 | xargs -0 cat -- >> "$p""$ofile"TMP; #quick cat all relevant files
        awk -F[$'\t'.] '{print $0"\t"$4}' "$p""$ofile"TMP | sort -t$'\t' -n -k1,1 -n -k7,7 | cut -d$'\t' -f1-6 > "$p""$ofile"; #use awk to get the genome index using two delimiters, it is added as a new, last column, then sort performed, first on blocklength, then on index
        rm "$p""$ofile"TMP;

        #test whether all were concatenated
        s1=$(find "$p" -maxdepth 1 -type f | grep -c "\.loc"); #count the number of files with "loc" in their name using grep. can't wc -l because too many arguments
        s2=$(wc -l "$p""$ofile" |  awk '{print $1}');
         if [ $s1 == $s2 ]; 
           then find "$p" -maxdepth 1 -name "*.loc*" -delete; #remove files when there are too many for rm
           else echo "s1="$s1", s2="$s2". Consolidated file "$ofile" no good. Quitting..." > err.txt;
           exit 1;
         fi;
      done; #$b
      
    #concatenate files for this Rgeo/Rpco
    statsfile="$p""+STATS."$t".R"$e"Tgen."$bmin"."$bmax".txt"$mx$de;
    echo blocklength$'\t'locus$'\t'sum$'\t'mean$'\t'sd$'\t'n > "$statsfile";
    >"$statsfile"TMP;
    #quick cat and sort (some out files from ceres are not sorted correctly by locus)
    cat "$p""+STATS.R"$e"Tgen.b"*".txt"$mx$de >> "$statsfile"TMP; #cat together the whole mess
    sortkey1=$(cut -d$'\t' -f2 "$statsfile"TMP | cut -d'.' -f3); #get a column containing only the locus index number, field 3 of column 2
    paste -d$'\t' "$statsfile"TMP <(echo "$sortkey1") | sort -t$'\t' -n -k1,1 -n -k7,7 | cut -d$'\t' -f1-6 > "$statsfile"TMP2; #add sort key to last column, sort by blocklength then locus index, remove last column
    cat "$statsfile"TMP2 >> "$statsfile";
    rm "$statsfile"TMP;
    rm "$statsfile"TMP2;

    #verify line counts are equal when $bmin != $bmax
    if [[ $dlc == "on" ]]; then
      #test whether all individual +STATS files have been added to the final concatenated file using number of lines
      s1=$(wc -l +STATS.R"$e"Tgen.b*.txt"$mx$de" | grep total | awk '{print $1}');
      s2=$(wc -l "$statsfile" | awk '{print $1}');
      s2=$(( $s2 - 1 )); #decount header
      if [ $s1 == $s2 ]; 
        then rm +STATS.R"$e"Tgen.b*.txt"$mx$de"; #remove files if number of lines in concat file is the same as the sum of all input files
        else echo "the number of lines ain't the same. something is 'crewed.";
      fi;
    else rm +STATS.R"$e"Tgen.b*.txt"$mx$de";
    fi;

  done; #$e
#**********BASHEND**********
#
#
#
#
#
#	add the proportion of NG,SS,NS sites from SorLocusWeights.txt to the +STATS summary file,
#	generating the +v2STATS file.  Requires the SorLocusWeights.txt file as input.  
#
#**********BASH**********
#use md5 on mac md5sum on linux
#modify weightsfile to include only the blocklengths analyzed
wtfhead=$(head -1 "$weightsfile");
wtfbod=$(awk -v blist="$blist" 'BEGIN{split(blist,t); for (i in t) vals[t[i]]} ($2 in vals)' "$weightsfile");
#below is slower, but more human-readable alternative
#wtfbod=$(for i in $blist;
#  do awk -F"\t" -v i=$i '$2 == i {print}' "$weightsfile"; #accumulate each blocklength needed
#  done;)
wtf=$(echo "$wtfhead"; echo "$wtfbod");

#echo "$wtf" > wtf.txt; 

#wtf=$(awk -F"\t" '$2 != "1" {print}' "$weightsfile"); #exclude blocklength=1 from the weightsfile for use with blocklength 2-200
for e in $elist;
  do statsfile="$p""+STATS."$t".R"$e"Tgen."$bmin"."$bmax".txt"$mx$de;
    s=$(cut -d$'\t' -f2 "$statsfile" | tail -n +2 | md5sum); #get md5 of locus ID column from the +STATS file
    w1=$(echo "$wtf" | cut -d$'\t' -f4 | sed 's/_/./g'); #get the locus ID less blocklength index column, replace _ with . for md5 check
    w2=$(echo "$wtf" | cut -d$'\t' -f3); #get the blocklength index column
    w=$(paste -d'.' <(echo "$w1") <(echo "$w2") | tail -n +2 | md5sum); #get md5 of locus ID column from the Weights file
    if [ "$s" == "$w" ]; #verify that files are in same order using a checksum of the unique locus ID column
      then boo=$(echo "$wtf" | sed 's/blocklength/blocklength2/g' | sed 's/locus/locus2/g'); #modify the header in the Weights file so that there are no redundant column names after combining
        paste -d$'\t' "$statsfile" <(echo "$boo") > "+v2STATS.Sorghum.R"$e"Tgen."$bmin"."$bmax".txt"$mx$de; #paste the statsfile and the weights file together
      else echo "md5s do not match, $e, aborting.";
        kill -INT $$; #terminate the script, return to the shell
    fi;
  done;
#**********BASHEND**********
#
#
#
#
#
#	Using the file +*SNPtoBlockMap.txt, calculate the sum of M+ enrichment value across blocklengths,
#	for each SNP, i.e. for each position in the genome.  These data will be extracted from the 
#	+v2STATS* files in several steps.  The first makes a file containing a table with the max 
#	enrichment value at each site for each block, *EnrichAcrossBlocks*. Step 1:
#
#**********BASHPARALLEL**********
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
        v2name="$p""+v2STATS."$v2p".R"$e"Tgen."$bmin"."$bmax".txt"$mx$de;
        
        echo "e=$e  b=$b";
        bout="$p""rr.R"$e"Tgen.b""$b"".tmp"$mx$de; #temporary outfile

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

#first parallel, distribute across envgeopco models and blocklengths
parallel --sshloginfile ~/machinesALLCPUSPEC --jobs 1 --env myp --env myqq --env mx --env de --env p --env v2p --env tname --env bmin --env bmax myp ::: $elist ::: $blist;
#parallel --jobs 40 --env myp --env mx --env de --env p --env v2p --env tname myp ::: $elist ::: $blist; #parallel command for head node only
#parallel --jobs 25 --env myp --env mx --env de --env p --env v2p --env tname myp ::: $elist ::: $blist; #parallel command for ceres
#parallel --env myp --env mx --env de --env p --env v2p --env tname myp ::: $elist ::: $blist; #parallel command for Mac

  
#reconstitute
echo "reconstituting..."
for e in $elist;
  do echo "e=$e";
    comm=""; #create a list of tmp files to paste together
    for b in $blist;
      do comm+="rr.R"$e"Tgen.b$b.tmp"$mx$de" ";
      done;
    pos=$(cut -d$'\t' -f4 "$tname" | sed 's/'^.*_'//g' | sed 's/locus/pos/g');
    h=$(cut -d$'\t' -f1-6 "$tname" | tr '\t' ' ');
    h2=$(paste -d' ' <(echo "$h") <(echo "$pos"));
    paste -d' ' <(echo "$h2") $comm > "+""$ts""EnrichAcrossBlocks.R"$e"Tgen."$bmin"."$bmax".txt"$mx$de;
    rm $comm; #clean up
  done;#e
  
#**********BASHEND**********
#	
#
#
#
#
#	Process each line (SNP position) of *EnrichAcrossBlocks*.txt to get sum across all blocklengths
#	of cumulative enrichment across all core sizes of max enrichment value from 10 reps of M+.
#	Produces a file called SUM*EnrichAcrossBlocks... Step 2:
#
#**********BASH**********
mysd2() {
       #calculate sum, mean, sd, n, mean,using Welford's algorithm and AWK
       #here, sum is the cumulative enrichment across all blocklengths, mean is the average
       #supply space delimited string of values as argument like: mysd2 $foo, where foo="1 2 3 4.5"
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
         #deal with situations where there is only one value
         if ( NF == 1 )
           { print sum " " sum/NF " -- " NF; }
         if ( NF != 1 )
           { var = S/(NF - 1);
           print sum " " sum/NF " " sqrt(var) " " NF; }
       }' <<< "$1"; # the first argument when calling mysd2 
}
export -f mysd2;

myp2() {
    l=$1

    pp=$(echo "$l" | cut -d' ' -f8-); #cut out the portion of the line with max enrichment values

    #skip loci that have missing blocks
    if [[ $pp != *"--"* ]]; then
      m=$(mysd2 "$pp"); #calculate sum, mean, sd, n
      chr=$(echo "$l" | cut -d' ' -f5);
      pos=$(echo "$l" | cut -d' ' -f7);
      i=$(echo "$l" | cut -d' ' -f1); #get genome index
      loc=$chr"."$pos"."$i;
      #echo "snp=$i";

      echo $i" "$chr" "$pos" "$loc" "$m > "$p"$i".tmp"$mx$de; #add new values to growing string, write to tmp file
    fi;
}
export -f myp2;

#parallel
for e in $elist;
do echo $e;
  #parallel command
  #for cluster
  tail -n +2 "+""$ts""EnrichAcrossBlocks.R"$e"Tgen."$bmin"."$bmax".txt"$mx$de | parallel  --sshloginfile ~/machinesALLCPUSPEC --jobs 1 --env myp2 --env mysd2 --env mx --env de --env p --pipe --round-robin /home/reevesp/bin/parallel -j96 --env myp2 --env mysd2 --env mx --env de --env p myp2;
  #for head node only
  #tail -n +2 "+""$ts""EnrichAcrossBlocks.R"$e"Tgen."$bmin"."$bmax".txt"$mx$de | parallel --jobs 1 --env myp2 --env mysd2 --env mx --env de --env p --pipe --round-robin /home/reevesp/bin/parallel -j96  --env myp2 --env mysd2 --env mx --env de --env p myp2;
  #for ceres N=26
  #tail -n +2 "+""$ts""EnrichAcrossBlocks.R"$e"Tgen."$bmin"."$bmax".txt"$mx$de | parallel  --jobs 26 --env myp2 --env mysd2 --env mx --env de --env p myp2;
  #for mac
  #tail -n +2 "+""$ts""EnrichAcrossBlocks.R"$e"Tgen."$bmin"."$bmax".txt"$mx$de | parallel  --env myp2 --env mysd2 --env mx --env de --env p myp2;
  
  #assemble results
  echo "concatenating...";
  ofile="+""$ts""SUMEnrichAcrossBlocks.R"$e"Tgen."$bmin"."$bmax".txt"$mx$de;
  o="genomeindex chr pos loc sum mean sd n";
  echo "$o" > "$ofile";
  
  find "$p"  -maxdepth 1 -name "*.tmp"$mx$de -print0 | xargs -0 cat | sort -n -t' ' -k1 >> "$ofile";
  
  #clean up
  echo "clean up...";
  find "$p" -name "*.tmp"$mx$de -print0 | xargs -0 rm;
done;
#**********ENDBASH**********
