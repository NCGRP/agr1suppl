Enriching collections at gene regions associated with GO terms
  Project work is saved at: /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/sagoAnalysis
  Project was computed using Linux on blip:compute-0-9



###Generate a corpus containing name and def lines for all GO terms###

#process go-basic.obo to generate complete GO corpus, Mac takes forever, use Linux
oboin="go-basic.obo";
csplit -n6 "$oboin" /^\\[Term\\]$/ {*}; #splits file on "[Term]", writes output files like xx000001
rm xx000000; #first file contains obo header stuff, remove
a=$(ls xx* | wc -l | awk '{print $1}'); #get total number of files generated
f="xx"$(printf "%.6d" "$a"); #create the file name for the last file, %.xd where x specifies number of characters
sed -n '/^\[Typedef\]$/q;p' "$f" > tmp.txt; #remove extra [Typedef] blocks at end of last file
mv tmp.txt "$f";

oftmp=tmp.corpus.txt;
of="$oboin".corpus.txt;
>"$oftmp";
>"$of";
for i in $(seq 1 $a); 
  do f="xx"$(printf "%.6d" "$i");
    echo "$f";
    id=$(grep ^"id: " "$f" | sed 's/id: //');
    name=$(grep ^"name: " "$f" | sed 's/name: //');
    namespace=$(grep ^"namespace: " "$f" | sed 's/namespace: //');
    def=$(grep ^"def: " "$f" | gsed 's/ \[.*\]//g' | gsed 's/^def: //g' | sed 's/"//g');
    echo "$id"@"$name"@"$namespace"@"$def" >> "$oftmp"; #write useful info for each [Term] to output, delimit using @, which doesn't exist in go-basic.obo or goslim_plant.obo
  done;

grep -v OBSOLETE "$oftmp" > "$of"; #remove terms marked as obsolete

awk -F@ '{print NF}' "$of" | sort -u; #verify that all lines have 4 columns

rm xx*; #clean up
rm "$oftmp";

###END Generate a corpus containing name and def lines of all GO terms###




###Generate species specific corpus of GO terms###

spp="Sor"; #"At", "Pop", "Sor"
cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/sagoAnalysis;
b=$(cut -d$'\t' -f5 "$spp"RawGAF.txt | sort -u); #extract all unique GO terms from the RawGAF file (At,6928; Pop,2841; Sor,2702)

#extract all lines in go-basic.obo.corpus.txt that match the species specific GO terms ($b)
oboout="$spp".GOcorpus.txt;
>"$oboout";
for i in $b;
  do echo "$i";
    grep ^"$i"@ go-basic.obo.corpus.txt >> "$oboout";
  done;

###END Generate species specific corpus of GO terms###


 


###Run latent semantic analysis###

#Notes: the TruncatedSVD(n_components= setting determines the point at which the ranking of
#GO terms will be stable. If n_components=2500 (i.e. > the 2230 terms in the corpus), then the first cycle will be stable and there
#is no need for progressive reduction of the corpus size to reach a stable ranking of
#GO terms with respect to the query.  So that is probably how it should be done and what follows
#immediately below marked by an ##EXPLORATORY## block is just that, exploratory.
#The final analysis is further down

##EXPLORATORY##
#Concatenate journal articles on Genetics, Climate Change, United States, Agriculture
genetics climate adaptation agriculture united states
cat Abberton2015.txt Lobell2014.txt Schlenker2009.txt > ConcatClimateChange.txt;
cat ConcatClimateChange.txt "SYR_AR5_FINAL_full_wcover.txt" "Climate Change and Agriculture in the United States_ Effects and.txt" > All.txt;
cat Lobell2014.txt Schlenker2009.txt "Climate Change and Agriculture in the United States_ Effects and.txt" > US.txt;

#setup
query1="ConcatClimateChange.txt";
query2="SYR_AR5_FINAL_full_wcover.txt";
query3="Climate Change and Agriculture in the United States_ Effects and.txt";
query4="All.txt";
query5="US.txt";
iconv -f utf-8 -t utf-8 -c "$query2" > tmp.tmp; #strip any non utf-8 characters from the file converted to plain text by adobe
mv tmp.tmp "$query2";
iconv -f utf-8 -t utf-8 -c "$query3" > tmp.tmp; #strip any non utf-8 characters from the file converted to plain text by adobe
mv tmp.tmp "$query3";

off="SorGOout";
kmax=5; #maximum number of random reps
nlc=30; #smallest number of lines allowed in output file (nlcutoff), also the number of lines included in summary statistics
trim=0.25; #proportion of output to retain during corpus reduction

j=4; #j indexes the query number
for q in "$query4" "$query5";
do i=1; #i indexes the reduction cycle
  #refresh the corpus and remove cellular component terms between queries
  corpus="Sor.GOcorpus.noCC.txt";
  grep -v "cellular_component" Sor.GOcorpus.txt > "$corpus"; 
  while(true);
    do k=1; #k indexes the number of reps before taking best
      sbest=1; #sbest is the mean cosine goodness of fit value for the rep, lower is better 
      while((k<=kmax));
        do true;
          ofg="$off.q$j.$i.r$k.txt"; #output file name for cycle-rep
          #sago.py corpusfile queryfile outputfile
          python sago.py "$corpus" "$q" "$ofg";

#          s=$(tail -n+2 "$ofg" | head -$nlc | cut -d$'\t' -f6 | awk -v nlc=$nlc '{sum+=$1}END{print sum/nlc}'); 
          w=$(wc -l "$ofg" | awk '{print $1}'); #number of lines in output file
          s=$(tail -n +2 "$ofg" | cut -d$'\t' -f6 | awk -v w=$w '{sum+=$1}END{print sum/w}');
          sha=$(tail -n +2 "$ofg" | head -$nlc | cut -d$'\t' -f1-5 | shasum); #calculate sha1sum to determine convergence on order 
          
          #test whether a new best rep has been found
          if (( $(bc <<< "$s <= $sbest") )); then
            sbest=$s;
            ofbest="$ofg";
          fi;
          
          echo "q$j reduce$i rep$k w$w $s $sbest $sha"; #print out the cycle-rep and the mean fit stat

          ((k++));
        done; #k random reps
        
      a=$(wc -l "$ofbest" | awk '{print $1}'); #number of lines in output
      nlx=$(echo "$a*$trim" | bc | xargs printf "%.*f\n" 0); #keep best $trim% of output (top half) to form next corpus, round using printf, last number (0) is precision
      if (( nlx < nlc )); then echo; break; fi; #if there are fewer than 30 lines in the output, exit loop
      #keep the best x% of the GO terms to form the corpus for the next round
      tail -n +2 "$ofbest" | head -$nlx | cut -d$'\t' -f2-5  | tr '\t' '@' > c$((i+1)).tmp;
      corpus="c$((i+1)).tmp";
      ((i++));
    done; #nlx > nlcutoff
  rm c*.tmp;
  ((j++));
done;
##END EXPLORATORY##


##FINAL LSA##

#Change sago.py so that TruncatedSVD(ncomponents) > number of GO terms in corpus.  Use 2500 in this case

#Concatenate journal articles on 'genetics climate adaptation agriculture united states'
cat ConcatClimateChange.txt "SYR_AR5_FINAL_full_wcover.txt" "Climate Change and Agriculture in the United States_ Effects and.txt" > All.txt;
cat Lobell2014.txt Schlenker2009.txt "Climate Change and Agriculture in the United States_ Effects and.txt" > US.txt;

#setup
query1="All.txt";
query2="US.txt";
iconv -f utf-8 -t utf-8 -c "$query2" > tmp.tmp; #strip any non utf-8 characters from the file converted to plain text by adobe
mv tmp.tmp "$query2";
iconv -f utf-8 -t utf-8 -c "$query3" > tmp.tmp; #strip any non utf-8 characters from the file converted to plain text by adobe
mv tmp.tmp "$query3";

j=1; #j indexes the query number
for q in "$query1" "$query2";
do echo "q$j";

  #refresh the corpus and remove cellular component terms between queries
  corpus="Sor.GOcorpus.noCC.txt";
  grep -v "cellular_component" Sor.GOcorpus.txt > "$corpus"; 

  ofg="SorGOout.q$j.txt"; #output file name for cycle-rep
  python sago.py "$corpus" "$q" "$ofg"; #sago.py corpusfile queryfile outputfile

  ((j++));
done;

##END FINAL LSA##


###END Run latent semantic analysis###








#The final LSA produces two files, SorGOout.q[12].txt. q1 has data from all input files, q2 from US only.  Use
#US only going forward (SorGOout.q2.txt).

#print out the top 20 most semantically related GO terms
head -21 SorGOout.q2.txt | cut -d$'\t' -f2,3;

			GO:0010114	response to red light	0.7674808314389545
			GO:0010018	far-red light signaling pathway	0.779099077556022
			GO:0009408	response to heat	0.8336608045318876
			GO:0015979	photosynthesis	0.8417371556147005
			GO:0009409	response to cold	0.8418193256108926
			GO:0031990	mRNA export from nucleus in response to heat stress	0.8439705081265143
			GO:0009640	photomorphogenesis	0.8446190468233951
			GO:0009651	response to salt stress	0.8474890582535507
			GO:0009555	pollen development	0.8488576602971163
			GO:0006979	response to oxidative stress	0.8500582705918619
			GO:0034599	cellular response to oxidative stress	0.8518270247163149
			GO:0009266	response to temperature stimulus	0.8577368206092867
			GO:0042538	hyperosmotic salinity response	0.8578576229149315
			GO:0009631	cold acclimation	0.8590780167526975
			GO:0071470	cellular response to osmotic stress	0.8591599987504168
			GO:0006338	chromatin remodeling	0.8598972189859403
			GO:0070483	detection of hypoxia	0.8616666540402179
			GO:0050826	response to freezing	0.8626040239137245
			GO:0048577	negative regulation of short-day photoperiodism, flowering	0.8626713873136767
			GO:0048364	root development	0.8629356897204199

#figure out how many genes we're talking about, some oddball nomenclature so remove those (SORBI,ndh,SB234M12) 
#scan across the top 50 to see if there is a better cutoff than 20, figure out if there are oddball names turning up, remove those in a second round
>SorLSAgenes.txt;
for i in $(head -51 SorGOout.q2.txt | tail -n +2 | cut -d$'\t' -f2 | tr "\n" " ");
  do grep "$i" SorRawGAF.txt | grep -v SORBI_ | grep -v ndh | grep -v SB234M12 | grep -v NIP2 | grep -v AMT2 | awk -F$'\t' '$3!=13{print $0}' >> SorLSAgenes.txt;
    echo "$i "$(wc -l SorLSAgenes.txt);
  done;

#20 seems reasonable, both in terms of the specific GO terms found (they don't all make some sense after that)
#and in terms of the asymptotic curve shape and in terms of the number of genes found.
#Using 50 instead because there are too few remaining loci at blocklength = 50 with first 20 terms
#There are 115 genetic loci implicated for the top 20 terms.
#There are 543 genetic loci implicated for the top 50 terms. (Actually 505 after excluding supernumerary "chromosomes" and duplicates)

#get the genome positions for these genes from Sb14genes.bed, save to a file called LSAGORegionsToPipeTop50.txt
glist=$(cut -d$'\t' -f3 SorLSAgenes.txt);
(for g in $glist;
  do echo -n "$g ";
    a=$(awk -F$'\t' -v g=$g '$4==g{print $0}' Sb14genes.bed);
    #test whether more than one line is returned
    if (( $(echo "$a" | wc -l) > 1 )); then 
      echo "more than one line returned for $g";
      break;
    fi;
    c=$(echo "$a" | cut -d$'\t' -f1 | cut -d_ -f2); #get the chromosome
    d=$(echo "$a" | cut -d$'\t' -f2); #get bp start
    e=$(echo "$a" | cut -d$'\t' -f3); #get bp end
    b=$(echo "$c"."$d":"$c"."$e"); #get string defining start and end, suitable for haplotypista input
    echo "$b";
  done;) | sort -t' ' -k2,2n > LSAGO.tmp;
grep -v 's' LSAGO.tmp | cut -d' ' -f2  | sort -u > LSAGORegionsToPipeTop50.txt; #remove all lines with lowercase 's', these are supernumerary chromosomes (i.e. not 1-10)
                                                                                #also remove all duplicates, which sometimes happens due to assignment to multiple GO terms
rm LSAGO.tmp;
#There are now 505 genomic segments. 
  



 ############BUILD DATASETS AND PROCESS############

#in the folder sagoAnalysis/LSAcliDeploymentTop[25]0 (the destination folder, $depath), have a set of folders labeled
#2,4,6,8,...50, each with the unmodified .dat and .var files output by various applescripts
#(see README Final Analysis2.rev expt 2 with hclust.txt), and a copy of m+1 for ceres.
#For current purposes you can just copy the .dat and .var and m+1 from the Final Analysis 3 deployment:
l="20"; #Top20 or Top50
cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/sagoAnalysis;
mkdir LSAcliDeploymentTop"$l";
cd LSAcliDeploymentTop"$l";
seq 2 2 50 | parallel mkdir {};
seq 2 2 50 | parallel cp "/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Sorghum/Experiment2/+revision\ with\ hclust/m+deployment/hclust/+forRgenruns/rarefaction/RgenDeploymentPaper3/"{}"/*[vdm][a+][rt1]" {};
cd ../;


#Modify var files and dat files to focus on loci in LSAGORegionsToPipeTop[25]0.txt
#setup
t="Sor";
export t;

#make working copies of the original var files
for b in $(seq 2 2 50);
  do cp ./LSAcliDeploymentTop"$l"/$b/"$t"SNP.b"$b".RgenTenv.var ./LSAcliDeploymentTop"$l"/$b/"$t"SNP.b"$b".RgenTenv.varORIG;
  done;

#turn all mstrat variable codings to ignore, 10115
for b in $(seq 2 2 50);
  do echo "$b";
    inf=./LSAcliDeploymentTop"$l"/$b/"$t"SNP.b"$b".RgenTenv.var; #define input file
    cat "$inf" | tr "\t" " " | sed 's/2 1 0 1 5/1 0 1 1 5/g' | sed 's/2 0 1 1 5/1 0 1 1 5/g' > tmp.txt;
    mv tmp.txt "$inf";
  done;

#make single space delimiter in dat file, set up as new file name, envgeopcolig, which will
#get the reference loci duplicated onto its end
cd LSAcliDeploymentTop"$l";
for i in $(seq 2 2 50);
  do echo "$i";
    sed 's/  */ /g' "$i/SorSNP.b"$i".envgeopco.dat" > "$i/SorSNP.b"$i".envgeopcolig.dat";
  done;



#identify snps that belong to 'climate' regions, add those as 21015 reference to end of var and dat files
cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/sagoAnalysis;

#or, on blip compute-0-9:
mkdir /scratch/reevesp/m+3/sagoAnalysis;
cd /scratch/reevesp/m+3/sagoAnalysis;
rsync -aP shrub@10.177.9.250:"/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/sagoAnalysis/LSAcliDeploymentTop"$l ./;
rsync -aP shrub@10.177.9.250:"/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/sagoAnalysis/LSAGORegionsToPipeTop"$l".txt" ./;


####BEGIN BASH####
mypp() {
     b=$1; #blocklength
     varfile="$wd/LSAcliDeploymentTop"$l"/$b/$t"SNP.b"$b".RgenTenv.var;
     rangefile="$wd/LSAGORegionsToPipeTop"$l".txt";  #define file containing genomic ranges to treat as reference loci
     datfile="$wd/LSAcliDeploymentTop"$l"/$b/"$t"SNP.b"$b".envgeopcolig.dat"; #the starting dat file
     
     #define coding
     referencea='2 1 0 1 5';
     target='2 0 1 1 5';
     ignore='1 0 1 1 5';

     #get the list of chromosomes that are non-relevant
     relchrs=$(cut -d. -f1 "$rangefile" | sort -nu); #range file chrs, reverse sort so 10 comes before 1, important for removal later
     varchrs=$(tail -n +4 "$varfile" | grep '\.' | cut -d. -f1 | sort -nu); #var file chrs
     nonrelchrs="$varchrs";
     for i in $relchrs;
       do nonrelchrs=$(echo "$nonrelchrs" | grep -v ^$i$);
       done;
     
     #load rangefile into memory
     rinmem=$(cat "$rangefile" | sed 's/:/./g');
     
     #refine a list of loci that might fall in the ranges covered by genes of interest
     locusstr=$(grep '\.' "$varfile" | cut -d$' ' -f1 | uniq); #read all genetic locus names in the var file into a space delimited string, works because only genetic loci have . character in them
     #remove all non-relevant chromosomes
     for i in $nonrelchrs;
       do locusstr=$(echo "$locusstr" | grep -v ^$i\.);
       done;
     #accumulate loci from the $locusstr that fall within the specified ranges
     locusstr2=$(for i in $rinmem;
                  do c=$(echo "$i" | cut -d. -f1); #get chr
                    l=$(echo "$i" | cut -d. -f2); #get low of range
                    h=$(echo "$i" | cut -d. -f4); #get high of range
                    echo "$locusstr" | awk -F'.' -v OFS='.' -v c="$c" -v l="$l" -v h="$h" '$1 == c && $2 >= l && $2 <= h {print $0}'; #add loci from locusstr that have same chromosome and fall within desired range
                  done;
                );

      #if any loci are found, write new var file and dat file
      #add 21015 reference coding and 'lig' prefix to $locusstr2
      if [[ $locusstr2 != "" ]];
      then
        
        #write new dat file
        rlist=$(for i in $locusstr2;
                 do grep -n ^"$i " "$varfile" | cut -d: -f1; #get the line numbers of reference loci from var file, these are column numbers in dat file
                 done;
               );
        for r in $rlist;
          do c=$(cut -d' ' -f$r "$datfile");
            allc=$(paste -d' ' <(echo "$allc") <(echo "$c")); #get all columns containing reference data pasted together
          done;
        allc=$(echo "$allc" | sed 's/^ //g'); #remove line leading space
        paste -d' ' "$datfile" <(echo "$allc") > "$datfile"2;

        #write new var file
        locusstr2=$(echo "$locusstr2" | sed 's/$/ '"$referencea"'/g' | sed 's/^/lig/g');
        cat "$varfile" <(echo "$locusstr2") > "$wd/LSAcliDeploymentTop"$l"/$b/$t"SNP.b"$b".Rlig.var; #write reference loci to end of new var file
      else
        cat "$varfile" > "$wd/LSAcliDeploymentTop"$l"/$b/$t"SNP.b"$b".Rlig.var; #write new var file with no changes from old
        cp "$datfile" "$datfile"2; #just copy datfile to new name, no changes required
      fi;

      echo "b$b";
}
export -f mypp;

t="Sor"; export t;
l="20"; export l;
wd=$(pwd); export wd;

seq 2 2 50 | parallel --env t --env l --env wd --env mypp mypp;

####END BASH####
	Takes 8 minutes on Shrub, ~5 on blip, for Top20.  At least 2.5 hours for Top50 on Shrub, 20 minutes on blip.

#from blip:
rsync -aP LSAcliDeploymentTop"$l" shrub@10.177.9.250:"/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/sagoAnalysis";



#calculate the number of LSA based 'climate' loci remaining for each blocklength
cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/LSAAnalysis/LSAcliDeploymentTop"$l";
g=$(for i in $(find . -name "*Rlig.var"); do echo -n "$i "; grep " 2 1 0 1 5" "$i" | wc -l; done;);
echo "$g" | sort -t/ -k3,3n | grep Top"$l";
		Result for Top20:
			./LSAcliDeploymentTop20/2/SorSNP.b2.Rlig.var       55
			./LSAcliDeploymentTop20/4/SorSNP.b4.Rlig.var       25
			./LSAcliDeploymentTop20/6/SorSNP.b6.Rlig.var       15
			./LSAcliDeploymentTop20/8/SorSNP.b8.Rlig.var       12
			./LSAcliDeploymentTop20/10/SorSNP.b10.Rlig.var        6
			./LSAcliDeploymentTop20/12/SorSNP.b12.Rlig.var        6
			./LSAcliDeploymentTop20/14/SorSNP.b14.Rlig.var        4
			./LSAcliDeploymentTop20/16/SorSNP.b16.Rlig.var        5
			./LSAcliDeploymentTop20/18/SorSNP.b18.Rlig.var        4
			./LSAcliDeploymentTop20/20/SorSNP.b20.Rlig.var        1
			./LSAcliDeploymentTop20/22/SorSNP.b22.Rlig.var        2
			./LSAcliDeploymentTop20/24/SorSNP.b24.Rlig.var        2
			./LSAcliDeploymentTop20/26/SorSNP.b26.Rlig.var        1
			./LSAcliDeploymentTop20/28/SorSNP.b28.Rlig.var        2
			./LSAcliDeploymentTop20/30/SorSNP.b30.Rlig.var        2
			./LSAcliDeploymentTop20/32/SorSNP.b32.Rlig.var        1
			./LSAcliDeploymentTop20/34/SorSNP.b34.Rlig.var        2
			./LSAcliDeploymentTop20/36/SorSNP.b36.Rlig.var        1
			./LSAcliDeploymentTop20/38/SorSNP.b38.Rlig.var        0
			./LSAcliDeploymentTop20/40/SorSNP.b40.Rlig.var        0
			./LSAcliDeploymentTop20/42/SorSNP.b42.Rlig.var        1
			./LSAcliDeploymentTop20/44/SorSNP.b44.Rlig.var        0
			./LSAcliDeploymentTop20/46/SorSNP.b46.Rlig.var        4
			./LSAcliDeploymentTop20/48/SorSNP.b48.Rlig.var        1
			./LSAcliDeploymentTop20/50/SorSNP.b50.Rlig.var        0
		Result for Top50:
			./LSAcliDeploymentTop50/2/SorSNP.b2.Rlig.var     1582
			./LSAcliDeploymentTop50/4/SorSNP.b4.Rlig.var      717
			./LSAcliDeploymentTop50/6/SorSNP.b6.Rlig.var      439
			./LSAcliDeploymentTop50/8/SorSNP.b8.Rlig.var      302
			./LSAcliDeploymentTop50/10/SorSNP.b10.Rlig.var      227
			./LSAcliDeploymentTop50/12/SorSNP.b12.Rlig.var      166
			./LSAcliDeploymentTop50/14/SorSNP.b14.Rlig.var      140
			./LSAcliDeploymentTop50/16/SorSNP.b16.Rlig.var      115
			./LSAcliDeploymentTop50/18/SorSNP.b18.Rlig.var      101
			./LSAcliDeploymentTop50/20/SorSNP.b20.Rlig.var       86
			./LSAcliDeploymentTop50/22/SorSNP.b22.Rlig.var       73
			./LSAcliDeploymentTop50/24/SorSNP.b24.Rlig.var       70
			./LSAcliDeploymentTop50/26/SorSNP.b26.Rlig.var       57
			./LSAcliDeploymentTop50/28/SorSNP.b28.Rlig.var       52
			./LSAcliDeploymentTop50/30/SorSNP.b30.Rlig.var       68
			./LSAcliDeploymentTop50/32/SorSNP.b32.Rlig.var       57
			./LSAcliDeploymentTop50/34/SorSNP.b34.Rlig.var       58
			./LSAcliDeploymentTop50/36/SorSNP.b36.Rlig.var       40
			./LSAcliDeploymentTop50/38/SorSNP.b38.Rlig.var       41
			./LSAcliDeploymentTop50/40/SorSNP.b40.Rlig.var       39
			./LSAcliDeploymentTop50/42/SorSNP.b42.Rlig.var       38
			./LSAcliDeploymentTop50/44/SorSNP.b44.Rlig.var       24
			./LSAcliDeploymentTop50/46/SorSNP.b46.Rlig.var       47
			./LSAcliDeploymentTop50/48/SorSNP.b48.Rlig.var       30
			./LSAcliDeploymentTop50/50/SorSNP.b50.Rlig.var       31

#clean up folders to contain just the files needed for execution
cd LSAcliDeploymentTop"$l";
for i in $(seq 2 2 50);
  do echo "$i";
    mv "$i"/SorSNP.b"$i".envgeopcolig.dat2 "$i"/SorSNP.b"$i".envgeopco.dat;
    rm "$i"/SorSNP.b"$i".envgeopcolig.dat;
    rm "$i/SorSNP.b"$i".RgenTenv.var" "$i/SorSNP.b"$i".RgenTenv.varORIG";
  done;

#remove null bytes for some reason polluting the var files
for i in $(seq 2 2 50);
  do echo $i;
    gsed -i 's/\x0//g' $i/*.var; #use gsed on Mac only, sed otherwise
  done;


#Run 9dSlicedBatchGeneratorRaref.sh to create slurm submit files. It will probably not be necessary to slice
#the data to fit in a 48 hour run time, since there are so many fewer loci to optimize,
#but let's just stick with that framework for consistency.
#For instructions on how to run 9dSlicedBatchGeneratorRaref.sh, see its header.
#In folder "/Volumes/J22/M+ Study/Analysis/Final analysis for 3rd paper/ForSlicedSearch" issue:

cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/ForSlicedSearch;
./9dSlicedBatchGeneratorRaref.sh;

#transfer to ceres and make some mods:
rsync -aP /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/sagoAnalysis/LSAcliDeploymentTop"$l"/* pat.reeves@login.scinet.science:"/home/pat.reeves/mplusrunsSor/lsa";
cd /home/pat.reeves/mplusrunsSor/lsa/+SorSlicedTools;
sed -i 's/short/short,medium,long,mem,mem768/' slurmArraySorb*Sliced.sh; #add more node possibilities to all slurm files
sed -i 's/--mem=0/--mem=124G/' slurmArraySorb*Sliced.sh;
sed -i 's:gcc/64/::g' slurmArraySorb*Sliced.sh;
sed -i 's:_libs::g' slurmArraySorb*Sliced.sh;
cd /home/pat.reeves/mplusrunsSor/lsa;
for i in $(seq 2 2 50);
  do echo $i;
    cp /home/pat.reeves/m+1/m+1 $i;
  done;
cd /home/pat.reeves/mplusrunsSor/lsa/+SorSlicedTools;
chgrp -Rv proj-patellifolia ~/mplusrunsSor/lsa; #chgrp for all folders and files so they can be written to by my or ann's account
chmod -Rv g+w ~/mplusrunsSor/lsa; #allow write on all objects for all members of group
cd /home/pat.reeves/mplusrunsSor/lsa/+SorSlicedTools;
grep ": sbatch" *Sliced.sh | cut -d: -f3 | sort -t'b' -k3,3n; #print all the sbatch commands

#or move to blip, from blip issue:
mkdir /home/reevesp/mplusrunsSor/lsa;
rsync -aP shrub@10.177.9.250:"/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/sagoAnalysis/LSAcliDeploymentTop"$l"/*" /home/reevesp/mplusrunsSor/lsa;

#convert ceres SLURM files to work on blip, add a set that maximize on raw allele account ("slurmArray*.L.sh")
cd /home/reevesp/mplusrunsSor/lsa/+SorSlicedTools;
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
cd /home/reevesp/mplusrunsSor/lsa;
for i in $(seq 2 2 50);
  do echo $i;
    cp /home/reevesp/Mplus/m+1 $i;
    cp /home/reevesp/Mplus/m+1L $i;
  done;

#make a copy of the 9d* execution scripts and modify for use with m+1L
cd /home/reevesp/mplusrunsSor/lsa;
for i in $(seq 2 2 50);
  do echo $i;
    ls $i/9d* | parallel "sed 's/m+1/m+1L/g' {} > {}.tmp";
    ls $i/9d*.tmp | cut -d'.' -f1-2 | parallel "mv {}.sh.tmp {}.L.sh";
    chmod a+x $i/*.sh;
  done;
 
#run on blip, takes ~17 hrs for Top20, ~28 hrs for Top50 
cd /home/reevesp/mplusrunsSor/lsa/+SorSlicedTools;
grep ": sbatch" *Sliced.sh | cut -d: -f3 | sort -t'b' -k3,3n; #print all the sbatch commands
grep ": sbatch" *Sliced.L.sh | cut -d: -f3 | sort -t'b' -k3,3n; #print all the sbatch commands
 
  
  
