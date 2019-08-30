	Expert filtering of GO terms found in Sorghum to include things related to lignin
	biosynthesis and maybe, more broadly, biofuel production.
	
	Working folder is /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/GOAnalysis.
	

############DETERMINE RELEVANT LOCI############

#extract all GO terms annotated in Sorghum genome
cut -d$'\t' -f5 SorRawGAF.txt | sort -u > AllGOTermsInSorghum.txt; #there are 2702 terms

###Generate species specific corpus of GO terms###
#go-basic.obo.corpus.txt comes from READMEsago.txt
#you don't have to do this, you can just get Sor.GOcorpus.txt from sagoAnalysis folder

spp="Sor"; #"At", "Pop", "Sor"
cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/GOAnalysis;
b=$(cut -d$'\t' -f5 "$spp"RawGAF.txt | sort -u); #extract all unique GO terms from the RawGAF file (At,6928; Pop,2841; Sor,2702)

#extract all lines in go-basic.obo.corpus.txt that match the species specific GO terms ($b), 2697/2702 are found
oboout="$spp".GOcorpus.txt;
>"$oboout";
for i in $b;
  do echo "$i";
    grep ^"$i"@ go-basic.obo.corpus.txt >> "$oboout";
  done;
###END Generate species specific corpus of GO terms###


####LIGNIN BLOCK#####

#figure out which GO terms are related to lignin using a simple text search
grep lignin Sor.GOcorpus.txt;
		Returns:
		GO:0005618@cell wall@cellular_component@The rigid or semi-rigid envelope lying outside the cell membrane of plant, fungal, most prokaryotic cells and some protozoan parasites, maintaining their shape and protecting them from osmotic lysis. In plants it is made of cellulose and, often, lignin; in fungi it is composed largely of polysaccharides; in bacteria it is composed of peptidoglycan; in protozoan parasites such as Giardia species, it's made of carbohydrates and proteins.
		GO:0009808@lignin metabolic process@biological_process@The chemical reactions and pathways involving lignins, a class of polymers of phenylpropanoid units.
		GO:0009809@lignin biosynthetic process@biological_process@The chemical reactions and pathways resulting in the formation of lignins, a class of polymers formed by the dehydrogenetive radical polymerization of various phenylpropanoid monomers.
		
#figure out how many genes we're talking about (just for GO:000980[89], ignoring GO:0005618, cell wall cellular component for the moment
grep 'GO:000980'[89] SorRawGAF.txt > LigninProcessGOgenes.txt;
#this list requires some manual processing, the SORBI_ genes are all duplicated with Sb genes, so remove those
grep 'GO:000980'[89] SorRawGAF.txt | grep -v SORBI_ > LigninProcessGOgenes.txt;

#this finds 15 genes, X means they are in Hannah's expert curated list of monolignols
#some of those that are not in Hannah's list are involved in lignin catabolism, which wouldn't be included since hers are lignin biosynthesis

for i in Sb01g039690 Sb03g038160 Sb09g022510 Sb09g022460 Sb03g028920 Sb01g041085 Sb03g037380 Sb03g039520 Sb03g039570 Sb01g017270 Sb02g002630 C4H Sb03g039530 Sb05g007210 Sb09g024210;
  do echo -n $i; 
    if grep -q $i HannahsGenes.txt; then echo " X"; else echo; fi;
  done;
		Sb01g039690
		Sb03g038160 X
		Sb09g022510
		Sb09g022460
		Sb03g028920
		Sb01g041085
		Sb03g037380 X
		Sb03g039520
		Sb03g039570
		Sb01g017270 X
		Sb02g002630 X
		C4H
		Sb03g039530
		Sb05g007210 X
		Sb09g024210 X

#get the genome positions for these genes from Sb14genes.bed, manually save to a file called LigninGORegionsToPipe.txt
glist=$(cut -d$'\t' -f3 LigninProcessGOgenes.txt);
for g in $glist;
  do echo -n "$g ";
    a=$(grep "$g" Sb14genes.bed);
    c=$(echo "$a" | cut -d$'\t' -f1 | cut -d_ -f2); #get the chromosome
    d=$(echo "$a" | cut -d$'\t' -f2); #get bp start
    e=$(echo "$a" | cut -d$'\t' -f3); #get bp end
    b=$(echo "$c"."$d":"$c"."$e"); #get string defining start and end, suitable for haplotypista input
    echo "$b";
  done;
		Result (total of 14 genes):
		Sb01g039690 1.63184944:1.63187783
		Sb03g038160 3.66089840:3.66091809
		Sb09g022510 9.52167275:9.52171247
		Sb09g022460 9.52116619:9.52118445
		Sb03g028920 3.57076630:3.57078499
		Sb01g041085 1.64371029:1.64371593
		Sb03g037380 3.65339179:3.65342378
		Sb03g039520 3.67197922:3.67201832
		Sb03g039570 3.67242056:3.67244587
		Sb01g017270 1.17757118:1.17761789
		Sb02g002630 2.2782642:2.2784313
		C4H .:.                           #this gene naming scheme not used in annotation Sb14genes.bed, so not found
		Sb03g039530 3.67205947:3.67208947
		Sb05g007210 5.12755888:5.12758128
		Sb09g024210 9.53754404:9.53758600 
 
#if we also include GO0005618 cell wall cellular component, we find a total of 124 genes:
#Decided, for now, not to use this larger set of genes
grep 'GO:000980'[89] SorRawGAF.txt | grep -v SORBI_ > AllLigninGOgenes.txt;
grep 'GO:0005618' SorRawGAF.txt | grep -v SORBI_ >> AllLigninGOgenes.txt;
#determine how many of HannahsGenes are in this larger set:
#turns out just 6, same as above, so adding cell wall cellular component doesn't pick up any
#more of the expert curated genes
for i in $(cut -d$'\t' -f3 AllLigninGOgenes.txt);
  do echo -n $i; 
    if grep -q $i HannahsGenes.txt; then echo " X"; else echo; fi;
  done;


#since probably not all cell wall cellular component genes are relevant to lignin biosynthesis
#try the "second order approach": find all GO terms corresponding to genes found for
#		GO:0005618@cell wall@cellular_component@The rigid or semi-rigid envelope lying outside the cell membrane of plant, fungal, most prokaryotic cells and some protozoan parasites, maintaining their shape and protecting them from osmotic lysis. In plants it is made of cellulose and, often, lignin; in fungi it is composed largely of polysaccharides; in bacteria it is composed of peptidoglycan; in protozoan parasites such as Giardia species, it's made of carbohydrates and proteins.
#then remove
#		GO:0009808@lignin metabolic process@biological_process@The chemical reactions and pathways involving lignins, a class of polymers of phenylpropanoid units.
#		GO:0009809@lignin biosynthetic process@biological_process@The chemical reactions and pathways resulting in the formation of lignins, a class of polymers formed by the dehydrogenetive radical polymerization of various phenylpropanoid monomers.
#then let Hannah decide if they are liable to be important.  This is essentially a filter on GO:0005618.

a=$(grep 'GO:0005618' SorGAF.txt | grep -v 'GO:0009808' | grep -v 'GO:0009809' | cut -d$'\t' -f2 | tr ';' '\n' | sort -u);
		Returns 12 unique terms (11 plus query)
		GO:0004252
		GO:0004565
		GO:0005618
		GO:0005773
		GO:0009506
		GO:0009834
		GO:0010497
		GO:0016759
		GO:0030244
		GO:0033609
		GO:0046564
		GO:2000280

for i in $a; do grep "$i" Sor.GOcorpus.txt; done; #show the description for the secondary GO terms
		GO:0004252@serine-type endopeptidase activity@molecular_function@Catalysis of the hydrolysis of internal, alpha-peptide bonds in a polypeptide chain by a catalytic mechanism that involves a catalytic triad consisting of a serine nucleophile that is activated by a proton relay involving an acidic residue (e.g. aspartate or glutamate) and a basic residue (usually histidine).
		GO:0004565@beta-galactosidase activity@molecular_function@Catalysis of the hydrolysis of terminal, non-reducing beta-D-galactose residues in beta-D-galactosides.
		GO:0005618@cell wall@cellular_component@The rigid or semi-rigid envelope lying outside the cell membrane of plant, fungal, most prokaryotic cells and some protozoan parasites, maintaining their shape and protecting them from osmotic lysis. In plants it is made of cellulose and, often, lignin; in fungi it is composed largely of polysaccharides; in bacteria it is composed of peptidoglycan; in protozoan parasites such as Giardia species, it's made of carbohydrates and proteins.
		GO:0005773@vacuole@cellular_component@A closed structure, found only in eukaryotic cells, that is completely surrounded by unit membrane and contains liquid material. Cells contain one or several vacuoles, that may have different functions from each other. Vacuoles have a diverse array of functions. They can act as a storage organelle for nutrients or waste products, as a degradative compartment, as a cost-effective way of increasing cell size, and as a homeostatic regulator controlling both turgor pressure and pH of the cytosol.
		GO:0009506@plasmodesma@cellular_component@A fine cytoplasmic channel, found in all higher plants, that connects the cytoplasm of one cell to that of an adjacent cell.
		GO:0009834@plant-type secondary cell wall biogenesis@biological_process@A cellular process that results in the biosynthesis of constituent macromolecules, assembly, and arrangement of constituent parts of inextensible cellulose- and pectin-containing cell walls that are formed between the plasma membrane and primary cell wall after cell expansion is complete. An example of this is found in Arabidopsis thaliana.
		GO:0010497@plasmodesmata-mediated intercellular transport@biological_process@The movement of substances between cells via plasmodesmata. Plasmodesmata is a fine cytoplasmic channel, found in all higher plants, that connects the cytoplasm of one cell to that of an adjacent cell.
		GO:0016759@cellulose synthase activity@molecular_function@Catalysis of the reaction: nucleoside-disphosphate-glucose + ((1,4)-beta-D-glucosyl)(n) = nucleoside-disphosphate + ((1,4)-beta-D-glucosyl)(n+1).
		GO:0030244@cellulose biosynthetic process@biological_process@The chemical reactions and pathways resulting in the formation of cellulose, a linear beta1-4 glucan of molecular mass 50-400 kDa with the pyranose units in the -4C1 conformation.
		GO:0033609@oxalate metabolic process@biological_process@The chemical reactions and pathways involving oxalate, the organic acid ethanedioate.
		GO:0046564@oxalate decarboxylase activity@molecular_function@Catalysis of the reaction: H(+) + oxalate = CO(2) + formate.
		GO:2000280@regulation of root development@biological_process@Any process that modulates the frequency, rate or extent of root development.

#identify unique patterns of GO term association
b=$(grep 'GO:0005618' SorGAF.txt | grep -v 'GO:0009808' | grep -v 'GO:0009809' | cut -d$'\t' -f2 | sort -u);
c=$(for i in $b; do echo "$i" | tr ';' '\n' | sort -u | tr '\n' ';'; echo; done;); #sort each line, remove redundancies
echo "$c" | sort -u; #remove any new redundancies
		All genes that include cell wall (GO:0005618), and exclude GO:0009808;GO:0009809
		can be categorized into 5 unique groups of GO terms:
		GO:0004252;GO:0005618;  
			+ GO:0004252@serine-type endopeptidase activity@molecular_function@Catalysis of the hydrolysis of internal, alpha-peptide bonds in a polypeptide chain by a catalytic mechanism that involves a catalytic triad consisting of a serine nucleophile that is activated by a proton relay involving an acidic residue (e.g. aspartate or glutamate) and a basic residue (usually histidine).
		GO:0004565;GO:0005618;GO:0005773;
			+ GO:0004565@beta-galactosidase activity@molecular_function@Catalysis of the hydrolysis of terminal, non-reducing beta-D-galactose residues in beta-D-galactosides.
			+ GO:0005773@vacuole@cellular_component@A closed structure, found only in eukaryotic cells, that is completely surrounded by unit membrane and contains liquid material. Cells contain one or several vacuoles, that may have different functions from each other. Vacuoles have a diverse array of functions. They can act as a storage organelle for nutrients or waste products, as a degradative compartment, as a cost-effective way of increasing cell size, and as a homeostatic regulator controlling both turgor pressure and pH of the cytosol.
		GO:0005618;GO:0009506;GO:0010497;GO:0033609;GO:0046564;GO:2000280;
			+ GO:0009506@plasmodesma@cellular_component@A fine cytoplasmic channel, found in all higher plants, that connects the cytoplasm of one cell to that of an adjacent cell.
			+ GO:0010497@plasmodesmata-mediated intercellular transport@biological_process@The movement of substances between cells via plasmodesmata. Plasmodesmata is a fine cytoplasmic channel, found in all higher plants, that connects the cytoplasm of one cell to that of an adjacent cell.
			+ GO:0033609@oxalate metabolic process@biological_process@The chemical reactions and pathways involving oxalate, the organic acid ethanedioate.
			+ GO:0046564@oxalate decarboxylase activity@molecular_function@Catalysis of the reaction: H(+) + oxalate = CO(2) + formate.
			+ GO:2000280@regulation of root development@biological_process@Any process that modulates the frequency, rate or extent of root development.
		GO:0005618;GO:0009834;GO:0016759;GO:0030244;
			+ GO:0009834@plant-type secondary cell wall biogenesis@biological_process@A cellular process that results in the biosynthesis of constituent macromolecules, assembly, and arrangement of constituent parts of inextensible cellulose- and pectin-containing cell walls that are formed between the plasma membrane and primary cell wall after cell expansion is complete. An example of this is found in Arabidopsis thaliana.
			+ GO:0016759@cellulose synthase activity@molecular_function@Catalysis of the reaction: nucleoside-disphosphate-glucose + ((1,4)-beta-D-glucosyl)(n) = nucleoside-disphosphate + ((1,4)-beta-D-glucosyl)(n+1).
			+ GO:0030244@cellulose biosynthetic process@biological_process@The chemical reactions and pathways resulting in the formation of cellulose, a linear beta1-4 glucan of molecular mass 50-400 kDa with the pyranose units in the -4C1 conformation.
		GO:0005618;GO:0033609;GO:0046564;
			+ GO:0033609@oxalate metabolic process@biological_process@The chemical reactions and pathways involving oxalate, the organic acid ethanedioate.
			+ GO:0046564@oxalate decarboxylase activity@molecular_function@Catalysis of the reaction: H(+) + oxalate = CO(2) + formate.


#identify all the GO combinations that correspond to Hannah's genes:
for i in $(sed 's/\.1//g' HannahsGenes.txt); 
  do grep "$i" SorGAF.txt;
  done;
		Returns 38/135:
		Sb01g017270	GO:0016709;GO:0016020;GO:0009809
		Sb04g031110	GO:0043565;GO:0030154;GO:0006357;GO:0001135;GO:0005634;GO:0000981;GO:0044212
		Sb07g003860	GO:0008757;GO:0008171;GO:0019438;GO:0005829
		Sb09g024210	GO:0016709;GO:0016020;GO:0009809
		Sb01g042900	GO:0008757;GO:0008171;GO:0019438;GO:0005829
		Sb02g002630	GO:0016709;GO:0016020;GO:0009809
		Sb02g035230	GO:0007165
		Sb03g034845	GO:0008757;GO:0008171;GO:0019438
		Sb03g037380	GO:0016709;GO:0016020;GO:0009809
		Sb03g038160	GO:0016020;GO:0044550;GO:0009808;GO:0016710
		Sb04g011090	GO:0008757;GO:0008171;GO:0019438
		Sb04g017460	GO:0016709;GO:0016020;GO:0044550
		Sb04g036890	GO:0008757;GO:0008171;GO:0019438
		Sb04g036900	GO:0008757;GO:0008171;GO:0019438
		Sb04g037820	GO:0008757;GO:0008171;GO:0019438;GO:0005829
		Sb05g003780	GO:0008757;GO:0008171;GO:0019438
		Sb05g007043	GO:0008757;GO:0008171;GO:0019438
		Sb05g007210	GO:0016709;GO:0016020;GO:0009809
		Sb05g008830	GO:0008757;GO:0008171;GO:0019438
		Sb05g010100	GO:0008757;GO:0008171;GO:0019438
		Sb05g026710	GO:0008757;GO:0008171;GO:0019438
		Sb05g026730	GO:0008757;GO:0008171;GO:0019438
		Sb06g000816	GO:0008757;GO:0008171;GO:0019438
		Sb07g004680	GO:0008757;GO:0008171;GO:0019438;GO:0005829
		Sb07g004690	GO:0008757;GO:0008171;GO:0019438;GO:0005829
		Sb07g004710	GO:0008757;GO:0008171;GO:0019438;GO:0005829
		Sb07g005970	GO:0008757;GO:0008171;GO:0019438
		Sb07g011460	GO:0008757;GO:0008171;GO:0019438
		Sb07g024270	GO:0008757;GO:0008171;GO:0019438
		Sb08g001245	GO:0008757;GO:0008171;GO:0019438
		Sb08g005125	GO:0008757;GO:0008171;GO:0019438
		Sb09g003740	GO:0008757;GO:0008171;GO:0019438
		Sb09g025510	GO:0008757;GO:0008171;GO:0019438
		Sb09g025530	GO:0008757;GO:0008171;GO:0019438
		Sb09g025540	GO:0008757;GO:0008171;GO:0019438
		Sb09g025550	GO:0008757;GO:0008171;GO:0019438
		Sb09g025570	GO:0008757;GO:0008171;GO:0019438
		Sb10g027640	GO:0008757;GO:0008171;GO:0019438
		
#condense to unique sets
d=$(for i in $(sed 's/\.1//g' HannahsGenes.txt); 
  do grep "$i" SorGAF.txt;
  done;
);
echo "$d" | cut -d$'\t' -f2 | sort | uniq -c; #count occurrence of unique GO sets
		Returns
		   1 GO:0007165
		     + GO:0007165@signal transduction@biological_process@The cellular process in which a signal is conveyed to trigger a change in the activity or state of a cell. Signal transduction begins with reception of a signal (e.g. a ligand binding to a receptor or receptor activation by a stimulus such as light), or for signal transduction in the absence of ligand, signal-withdrawal or the activity of a constitutively active receptor. Signal transduction ends with regulation of a downstream cellular process, e.g. regulation of transcription or regulation of a metabolic process. Signal transduction covers signaling from receptors located on the surface of the cell and signaling via molecules located within the cell. For signaling between cells, signal transduction is restricted to events at and within the receiving cell.
		  23 GO:0008757;GO:0008171;GO:0019438
		     + GO:0008757@S-adenosylmethionine-dependent methyltransferase activity@molecular_function@Catalysis of the transfer of a methyl group from S-adenosyl-L-methionine to a substrate.
		     + GO:0008171@O-methyltransferase activity@molecular_function@Catalysis of the transfer of a methyl group to the oxygen atom of an acceptor molecule.
		     + GO:0019438@aromatic compound biosynthetic process@biological_process@The chemical reactions and pathways resulting in the formation of aromatic compounds, any substance containing an aromatic carbon ring.
		   6 GO:0008757;GO:0008171;GO:0019438;GO:0005829
		     + GO:0008757@S-adenosylmethionine-dependent methyltransferase activity@molecular_function@Catalysis of the transfer of a methyl group from S-adenosyl-L-methionine to a substrate.
		     + GO:0008171@O-methyltransferase activity@molecular_function@Catalysis of the transfer of a methyl group to the oxygen atom of an acceptor molecule.
		     + GO:0019438@aromatic compound biosynthetic process@biological_process@The chemical reactions and pathways resulting in the formation of aromatic compounds, any substance containing an aromatic carbon ring.
		     + GO:0005829@cytosol@cellular_component@The part of the cytoplasm that does not contain organelles but which does contain other particulate matter, such as protein complexes.
		   1 GO:0016020;GO:0044550;GO:0009808;GO:0016710
		     + 
		     + 
		     + 
		     + 
		   5 GO:0016709;GO:0016020;GO:0009809
		     + 
		     + 
		     + 
		   1 GO:0016709;GO:0016020;GO:0044550
		     + 
		     + 
		     + 
		   1 GO:0043565;GO:0030154;GO:0006357;GO:0001135;GO:0005634;GO:0000981;GO:0044212
		     + 
		     + 
		     + 
		     + 
		     + 
		     + 
		     + 
####END LIGNIN BLOCK#####



####DEFENSE BLOCK#####
#figure out which GO terms are related to defense using a simple text search
grep defense Sor.GOcorpus.txt;
		Returns:
		GO:0006952@defense response@biological_process@Reactions, triggered in response to the presence of a foreign body or the occurrence of an injury, which result in restriction of damage to the organism attacked or prevention/recovery from the infection caused by the attack.
		GO:0031347@regulation of defense response@biological_process@Any process that modulates the frequency, rate or extent of a defense response.
		GO:0045087@innate immune response@biological_process@Innate immune responses are defense responses mediated by germline encoded components that directly recognize components of potential pathogens.
		GO:0098542@defense response to other organism@biological_process@Reactions triggered in response to the presence of another organism that act to protect the cell or organism from damage caused by that organism.
		
#figure out how many genes we're talking about, some oddball nomenclature so remove those (SORBI,15A,15B,15C,18 
>DefenseGOgenes.txt;
for i in GO:0006952 GO:0031347 GO:0045087 GO:0098542;
  do grep "$i" SorRawGAF.txt | grep -v SORBI_ | grep -v $'\t'15A$'\t' | grep -v $'\t'15B$'\t' | grep -v $'\t'15C$'\t' | grep -v $'\t'18$'\t' >> DefenseGOgenes.txt;
  done;

#this finds 187 genes
#get the genome positions for these genes from Sb14genes.bed, save to a file called DefenseGORegionsToPipe.txt
glist=$(cut -d$'\t' -f3 DefenseGOgenes.txt);
(for g in $glist;
  do echo -n "$g ";
#    a=$(grep "$g" Sb14genes.bed);
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
  done;) | sort -t' ' -k2,2n > DefenseGORegionsToPipe.txt;



############BUILD DATASETS AND PROCESS############

#in the folder GOAnalysis/GOdefDeployment (the destination folder, $depath), have a set of folders labeled
#2,4,6,8,...50, each with the unmodified .dat and .var files output by various applescripts
#(see README Final Analysis2.rev expt 2 with hclust.txt), and a copy of m+1 for ceres.
#For current purposes you can just copy the .dat and .var and m+1 from the Final Analysis 3 deployment:
mkdir GOdefDeployment;
cd GOdefDeployment;
seq 2 2 50 | parallel mkdir {};
seq 2 2 50 | parallel cp "/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Sorghum/Experiment2/+revision\ with\ hclust/m+deployment/hclust/+forRgenruns/rarefaction/RgenDeployment/"{}"/*[vdm][a+][rt1]" {};


#Modify var files and dat files to focus on loci in DefenseGORegionsToPipe.txt
#setup
t="Sor";
export t;

#make working copies of the original var files
for b in $(seq 2 2 50);
  do cp ./GOdefDeployment/$b/"$t"SNP.b"$b".RgenTenv.var ./GOdefDeployment/$b/"$t"SNP.b"$b".RgenTenv.varORIG;
  done;

#turn all mstrat variable codings to ignore, 10115
for b in $(seq 2 2 50);
  do echo "$b";
    inf=./GOdefDeployment/$b/"$t"SNP.b"$b".RgenTenv.var; #define input file
    cat "$inf" | tr "\t" " " | sed 's/2 1 0 1 5/1 0 1 1 5/g' | sed 's/2 0 1 1 5/1 0 1 1 5/g' > tmp.txt;
    mv tmp.txt "$inf";
  done;

#make single space delimiter in dat file, set up as new file name, envgeopcolig, which will
#get the reference loci duplicated onto its end
cd GOdefDeployment;
for i in $(seq 2 2 50);
  do echo "$i";
    sed 's/  */ /g' "$i/SorSNP.b"$i".envgeopco.dat" > "$i/SorSNP.b"$i".envgeopcolig.dat";
  done;



#identify snps that belong to 'defense' regions, add those as 21015 reference to end of var and dat files
cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/GOAnalysis;

####BEGIN BASH####
mypp() {
     b=$1; #blocklength
     varfile="$wd/GOdefDeployment/$b/$t"SNP.b"$b".RgenTenv.var;
     rangefile="$wd/DefenseGORegionsToPipe.txt";  #define file containing genomic ranges to treat as reference loci
     datfile="$wd/GOdefDeployment/$b/"$t"SNP.b"$b".envgeopcolig.dat"; #the starting dat file
     
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
        cat "$varfile" <(echo "$locusstr2") > "$wd/GOdefDeployment/$b/$t"SNP.b"$b".Rlig.var; #write reference loci to end of new var file
      else
        cat "$varfile" > "$wd/GOdefDeployment/$b/$t"SNP.b"$b".Rlig.var; #write new var file with no changes from old
        cp "$datfile" "$datfile"2; #just copy datfile to new name, no changes required
      fi;

      echo "b$b";
}
export -f mypp;

t="Sor";
export t;
wd=$(pwd);
export wd

seq 2 2 50 | parallel --env t --env wd --env mypp mypp;

####END BASH####
	Takes 8 minutes on Shrub, ~5 on blip




#calculate the number of GO based 'defense' loci remaining for each blocklength
cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/GOAnalysis/GOdefDeployment;
g=$(for i in $(find . -name "*Rlig.var"); do echo -n "$i "; grep " 2 1 0 1 5" "$i" | wc -l; done;);
echo "$g" | sort -t/ -k2,2n;
		Result:
		./2/SorSNP.b2.Rlig.var      101
		./4/SorSNP.b4.Rlig.var       44
		./6/SorSNP.b6.Rlig.var       24
		./8/SorSNP.b8.Rlig.var       18
		./10/SorSNP.b10.Rlig.var       10
		./12/SorSNP.b12.Rlig.var        8
		./14/SorSNP.b14.Rlig.var        8
		./16/SorSNP.b16.Rlig.var        8
		./18/SorSNP.b18.Rlig.var        8
		./20/SorSNP.b20.Rlig.var        1
		./22/SorSNP.b22.Rlig.var        6
		./24/SorSNP.b24.Rlig.var        2
		./26/SorSNP.b26.Rlig.var        3
		./28/SorSNP.b28.Rlig.var        3
		./30/SorSNP.b30.Rlig.var        2
		./32/SorSNP.b32.Rlig.var        1
		./34/SorSNP.b34.Rlig.var        2
		./36/SorSNP.b36.Rlig.var        0
		./38/SorSNP.b38.Rlig.var        2
		./40/SorSNP.b40.Rlig.var        0
		./42/SorSNP.b42.Rlig.var        1
		./44/SorSNP.b44.Rlig.var        2
		./46/SorSNP.b46.Rlig.var        1
		./48/SorSNP.b48.Rlig.var        0
		./50/SorSNP.b50.Rlig.var        3

#clean up folders to contain just the files needed for execution
cd GOdefDeployment;
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

#move to blip, from blip issue:
mkdir /home/reevesp/mplusrunsSor/go;
rsync -aP shrub@10.177.9.250:"/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/GOAnalysis/GOdefDeployment/*" /home/reevesp/mplusrunsSor/go;

#convert ceres SLURM files to work on blip
cd /home/reevesp/mplusrunsSor/go/+SorSlicedTools;
for i in $(seq 2 2 50);
  do echo $i;
    sed -i 's/-p short/-p CLUSTER/' slurmArraySorb"$i"Sliced.sh;
    sed -i 's/-n 40/-n 24/' slurmArraySorb"$i"Sliced.sh;
    sed -i 's/--mem=0/--mem=44G/' slurmArraySorb"$i"Sliced.sh;
    sed -i 's/-t 48:00:00/-t 0/' slurmArraySorb"$i"Sliced.sh;
    sed -i 's/^module/#module/g' slurmArraySorb"$i"Sliced.sh;
    #sed -i 's:#rm -fr $TMPDIR:rm -fr $TMPDIR/*:' slurmArraySorb"$i"Sliced.sh;
  done;
  
#reduce to -n 12 for b2 runs on blip to not exceed memory
sed -i 's/-n 24/-n 12/' slurmArraySorb2Sliced.sh;

#get the proper m+1
cd /home/reevesp/mplusrunsSor/go;
for i in $(seq 2 2 50);
  do echo $i;
    cp /home/reevesp/Mplus/m+1 $i;
  done;

#run on blip, takes ~17 hrs
cd /home/reevesp/mplusrunsSor/go/+SorSlicedTools;
grep ": sbatch" *.sh | cut -d: -f3 | sort -t'b' -k3,3n; #print all the sbatch commands


#####END DEFENSE BLOCK#####







#####OTHER CLASSES BLOCK#####
grep respons Sor.GOcorpus.txt > respons1.txt; #79 terms
awk -F'@' '$2~/respons/{print $0}' Sor.GOcorpus.txt > respons2.txt; #52 terms
awk -F'@' '$2~/regulat/{print $0}' Sor.GOcorpus.txt > regulat.txt; #161 terms
grep photosynth Sor.GOcorpus.txt > photosynth.txt; #16 terms
grep seed Sor.GOcorpus.txt > seed.txt; #8 terms

for i in respons1 respons2 regulat photosynth seed;
  do echo "$i";
    off="$i"GOgenes.txt; 
    >"$off";#init outfile 
    aa=$(cut -d'@' -f1 "$i".txt | tr "\n" " "); #get the list of GO terms
    for a in $aa;
      do echo " $a";
        grep "$a" SorRawGAF.txt | grep -v SORBI_ >> "$off"
      done;
  done;

wc -l *GOgenes.txt; #how many genes per set of GO terms? --this is preliminary, prior to removing redundancies in the annotation. see below for real number
     124 AllLigninGOgenes.txt
     187 DefenseGOgenes.txt
      15 LigninProcessGOgenes.txt
     170 photosynthGOgenes.txt
    1447 regulatGOgenes.txt
     699 respons1GOgenes.txt
     524 respons2GOgenes.txt
      27 seedGOgenes.txt
#continue using regulat, photosynth and respons2

#get the genome positions for these genes from Sb14genes.bed, save to a file called *GORegionsToPipe.txt
for i in regulat photosynth respons2;
  do glist=$(cut -d$'\t' -f3 "$i"GOgenes.txt | grep ^Sb | sort -u); #sort unique to limit to unique 'Sb...' gene names
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
      done;) | gsort -t' ' -V -k2,2 -u | grep ^Sb[01][0-9]g > "$i"GORegionsToPipe.txt; #sort like version, remove non chr1-10 Sb gene names
    cut -d' ' -f2- "$i"GORegionsToPipe.txt > "$i"GORegionsToPipeHaplotypista.txt;
  done;

wc -l *GORegionsToPipeHaplotypista.txt; #how many unique genes on chr1-10 per set of GO terms?
		101  photosynthGORegionsToPipeHaplotypista.txt
	   1040  regulatGORegionsToPipeHaplotypista.txt
		455  respons2GORegionsToPipeHaplotypista.txt


  ##BUILD DATASETS AND PROCESS##

#apply below code to regulat photosynth respons2:
#in the folder GOAnalysis/GO*Deployment (the destination folder, $depath), have a set of folders labeled
#2,4,6,8,...50, each with the unmodified .dat and .var files output by various applescripts
#(see README Final Analysis2.rev expt 2 with hclust.txt), and a copy of m+1 for ceres.
#For current purposes you can just copy the .dat and .var and m+1 from the Final Analysis 3 deployment:
here=$(pwd);
for r in regulat photosynth respons2;
  do mkdir GO"$r"Deployment;
    cd GO"$r"Deployment;
    seq 2 2 50 | parallel mkdir {};
    seq 2 2 50 | parallel cp "/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 2nd\ paper/Sorghum/Experiment2/+revision\ with\ hclust/m+deployment/hclust/+forRgenruns/rarefaction/RgenDeploymentPaper3/"{}"/*[vdm][a+][rt1]" {};
	cd "$here";
  done;

#Modify var files and dat files to focus on loci in *GORegionsToPipe.txt
#setup
t="Sor";

#make working copies of the original var files
for r in regulat photosynth respons2;
  do for b in $(seq 2 2 50);
    do cp ./GO"$r"Deployment/$b/"$t"SNP.b"$b".RgenTenv.var ./GO"$r"Deployment/$b/"$t"SNP.b"$b".RgenTenv.varORIG;
    done;
  done;

#turn all mstrat variable codings to ignore, 10115
for r in regulat photosynth respons2;
  do for b in $(seq 2 2 50);
    do echo "$b";
      inf=./GO"$r"Deployment/$b/"$t"SNP.b"$b".RgenTenv.var; #define input file
      cat "$inf" | tr "\t" " " | sed 's/2 1 0 1 5/1 0 1 1 5/g' | sed 's/2 0 1 1 5/1 0 1 1 5/g' > tmp.txt;
      mv tmp.txt "$inf";
    done;
  done;

#make single space delimiter in dat file, set up as new file name, envgeopcolig, which will
#get the reference loci duplicated onto its end
here=$(pwd);
for r in regulat photosynth respons2;
  do cd GO"$r"Deployment;
    for i in $(seq 2 2 50);
      do echo "$i";
        sed 's/  */ /g' "$i/SorSNP.b"$i".envgeopco.dat" > "$i/SorSNP.b"$i".envgeopcolig.dat";
      done;
    cd "$here";
  done;



#identify snps that belong to 'regulat photosynth respons2' regions, add those as 21015 reference to end of var and dat files
cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/GOAnalysis;

####BEGIN BASH####
mypp() {
     b=$1; #blocklength
     st=$2; #search term (regulat photosynth respons2)
     varfile="$wd/GO"$st"Deployment/$b/$t"SNP.b"$b".RgenTenv.var;
     rangefile="$wd/"$st"GORegionsToPipe.txt";  #define file containing genomic ranges to treat as reference loci
     datfile="$wd/GO"$st"Deployment/$b/"$t"SNP.b"$b".envgeopcolig.dat"; #the starting dat file
     
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
        cat "$varfile" <(echo "$locusstr2") > "$wd/GO"$st"Deployment/$b/$t"SNP.b"$b".Rlig.var; #write reference loci to end of new var file
      else
        cat "$varfile" > "$wd/GO"$st"Deployment/$b/$t"SNP.b"$b".Rlig.var; #write new var file with no changes from old
        cp "$datfile" "$datfile"2; #just copy datfile to new name, no changes required
      fi;

      echo "$st b$b";
}
export -f mypp;

#perform analysis on blip compute-0-9
cd /scratch/reevesp/M+3/GOAnalysis; # here you have the GO*Deployment folders and the 
wd=$(pwd);
t="Sor";
export wd t;

parallel --env t --env wd --env mypp mypp ::: $(seq 2 2 50) ::: regulat photosynth respons2;

#for st in regulat photosynth respons2;
#  do export st;
#    echo "$st";
#    seq 2 2 50 | parallel --env t --env wd --env st --env mypp mypp;
#  done;

####END BASH####
	Takes 8 minutes on Shrub, ~5 on blip for one data set




#calculate the number of GO based 'regulat photosynth respons2' loci remaining for each blocklength
cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/GOAnalysis/GOdefDeployment;
for st in regulat photosynth respons2;
  do cd /scratch/reevesp/M+3/GOAnalysis/GO"$st"Deployment;
    g=$(for i in $(find . -name "*Rlig.var"); do echo -n "$i "; grep " 2 1 0 1 5" "$i" | wc -l; done;);
    echo "$g" | sort -t/ -k2,2n | sed 's/^/'"$st "'/g';
  done;
		Result:
		regulat ./2/SorSNP.b2.Rlig.var 804
		regulat ./4/SorSNP.b4.Rlig.var 370
		regulat ./6/SorSNP.b6.Rlig.var 225
		regulat ./8/SorSNP.b8.Rlig.var 182
		regulat ./10/SorSNP.b10.Rlig.var 126
		regulat ./12/SorSNP.b12.Rlig.var 86
		regulat ./14/SorSNP.b14.Rlig.var 77
		regulat ./16/SorSNP.b16.Rlig.var 68
		regulat ./18/SorSNP.b18.Rlig.var 60
		regulat ./20/SorSNP.b20.Rlig.var 72
		regulat ./22/SorSNP.b22.Rlig.var 58
		regulat ./24/SorSNP.b24.Rlig.var 46
		regulat ./26/SorSNP.b26.Rlig.var 30
		regulat ./28/SorSNP.b28.Rlig.var 32
		regulat ./30/SorSNP.b30.Rlig.var 36
		regulat ./32/SorSNP.b32.Rlig.var 24
		regulat ./34/SorSNP.b34.Rlig.var 25
		regulat ./36/SorSNP.b36.Rlig.var 14
		regulat ./38/SorSNP.b38.Rlig.var 25
		regulat ./40/SorSNP.b40.Rlig.var 21
		regulat ./42/SorSNP.b42.Rlig.var 29
		regulat ./44/SorSNP.b44.Rlig.var 20
		regulat ./46/SorSNP.b46.Rlig.var 20
		regulat ./48/SorSNP.b48.Rlig.var 20
		regulat ./50/SorSNP.b50.Rlig.var 26
		photosynth ./2/SorSNP.b2.Rlig.var 110
		photosynth ./4/SorSNP.b4.Rlig.var 56
		photosynth ./6/SorSNP.b6.Rlig.var 30
		photosynth ./8/SorSNP.b8.Rlig.var 18
		photosynth ./10/SorSNP.b10.Rlig.var 16
		photosynth ./12/SorSNP.b12.Rlig.var 11
		photosynth ./14/SorSNP.b14.Rlig.var 10
		photosynth ./16/SorSNP.b16.Rlig.var 8
		photosynth ./18/SorSNP.b18.Rlig.var 7
		photosynth ./20/SorSNP.b20.Rlig.var 9
		photosynth ./22/SorSNP.b22.Rlig.var 4
		photosynth ./24/SorSNP.b24.Rlig.var 6
		photosynth ./26/SorSNP.b26.Rlig.var 6
		photosynth ./28/SorSNP.b28.Rlig.var 3
		photosynth ./30/SorSNP.b30.Rlig.var 5
		photosynth ./32/SorSNP.b32.Rlig.var 4
		photosynth ./34/SorSNP.b34.Rlig.var 4
		photosynth ./36/SorSNP.b36.Rlig.var 4
		photosynth ./38/SorSNP.b38.Rlig.var 2
		photosynth ./40/SorSNP.b40.Rlig.var 4
		photosynth ./42/SorSNP.b42.Rlig.var 1
		photosynth ./44/SorSNP.b44.Rlig.var 4
		photosynth ./46/SorSNP.b46.Rlig.var 1
		photosynth ./48/SorSNP.b48.Rlig.var 1
		photosynth ./50/SorSNP.b50.Rlig.var 0
		respons2 ./2/SorSNP.b2.Rlig.var 227
		respons2 ./4/SorSNP.b4.Rlig.var 107
		respons2 ./6/SorSNP.b6.Rlig.var 60
		respons2 ./8/SorSNP.b8.Rlig.var 43
		respons2 ./10/SorSNP.b10.Rlig.var 27
		respons2 ./12/SorSNP.b12.Rlig.var 19
		respons2 ./14/SorSNP.b14.Rlig.var 15
		respons2 ./16/SorSNP.b16.Rlig.var 16
		respons2 ./18/SorSNP.b18.Rlig.var 16
		respons2 ./20/SorSNP.b20.Rlig.var 8
		respons2 ./22/SorSNP.b22.Rlig.var 13
		respons2 ./24/SorSNP.b24.Rlig.var 9
		respons2 ./26/SorSNP.b26.Rlig.var 9
		respons2 ./28/SorSNP.b28.Rlig.var 8
		respons2 ./30/SorSNP.b30.Rlig.var 11
		respons2 ./32/SorSNP.b32.Rlig.var 8
		respons2 ./34/SorSNP.b34.Rlig.var 13
		respons2 ./36/SorSNP.b36.Rlig.var 4
		respons2 ./38/SorSNP.b38.Rlig.var 6
		respons2 ./40/SorSNP.b40.Rlig.var 4
		respons2 ./42/SorSNP.b42.Rlig.var 3
		respons2 ./44/SorSNP.b44.Rlig.var 6
		respons2 ./46/SorSNP.b46.Rlig.var 5
		respons2 ./48/SorSNP.b48.Rlig.var 4
		respons2 ./50/SorSNP.b50.Rlig.var 8

#return the files to shrub
rsync -aP GO*Deployment shrub@10.177.9.250:"/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/GOAnalysis";

#perform some clean up operations on files
for st in regulat photosynth respons2;
  do cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/GOAnalysis/GO"$st"Deployment;

    #clean up folders to contain just the files needed for execution
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
  done;


#Run 9dSlicedBatchGeneratorRaref.sh to create slurm submit files. It will probably not be necessary to slice
#the data to fit in a 48 hour run time, since there are so many fewer loci to optimize,
#but let's just stick with that framework for consistency.
#For instructions on how to run 9dSlicedBatchGeneratorRaref.sh, see its header.
#You will need to change depath= for each 'regulat photosynth respons2'
#In folder "/Volumes/J22/M+ Study/Analysis/Final analysis for 3rd paper/ForSlicedSearch" issue:

cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/ForSlicedSearch;
./9dSlicedBatchGeneratorRaref.sh;

#transfer to ceres and make some mods:
rsync -aP /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/GOAnalysis/GOregulatDeployment/* pat.reeves@login.scinet.science:"/home/pat.reeves/mplusrunsSor/go";
cd /home/pat.reeves/mplusrunsSor/go/+SorSlicedTools;
sed -i 's/short/short,medium,long,mem,mem768/' slurmArraySorb*Sliced.sh; #add more node possibilities to all slurm files
sed -i 's/--mem=0/--mem=124G/' slurmArraySorb*Sliced.sh;
sed -i 's:gcc/64/::g' slurmArraySorb*Sliced.sh;
sed -i 's:_libs::g' slurmArraySorb*Sliced.sh;
cd /home/pat.reeves/mplusrunsSor/go;
for i in $(seq 2 2 50);
  do echo $i;
    cp /home/pat.reeves/m+1/m+1 $i;
  done;
cd /home/pat.reeves/mplusrunsSor/go/+SorSlicedTools;
chgrp -Rv proj-patellifolia ~/mplusrunsSor/go; #chgrp for all folders and files so they can be written to by my or ann's account
chmod -Rv g+w ~/mplusrunsSor/go; #allow write on all objects for all members of group
cd /home/pat.reeves/mplusrunsSor/go/+SorSlicedTools;
grep ": sbatch" *Sliced.sh | cut -d: -f3 | sort -t'b' -k3,3n; #print all the sbatch commands


#or move to blip, from blip issue:
mkdir /home/reevesp/mplusrunsSor/go;
rsync -aP shrub@10.177.9.250:"/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/GOAnalysis/GOregulatDeployment/*" /home/reevesp/mplusrunsSor/go;

#convert ceres SLURM files to work on blip, add a set that maximize on raw allele account ("slurmArray*.L.sh")
cd /home/reevesp/mplusrunsSor/go/+SorSlicedTools;
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
cd /home/reevesp/mplusrunsSor/go;
for i in $(seq 2 2 50);
  do echo $i;
    cp /home/reevesp/Mplus/m+1 $i;
    cp /home/reevesp/Mplus/m+1L $i;
  done;

#make a copy of the 9d* execution scripts and modify for use with m+1L
cd /home/reevesp/mplusrunsSor/go;
for i in $(seq 2 2 50);
  do echo $i;
    ls $i/9d* | parallel "sed 's/m+1/m+1L/g' {} > {}.tmp";
    ls $i/9d*.tmp | cut -d'.' -f1-2 | parallel "mv {}.sh.tmp {}.L.sh";
    chmod a+x $i/*.sh;
  done;
  
#if necessary change reps=1 to reps=10
cd /home/reevesp/mplusrunsSor/go;
for i in $(seq 2 2 50);
  do echo $i;
    find $i -name "*.sh" -exec sed -i 's/reps=1;/reps=10;/g' {} \;
  done;
 
#run on blip, takes ~17 hrs with one rep
cd /home/reevesp/mplusrunsSor/go/+SorSlicedTools;
grep ": sbatch" *Sliced.sh | cut -d: -f3 | sort -t'b' -k3,3n; #print all the sbatch commands
grep ": sbatch" *Sliced.L.sh | cut -d: -f3 | sort -t'b' -k3,3n; #print all the sbatch commands



#####END OTHER CLASSES BLOCK#####

YOU ARE HERE



#in a screen, run the 9dSlicedSumMonitor.sh by cut and paste

#you can verify that the new SUM files differ from the old ones by the addition of the reference loci, which
#are now also treated as target.  This will allow the +v2STATS file to be calculated correctly later.

for i in $(seq 2 2 50);
  do echo $i;
    diff <(cut -d$'\t' -f2 ./$i/SUM.b"$i".RligTgen.txt) <(cut -d$'\t' -f2 ../monolignolOLD/work/SUM.b"$i".RligTgen.txt)  | grep '<' | sort -u;
  done;




###ANALYSIS OF GO-SELECTED "LIGNIN" LOCI###


***SORGHUM***
	Calculate mean M+ enrichment across core sizes, for each locus, for each blocklength. This uses the
	SUMfiles as input and produces the +STATS files.  Do this on the cluster.

	You can only go as far as b=32 for GO-selected loci because there were no available
	reference loci for 34,36,42,48.
	

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
bmax=32; #max blocklength
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
        > "$p""$ofile"TMP;
       
        #fast cat
        find "$p" -maxdepth 1 -name "+STATS.R"$e"Tgen.loc*.b"$b".txt" -print0 | xargs -0 cat -- >> "$p""$ofile"TMP; #quick cat all relevant files
        awk -F[$'\t'.] '{print $0"\t"$4}' "$p""$ofile"TMP | sort -t$'\t' -n -k1,1 -n -k7,7 | cut -d$'\t' -f1-6 > "$p""$ofile"; #use awk to get the genome index using two delimiters, it is added as a new, last column, then sort performed, first on blocklength, then on index
        rm "$p""$ofile"TMP;

         #slow cat
#        for l in $uloc;
#          do echo "concatenating loc file $l to +STATS.R"$e"Tgen.b"$b".txt ...";
#            cat "$p""+STATS.R"$e"Tgen.loc"$l".b"$b".txt" >> "$p""$ofile";
#          done;

        #test whether all were concatenated
        s1=$(find "$p" -maxdepth 1 -type f | grep -c "loc"); #count the number of files with "loc" in their name using grep. can't wc -l because too many arguments
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
      then echo; # rm +STATS.R"$e"Tgen.b*.txt; #remove files if number of lines in concat file is the same as the sum of all input files
      else echo "the number of lines ain't the same. something is 'crewed.";
    fi;

  done; #$e
**********BASHEND**********
	Takes ~20 minutes on blip, 9 nodes. 

	****WATCH FOR ERROR: Signal SIGCHLD received, but no signal handler set.**** 




	add the proportion of NG,SS,NS sites from SorLocusWeights.txt to the +STATS summary file,
	generating the +v2STATS file.  Requires the SorLocusWeights.txt file as input.  

**********BASH**********
#use md5 on mac md5sum on linux
#modify weightsfile to include only the blocklengths analyzed
weightsfile="SorLocusWeights.txt";
bmin=2; #min blocklength
bmax=32; #max blocklength
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

#fix some file names
mv +STATS.Sorghum.RligTgen.txt +STATS.Sorghum.RgoTgen.txt;
mv +v2STATS.Sorghum.RligTgen.txt +v2STATS.Sorghum.RgoTgen.txt;



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
maxbl=32; #50,200
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
elist="go";
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
  #find "$path" -name "*.tmp" -print0 | xargs -0 rm;
done;
**********ENDBASH**********
	Takes ~10 minutes, about 7 hours on ceres N=26


#fix some file names if necessary
mv +SorEnrichAcrossBlocks.RligTgen.txt +SorEnrichAcrossBlocks.RgoTgen.txt;
mv +SorSUMEnrichAcrossBlocks.RligTgen.txt +SorSUMEnrichAcrossBlocks.RgoTgen.txt;



	Use R to make some pretty plots of significant enrichment across the genome.

**********BEGINR**********
#install.packages("fitdistrplus")
#install.packages("plotrix")
options(error = recover)
rm(list=ls()) 
setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/GOAnalysis/GOAnalysis/+results/rarefaction/+analysis")
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

			i="go"
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






















