Plot Rgeo and Rpco on top of Rgen to show how much better gen is.
	
	First, run READMEConsolidatedAnalysis.txt against Rgen Rgeo Rpco Rlig Rgo Rlsa over the
	range 2..50.

cd /share/space/reevesp/ConsolidatedR; #on blip, contains SUM* files from M+ runs
./READMEConsolidatedAnalysis.txt 2 50;
	
	
	
	Make publication quality pretty plots of significant enrichment across the genome.
	All data sets must have the same range of blocklengths considered
	Overlay only
**********BEGINR**********
#install.packages("fitdistrplus")
#install.packages("plotrix")
options(error = recover)
rm(list=ls()) 
setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/ComparisonFigs")
library(fitdistrplus)
library(plotrix)


par(mfrow=c(3,1))

usezero="no"            #this sets whether to include significant sites only when they are also > or < 0
                        #to use the zero criterion, say "yes", to just use the ecdf around the mean, say "no"
dirlist=c("Arabidopsis","Populus","Sorghum") # c("Arabidopsis","Populus","Sorghum")
tlist=c("At","Pop","Sor") # c("At","Pop","Sor")
bmin=2
bmax=50

#summary stats output file
summaryoutputfile=paste("SummaryStatsForOverlay.txt")
write("Taxon	Ref	mean	sd	n", file=summaryoutputfile, append=FALSE, sep="\t")

for (k in 1:length(dirlist))
{
	wd=dirlist[k]
	t=tlist[k]
	cat(t,"\n",sep="")

	ggen <- read.table(paste("+",t,"SUMEnrichAcrossBlocks.RgenTgen.",bmin,".",bmax,".txt", sep=""), header=TRUE, sep=" ")
	ggen$index <- 1:nrow(ggen) #add a column with the sequential index
	
	#add a column with the genomic position for a single linearized genome (across chromosomes)
	#iterate across chromosomes, using the last bp position as the beginning to add to for the new chromosome
	dd=0
	cc=vector() #initialize vector to hold genome position
	for (i in 1:max(ggen$chr))
	{
		cat("chr",i,"\n",sep="")
		tmp=ggen[ggen$chr==i,] #subset for curr chr
		cc=c(cc,(dd+tmp$pos))
		dd=dd+max(tmp$pos) #before looping, set the new constant to add as the largest value in the last chromosome
	}

	ggen$genindex=cc #add a column with the sequential index



	ggeo <- read.table(paste("+",t,"SUMEnrichAcrossBlocks.RgeoTgen.",bmin,".",bmax,".txt", sep=""), header=TRUE, sep=" ")
	ggeo$index <- 1:nrow(ggeo) #add a column with the sequential index
	ggeo$genindex=cc #add a column with the sequential index

	gpco <- read.table(paste("+",t,"SUMEnrichAcrossBlocks.RpcoTgen.",bmin,".",bmax,".txt", sep=""), header=TRUE, sep=" ")
	gpco$index <- 1:nrow(gpco) #add a column with the sequential index
	gpco$genindex=cc #add a column with the sequential index
	
#	glig <- read.table(paste("+",t,"SUMEnrichAcrossBlocks.RligTgen.",bmin,".",bmax,".txt", sep=""), header=TRUE, sep=" ")
#	glig$index <- 1:nrow(glig) #add a column with the sequential index
#	glig$genindex=cc #add a column with the sequential index
#
#	ggo <- read.table(paste("+",t,"SUMEnrichAcrossBlocks.RgoTgen.",bmin,".",bmax,".txt", sep=""), header=TRUE, sep=" ")
#	ggo$index <- 1:nrow(ggo) #add a column with the sequential index
#	ggo$genindex=cc #add a column with the sequential index

#	glsa <- read.table(paste("+",t,"SUMEnrichAcrossBlocks.RlsaTgen.",bmin,".",bmax,".txt", sep=""), header=TRUE, sep=" ")
#	glsa$index <- 1:nrow(ggo) #add a column with the sequential index
#	glsa$genindex=cc #add a column with the sequential index

	cat("plotting..\n")
	
	#calculate the ends of chromosomes
	chrbreaks=vector()
	for (i in 1:max(ggen$chr))
	{
		tmp=ggen[ggen$chr==i,] #subset for curr chr
		q=subset(tmp,pos==max(tmp$pos),genindex) #find the genindex of the largest pos, this is the last SNP in the chromosome
		chrbreaks=c(chrbreaks, q$genindex[1]) #convert vector returned by subset to value, add to list of breakpoints
		dd=dd+max(tmp$pos) #before looping, set the new constant to add as the largest value in the last chromosome
	}
	
	#find minimum and maximum to define y-axis limits
	mmum=min(c(ggen$sum,ggeo$sum,gpco$sum))
	mmax=max(c(ggen$sum,ggeo$sum,gpco$sum))
	
	#generally, you plot so that opacity is 1/(number files underneath +1). In this case we mess with base opacity, which helps to deal with very different numbers of SNPs between species
	op=30000/length(ggen$genindex) 	#base level of alpha-opacity calculated as 50K divided by the number of data points to achieve approximate equality of intensity among species
	
	#make output file
	png(file=paste(t,"EnrichOverlay", ".png", sep=""),width=2000,height=1000,res=300,type="quartz",bg="transparent")
	plot(ggen$genindex, ggen$sum, type="p", pch=16, col=rgb(red=182/255, green=109/255, blue=255/255, alpha=(op/1)),cex=0.25,ann=FALSE,xaxt='n',ylim=c(mmum,mmax)) #gen=purple (8) http://mkweb.bcgsc.ca/biovis2012/color-blindness-palette.png

	#plot set 1
	points(ggen$genindex, y=ggeo$sum, type="p", pch=16, col=rgb(red=255/255, green=255/255, blue=109/255, alpha=(op/2)),cex=0.25,ann=FALSE,xaxt='n') #geo=yellow (15)
	points(ggen$genindex, y=gpco$sum, type="p", pch=16, col=rgb(red=0/255, green=146/255, blue=146/255, alpha=(op/3)),cex=0.25,ann=FALSE,xaxt='n') #pco=blue (3)

	#plot set 2
#	points(ggen$genindex, y=glig$sum, type="p", pch=16, col=rgb(red=36/255, green=255/255, blue=36/255, alpha=(op/2)),cex=0.25,ann=FALSE,xaxt='n') #lig=green (14)
#	points(ggo$genindex, y=ggo$sum, type="p", pch=16, col=rgb(red=219/255, green=109/255, blue=0/255, alpha=(op/3)),cex=0.25,ann=FALSE,xaxt='n') #go=orange (13)
#2	points(ggen$genindex, y=gpco$sum, type="p", pch=16, col=rgb(red=0/255, green=146/255, blue=146/255, alpha=(op/2)),cex=0.25,ann=FALSE,xaxt='n') #pco=blue (3)
#	points(ggen$genindex, y=gpco$sum, type="p", pch=16, col=rgb(red=219/255, green=109/255, blue=0/255, alpha=(op/2)),cex=0.25,ann=FALSE,xaxt='n') #pco=orange(13)
#	points(ggen$genindex, y=ggeo$sum, type="p", pch=16, col=rgb(red=0/255, green=146/255, blue=146/255, alpha=(op/3)),cex=0.25,ann=FALSE,xaxt='n') #geo=blue (3)
#	points(ggen$genindex, y=ggeo$sum, type="p", pch=16, col=rgb(red=219/255, green=109/255, blue=0/255, alpha=(op/3)),cex=0.25,ann=FALSE,xaxt='n') #geo=orange (13)
#	points(ggen$genindex, y=ggeo$sum, type="p", pch=16, col=rgb(red=36/255, green=255/255, blue=36/255, alpha=(op/3)),cex=0.25,ann=FALSE,xaxt='n') #geo=green (14)
#2	points(ggen$genindex, y=ggeo$sum, type="p", pch=16, col=rgb(red=255/255, green=255/255, blue=109/255, alpha=(op/3)),cex=0.25,ann=FALSE,xaxt='n') #geo=yellow (15)
	abline(v=chrbreaks[-length(chrbreaks)],col="grey59") #put lines at each chromosome end, except the last chromosome
	dev.off()
	
	#print some summary stats to an output file
	write(paste(t,"gen",mean(ggen$sum),sd(ggen$sum),length(ggen$sum),sep="\t"), file=summaryoutputfile, append=TRUE, sep="\t")
	write(paste(t,"geo",mean(ggeo$sum),sd(ggeo$sum),length(ggeo$sum),sep="\t"), file=summaryoutputfile, append=TRUE, sep="\t")
	write(paste(t,"pco",mean(gpco$sum),sd(gpco$sum),length(gpco$sum),sep="\t"), file=summaryoutputfile, append=TRUE, sep="\t")
	
	ee=sum(ggen$sum-ggeo$sum > 0 )/nrow(ggen) #number of sites where ggen is better than ggeo
	cc=sum(ggen$sum-gpco$sum > 0 )/nrow(ggen) #number of sites where ggen is better than gpco
	write(paste(t,"ggen>ggeo=",ee,"ggen>gpco=",cc,sep="\t"), file=summaryoutputfile, append=TRUE, sep="\t")
}

#Calculate pairwise t-tests using R
#####BEGINR#####
#install.packages("BSDA")
options(error = recover)
rm(list=ls()) 
setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/ComparisonFigs")
library(BSDA)
library(ggplot2)

d = read.csv(file="SummaryStatsForOverlay.txt", header=TRUE, sep="\t")

#perform pairwise t-tests by taxon
#when tsum.test prints out "p-value < 2.2e-16", that is just a .Machine limitation
#on $double.eps (smallest positive floating-point number)
for (i in levels(d$Taxon))
{
    b=d[d$Taxon==i,] #subset by factor (Taxon name)
    #gen x geo
    mx=b[b$Ref=="gen",'mean'] #mean of Rgen
    sx=b[b$Ref=="gen",'sd']
    nx=b[b$Ref=="gen",'n']
    my=b[b$Ref=="geo",'mean']
    sy=b[b$Ref=="geo",'sd']
    ny=b[b$Ref=="geo",'n']
    t=tsum.test(mean.x=mx, s.x=sx, n.x=nx, mean.y=my, s.y=sy, n.y=ny, alternative="two.sided", var.equal=FALSE)
    
    cat("i","gen x geo",t$p.value,sep=" ")
    print(t)
    cat("\n")

    
    #gen x pco
    mx=b[b$Ref=="gen",'mean'] #mean of Rgen
    sx=b[b$Ref=="gen",'sd']
    nx=b[b$Ref=="gen",'n']
    my=b[b$Ref=="pco",'mean']
    sy=b[b$Ref=="pco",'sd']
    ny=b[b$Ref=="pco",'n']
    t=tsum.test(mean.x=mx, s.x=sx, n.x=nx, mean.y=my, s.y=sy, n.y=ny, alternative="two.sided", var.equal=FALSE)
    
    cat("i","gen x pco",t$p.value,sep=" ")
    print(t)
    cat("\n")
    
    #pco x geo
    mx=b[b$Ref=="pco",'mean'] #mean of Rgen
    sx=b[b$Ref=="pco",'sd']
    nx=b[b$Ref=="pco",'n']
    my=b[b$Ref=="geo",'mean']
    sy=b[b$Ref=="geo",'sd']
    ny=b[b$Ref=="geo",'n']
    t=tsum.test(mean.x=mx, s.x=sx, n.x=nx, mean.y=my, s.y=sy, n.y=ny, alternative="two.sided", var.equal=FALSE)
    
    cat("i","pco x geo",t$p.value,sep=" ")
    print(t)
    cat("\n")

}





for (i in c("2_2","2_4","2_6","2_8","2_10","2_12","2_14","2_16","2_18","2_20","2_22","2_24","2_26","2_28","2_30","2_32","2_34","2_36","2_38","2_40","2_42","2_44","2_46","2_48","2_50"))
{
  a=d[d$max.block.length==i,] #subset by maxblocklength
  for (j in c("lig.Lig","lsa.LSA","regulat.reg"))
  {
    b=a[a$modl==j,] #subset by dataset
    mx=b[b$ref=="tomarkgen$sum",'mean'] #mean of Rgen
    sx=b[b$ref=="tomarkgen$sum",'sd']
    nx=b[b$ref=="tomarkgen$sum",'n']
    my=b[b$ref=="tomarkmodl$sum",'mean']
    sy=b[b$ref=="tomarkmodl$sum",'sd']
    ny=b[b$ref=="tomarkmodl$sum",'n']
    t=tsum.test(mean.x=mx, s.x=sx, n.x=nx, mean.y=my, s.y=sy, n.y=ny, alternative="two.sided", var.equal=FALSE)
    cat(i,j,t$p.value,sep=" ")
    cat("\n")
  }
}


#####ENDR#####

























	Overlay and highlight specific positions
**********BEGINR**********
#install.packages("fitdistrplus")
#install.packages("plotrix")
options(error = recover)
rm(list=ls()) 
library(fitdistrplus)
library(plotrix)
library(shape)
library(ggplot2)


#repeat script over a bunch of different max block lengths
for (b in c("2","4","6","8","10","12","14","16","18","20","22","24","26","28","30","32","34","36","38","40","42","44","46","48","50")) 
{
	tbl=50 #terminal block length, usually 50, but may depend on range considered
	
	#setwd(paste("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/ComparisonFigs/AnalysesForRanges/2_",b,"/bestEnrich",sep=""))
	#setwd(paste("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/ComparisonFigs/AnalysesFor1Length/",b,".",b,"/bestEnrich",sep=""))
	#setwd(paste("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/ComparisonFigs/AnalysesFor1Length/TestLigBlip/",b,".",b,"/bestEnrich",sep=""))
	setwd(paste("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/ComparisonFigs/AnalysesForRanges/TestGoLigCeresRegulatBlip/",2,"_",b,"/bestEnrich",sep=""))
	par(mfrow=c(3,1))

	usezero="no"            #this sets whether to include significant sites only when they are also > or < 0
							#to use the zero criterion, say "yes", to just use the ecdf around the mean, say "no"
	dirlist=c("Sorghum")    # c("Arabidopsis","Populus","Sorghum")
	tlist=c("Sor")          # c("At","Pop","Sor")

	de="" #estimator of allelic diversity at target loci, use nothing ("") for m+ estimator, use "ALLELECNT" for allele counts

	#rfilelist=c("regulatGORegionsToPipeHaplotypista.txt", "LigninRegionsToPipeHaplotypista.txt", "LSA50GORegionsToPipeHaplotypista.txt") #file containing genomic regions to highlight, in haplotypista format
	#rfilelist=c("photosynthGORegionsToPipeHaplotypista.txt", "respons2GORegionsToPipeHaplotypista.txt", "regulatGORegionsToPipeHaplotypista.txt", "LigninRegionsToPipeHaplotypista.txt", "LSA50GORegionsToPipeHaplotypista.txt") #file containing genomic regions to highlight, in haplotypista format
	rfilelist=c("LigninRegionsToPipeHaplotypista.txt","regulatGORegionsToPipeHaplotypista.txt", "LSA50GORegionsToPipeHaplotypista.txt") #file containing genomic regions to highlight, in haplotypista format

	bmin=2 #use when considering a range of block lengths
	#bmin=as.numeric(b) #use when only considering one block length
	bmax=as.numeric(b)

	for (k in 1:length(dirlist))
	{
		wd=dirlist[k]
		t=tlist[k]
		cat(t,"\n",sep="")
		ofdall=data.frame() #initialize a data.frame to receive summary stats

		#create table for gen data
		ggen <- read.table(paste("+",t,"SUMEnrichAcrossBlocks.RgenTgen.",bmin,".",bmax,".txt",de, sep=""), header=TRUE, sep=" ")
		ggen$index <- 1:nrow(ggen) #add a column with the sequential index

		#add a column with the genomic position for a single linearized genome (across chromosomes)
		#iterate across chromosomes, using the last bp position as the beginning to add to for the new chromosome
		dd=0
		cc=vector() #initialize vector to hold genome position
		for (i in 1:max(ggen$chr))
		{
			cat("chr",i,"\n",sep="")
			tmp=ggen[ggen$chr==i,] #subset for curr chr
			cc=c(cc,(dd+tmp$pos))
			dd=dd+max(tmp$pos) #before looping, set the new constant to add as the largest value in the last chromosome
		}

		ggen$genindex=cc #add a column with the genome index
		ggen$cex=rep("0.25",length(ggen$genindex)) #add a column with plotted marker size

		#create tables for geo, pco data
		ggeo <- read.table(paste("+",t,"SUMEnrichAcrossBlocks.RgeoTgen.",bmin,".",bmax,".txt",de, sep=""), header=TRUE, sep=" ")
		ggeo$index <- 1:nrow(ggeo) #add a column with the sequential index
		ggeo$genindex=cc #add a column with the genome index
		ggeo$cex=rep("0.25",length(ggen$genindex)) #add a column with plotted marker size

		gpco <- read.table(paste("+",t,"SUMEnrichAcrossBlocks.RpcoTgen.",bmin,".",bmax,".txt",de, sep=""), header=TRUE, sep=" ")
		gpco$index <- 1:nrow(gpco) #add a column with the sequential index
		gpco$genindex=cc #add a column with the genome index
		gpco$cex=rep("0.25",length(ggen$genindex)) #add a column with plotted marker size

		#repeat over reference loci and Rlig Rgo Rlsa
		for (rfile in rfilelist)
		{
			#set the relevant models to explore (models = source of reference data)
			if (rfile=="LigninRegionsToPipeHaplotypista.txt") modllist=c("gen","lig")
			if (rfile=="DefenseGORegionsToPipeHaplotypista.txt") modllist=c("gen","go")
			if (rfile=="LSA50GORegionsToPipeHaplotypista.txt") modllist=c("gen","lsa")
			if (rfile=="LSA20GORegionsToPipeHaplotypista.txt") modllist=c("gen","lsa20")
			if (rfile=="regulatGORegionsToPipeHaplotypista.txt") modllist=c("gen","regulat")
			if (rfile=="respons2GORegionsToPipeHaplotypista.txt") modllist=c("gen","respons2")
			if (rfile=="photosynthGORegionsToPipeHaplotypista.txt") modllist=c("gen","photosynth")


			hasRun=0 #switch to make Rgeo Rgen Rpco marking routine run just once per modllist
			for (modl in modllist)
			{
				cat("bmax=",bmax,"\n",sep="")
				cat(modl, ".", substr(rfile,1,5), "\n", sep="")
		
				#create tables for lig, go, lsa data
				gmodl <- read.table(paste("+",t,"SUMEnrichAcrossBlocks.R",modl,"Tgen.",bmin,".",bmax,".txt",de, sep=""), header=TRUE, sep=" ")
				gmodl$index <- 1:nrow(gmodl) #add a column with the sequential index
				gmodl$genindex=cc #add a column with the genome index
				gmodl$cex=rep("0.25",length(ggen$genindex)) #add a column with plotted marker size

				#extract desired genomic ranges, these will be overlaid with a different size and possibly color
				regions <- read.table(text = gsub("\\.", ":", readLines(rfile)),sep=":",header=FALSE)
				colnames(regions)=c("chr", "startpos", "endchr", "endpos")
				regions=regions[order(regions$chr,regions$startpos),] #sort
	
				datalist=list()
				tomarkmodl=gmodl[0,] #create an empty table with same colnames as gmodl
				cat("nloci ",nrow(regions),"\n",sep="")
				for (i in 1:nrow(regions))
				{
					datalist[[i]]=subset(gmodl, (chr==regions$chr[i] & pos>=regions$startpos[i] & pos<=regions$endpos[i]))
				}
				tomarkmodl=dplyr::bind_rows(datalist) # or use: do.call(rbind,datalist)
				tomarkmodl$cex=rep("1",length(tomarkmodl$cex)) #modify column to enlarge marker size


				#extract values at loci in the desired genomic range for Rgen, Rgeo, Rpco
				cat("marking Rgen Rgeo Rpco..\n")
				tomarkgen = subset(ggen, subset = genomeindex %in% tomarkmodl$genomeindex)
				tomarkgen$cex=rep("1",length(tomarkgen$cex)) #modify column to enlarge marker size
				tomarkgeo = subset(ggeo, subset = genomeindex %in% tomarkmodl$genomeindex)
				tomarkgeo$cex=rep("1",length(tomarkgeo$cex)) #modify column to enlarge marker size
				tomarkpco = subset(gpco, subset = genomeindex %in% tomarkmodl$genomeindex)
				tomarkpco$cex=rep("1",length(tomarkpco$cex)) #modify column to enlarge marker size




				cat("plotting..\n")
	
				#calculate the ends of chromosomes
				chrbreaks=vector()
				for (i in 1:max(ggen$chr))
				{
					tmp=ggen[ggen$chr==i,] #subset for curr chr
					q=subset(tmp,pos==max(tmp$pos),genindex) #find the genindex of the largest pos, this is the last SNP in the chromosome
					chrbreaks=c(chrbreaks, q$genindex[1]) #convert vector returned by subset to value, add to list of breakpoints
					dd=dd+max(tmp$pos) #before looping, set the new constant to add as the largest value in the last chromosome
				}
	
				#generally, you plot so that opacity is 1/(number files underneath +1). In this case we mess with base opacity, which helps to deal with very different numbers of SNPs between species
				op=30000/length(ggen$genindex) 	#base level of alpha-opacity calculated as 50K divided by the number of data points to achieve approximate equality of intensity among species
	
				#make output file
				png(file=paste(t, modl, ".", substr(rfile,1,5), ".EnrichHighlight.",bmin,".",bmax,de,".png", sep=""),width=2000,height=1000,res=300,type="quartz",bg="transparent")
				plot(ggen$genindex, ggen$sum, type="p", pch=16, col=rgb(red=182/255, green=109/255, blue=255/255, alpha=(op/1)),cex=0.25,ann=FALSE,xaxt='n',main=modl) #gen=purple (8) http://mkweb.bcgsc.ca/biovis2012/color-blindness-palette.png

				#plot set 1, scatter
				points(ggeo$genindex, y=ggeo$sum, type="p", pch=16, col=rgb(red=255/255, green=255/255, blue=109/255, alpha=(op/2)),cex=0.25,ann=FALSE,xaxt='n') #geo=yellow (15)
				points(gpco$genindex, y=gpco$sum, type="p", pch=16, col=rgb(red=0/255, green=146/255, blue=146/255, alpha=(op/3)),cex=0.25,ann=FALSE,xaxt='n') #pco=blue (3)
				points(gmodl$genindex, y=gmodl$sum, type="p", pch=16, col=rgb(red=219/255, green=109/255, blue=0/255, alpha=(op/4)),cex=0.25,ann=FALSE,xaxt='n') #modl=orange (13)

				#plot set 1, enlarged/marked points
				#toggle below to plot white x's at enrichment level for targeted genes using Rgen
				#points(tomarkgen$genindex, y=tomarkgen$sum, type="p", pch=4, col="white",cex=0.25,ann=FALSE,xaxt='n',main=modl) #tomarkgen=white 'x' (8) http://mkweb.bcgsc.ca/biovis2012/color-blindness-palette.png
				points(tomarkmodl$genindex, y=tomarkmodl$sum, type="p", pch=16, col=rgb(red=0/255, green=0/255, blue=0/255, alpha=(0.5)),cex=0.25,ann=FALSE,xaxt='n') #tomarkmodl=black (1)

				Arrows(x=max(ggen$genindex)+1,y=mean(ggen$sum),x1=max(ggen$genindex),y1=mean(ggen$sum), col=rgb(red=182/255, green=109/255, blue=255/255), lwd=0.5, arr.adj=1, arr.length=0.2) #put a purple arrowhead at mean of RgenTgen enrichment, right side
				Arrows(x=min(ggeo$genindex)-1,y=mean(ggeo$sum),x1=min(ggeo$genindex),y1=mean(ggeo$sum), col=rgb(red=255/255, green=255/255, blue=109/255), lwd=0.5, arr.adj=1, arr.length=0.2) #put a yellow arrowhead at mean of RgeoTgen enrichment, left side
				Arrows(x=max(gpco$genindex)+1,y=mean(gpco$sum),x1=max(gpco$genindex),y1=mean(gpco$sum), col=rgb(red=0/255, green=146/255, blue=146/255), lwd=0.5, arr.adj=1, arr.length=0.2) #put a blue arrowhead at mean of RpcoTgen enrichment, right side
				Arrows(x=min(gmodl$genindex)-1,y=mean(gmodl$sum),x1=min(gmodl$genindex),y1=mean(gmodl$sum), col=rgb(red=219/255, green=109/255, blue=0/255), lwd=0.5, arr.adj=1, arr.length=0.2) #put an orange arrowhead at mean of RmodlTgen enrichment, left side
				abline(v=chrbreaks[-length(chrbreaks)],col="grey59",lwd=0.5) #put lines at each chromosome end, except the last chromosome
				abline(h=0,col="black",lty="solid",lwd=0.25) #put line at zero, the random expectation

				abline(h=mean(tomarkgen$sum),col=rgb(red=182/255, green=109/255, blue=255/255),lty="32",lwd=0.75) #put a purple dotted line at mean of sites in LigninRegionsToPipe.txt enriched using Rgen
				abline(h=mean(tomarkgeo$sum),col=rgb(red=255/255, green=255/255, blue=109/255),lty="32",lwd=0.75) #put a yellow dotted line at mean of sites in LigninRegionsToPipe.txt enriched using Rgeo
				abline(h=mean(tomarkpco$sum),col=rgb(red=0/255, green=146/255, blue=146/255),lty="32",lwd=0.75) #put a blue dotted line at mean of sites in LigninRegionsToPipe.txt enriched using Rpco
				abline(h=mean(tomarkmodl$sum),col="black",lty="32",lwd=0.75) #put a black dotted line at mean of sites in LigninRegionsToPipe.txt enriched using R*modl*

				dev.off()

				#print out some summary stats
				cat(modl, ".", substr(rfile,1,5), " mean(ggen$sum) ", mean(ggen$sum),"\n",sep="")
				cat(modl, ".", substr(rfile,1,5), " mean(tomarkgen$sum) ", mean(tomarkgen$sum),"\n",sep="")
				cat(modl, ".", substr(rfile,1,5), " mean(gmodl$sum) ", mean(gmodl$sum),"\n",sep="")
				cat(modl, ".", substr(rfile,1,5), " mean(tomarkmodl$sum) ", mean(tomarkmodl$sum),"\n",sep="")
				cat(modl, ".", substr(rfile,1,5), " mean(ggeo$sum) ", mean(ggeo$sum),"\n",sep="")
				cat(modl, ".", substr(rfile,1,5), " mean(tomarkgeo$sum) ", mean(tomarkgeo$sum),"\n",sep="")
				cat(modl, ".", substr(rfile,1,5), " mean(gpco$sum) ", mean(gpco$sum),"\n",sep="")
				cat(modl, ".", substr(rfile,1,5), " mean(tomarkpco$sum) ", mean(tomarkpco$sum),"\n",sep="")
			
				#save some summary stats
				mz=paste(modl, ".", substr(rfile,1,5),sep="") #model x ref id string
				ofd=data.frame(list(mz,"ggen$sum",mean(ggen$sum),sd(ggen$sum),length(ggen$sum)),stringsAsFactors=FALSE)
				names(ofd)=c("modl","ref","mean","sd","n")
				ofd[nrow(ofd) + 1,] = list(mz,"tomarkgen$sum",mean(tomarkgen$sum),sd(tomarkgen$sum),length(tomarkgen$sum))
				ofd[nrow(ofd) + 1,] = list(mz,"gmodl$sum",mean(gmodl$sum),sd(gmodl$sum),length(gmodl$sum))
				ofd[nrow(ofd) + 1,] = list(mz,"tomarkmodl$sum",mean(tomarkmodl$sum),sd(tomarkmodl$sum),length(tomarkmodl$sum))
				ofd[nrow(ofd) + 1,] = list(mz,"ggeo$sum",mean(ggeo$sum),sd(ggeo$sum),length(ggeo$sum))
				ofd[nrow(ofd) + 1,] = list(mz,"tomarkgeo$sum",mean(tomarkgeo$sum),sd(tomarkgeo$sum),length(tomarkgeo$sum))
				ofd[nrow(ofd) + 1,] = list(mz,"gpco$sum",mean(gpco$sum),sd(gpco$sum),length(gpco$sum))
				ofd[nrow(ofd) + 1,] = list(mz,"tomarkpco$sum",mean(tomarkpco$sum),sd(tomarkpco$sum),length(tomarkpco$sum))
				ofdall=dplyr::bind_rows(ofdall,ofd)

			} #modl in modllist
		} #rfile in rfilelist

		#write out summary stats
		summaryoutputfile=paste(t,"SummaryStats",bmin,"_",bmax,de,"_",format(Sys.time(), "%m.%d.%Y_%H.%M.%S"),".to",tbl,".txt",sep="")
		ofdall2 = data.frame(lapply(ofdall, trimws), stringsAsFactors = FALSE) #remove whitespace from entries
		ofdall2$mean=as.numeric(ofdall2$mean) #restore numeric type
		ofdall2$sd=as.numeric(ofdall2$sd) #restore numeric type
		write.csv(format(ofdall2,digits=4,nsmall=4), file=summaryoutputfile, row.names=FALSE, quote=FALSE)

	} #k in dirlist

	#continue the same Rscript by summarizing the summary file with plots of mean enrichment +- 1 sd
	d = read.csv(file=summaryoutputfile, header=TRUE, sep=",")
	d2 = subset(d, modl %in% c("lig.Ligni","lsa.LSA50","regulat.regul")) #remove the Rgen controls
	d3 = subset(d2, ref %in% c("tomarkmodl$sum","tomarkgen$sum","tomarkgeo$sum","tomarkpco$sum")) #remove ggen$sum,ggeo$sum,gpco$sum since they are the same across all
	d3$ref = factor(d3$ref, levels = c("tomarkgen$sum","tomarkmodl$sum","tomarkgeo$sum","tomarkpco$sum"))

	#make output file
	png(file=paste(t,"SummaryStats",bmin,"_",bmax,de,".to",tbl,".png",sep=""),width=2000,height=1000,res=300,type="quartz",bg="transparent")
	ggplot(d3, aes(x=modl, y=mean, fill=ref)) + 
		geom_bar(position=position_dodge(),stat="identity") +
		geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
					  width=0,
					  position=position_dodge(0.9)) +
		labs(title=paste(bmin,"_",bmax,sep=""))
	dev.off()

	#plot after normalizing to account for different number of block lengths in, say 2_8 vs 2_50
	nrm=bmax/2 #number of block lengths considered for summary
	d3$meannrm=d3$mean/nrm

	png(file=paste(t,"SummaryStats",bmin,"_",bmax,de,".to",tbl,"NORM.png",sep=""),width=2000,height=1000,res=300,type="quartz",bg="transparent")
	ggplot(d3, aes(x=modl, y=meannrm, fill=ref)) + 
		geom_bar(position=position_dodge(),stat="identity") +
		labs(title=paste(bmin,"_",bmax,"NORM", sep=""))
	dev.off()
}

#####ENDR#####

#use some bash to get the data that is useful out of the summary files
typ="AnalysesForRanges"; #AnalysesForRanges AnalysesFor1Length, specifies type of analysis in output file

#AnalysesForRanges
#cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/ComparisonFigs/AnalysesForRanges;
cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/ComparisonFigs/AnalysesForRanges/TestGoLigCeresRegulatBlip;
#cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/ComparisonFigs/AnalysesForRanges/TestLigBlip;
a=$(find 2_*/bestEnrich -name "SorSummaryStats2_*_*.*.*_*.*.*.txt" | sort -t_ -k2,2n); #for AnalysesForRanges

#AnalysesFor1Length
#cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/ComparisonFigs/AnalysesFor1Length;
#cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/ComparisonFigs/AnalysesFor1Length/TestLigCeres;
#cd /Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/ComparisonFigs/AnalysesFor1Length/TestLigBlip;
#a=$(find */bestEnrich -name "SorSummaryStats*_*_*.*.*_*.*.*.txt" | sort -t_ -k2,2n);


#step 1: Enrichment at loci of interest
#extract useful data to chart, e.g. effect of Rgen/Rgeo/Rpco/Rlig on enrichment at lig sites. These are the 'tomark' lines.
#normalize relative to max block length (i.e. the number of block lengths considered for a given max block length, ie 1,2,3,12,25)
b=$(for i in $a;
  do bb=$(echo "$i" | rev | cut -d'/' -f1 | rev | cut -d'_' -f2); #extract max block length (2,4,8,24,50)
    if [[ "$typ" == "AnalysesForRanges" ]]; then
      cc=$(echo "$bb/2" | bc); #calculate number of block lengths for max block length
      grep -v "^gen" "$i" | grep -v "modl,ref,mean,sd,n" | grep "tomark" | sed 's/^/2_'$bb',/g' | awk -F, -v cc=$cc '{print $0 "," ($4/cc)}'; #remove 'gen' lines, add a header row containing the max block length, calcluate a 7th row (mean normalized to maxbl)
    else
      grep -v "^gen" "$i" | grep -v "modl,ref,mean,sd,n" | grep "tomark" | sed 's/^/2_'$bb',/g' | awk -F, '{print $0 "," $4}'; #remove 'gen' lines, add a header row containing the max block length, calcluate a 7th row (repeat of mean, not normalized because only one block length is used)
    fi;
  done;)
b=$(echo "$b" | sed 's/ //g'); #remove randomly distributed spaces left over from R

#baseline the data to Rgen
f=$(for ii in $(seq 2 2 50);
  do  i="2_"$ii;
    c=$(echo "$b" | grep "^$i,"); #get data for each max block size
    for j in regulat lig lsa;
      do d=$(echo "$c" | grep ",$j\."); #get data for each subset of reference genes
         e=$(echo "$d" | awk -F, '$3~/tomarkgen/{print $7}'); #get the mean normalized enrichment value for Rgen, which will be used as the baseline value
         echo "$d" | awk -F, -v e=$e '{print $0 "," ($7-e)}'; #print the baselined column as $7
      done;
  done;)

f=$(echo "$f" | sort -t_ -k2,2n | sort -t, --stable -k3,3 -k2,2 -k1,1n); #sort data sensibly
#echo "$f" > "/Volumes/J22/M+ Study/Analysis/Final analysis for 3rd paper/ComparisonFigs/+Summary/EffectOfMaxBlockLengthToMark."$typ".txt";
echo "$f" > "/Volumes/J22/M+ Study/Analysis/Final analysis for 3rd paper/ComparisonFigs/+Summary/repeated test analyses/EffectOfMaxBlockLengthToMark."$typ".TestGoLigCeresRegulatBlip.txt";


#step 2: Enrichment across genome from using loci of interest as reference
#extract useful data to chart gmodl$sum vs ggen$sum over varying maxblock lengths.  This is enrichment across whole genome
#caused by using either the whole genome or the loci of interest as the Reference
#normalize relative to max block length (or, the number of block lengths considered for a given max block length, ie 1,2,3,12,25)
b=$(for i in $a;
  do bb=$(echo "$i" | rev | cut -d'/' -f1 | rev | cut -d'_' -f2); #extract max block length (2,4,6,...,24,50)
    if [[ "$typ" == "AnalysesForRanges" ]]; then
      cc=$(echo "$bb/2" | bc); #calculate number of block lengths for max block length
      grep -v "^gen" "$i" | grep -v "modl,ref,mean,sd,n" | grep -v "tomark" | sed 's/^/2_'$bb',/g' | awk -F, -v cc=$cc '{print $0 "," ($4/cc)}'; #remove 'gen' lines, add a header row containing the max block length, calcluate a 7th row (mean normalized to maxbl)
    else
      grep -v "^gen" "$i" | grep -v "modl,ref,mean,sd,n" | grep -v "tomark" | sed 's/^/2_'$bb',/g' | awk -F, '{print $0 "," $4}'; #remove 'gen' lines, add a header row containing the max block length, calcluate a 7th row (repeat of mean, not normalized because only one block length is used)
    fi;
  done;)
b=$(echo "$b" | sed 's/ //g'); #remove randomly distributed spaces left over from R

#baseline the data to Rgen
f=$(for ii in $(seq 2 2 50);
  do i="2_"$ii;
    c=$(echo "$b" | grep "^$i,"); #get data for each max block size
    for j in regulat lig lsa;
      do d=$(echo "$c" | grep ",$j\."); #get data for each subset of reference genes
         e=$(echo "$d" | awk -F, '$3~/ggen/{print $7}'); #get the mean normalized enrichment value for Rgen, which will be used as the baseline value
         echo "$d" | awk -F, -v e=$e '{print $0 "," ($7-e)}'; #print the baselined column as $7
      done;
  done;)

f=$(echo "$f" | sort -t_ -k2,2n | sort -t, --stable -k3,3 -k2,2 -k1,1n); #sort data sensibly
#echo "$f" > "/Volumes/J22/M+ Study/Analysis/Final analysis for 3rd paper/ComparisonFigs/+Summary/EffectOfMaxBlockLengthGgenVGmodl."$typ".txt";
echo "$f" > "/Volumes/J22/M+ Study/Analysis/Final analysis for 3rd paper/ComparisonFigs/+Summary/repeated test analyses/EffectOfMaxBlockLengthGgenVGmodl."$typ".TestGoLigCeresRegulatBlip.txt";


#Plot these in Excel in a file called /Volumes/J22/M+ Study/Analysis/Final analysis for 3rd paper/ComparisonFigs/+Summary/EffectOfMaxBlockLength.xlsx

#Calculate pairwise t-tests using R for test of enrichment at target loci using various reference data
#For example, how well does whole genome optimization predict enrichment at target gene sets?
#Or, how well does lignin biosynthetic gene set optimization predict enrichment at lignin biosynthetic genes?
#Or, how well does environmental data optimization predict enrichment at target gene sets?
#####BEGINR#####
#install.packages("BSDA")
options(error = recover)
rm(list=ls()) 
#setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/ComparisonFigs/+Summary/")
setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/ComparisonFigs/+Summary/repeated\ test\ analyses")
library(BSDA)
library(ggplot2)

atype="AnalysesForRanges" #"AnalysesForRanges", "AnalysesFor1Length"
d = read.csv(file=paste("EffectOfMaxBlockLengthToMark.",atype,".TestLigCeres.txt",sep=""), header=FALSE, sep=",",col.names=c("max block length","modl","ref","mean","sd","n","meannrm","relativetoRgen"))
op=setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("blocklength", "modl", "p"))


for (i in c("2_2","2_4","2_6","2_8","2_10","2_12","2_14","2_16","2_18","2_20","2_22","2_24","2_26","2_28","2_30","2_32","2_34","2_36","2_38","2_40","2_42","2_44","2_46","2_48","2_50"))
{
  a=d[d$max.block.length==i,] #subset by maxblocklength
  for (j in c("lig.Ligni","lsa.LSA50","regulat.regul"))
  {
    b=a[a$modl==j,] #subset by dataset
    mx=b[b$ref=="tomarkgen$sum",'mean'] #mean of Rgen
    sx=b[b$ref=="tomarkgen$sum",'sd']
    nx=b[b$ref=="tomarkgen$sum",'n']
    my=b[b$ref=="tomarkmodl$sum",'mean']
    sy=b[b$ref=="tomarkmodl$sum",'sd']
    ny=b[b$ref=="tomarkmodl$sum",'n']
    t=tsum.test(mean.x=mx, s.x=sx, n.x=nx, mean.y=my, s.y=sy, n.y=ny, alternative="two.sided", var.equal=FALSE)
    op[nrow(op) + 1,]=list(i,j,t$p.value) #add pvalue to output data frame

    cat(i,j,t$p.value,sep=" ")
    cat("\n")
  }
}

pholm=p.adjust(op$p, method = "holm", n = length(op$p))
op[,"pholm"]=pholm
print(op)

#####ENDR#####



#Calculate pairwise t-tests using R for test of enrichment at whole genome using target loci as reference
#For example how well does climate adaptation gene set optimization predict enrichment across the whole genome?
#Target gene set optimization results in significantly less haplotype capture across the
#whole genome than whole genome optimization, p-holm <<0.001 for all tests.
#See Excel worksheet EffectOfMaxBlockLength for table of p-values.
#####BEGINR#####
#install.packages("BSDA")
options(error = recover)
rm(list=ls()) 
#setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/ComparisonFigs/+Summary/")
setwd("/Volumes/J22/M+\ Study/Analysis/Final\ analysis\ for\ 3rd\ paper/ComparisonFigs/+Summary/repeated\ test\ analyses")
library(BSDA)
library(ggplot2)

atype="AnalysesForRanges" #"AnalysesForRanges", "AnalysesFor1Length"
d = read.csv(file=paste("EffectOfMaxBlockLengthGgenVGmodl.",atype,".TestLigCeres.txt",sep=""), header=FALSE, sep=",",col.names=c("max block length","modl","ref","mean","sd","n","meannrm","relativetoRgen"))

###compare 3 target sets with gen##
cat("Rmodl vs ggen\n")
op=NULL #initialize
op=setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("blocklength", "modl", "mean.x", "mean.y", "p"))

for (i in c("2_2","2_4","2_6","2_8","2_10","2_12","2_14","2_16","2_18","2_20","2_22","2_24","2_26","2_28","2_30","2_32","2_34","2_36","2_38","2_40","2_42","2_44","2_46","2_48","2_50"))
{
  a=d[d$max.block.length==i,] #subset by maxblocklength
  for (j in c("lig.Ligni","lsa.LSA50","regulat.regul"))
  {
    b=a[a$modl==j,] #subset by dataset
    mx=b[b$ref=="ggen$sum",'mean'] #mean of Tgen
    sx=b[b$ref=="ggen$sum",'sd']
    nx=b[b$ref=="ggen$sum",'n']
    my=b[b$ref=="gmodl$sum",'mean'] #mean of Tmodl
    sy=b[b$ref=="gmodl$sum",'sd']
    ny=b[b$ref=="gmodl$sum",'n']
    t=tsum.test(mean.x=mx, s.x=sx, n.x=nx, mean.y=my, s.y=sy, n.y=ny, alternative="two.sided", var.equal=FALSE)
    op[nrow(op) + 1,]=list(i,j,mx,my,t$p.value) #add pvalue to output data frame

    cat(i,j,mx,my,t$p.value,sep=" ")
    cat("\n")
  }
}

pholm=p.adjust(op$p, method = "holm", n = length(op$p))
op[,"pholm"]=pholm
print(op)


###compare 3 target sets with geo##
#Does optimization of target gene set diversity capture whole genome diversity better than
#geographic data? Yes, p-holm <<0.001 for all tests.
cat("RmodlTgen vs ggeo")
op=NULL #initialize
op=setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("blocklength", "modl", "mean.x", "mean.y", "p"))

for (i in c("2_2","2_4","2_6","2_8","2_10","2_12","2_14","2_16","2_18","2_20","2_22","2_24","2_26","2_28","2_30","2_32","2_34","2_36","2_38","2_40","2_42","2_44","2_46","2_48","2_50"))
{
  a=d[d$max.block.length==i,] #subset by maxblocklength
  for (j in c("lig.Ligni","lsa.LSA50","regulat.regul"))
  {
    b=a[a$modl==j,] #subset by dataset
    mx=b[b$ref=="ggeo$sum",'mean'] #mean of Tgeo
    sx=b[b$ref=="ggeo$sum",'sd']
    nx=b[b$ref=="ggeo$sum",'n']
    my=b[b$ref=="gmodl$sum",'mean'] #mean of Tmodl
    sy=b[b$ref=="gmodl$sum",'sd']
    ny=b[b$ref=="gmodl$sum",'n']
    t=tsum.test(mean.x=mx, s.x=sx, n.x=nx, mean.y=my, s.y=sy, n.y=ny, alternative="two.sided", var.equal=FALSE)
    op[nrow(op) + 1,]=list(i,j,mx,my,t$p.value) #add pvalue to output data frame

    cat(i,j,mx,my,t$p.value,sep=" ")
    cat("\n")
  }
}

pholm=p.adjust(op$p, method = "holm", n = length(op$p))
op[,"pholm"]=pholm
print(op)


###compare 3 target sets with pco##
#Does optimization of target gene set diversity capture whole genome diversity better than
#environmental data? Yes, p-holm <<0.001 for all tests.
cat("RmodlTgen vs gpco")
op=NULL #initialize
op=setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("blocklength", "modl", "mean.x", "mean.y", "p"))

for (i in c("2_2","2_4","2_6","2_8","2_10","2_12","2_14","2_16","2_18","2_20","2_22","2_24","2_26","2_28","2_30","2_32","2_34","2_36","2_38","2_40","2_42","2_44","2_46","2_48","2_50"))
{
  a=d[d$max.block.length==i,] #subset by maxblocklength
  for (j in c("lig.Ligni","lsa.LSA50","regulat.regul"))
  {
    b=a[a$modl==j,] #subset by dataset
    mx=b[b$ref=="gpco$sum",'mean'] #mean of Tpco
    sx=b[b$ref=="gpco$sum",'sd']
    nx=b[b$ref=="gpco$sum",'n']
    my=b[b$ref=="gmodl$sum",'mean'] #mean of Tmodl
    sy=b[b$ref=="gmodl$sum",'sd']
    ny=b[b$ref=="gmodl$sum",'n']
    t=tsum.test(mean.x=mx, s.x=sx, n.x=nx, mean.y=my, s.y=sy, n.y=ny, alternative="two.sided", var.equal=FALSE)
    op[nrow(op) + 1,]=list(i,j,mx,my,t$p.value) #add pvalue to output data frame

    cat(i,j,mx,my,t$p.value,sep=" ")
    cat("\n")
  }
}

pholm=p.adjust(op$p, method = "holm", n = length(op$p))
op[,"pholm"]=pholm
print(op)

#####ENDR#####




























