#! /is1/data/users/plapie01/software/R/R-3.0.1/bin/Rscript
#To use: Rscript Cloaked-bear.R countFile.txt annotation.gff OuputDirectoryName
#Output directory will be created if it does not exist.
#wILL OVERWRITE PLOTS AND FILES IF NEEDED.
#Will run with 4 threads by default. This can be changed by modifing the numThreads value
#This analysis is based on http://cgrlucb.wikispaces.com/edgeR+spring2013
#Most parts are a copy of what can be found on the website...
#There is no guarantee that this code will provide good results.

#To use: change the TRUE and FALSE switch to obtain the desired result...
#Add a few small lines of code where specified (between the ! lines).
#######################################################################################
##                                                                                   ##
##      YOU ARE RESPONSIBLE OF YOUR ANALYSIS. THIS CODE IS ONLY A TEMPLATE.          ##
##                                                                                   ##
#######################################################################################



library(Rsubread)
library(RColorBrewer)
library(edgeR)
library(gplots);

#Using more thread is not recommended...
numThreads=1;
numOfGene=500; #for heatMap
colors=brewer.pal(9, "Set1");

args=commandArgs(trailingOnly=TRUE)
outputDirectory=args[length(args)];
sampleInfoFile=args[(length(args)-2)]
countFile=args[length(args)-1];

#print(annotationFile);
dir.create(outputDirectory);

#Will write a Log file and print message to screen. Should be in the outputDir...
#sink(file="SessionLog.txt", append=FALSE, split=TRUE);

#lets start to work with some data!
#count=featureCounts(files=fileList,nthreads=numThreads, annot=annotationFile, isGTFAnnotationFile=TRUE, GTF.attrType="gene_name", isPairedEnd=TRUE, );
#countTable=count$count;
gene=read.table(countFile);

##############################
# select the needed column
# i.e: countTable=countLevel[,1:10];
# if you are not sure, try to do it interactively and copy your code here!
#New: comment the one you need!
#################################################
#gene$sample_AGM1.21.MED.sam=NULL
#gene$sample_AGM1.22.LPS.sam=NULL
#gene$sample_AGM1.23.CLO.sam=NULL
#gene$sample_AGM2.26.MED.sam=NULL
#gene$sample_AGM2.27.LPS.sam=NULL
#gene$sample_AGM2.28.CLO.sam=NULL
#gene$sample_AGM3.31.MED.sam=NULL
#gene$sample_AGM3.32.LPS.sam=NULL
#gene$sample_AGM3.33.CLO.sam=NULL
#gene$sample_AGM5.CA1I1.CLO.sam=NULL
#gene$sample_AGM5.CA1I1.LPS.sam=NULL
#gene$sample_AGM5.CA1I1.MED.sam=NULL
#gene$sample_AGM6.CA1I2.CLO.sam =NULL
#gene$sample_AGM6.CA1I2.LPS.sam=NULL
#gene$sample_AGM6.CA1I2.MED.sam=NULL
#gene$sample_CH1.I4.CLO.sam=NULL
#gene$sample_CH1.I4.LPS.sam=NULL
#gene$sample_CH1.I4.MED.sam =NULL
#gene$sample_CH2.I5.CLO.sam=NULL
#gene$sample_CH2.I5.LPS.sam=NULL
#gene$sample_CH2.I5.MED.sam=NULL
#gene$sample_IND1.1.MED.sam=NULL
#gene$sample_IND1.2.LPS.sam=NULL
#gene$sample_IND1.3.CLO.sam=NULL
#gene$sample_IND2.6.MED.sam=NULL
#gene$sample_IND2.7.LPS.sam=NULL
#gene$sample_IND2.8.CLO.sam=NULL
#gene$sample_IND3.M1I1.CLO.sam=NULL
#gene$sample_IND3.M1I1.LPS.sam=NULL
#gene$sample_IND3.M1I1.MED.sam=NULL
#gene$sample_IND4.16.MED.sam=NULL
#gene$sample_IND4.17.LPS.sam=NULL
#gene$sample_IND4.18.CLO.sam =NULL
#gene$sample_IND4.M1I3.CLO.sam=NULL
#gene$sample_IND4.M1I3.LPS.sam =NULL
#gene$sample_IND4.M1I3.MED.sam =NULL
#gene$sample_RM3.36.MED.sam=NULL
#gene$sample_RM3.37.LPS.sam=NULL
#gene$sample_RM3.38.CLO.sam =NULL
#gene$sample_RM4.41.MED.sam =NULL
#gene$sample_RM4.42.LPS.sam =NULL
#gene$sample_RM4.43.CLO.sam=NULL





##################################################################
countTable=gene;



means=rowMeans(countTable);
filter=means>=10;

#set to true to print a table showing the number of exon with an expression level >10
if(FALSE){
	table(filter);
}

geneCounts=countTable[filter,];
#Now that files are read, we can change the wd.
#setwd(outputDirectory);

##Uncomment and adapt if graph is wanted. A new file describing sampleprep and condition could be nedded...
if(FALSE){
	totCounts <- colSums(geneLevelCounts)
	##laneInfo is a table containing the samples description.Change the ...
	laneInfo=...
	barplot(totCounts, las=2, col=colors[laneInfo[,2]])
	barplot(totCounts, las=2, col=colors[laneInfo[,4]])
}


#BoxPlot of difference in distribution between each experiment.Set to true if wanted...
if(FALSE){
	png("BoxPlot_DistributionPerExperiment.png");
	#Need to add flafla to graph
	boxplot(log2(geneCounts+1), las=2);
	dev.off();
}

#The next section use EdgeR. Documentation: http://www.bioconductor.org/packages/2.12/bioc/manuals/edgeR/man/edgeR.pdf
#There is no best way to analyse RNA-seq because every experiment is different.
# I recommand to carefully look at the code bellow.


#######################################################################################
##This section is for a two component analysis (2 groups, 2 samples or 2 conditions) ##
#######################################################################################
#Set to true if needed. Set to false otherwise.
#You must create a table containing only the necessary samples and a group variable
#the group variable should look like this for 6 samples splitted in 2 groups(Del and Gly)
# group
#[1] Del Del Del Gly Gly Gly
#dataToCreateDGE is a selection of column from the original table.
# CHange the ... below
if(FALSE){
	group=c("MED","MED","MED","MED","MED","CLO","CLO","CLO","CLO","CLO");
	dataToCreateDGE=geneCounts;

	cds=DGEList(dataToCreateDGE, group=group);
	#If true will print the number of having 0 counts across all samples.
	if(TRUE){
		print("0 counts across all samples");
		sum(cds$all.zeros)
	}
	#The normalisation method use, possibility: "TMM", "RLE", "upperquartile", "non"
	normalisationMethod="TMM"
	cds=calcNormFactors(cds, method=normalisationMethod);
	
	#pseudo-count creation for cluster analysis. Needed for heat-map...
	if(TRUE){
		scale=cds$samples$lib.size*cds$samples$norm.factors;
		normCounts=round(t(t(dataToCreateDGE)/scale)*mean(scale));
		#BoxPLot of the normalized pseudo-count
		if(TRUE){
			png("BoxPlot_DistributionNormlizedPseudoCounts.png");
			boxplot(log2(normCounts+1), las=2); #could add color
			#need to add flafla to graph
			dev.off();
		}
	}
	#An MDS plot to measures the similarity of the samples. Usefull for QC and sample visualization.
	#Does not work...
	if(FALSE){
		png("MDSplot_ForCountData.png");
		plotMDS(cds, main="MDS Plot for normalized count data", labels=colnames(cds$counts));
		dev.off();
		#flafla...
	}
	#We need to estimate the common dispersion: this assumes that all the genes have the same dispersion.
	#For more explication see http://cgrlucb.wikispaces.com/edgeR+spring2013
	cds=estimateCommonDisp(cds, verbose=TRUE);
	cds=estimateTagwiseDisp(cds);
	
	#To plot the biological coefficient of variation. Can be set to FALSE
	if(TRUE){
		png("Plot_BiologicalCoefficientOfVariation.png");
		plotBCV(cds);
		dev.off();
	}
	#To plot mean-variance relation. 
	if(TRUE){
		png("Plot_Mean-Variance.png");
		#First line is for data, second is for format. You will need to do some adjustment on the first line...
		plotMeanVar(cds, show.raw.vars=TRUE, show.tagwise.vars=TRUE, show.binned.common.disp.vars=FALSE, show.ave.raw.vars=FALSE, NBline = TRUE , nbins = 100 , 
		pch = 16 , xlab ="Mean Expression (Log10 Scale)" , ylab = "Variance (Log10 Scale)" , main = "Mean-Variance Plot" );
		dev.off();
	}
	#I'n not sure if unique is the best way to solve the problem...
	et <- exactTest(cds, pair=unique(group))
	topTags(et)
	top <- topTags(et, n=nrow(cds$counts))$table
	head(top)
	de <- rownames(top[top$FDR<0.1,])
	print("Length DE");
	length(de)
	head(de)
	if(TRUE){
		png("Histogram_distributionPvalue.png");
		hist(top$PValue, breaks=20)
		dev.off()
	}
	if(FALSE){
		png("PlotSmear_mean-diffrence.png");
		plotSmear(cds , de.tags=de)
		abline(h=c(-2, 2), col=colors[2])
	}
	if(TRUE){
		png("VolcanoPlot_logFoldChangesAndPvalues.png");
		plot(top$logFC, -log10(top$PValue), pch=20, cex=.5, ylab="-log10(p-value)", xlab="logFC", col=as.numeric(rownames(top) %in% de)+1)
		abline(v=c(-2, 2), col=colors[2])
		dev.off()
	}
	if(TRUE){
		if(length(de)<numOfGene){
			numOfGene=length(de);
		}
		pdf("HeatMap_topGenes.pdf", height=15);
		#fpar(font.axis="Courier");
		#The heat map really need optimisation...
		#colors <- colorRampPalette(c("white","green","green4","violet","purple"))(100)
		colors <- colorRampPalette(c("blue", "green", "yellow", "orange", "red"))(100)
		heatmap.2(log2(normCounts[de[1:numOfGene],]+1), col=colors, mar=c(20,20), scale="none",trace="none", cexRow=0.4 ,cexCol=0.6,  main="Heat-map Log2")
		#heatmap(log(normCounts[de[1:numOfGene],]+1), col=colors, margins=c(15,20), scale="none", 
		dev.off();
	}
	if(TRUE){
		write.table(top, file="two-class-results.txt", sep='\t', quote=FALSE)
	}
}

#############################################################################
#3 components comparison using the Generalized Linear Model (GLM) approach. #
#############################################################################
if(TRUE){
	sampleInfo=read.delim(sampleInfoFile);
	rownames(sampleInfo)=sampleInfo$Sample;

	#Change group
	group=sampleInfo[,2];
	#For GLM comparison based on monkey type
	setwd(outputDirectory);
	if(TRUE){
		design=model.matrix(~sampleInfo$TypeSinge)
	}
	#For GLM comparison based on treatement
	if(FALSE){
		design=model.matrix(~sampleInfo$Traitement)
	}
	#For GLM comparison one to one considering the other... Not sure to get this.
	if(FALSE){
		typeSinge=factor(paste(sampleInfo$TypeSinge,sep="_"));
		singe=factor(sampleInfo$TypeSinge); #Based on Charles date, we are blocking on the monkeys.
		traitement=factor(paste(sampleInfo$TypeSinge,sep="_"));
	}
	
	cds= DGEList(geneCounts, group=group);
	
	design=model.matrix(~group)
	#Ste to TRUE to see the design matrix.
	if(FALSE){
		design
	}
	#The normalisation method use, possibility: "TMM", "RLE", "upperquartile", "non"
        normalisationMethod="TMM"
	cds=calcNormFactors(cds, method=normalisationMethod);
	
	#For pseudo-counts. Needed for heat-map
	if(TRUE){
		scale=cds$samples$lib.size*cds$samples$norm.factors
		#normCounts round(t(t(counts)/scale)*mean(scale))
		normCounts=round(t(t(geneCounts)/scale)*mean(scale))
		#Boxplot- distribution of pseudo counts. Not sure if log or log2
		if(TRUE){
			png("Boxplot_GLM_DistributionNormlizedPseudoCounts.png")
			boxplot(log2(normCounts+1), las=2, col=colors[geneCounts$samples$group])
			#flafla...
			dev.off();
		}
	}
	#MDS plot->similarity of the samples. For QC and samples vizualisation
	if(FALSE){
		png("MDSplot_GLM_PlotForCOuntData.png");
		plotMDS(cds,main-"MDS Plot for Count Data", labels=colnames(cds$counts));
		dev.off();
	}
	
	#dispersion estimation
	cds=estimateGLMCommonDisp(cds, design, verbose=TRUE);
	cds=estimateGLMTagwiseDisp(cds, design);
	#Plot biological coeefficient.
	if(TRUE){
		png("Plot_GLM_BiologicalCoefficientVariation.png");
		plotBCV(cds);
		dev.off();
	}
	#Mean-variance relation plot. You might prefer to use the next plot(dglmStdResid) instead.
	#Search for dglmStdResid in edgeR doc for more explication.
	if(FALSE){
		png("Plot_GLM_Mean-variance.png");
		meanVar=plotMeanVar(cds, show.raw.vars=TRUE, show.tagwise.vars=TRUE, show.binned.common.diusp.vars=FALSE, show.ave.raw.vars=FALSE,
		NBline=TRUE, nbins=100, pch=16, xlab="Mean Expression (log10)", ylab="Variance (log10)", main="Mean-Variance Plot");
		dev.off();
	}
	#A probably more suitable mean-variance relation plot for 3 component analysis.
	#Will need deep testing...and modification
	if(FALSE){
		png("Plot_GLM_StdResidue.png");
		dglmStdResid(cds$counts, design, ...);
		dev.off();
	}
	
	#Dif. exp. You might want to change the coef value.
	fit= glmFit(cds, design);
	lrt= glmLRT(fit, coef=2:3);
	top=topTags(lrt, n=nrow(cds$counts))$table;
	
	de=rownames(top[top$FDR<0.1,]);
	#Histogram of pvalue
	if(TRUE){
		png("Histogram_GLM_Pvalue.png")
		hist(top$PValue, breaks=20);
		dev.off();
	}
	#I could add plotSmear...

	#heat-map
        if(TRUE){
                if(length(de)<numOfGene){
                        numOfGene=length(de);
                }
                pdf("HeatMap_topGenes.pdf", height=35);
                colors <- colorRampPalette(c("blue", "green", "yellow", "orange", "red"))(100)
                heatmap.2(log2(normCounts[de[1:numOfGene],]+1), col=colors, mar=c(20,20), scale="none",trace="none", cexRow=0.4 ,cexCol=0.4,  main="Heat-map Log2")
                dev.off();
        }
	#Output results table...
	if(TRUE){
		write.table(top, file="Three-component-analysis.txt", sep='\t', quote=FALSE);
	}
}

#Will print session info...
if(TRUE){
	print("This analysis was performed on:")
	sessionInfo();
}




