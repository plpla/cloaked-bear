#! /is1/data/users/plapie01/software/R/R-3.0.1/bin/Rscript
#To use: Rscript Cloaked-bear.R file1 file2 file3 ... annotation.gff OuputDirectoryName
#Output directory will be created if it does not exist.
#wILL OVERWRITE PLOTS AND FILES IF NEEDED.
#Will run with 4 threads by default. This can be changed by modifing the numThreads value
#This analysis is based on http://cgrlucb.wikispaces.com/edgeR+spring2013
#Most parts are a copy of what can be found on the website...
#There is no guarantee that this code will provide good results.

#To use: change the TRUE and FALSE switch to obtain the desired result...
#Add a few small lines of code where specified (between the ! lines).
#######################################################################################
##										     ##
##	YOU ARE RESPONSIBLE OF YOUR ANALYSIS. THIS CODE IS ONLY A TEMPLATE.	     ##
##										     ##
#######################################################################################



library(Rsubread)
library(RColorBrewer)
library(edgeR)
library(gplots);

numThreads=1;
numOfGene=500; #for heatMap
colors=brewer.pal(9, "Set1");

args=commandArgs(trailingOnly=TRUE)
print(args)
outputDirectory=args[length(args)];
annotationFile=args[(length(args)-1)]
fileList=args[0:(length(args)-2)];
print(annotationFile);
dir.create(outputDirectory);

#Will write a Log file and print message to screen.
sink(file="SessionLog_AGM_CLO.1.txt", append=FALSE, split=TRUE);

#lets start to work with some data!
count=featureCounts(files=fileList,nthreads=numThreads, annot=annotationFile, isGTFAnnotationFile=TRUE, GTF.attrType="gene_name", isPairedEnd=TRUE, );
countTable=count$count;
means=rowMeans(countTable);
filter=means>=10;

#set to true to print a table showing the number of exon with an expression level >10
if(FALSE){
	table(filter);
}

geneCounts=countTable[filter,];
#Now that files are read, we can change the wd.
setwd(outputDirectory);

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
if(TRUE){
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
if(FALSE){
	#Group is a variable containing caracteristic to group samples in 3 groups. i.e. for 14 samples
	#group
	# [1] YPD YPD YPD YPD YPD YPD YPD YPD Del Del Del Gly Gly Gly
	#Levels: Del Gly YPD
	#Change the ...
	group=...
	cds= DGEList(geneCounts, group=group);
	
	design=model.matrix(~group)
	#Ste to TRUE to see the design matrix.
	if(FALSE){
		design
	}
	#The normalisation method use, possibility: "TMM", "RLE", "upperquartile", "non"
        normalisationMethod="TMM"
	cds=calcNormFactors(geneCounts, method=normalisationMethod);
	
	#For pseudo-counts. Needed for heat-map
	if(TRUE){
		scale=geneCounts$samples$lib.size*geneCounts$samples$norm.factors
		normCounts=round(t(t(counts)/scale)*mean(scale))
		#Boxplot- distribution of pseudo counts. Not sure if log or log2
		if(FALSE){
			png("Boxplot_GLM_DistributionNormlizedPseudoCounts.png")
			boxplot(log2(normCounts+1), las=2, col=colors[geneCounts$samples$group])
			#flafla...
			dev.off();
		}
	}
	#MDS plot->similarity of the samples. For QC and samples vizualisation
	if(FALSE){
		png("MDSplot_GLM_PlotForCOuntData.png");
		plotMDS(geneCounts,main-"MDS Plot for Count Data", labels=colnames(cds$counts));
		dev.off();
	}
	
	#dispersion estimation
	geneCounts=estimateGLMCommonDisp(geneCounts, design, verbose=TRUE);
	geneCounts=estimateGLMTagwiseDisp(geneCounts, design);
	#Plot biological coeefficient.
	if(FALSE){
		png("Plot_GLM_BiologicalCoefficientVariation.png");
		plotBCV(geneCounts);
		dev.off();
	}
	#Mean-variance relation plot. You might prefer to use the next plot(dglmStdResid) instead.
	#Search for dglmStdResid in edgeR doc for more explication.
	if(FALSE){
		png("Plot_GLM_Mean-variance.png");
		meanVar=plotMeanVar(geneCounts, show.raw.vars=TRUE, show.tagwise.vars=TRUE, show.binned.common.diusp.vars=FALSE, show.ave.raw.vars=FALSE,
		NBline=TRUE, nbins=100, pch=16, xlab="Mean Expression (log10)", ylab="Variance (log10)", main="Mean-Variance Plot");
		dev.off();
	}
	#A probably more suitable mean-variance relation plot for 3 component analysis.
	#Will need deep testing...and modification
	if(FALSE){
		png("Plot_GLM_StdResidue");
		dglmStdResid(geneCounts$counts, design, ...);
		dev.off();
	}
	
	#Dif. exp. You might want to change the coef value.
	fit= glmFit(geneCounts, design);
	lrt= glmLRT(fit, coef=2:3);
	top=topTags(lrt, n=nrow(geneCounts$counts))$table;
	
	de=rownames(top[top$FDR<0.1,]);
	#Histogram of pvalue
	if(TRUE){
		png("Histogram_GLM_Pvalue")
		hist(top$PValue, breaks=20);
		dev.off();
	}
	#I could add plotSmear...

	#heat-map
	if(TRUE){
		png("HeatMap_GLM_HeatMap");
		heatmap(log(normCounts[de,]+1), ColSideColor=colors[group]);
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

sink(file=NULL);



