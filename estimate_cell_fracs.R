#!/usr/bin/env Rscript
# this program is to estimate cell fractions using Houseman algorithm on DNA methylation data.

source("../my_estimateCellCounts.R") # load the file modified from
# minfi package, and loading must be before the definitions of
# functions

# functions here suppress all text output except plots
estimate_cell_fracs<-function(rgSet, normType, funcName="my_estimateCellCounts", ...)
{
	if(any(grepl(normType, "function", ignore=T)))
	{
		normType<-"function";
		processWay<-"preprocessFunnorm"
	}else
	{
		if(any(grepl(normType, "quantile", ignore=T)))
		{
				normType<-"quantile";
				processWay<-"preprocessQuantile"
		}else
		{
			stop("Unknown normalization type\n")
		}
	}
	estimateFracs<-get(funcName)
	results<-estimateFracs(rgSet, processMethod=processWay, ...)
	cellFracs=results;
	if(is(results,"list")) {cellFracs=results$counts}
	cellFracs<-as.data.frame(cellFracs)
	cellFracs$normalizeMethod<-normType
	cellFracs$arrayId<-rownames(cellFracs)
	if(is(results,"list")) {results$counts=cellFracs} else {results=cellFracs}
	return(invisible(results))
}

qc_plots<-function(Mset,f,...)
{
	if(is(Mset,"RGChannelSet")) {Mset<-preprocessRaw(Mset)}
	
	qc<-getQC(Mset)
	png(f,...) # set the canvas
	par(mar=newmar)
	par(mgp=newmgp)
	layout(matrix(1:2,nc=1))
	plotQC(qc)
	#densityPlot(Mset, sampGroups=pData(Mset)$sex)
	densityBeanPlot(Mset, sampGroups=pData(Mset)$sex)
	#controlStripPlot(RGset, controls="BISULFITE CONVERSION I")
	#controlStripPlot(RGset, controls="BISULFITE CONVERSION II")
	dev.off()
	return(invisible(Mset))
}

make_mean_plot<-function(mp,f,...)
{
	png(f,...)
	par(mar=newmar)
	par(mgp=newmgp)
	mean_plot(mp$smeans, mp$sampleMeans, mp$sampleType) # call the function in the source file
	dev.off()
}

scatter_plot_frac<-function(ctNames,d) # cell type names and the data frame
{
	nc<-3
	nr<-ceiling(length(ctNames)/nc)
	layout(matrix(seq_len(nc*nr), nc=nc,byrow=T))
	cat("  \n"); # prevent figure float
	for(ct in ctNames)
	{
		tmp.pred<-paste(ct,"pred",sep=".")
		with(d, plot(get(ct),get(tmp.pred), xlab="Observed", ylab="Predicted", cex.main=0.9, cex.axis=0.7, cex.lab=0.8, cex=0.6, main=ct))
		abline(0,1,lwd=0.8,col="blue")
		tmp.test<-with(d, cor.test(get(ct),get(tmp.pred)))
		tmp.R<-formatC(tmp.test$est, digits=2, format="g")
		tmp.P<-formatC(tmp.test$p.val, digits=2, format="g")
		tmp.ranges<-par("usr")
		text(mean(tmp.ranges[1:2]), tmp.ranges[4], labels=bquote(atop(R==.(tmp.R),italic(P)==.(tmp.P))), cex=0.9, pos=1)
		#readline("Enter to continue")
	}
	cat("  \n")
	cat("The solid line in each plot marks the spots when prediction matches observation exactly\n\n")
}

process_intensities<-function(datFile, infoDat, fieldSep="\t", pCutoff=NA)
{
	# uncompress the data beforehand, as fread failed due to memory limit
	if(grepl("\\.gz",datFile, ignore.case=T))
	{
		tmpFile<-paste("tmp_intensity", Sys.getpid(),"tsv", sep=".")
		system(paste("zcat",datFile,">",tmpFile), intern=F, wait=T)
		datFile<-tmpFile
	}
	intensity<-fread(datFile, sep=fieldSep)
	# get sample names from the data
	unmethCols<-grep("Unmethylate", colnames(intensity),ignore=T, value=T)
	methCols<-sub("Unmeth","Meth",unmethCols, ignore.case=T)
	stopifnot(length(unmethCols) > 0)
	sampleTitles<-sub("( .*|methylated)$","",methCols, ignore.case=T)
	## order the infoDat to make sure the sample order correspond those in intensity file
	if("sampleOrder" %in% colnames(infoDat)) { infoDat<-infoDat[order(infoDat$sampleOrder),] }
	infoDat$sampleTitle<-sampleTitles
	failedCpGs<-NULL
	if(!is.na(pCutoff)) {
		pvalCols<-grep("Pval", colnames(intensity),ignore=T,value=T)
		if(length(pvalCols) == length(unmethCols)) {
			failedRows<-apply(intensity[,..pvalCols],1, function(x) any(x > pCutoff)) # detection P value is too large
			#intensity<-intensity[!failedRows,] # good rows only
			failedCpGs<-intensity[[1]][failedRows] # mark failed rows only
			}else {
				warning("No detection P values are found, or P value
						columns don't match methylation columns")
			}
		}

# also get the reference set
## require(FlowSorted.Blood.450k)
## refRGset<-get("FlowSorted.Blood.450k")
## refMset<-preprocessRaw(refRGset)

	# construct the mSet object from real data
	sampleNames<-levels(infoDat$Sample_Name)[infoDat$Sample_Name]
	rownames(infoDat)<-sampleNames
	tmp.phenoData<-as(infoDat, "AnnotatedDataFrame")
	tmp.meth<-as.matrix(setnames(intensity[,..methCols],sampleNames));
	tmp.unmeth<-as.matrix(setnames(intensity[,..unmethCols],sampleNames));
	# storage.mode(tmp.meth)<-storage.mode(getMeth(refMset))
	# storage.mode(tmp.unmeth)<-storage.mode(getMeth(refMset))
	rownames(tmp.meth)<-intensity[[1]] # data.table is also a list
	rownames(tmp.unmeth)<-intensity[[1]]
	rm(intensity)
	mSet<-MethylSet(tmp.meth, tmp.unmeth, phenoData=tmp.phenoData)
	rm(tmp.meth, tmp.unmeth)
	#featureNames(mSet)<-intensity[,ID_REF]
	annotation(mSet)<-c(array="IlluminaHumanMethylation450k", annotation="ilmn12.hg19")
	if(file.exists(tmpFile)) { file.remove(tmpFile) } # remove temporary file
	return(list(mSet, failedCpGs))
}

parse_params<-function(f)
{
	requiredParams<-c("dataFile","infoFile","fieldSep","outFile","normType")
	dat<-read.delim(f, head=F, stringsAsFactors=F)
	params<-list()
	tmp<-lapply(seq_len(nrow(dat)), function(i) {r<-unlist(dat[i,]); params[[tolower(r[1])]]<<-r[2]} )
	if(!all(tolower(requiredParams) %in% names(params))) {
		stop("some of the following mandatory parameters are missing:\n",
		paste(requiredParams,collapse="\n"),"\n")
		}
	if(is.null(params$idat)) { params$idat<-F } else { params$idat<-as.logical(toupper(params$idat)) }
	if(!is.null(params$fieldsep) && params$fieldsep=="\\t") {
		params$fieldsep<-"\t" }
	if(!is.null(params$infosep)) {
		if(params$infosep=="\\t") { params$infosep<-"\t" } 
	} else { 
		params$infosep<-","	 # comma is default	
	}
	if(is.null(params$savemset)) {
		params$savemset<-F
	}else{
		params$savemset<-as.logical(params$savemset)
	}
	return(params)
}

usage<-function()
{
	message("
Usage: estimate_cell_fracs.R <input-file>

This program estimates cell fractions in blood samples based on
DNA methylation data measured on Illumina HM450K array.

The data files and settings can be specified in the <input-file>
in the following format, using tab as field separator:

param1	val1
param2	val2
...     ...

The available params are as follows:

dataFile:
a file of probe intensities (can be gzipped) or a folder
containing *.idat files.

idat:
a logical option. If T is given, the data are *.idat files,
otherwise probe intensities are assumed. Default is F for false.

fieldSep:
the field separator for the probe intensity file, if applicable.

normType:
the normalization type for the data. Either 'quantile' or 'function'
normalization.

savemSet:
a logical option. If T is given, then the nonnormalized mSet is saved.
Default is F, so no saving. If true, failed CpGs based on detection P
values *are* also saved in a separate file.

infoFile:
a file containing the samples' information, one row
per sample. One column 'sampleOrder' can be used to match the samples
in this file and those in the 'dataFile'. The field separator is
comma in default, but changeable by the next option.

infoSep:
the field separator for the 'infoFile'. Default is comma.

outFile:
the basename for output filenames. The outputs include both text
files and figures.

figWid:
the width of each figure, in the unit of inch. Default is 5.

figHei:
the height of each figure. In default, it is calculated accordingly,
which is ussually good.

	");
}


# analysis starts here
configFile<-commandArgs(T)
if(length(configFile) < 1) # no inputs
{
	usage();
	q("no")
}
library(data.table)
library(minfi)
message("# reading data")
param<-parse_params(configFile)
figWid<-5
figHei<-figWid
ppi<-300; # fig resolution
if(!is.null(param$figwid)) {figWid<-as.numeric(param$figwid)}
if(!is.null(param$fighei)) {figHei<-as.numeric(param$figHei)}
info<-read.delim(param$infofile, sep=param$infosep) # sample information
names(info)[grepl("gender",names(info),ignore=T)]<-"sex"; # make colnames consistent
processedData<-NULL
if(param$idat) # idat data
{
	stop("not implemented iDAT data yet")
}else # signal intensity
{
	# return mSet and failed CpGs
	processedData<-process_intensities(param$datafile, info, param$fieldsep, pCutoff=0.05)
}

mSet<-processedData[[1]]
failedCpGs<-processedData[[2]]
rm(processedData)
gc()

if(param$savemset)
{
	msetFile<-paste(param$outfile,"mset", "rda", sep=".")
	failedCpGFile<-paste(param$outfile,"failedCpGs", "csv", sep=".")
	message("# Saving non-normalized and FULL mSet to ", msetFile)
	message("# and failed CpGs to ",failedCpGFile)
	# mSet is not normalized yet
	if(is.na(preprocessMethod(mSet)["rg.norm"]))
	{
		# without this, the old version minfi reports errors
		mSet@preprocessMethod <- 
			c(rg.norm = "Raw (no normalization or bg correction)",
			  minfi = as.character(packageVersion("minfi")))
	}
	save(mSet, file=msetFile)
	write.csv(as.data.frame(failedCpGs),file=failedCpGFile,quote=F,row.names=F,col.names=F)
}

# from here, the mSet contains only good CpGs
if(!is.null(failedCpGs))
{
	mSet<-mSet[!is.element(rownames(mSet),failedCpGs),]
}

message("# quality control plots")
qcPlotFile<-paste(param$outfile,"qc", "png", sep=".")
qc_plots(mSet,qcPlotFile,width=figWid,height=figHei*2,units="in",res=ppi)

message("# predict sex")
# we need check the data are fine for sex prediction
gMset<-mapToGenome(mSet)
CN<-getCN(gMset)
xIndex <- which(seqnames(gMset) == "chrX")
yIndex <- which(seqnames(gMset) == "chrY")
xMed <- matrixStats::colMedians(CN, rows = xIndex, na.rm = TRUE)
yMed <- matrixStats::colMedians(CN, rows = yIndex, na.rm = TRUE)
dd<-yMed-xMed
if(all(is.finite(dd)))
{
	predictedSex<-getSex(gMset, cutoff=-2)
	info$sex.pred<-predictedSex$predictedSex
}else
{
	message("Infinite medians found, skipping sex prediction")
}
rm(gMset, CN)
gc()

message("# estimate cell fractions")
cellPropRes<-estimate_cell_fracs(mSet, normType=param$normtype, compositeCellType="Blood", probeSelect="both", withPlot=T)
cellProps<-cellPropRes$counts
meanPlot<-cellPropRes$plot
cellTypeCols<-1:(ncol(cellProps)-2)
colnames(cellProps)[cellTypeCols]<-paste(colnames(cellProps)[cellTypeCols], "pred",sep=".")

# make plots here
message("# make mean-plots for cell type deconvolution (whole blood samples are expected to be within the vertical range of purified samples in the plot)")
meanPlotFile<-paste(param$outfile,"mplot", "png", sep=".")
make_mean_plot(meanPlot,meanPlotFile,res=ppi,units="in",width=figWid,height=figHei)


message("# write results to output")
info<-cbind(info,cellProps) # the same sample order is assumed
if(!all(with(info, Sample_Name==arrayId))) { warning("The samples from infoFile and dataFile may be mismatched") }
resultFile<-paste(param$outfile,"_blood_cell_type_deconvolution.csv",sep="")
write.table(info, file=resultFile, quote=F, row=F, sep=",")

message("** Job is done **")


