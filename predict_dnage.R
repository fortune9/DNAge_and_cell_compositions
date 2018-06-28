#!/usr/bin/env Rscript

# this is the program to compute DNAge using Horvath or Hannum model

## functions
fill_NA<-function(val,replacement) # the two vectors must be in the same order
{
	stopifnot(length(val)==length(replacement))
	naIndex<-which(is.na(val))
	val[naIndex]<-replacement[naIndex]
	return(val)
}

# read the model data for Horvath clock
read_horvath<-function(f)
{
	model<-read.delim(f, skip=2, stringsAsFactors=F)
	# select needed columns only
	selected<-c("CpGmarker","CoefficientTraining","CoefficientTrainingShrunk","medianByCpG")
	model<-model[selected]
	names(model)<-c("cpg","coef","coefshrunk","medianByCpG")
	intercept<-unlist(model[1,c("coef","coefshrunk")])
	model<-model[-1,]
	return(list(dataFrame=model,intercept=intercept))
}

# read model data for Hannum clock
read_hannum<-function(f)
{
	model<-read.delim(f, stringsAsFactors=F)
	return(model)
}

# the help functions to convert DNAge and computed value in Horvath model,
# as the model was trained in a non-linear scale.
raw_to_dnage<-function(rawAge,adult_age=20)
{
	if(rawAge < 0) { return((1+adult_age)*exp(rawAge)-1) }
	else { return((1+adult_age)*rawAge+adult_age) }
}

dnage_to_raw<-function(dnage, adult_age=20)
{
	rawAge<-(dnage+1)/(1+adult_age)
	if(rawAge <= 1) {rawAge<-log(rawAge)} else { rawAge<-rawAge -1 }
	return(rawAge)
}

# select the beta rows according to the clock CpGs, and fill NAs if necessary
select_rows<-function(betas,clockCpGs)
{
	missingCpGs<-setdiff(clockCpGs,rownames(betas))
	if(length(missingCpGs) < 1 ) { # no missing CpGs
		return(betas[clockCpGs,]) 
	}
	# otherwise fill the missing CpGs with NAs
	warning(length(missingCpGs), " necessary CpGs are missing in the data\n")
	missingRows<-matrix(nr=length(missingCpGs),nc=ncol(betas),dimnames=list(missingCpGs,colnames(betas)))
	selectedBetas<-betas[intersect(clockCpGs,rownames(betas)),]
	selectedBetas<-rbind(selectedBetas,missingRows)
	selectedBetas<-selectedBetas[clockCpGs,] # reorder rows, VERY important

	return(selectedBetas)
}

run_horvath_model<-function(betas,coefField="coef")
{
	dat<-clockModel$dataFrame[!is.na(clockModel$dataFrame[,coefField]),] # remove NA rows
	selectedBetas<-select_rows(betas,dat[,"cpg"])
	# fill NA values in input data
	selectedBetas<-apply(selectedBetas,2,fill_NA,replacement=dat[,"medianByCpG"])
	stopifnot(all(rownames(selectedBetas)==dat[,"cpg"])) # the cpgs
	# must match
	rawAge<-crossprod(selectedBetas, dat[,coefField])+clockModel$intercept[coefField]
	dnage<-apply(rawAge,1, raw_to_dnage)
	return(dnage)
}

run_hannum_model<-function(betas)
{
	selectedBetas<-select_rows(betas,clockModel[,"cpg"]) # construct a matrix with NAs for missing CpGs
	stopifnot(all(rownames(selectedBetas)==clockModel[,"cpg"])) # the cpgs
	# must match
	selectedBetas<-apply(selectedBetas, 2, fill_NA, replacement=clockModel[,"medianBeta"])
	dnage<-crossprod(selectedBetas, clockModel[,"coef"])
	return(dnage[,1])
}

parse_params<-function(f)
{
     requiredParams<-c("betaFile","horvathFile","hannumFile","outFile","fieldSep","type")
     dat<-read.delim(f, head=F, stringsAsFactors=F, blank.lines.skip=T, comment.char = "#")
     params<-list()
     tmp<-lapply(seq_len(nrow(dat)), function(i) {r<-unlist(dat[i,]); params[[tolower(r[1])]]<<-r[2]} )
	 if(!is.null(params$fieldsep) && params$fieldsep=="\\t") { params$fieldsep<-"\t" }
	 if(is.null(params$fieldsep)) { params$fieldsep<-"," }
	 if(is.null(params$type)) { stop("The parameter 'type' is required") } else { params$type<-tolower(params$type) }
	 return(params)
}

usage<-function()
{
	message("
	Usage: predict_dnage.R <config-file>
	
	This program predicts DNAge based on DNA methylation data using
	the Horvath and Hannum clock models.
	
	The <config-file> contains the input parameters, with 2 columns
	(separated by a tab): 
	the first column is the parameter names and the second presents
	the values. Available parameters are as follows:
	
	type: which clock to run, Horvath or Hannum.
	
	betaFile:  the file containing beta values, one sample per column,
	and the first line should contain one field fewer than the rest of
	lines. From the second line, the first field are the CpG names. 
	The file can be gzipped.
	
	fieldSep:  the field separator for the betaFile. Default is comma.
	
	modelFile:  the Horvath or Hannum clock model file. Example files
	are at github
	
	outFile:  the file to store results.
	
	");
}

inputs<-commandArgs(T)

if(length(inputs) < 1)
{
	usage()
	q("no")
}

params<-parse_params(inputs[1])

if(params$type == "hannum" || params$type == "horvath")
{
	read_method<-get(paste("read", params$type, sep="_"))
	run_method<-get(paste("run", params$type, "model", sep="_"))
}else
{
	stop("Unknown clock type: ", params$type)
}
# setup

message("# Reading beta values")
#betaFile<-"GSE67751.beta.csv.gz"
betas<-NULL
if(grepl("\\.gz$",params$betafile,ignore.case=T))
{
	betas<-read.table(gzfile(params$betafile),sep=params$fieldsep)
}else
{
	betas<-read.table(params$betafile, sep=params$fieldsep)
}

message("# Computing DNAge")

# prepare the models
# horvathModel<-read_horvath(params$horvathfile)
# hannumModel<-read_hannum(params$hannumfile)
clockModel<-read_method(params$modelfile)

# dnage.horvath<-run_horvath_model(betas)
# dnage.hannum<-run_hannum_model(betas)
dnage<-run_method(betas)

results<-data.frame(names(dnage),dnage)
colnames(results)<-c("sample", paste("dnage",params$type,sep="."))
write.table(results, file=params$outfile, quote=F, sep=",", row=F)

message("Job done")
