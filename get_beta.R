#!/usr/bin/env Rscript

inputs<-commandArgs(T)

if(length(inputs) < 2)
{
	message("
	Usage: get_beta.R <mset-file> <outfile-base>

	This programs reads an mSet R object and output gzipped files
	of beta values, in both raw and quantile-normalized form.

	e.g.: get_beta.R GSE40279.mset.rda GSE40279.beta

	Note: the loaded R object of mSet must be named as 'mSet'.

			")
	q("no")
}


library(minfi)
outbase<-inputs[2]
load(inputs[1])

# raw beta
message("Getting raw beta\n")
betas<-getBeta(mSet)
outfile<-paste(outbase,"raw.csv.gz", sep=".")
zz<-gzfile(outfile,"w")
write.csv(format(betas,digits=4, scientific=F),file=zz,quote=F)
close(zz)

# quantile normalized
message("Getting quantile-normalized beta\n")
betas<-getBeta(preprocessQuantile(mSet))
outfile<-paste(outbase,"quant.csv.gz", sep=".")
zz<-gzfile(outfile,"w")
write.csv(format(betas,digits=4, scientific=F),file=zz,quote=F)
close(zz)

message("Job done")

