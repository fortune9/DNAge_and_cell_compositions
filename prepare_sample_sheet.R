#!/bin/env Rscript

# return terms with first-letter capitalized
# will be regarded as cell types and may be normalized by
# division with 100.
get_cell_type_name<-function(s)
{
	stopifnot(length(s)==1)
	s<-tolower(s)
	if(grepl("count",s)) { return(s) }
	if(grepl("ratio",s)) { return(s) }
	if(grepl("nadir",s)) { return(s) }
	if(grepl("low",s)) { return(s) }
	if(grepl("naive",s) || grepl("naïve",s))
	{
		if(grepl("cd4",s)) { return("cd4.naive") }
		if(grepl("cd8",s)) { return("cd8.naive") }
		return(s)
	}
	if(grepl("exhaust",s))
	{
		if(grepl("cd4",s)) { return("cd4.exhausted") }
		if(grepl("cd8",s)) { return("cd8.exhausted") }
		return(s)
	}
	if(grepl("cd8pcd28ncd45ran",s))
	{
		return("cd8pcd28ncd45ran")
	}
	
	if(grepl("cd3",s))
	{
		return("CD3T")
	}
	if(grepl("cd4",s))
	{
		return("CD4T")
	}
	if(grepl("cd8",s))
	{
		return("CD8T")
	}
	if(any(grepl("t cell",s)))
	{
		return("Tcell")
	}
	if(any(grepl("b cell",s)))
	{
		return("Bcell")
	}
	if(grepl("killer",s))
	{
		return("NK")
	}
	if(grepl("lymph",s)) # lymphocytes consist of T, B, and NK cells
	{
		return("Lymph")
	}

	if(grepl("mono",s))
	{
		return("Mono")
	}

	if(grepl("neut",s))
	{
		return("Neutro")
	}
	if(grepl("eosi",s))
	{
		return("Eosino")
	}
	if(grepl("baso",s))
	{
		return("Baso")
	}
	# granulocytes are the sum of the above three sub-types
	if(grepl("gran",s))
	{
		return("Gran")
	}

	# otherwise return the attribute value
	return(get_prefix(s))
}

get_prefix<-function(s, sep=":")
{
	s<-tolower(s)
	prefix<-sapply(strsplit(s,sep), function(x) x[1])
	prefix<-sub("^ +","",prefix)
	prefix<-sub(" +$","",prefix)
	prefix<-gsub(" +","_", prefix)
	return(prefix)
}

trim_blanks<-function(s)
{
	s<-sub("^ +","",s)
	s<-sub(" +$","",s)
	return(s)
}

# a function to process columns with different attributes mixed
process_mixed_columns<-function(s, sep=":")
{
	#s<-tolower(s)
	if(is.factor(s)) { s<-levels(s)[s] }
	prefix<-sapply(strsplit(s,sep), function(x) x[1])
	val<-sapply(strsplit(s,sep), function(x) x[2])
	val<-trim_blanks(val)
	prefix<-gsub(" +","_", trim_blanks(prefix))
	# the prefix may contain more than one value, mixing different
	# attributes to one column, so let's sep them into different
	# columns
	columnNames<-unique(prefix);
	columnNames<-columnNames[!grepl("^\\s*$",columnNames) &
							 !is.na(columnNames)] #remove empty prefix
	if(length(columnNames) > 1) { 
		message("more than one attribute [",
				paste(columnNames,collapse=","),
				"] found in one column") }
	vals<-lapply(columnNames, function(x) {
		   y<-rep(NA,length(val));
		   y[prefix %in% x]<-val[prefix %in% x];
		   return(y)
			}
		   )
	names(vals)<-columnNames
	vals<-as.data.frame(do.call(cbind, vals))
	return(vals)
}

validCellTypes<-c("Tcell","Bcell","NK","Gran","Mono","CD4T","CD8T",
				  "CD3T", "Baso","Neutro","Eosino","Lymph")

inputs<-commandArgs(T);
if(length(inputs) < 1)
{
	message("
	Usage: prepare_sample_sheet.R <*_series_matrix.txt.gz>
	[<iDAT-sample-sheet> <frac-already>]

	This program reads sample information from file
	'iDAT-sample-sheet' and parse sample information from
	'*_series_matrix.txt.gz' and merge together.

	'iDAT-sample-sheet' need be in csv format and optional. If it is
	not provided, or givn as 'NA', then only parsed information from the
	'series_matrix' file is output.

	'frac-already' is a logical value: 'T' or 'F', indicating whether
	the cell fractions are already in fraction scale [0,1] so they
	will not be further divided by 100.
	")

	q("no", status=1);
}

smatFile<-inputs[1];
sampleSheetFile<-inputs[2];
fracAlready<-as.logical(inputs[3]);
if(toupper(sampleSheetFile) == 'NA')
{
	sampleSheetFile<-NA
}
if(is.na(fracAlready))
{ fracAlready<-F }

# parse information from series matrix file and store them into a
# data.frame
tmpInfo<-read.delim(pipe(paste("zcat", smatFile, " | grep -E '^!(Sample_geo_accession|Sample_title|Sample_type|Sample_source_name|Sample_description|Sample_characteristics_ch)'")), head=F)
tmp1<-as.data.frame(t(tmpInfo[-1])) # the real data
names(tmp1)<-tmpInfo[,1]
dataCols<-grepl('Sample_characteristics_ch', names(tmp1)) # logicals

# now separate the prefix from the values in data columns,
# including cell fractions and others
dat<-tmp1[!dataCols]
tmp<-lapply(tmp1[dataCols], function(v) {
			processedCol<-process_mixed_columns(v);
			dat<<-cbind(dat,processedCol)
		   })
# rename column names
tmpNames<-sapply(colnames(dat), get_cell_type_name)
colnames(dat)<-tmpNames;

# now process cell-fraction columns which start with capital letetrs
## true cell types start with capital letters
cellTypeCols<-tmpNames[grepl('^[A-Z]',tmpNames)]

# get the numeric values of cell fractions,
# could be other information than cell type information
if(length(cellTypeCols)>0)
{
	tmpCellFrac<-apply(dat[cellTypeCols],2,function(x)
					   {if(is.factor(x)) {x<-levels(x)[x]}; 
					   return(as.numeric(x)) } )
	if(!fracAlready) { tmpCellFrac<-tmpCellFrac/100; }
	dat[cellTypeCols]<-tmpCellFrac
}

names(dat)[grep('geo_accession',names(dat))]<-'Sample_Name'
# combine the information from prepared iDAT files
if(!is.na(sampleSheetFile))
{
	sampleSheet<-read.csv(sampleSheetFile)
	tmp2<-merge(sampleSheet, dat, by="Sample_Name", all=T)
	dat<-tmp2
}
outFile<-"sample_sheet.tsv"
write.table(dat, file=outFile, quote=F, row=F, col.names=T, sep="\t")

message("Output is stored in ", outFile, "\n")

