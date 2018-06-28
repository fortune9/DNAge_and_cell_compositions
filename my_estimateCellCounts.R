library("genefilter")
library("quadprog")

my_estimateCellCounts<-function(Dset,...) # an interface function to both RGset and Mset data
{
	if(is(Dset, "RGChannelSet")) { return(my_estimateCellCounts_RG(Dset,...))}
	if(is(Dset, "MethylSet")) { return(my_estimateCellCounts_M(Dset,...))}
	stop(sprintf("Unknown class type %s for input\n", class(Dset)))
}

# the argument probeNum specify the half number of total probes for each cell type
my_estimateCellCounts_M <- function (mSet, compositeCellType = "Blood", processMethod = "auto", probeSelect = "auto", probeNum=50,
                                cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Gran"),
                                referencePlatform = c("IlluminaHumanMethylation450k", "IlluminaHumanMethylationEPIC", "IlluminaHumanMethylation27k"),
                                returnAll = FALSE, withPlot = FALSE, verbose=TRUE, ...) {
    
    #.isRGOrStop(rgSet)
    mSet <- as(mSet, "MethylSet")
    referencePlatform <- match.arg(referencePlatform)
    mPlatform <- sub("IlluminaHumanMethylation", "", annotation(mSet)[which(names(annotation(mSet))=="array")])
    platform <- sub("IlluminaHumanMethylation", "", referencePlatform) # reference platform
    if((compositeCellType == "CordBlood") && (!"nRBC" %in% cellTypes))
        message("[estimateCellCounts_M] Consider including 'nRBC' in argument 'cellTypes' for cord blood estimation.\n")   
    referencePkg <- sprintf("FlowSorted.%s.%s", compositeCellType, platform)
    subverbose <- max(as.integer(verbose) - 1L, 0L)
    if(!require(referencePkg, character.only = TRUE))
        stop(sprintf("Could not find reference data package for compositeCellType '%s' and referencePlatform '%s' (inferred package name is '%s')",
                     compositeCellType, platform, referencePkg))
    data(list = referencePkg) 
    referenceRGset <- get(referencePkg)
    if(mPlatform != platform) {
        mSet <- convertArray(mSet, outType = referencePlatform, verbose = subverbose)
    }
    if(! "CellType" %in% names(pData(referenceRGset)))
        stop(sprintf("the reference sorted dataset (in this case '%s') needs to have a phenoData column called 'CellType'"),
             names(referencePkg))
    if(sum(colnames(mSet) %in% colnames(referenceRGset)) > 0)
        stop("the sample/column names in the user set must not be in the reference data ")
    if(!all(cellTypes %in% referenceRGset$CellType))
        stop(sprintf("all elements of argument 'cellTypes' needs to be part of the reference phenoData columns 'CellType' (containg the following elements: '%s')",
                     paste(unique(referenceRGset$cellType), collapse = "', '")))
    if(length(unique(cellTypes)) < 2)
        stop("At least 2 cell types must be provided.")
    #if ((processMethod == "auto") && (compositeCellType %in% c("Blood", "DLPFC")))
        processMethod <- "preprocessQuantile" # for mSet, it can only be this method
    #if ((processMethod == "auto") && (!compositeCellType %in% c("Blood", "DLPFC")))
    #    processMethod <- "preprocessNoob"
    processMethod <- get(processMethod)
    if ((probeSelect == "auto") && (compositeCellType == "CordBlood")){
        probeSelect <- "any"} 
    if ((probeSelect == "auto") && (compositeCellType != "CordBlood")){
        probeSelect <- "both"}
    
    if(verbose) message("[estimateCellCounts_M] Combining user data with reference (flow sorted) data.\n")
	## here need convert reference dataset into MethylSet first and then combine
	refMset<-preprocessRaw(referenceRGset)
	
    newpd <- DataFrame(sampleNames = c(colnames(mSet), colnames(referenceRGset)),
                       studyIndex = rep(c("user", "reference"),
                                        times = c(ncol(mSet), ncol(referenceRGset))))
                       #stringsAsFactors = FALSE)
    referencePd <- pData(referenceRGset)
	userPd<-pData(mSet)
    combinedMset <- combineArrays(mSet, refMset, outType = "IlluminaHumanMethylation450k")
    pData(combinedMset) <- newpd
    colnames(combinedMset) <- newpd$sampleNames
	slot(combinedMset, "preprocessMethod")<-preprocessMethod(refMset) # this is needed for normalization
    rm(referenceRGset, refMset, mSet)
    
    if(verbose) message("[estimateCellCounts_M] Processing user and reference data together.\n")
    if (compositeCellType == "CordBlood"){
        ## Here Shan wants to discard probes that they have decided shouldn't be used, for example multi-mapping probes
        ## This is done by only using probes with names in the comptable.
        ## This is kind of ugly, and dataset dependent.
        combinedMset <- processMethod(combinedRGset, verbose=subverbose)
        compTable <- get(paste0(referencePkg, ".compTable"))
        combinedMset <- combinedMset[which(rownames(combinedMset) %in% rownames(compTable)),]
    } else {
        combinedMset <- processMethod(combinedMset) 
    }
    # rm(combinedRGset)
    
    ## Extracts normalized reference data and user data
    referenceMset <- combinedMset[, combinedMset$studyIndex == "reference"]
    pData(referenceMset) <- as(referencePd, "DataFrame")
    mSet <- combinedMset[, combinedMset$studyIndex == "user"]
    pData(mSet) <- as(userPd, "DataFrame")
    rm(combinedMset)
    
	## train model with reference data 
    if(verbose) message("[estimateCellCounts_M] Picking probes for composition estimation.\n")
    compData <- pickCompProbes(referenceMset, cellTypes = cellTypes, compositeCellType = compositeCellType, probeSelect = probeSelect, numProbes=probeNum)
    coefs <- compData$coefEsts
    if(!returnAll) { rm(referenceMset) }
    
	## estimate cell fractions on the real data
    if(verbose) message("[estimateCellCounts_M] Estimating composition.\n")
    counts <- projectCellType(getBeta(mSet)[rownames(coefs), ], coefs)
    rownames(counts) <- colnames(mSet)
    
    if (returnAll || withPlot) {
		if(verbose) message("[estimateCellCounts_RG] Generating mean_plot parameters.\n")
		smeans <- compData$sampleMeans
		smeans <- smeans[order(names(smeans))]
		coefs <- compData$coefEsts
		sampleMeans <- colMeans(getBeta(mSet)[rownames(coefs), ])

		plotParams<-list(smeans=smeans, sampleMeans=sampleMeans, sampleType=compositeCellType);
		if(returnAll) { list(counts = counts, compTable = compData$compTable,normalizedData = mSet, refData=referenceMset, plot=plotParams, coefs=coefs)  # returnAll has higher priority
		}else         { list(counts=counts, plot=plotParams) }
    } else {
        counts
    }
}

my_estimateCellCounts_RG <- function (rgSet, compositeCellType = "Blood", processMethod = "auto", probeSelect = "auto", probeNum=50,
                                cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Gran"),
                                referencePlatform = c("IlluminaHumanMethylation450k", "IlluminaHumanMethylationEPIC", "IlluminaHumanMethylation27k"),
                                returnAll = FALSE, meanPlot = FALSE, verbose=TRUE, withPlot=FALSE, ...) {
    
    #.isRGOrStop(rgSet)
    rgSet <- as(rgSet, "RGChannelSet")
    referencePlatform <- match.arg(referencePlatform)
    rgPlatform <- sub("IlluminaHumanMethylation", "", annotation(rgSet)[which(names(annotation(rgSet))=="array")])
    platform <- sub("IlluminaHumanMethylation", "", referencePlatform)
    if((compositeCellType == "CordBlood") && (!"nRBC" %in% cellTypes))
        message("[estimateCellCounts_RG] Consider including 'nRBC' in argument 'cellTypes' for cord blood estimation.\n")   
    referencePkg <- sprintf("FlowSorted.%s.%s", compositeCellType, platform)
    subverbose <- max(as.integer(verbose) - 1L, 0L)
    if(!require(referencePkg, character.only = TRUE))
        stop(sprintf("Could not find reference data package for compositeCellType '%s' and referencePlatform '%s' (inferred package name is '%s')",
                     compositeCellType, platform, referencePkg))
    data(list = referencePkg) 
    referenceRGset <- get(referencePkg)
    if(rgPlatform != platform) {
        rgSet <- convertArray(rgSet, outType = referencePlatform, verbose = subverbose)
    }
    if(! "CellType" %in% names(pData(referenceRGset)))
        stop(sprintf("the reference sorted dataset (in this case '%s') needs to have a phenoData column called 'CellType'"),
             names(referencePkg))
    if(sum(colnames(rgSet) %in% colnames(referenceRGset)) > 0)
        stop("the sample/column names in the user set must not be in the reference data ")
    if(!all(cellTypes %in% referenceRGset$CellType))
        stop(sprintf("all elements of argument 'cellTypes' needs to be part of the reference phenoData columns 'CellType' (containg the following elements: '%s')",
                     paste(unique(referenceRGset$cellType), collapse = "', '")))
    if(length(unique(cellTypes)) < 2)
        stop("At least 2 cell types must be provided.")
    if ((processMethod == "auto") && (compositeCellType %in% c("Blood", "DLPFC")))
        processMethod <- "preprocessQuantile"
    if ((processMethod == "auto") && (!compositeCellType %in% c("Blood", "DLPFC")))
        processMethod <- "preprocessNoob"
    processMethod <- get(processMethod)
    if ((probeSelect == "auto") && (compositeCellType == "CordBlood")){
        probeSelect <- "any"} 
    if ((probeSelect == "auto") && (compositeCellType != "CordBlood")){
        probeSelect <- "both"}
    
    if(verbose) message("[estimateCellCounts_RG] Combining user data with reference (flow sorted) data.\n")
    newpd <- DataFrame(sampleNames = c(colnames(rgSet), colnames(referenceRGset)),
                       studyIndex = rep(c("user", "reference"),
                                        times = c(ncol(rgSet), ncol(referenceRGset))))
                       #stringsAsFactors = FALSE)
    referencePd <- pData(referenceRGset)
    combinedRGset <- combineArrays(rgSet, referenceRGset, outType = "IlluminaHumanMethylation450k")
    pData(combinedRGset) <- newpd
    colnames(combinedRGset) <- newpd$sampleNames
    rm(referenceMset)
    
    if(verbose) message("[estimateCellCounts_RG] Processing user and reference data together.\n")
    if (compositeCellType == "CordBlood"){
        ## Here Shan wants to discard probes that they have decided shouldn't be used, for example multi-mapping probes
        ## This is done by only using probes with names in the comptable.
        ## This is kind of ugly, and dataset dependent.
        combinedMset <- processMethod(combinedRGset, verbose=subverbose)
        compTable <- get(paste0(referencePkg, ".compTable"))
        combinedMset <- combinedMset[which(rownames(combinedMset) %in% rownames(compTable)),]
    } else {
        combinedMset <- processMethod(combinedRGset) 
    }
    rm(combinedRGset)
    
	#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	# start from here, the analyses are the same between RGChannelSet and MethylSet
	#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    ## Extracts normalized reference data 
    referenceMset <- combinedMset[, combinedMset$studyIndex == "reference"]
    pData(referenceMset) <- as(referencePd, "DataFrame")
    mSet <- combinedMset[, combinedMset$studyIndex == "user"]
    pData(mSet) <- as(pData(rgSet), "DataFrame")
    rm(combinedMset)
    
    if(verbose) message("[estimateCellCounts_RG] Picking probes for composition estimation.\n")
    compData <- pickCompProbes(referenceMset, cellTypes = cellTypes, compositeCellType = compositeCellType, probeSelect = probeSelect, numProbes=probeNum)
    coefs <- compData$coefEsts
    if(!returnAll) { rm(referenceMset) }
    
    if(verbose) message("[estimateCellCounts_RG] Estimating composition.\n")
    counts <- projectCellType(getBeta(mSet)[rownames(coefs), ], coefs)
    rownames(counts) <- colnames(rgSet)
    
	# calculate mean_plot needed parameters
    if (returnAll || withPlot) {
		if(verbose) message("[estimateCellCounts_RG] Generating mean_plot parameters.\n")
		smeans <- compData$sampleMeans
		smeans <- smeans[order(names(smeans))]
		coefs <- compData$coefEsts
		sampleMeans <- colMeans(getBeta(mSet)[rownames(coefs), ])

		plotParams<-list(smeans=smeans, sampleMeans=sampleMeans, sampleType=compositeCellType);
		if(returnAll) { list(counts = counts, compTable = compData$compTable,normalizedData = mSet, refData=referenceMset, plot=plotParams, coefs=coefs)  # returnAll has higher priority
		}else         { list(counts=counts, plot=plotParams) }
    } else {
        counts
    }
}

pickCompProbes <- function(mSet, cellTypes = NULL, numProbes = 50, compositeCellType = compositeCellType, probeSelect = probeSelect) {
    splitit <- function(x) {
        split(seq(along=x), x)
    }
    p <- getBeta(mSet)
    pd <- as.data.frame(pData(mSet))
    if(!is.null(cellTypes)) {
        if(!all(cellTypes %in% pd$CellType))
            stop("elements of argument 'cellTypes' is not part of 'mSet$CellType'")
        keep <- which(pd$CellType %in% cellTypes)
        pd <- pd[keep,]
        p <- p[,keep]
    }
    ## make cell type a factor 
    pd$CellType <- factor(pd$CellType, levels = cellTypes)
    ffComp <- rowFtests(p, pd$CellType)
    prof <- sapply(splitit(pd$CellType), function(i) rowMeans(p[,i]))
    r <- matrixStats::rowRanges(p)
    compTable <- cbind(ffComp, prof, r, abs(r[,1] - r[,2]))
    names(compTable)[1] <- "Fstat"
    names(compTable)[c(-2,-1,0) + ncol(compTable)] <- c("low", "high", "range") 
    tIndexes <- splitit(pd$CellType)
	# compare beta values between each cell type and the others
    tstatList <- lapply(tIndexes, function(i) {
        x <- rep(0,ncol(p))
        x[i] <- 1
        return(rowttests(p, factor(x)))
    })
    # selecting the cell-type specific probes
    if (probeSelect == "any"){
        probeList <- lapply(tstatList, function(x) {
            y <- x[x[,"p.value"] < 1e-8,]
            yAny <- y[order(abs(y[,"dm"]), decreasing=TRUE),]      
            c(rownames(yAny)[1:(numProbes*2)])
        })
    } else {
        probeList <- lapply(tstatList, function(x) {
            y <- x[x[,"p.value"] < 1e-8,]
            yUp <- y[order(y[,"dm"], decreasing=TRUE),]
            yDown <- y[order(y[,"dm"], decreasing=FALSE),]
            c(rownames(yUp)[1:numProbes], rownames(yDown)[1:numProbes])
        })
    }
    
    trainingProbes <- unique(unlist(probeList))
    p <- p[trainingProbes,]
    
    pMeans <- colMeans(p)
    names(pMeans) <- pd$CellType
	
    form <- as.formula(sprintf("y ~ %s - 1", paste(levels(pd$CellType), collapse="+")))
    phenoDF <- as.data.frame(model.matrix(~pd$CellType-1))
    colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
	# essentially, the fitting below is get the meam beta value for each CpG in each cell type
    if(ncol(phenoDF) == 2) { # two group solution
        X <- as.matrix(phenoDF)
        coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
    } else { # > 2 group solution
        tmp <- validationCellType(Y = p, pheno = phenoDF, modelFix = form)
        coefEsts <- tmp$coefEsts
    }
    
    out <- list(coefEsts = coefEsts, compTable = compTable,
                sampleMeans = pMeans, csProbes = probeList) # now cell-type specific probes are also returned
    return(out)
}

projectCellType <- function(Y, coefCellType, contrastCellType=NULL, nonnegative=TRUE, lessThanOne=FALSE){ 
    if(is.null(contrastCellType))
        Xmat <- coefCellType
    else
        Xmat <- tcrossprod(coefCellType, contrastCellType) 
    
    nCol <- dim(Xmat)[2]
    if(nCol == 2) { # only 2 cell types
        Dmat <- crossprod(Xmat)
        mixCoef <- t(apply(Y, 2, function(x) { solve(Dmat, crossprod(Xmat, x)) }))
        colnames(mixCoef) <- colnames(Xmat)
        return(mixCoef)
    } else { # more than 2 cell types
        nSubj <- dim(Y)[2]
        
        mixCoef <- matrix(0, nSubj, nCol)
        rownames(mixCoef) <- colnames(Y) # samples
        colnames(mixCoef) <- colnames(Xmat) # cell types
        
        if(nonnegative){
            if(lessThanOne) {
                Amat <- cbind(rep(-1, nCol), diag(nCol))
                b0vec <- c(-1, rep(0, nCol))
            } else {
                Amat <- diag(nCol)
                b0vec <- rep(0, nCol)
            }
            for(i in 1:nSubj) { # loop for each subject
                obs <- which(!is.na(Y[,i])) 
                Dmat <- crossprod(Xmat[obs,])
                mixCoef[i,] <- solve.QP(Dmat, crossprod(Xmat[obs,], Y[obs,i]), Amat, b0vec)$sol # there is a constraint here, so QP is needed
            }
        } else {
            for(i in 1:nSubj) {
                obs <- which(!is.na(Y[,i])) 
                Dmat <- crossprod(Xmat[obs,])
                mixCoef[i,] <- solve(Dmat, t(Xmat[obs,]) %*% Y[obs,i]) # the vector can be negative, so simple solve() method is used.
            }
        }
        return(mixCoef)
    }
}

validationCellType <- function(Y, pheno, modelFix, modelBatch=NULL,
                               L.forFstat = NULL, verbose = FALSE){
    N <- dim(pheno)[1]
    pheno$y <- rep(0, N)
    xTest <- model.matrix(modelFix, pheno)
    sizeModel <- dim(xTest)[2]
    M <- dim(Y)[1]
    
    if(is.null(L.forFstat)) {
        L.forFstat <- diag(sizeModel)[-1,] # All non-intercept coefficients
        colnames(L.forFstat) <- colnames(xTest) 
        rownames(L.forFstat) <- colnames(xTest)[-1] 
    }

    ## Initialize various containers
    sigmaResid <- sigmaIcept <- nObserved <- nClusters <- Fstat <- rep(NA, M)
    coefEsts <- matrix(NA, M, sizeModel)
    coefVcovs <- list()

    if(verbose)
        cat("[validationCellType] ")
    for(j in 1:M) { # For each CpG
        ## Remove missing methylation values
        ii <- !is.na(Y[j,])
        nObserved[j] <- sum(ii)
        pheno$y <- Y[j,]
        
        if(j%%round(M/10)==0 && verbose)
            cat(".") # Report progress
        
        try({ # Try to fit a mixed model to adjust for plate
            if(!is.null(modelBatch)) {
                fit <- try(lme(modelFix, random=modelBatch, data=pheno[ii,]))
                OLS <- inherits(fit,"try-error") # If LME can't be fit, just use OLS
            } else
                OLS <- TRUE

            if(OLS) {
                fit <- lm(modelFix, data=pheno[ii,])
                fitCoef <- fit$coef
                sigmaResid[j] <- summary(fit)$sigma
                sigmaIcept[j] <- 0
                nClusters[j] <- 0
            } else { 
                fitCoef <- fit$coef$fixed
                sigmaResid[j] <- fit$sigma
                sigmaIcept[j] <- sqrt(getVarCov(fit)[1])
                nClusters[j] <- length(fit$coef$random[[1]])
            }
            coefEsts[j,] <- fitCoef
            coefVcovs[[j]] <- vcov(fit)
            
            useCoef <- L.forFstat %*% fitCoef
            useV <- L.forFstat %*% coefVcovs[[j]] %*% t(L.forFstat)
            Fstat[j] <- (t(useCoef) %*% solve(useV, useCoef))/sizeModel
        })
    }
    if(verbose)
        cat(" done\n")
    ## Name the rows so that they can be easily matched to the target data set
    rownames(coefEsts) <- rownames(Y)
    colnames(coefEsts) <- names(fitCoef)
    degFree <- nObserved - nClusters - sizeModel + 1

    ## Get P values corresponding to F statistics
    Pval <- 1-pf(Fstat, sizeModel, degFree)
    
    out <- list(coefEsts=coefEsts, coefVcovs=coefVcovs, modelFix=modelFix, modelBatch=modelBatch,
                sigmaIcept=sigmaIcept, sigmaResid=sigmaResid, L.forFstat=L.forFstat, Pval=Pval,
                orderFstat=order(-Fstat), Fstat=Fstat, nClusters=nClusters, nObserved=nObserved,
                degFree=degFree)
    
    out
}

# my functions
my_cor<-function(x,y,asString=F,minSize=5,...)
{
	if(length(x) < minSize) {res<-as.numeric(c(NA,NA));names(res)<-c("r","P"); return(ifelse(asString, "atop(rho==NA, italic(P)==NA)", c(res)))}
	r_test<-cor.test(x,y,...)
	if(!asString) {
	res<-c(r_test$est, r_test$p.val); names(res)<-c("r","P"); 
	return(res)
	}
	#values<-formatC(c(r_test$est, r_test$p.val), digits=2, format="g")
	lab<-sprintf("atop(rho==%.2g, italic(P)==%.3g)", r_test$est, r_test$p.val)
	return(lab)
	#return(as.character(parse(text=lab)))
}

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

qc_plots<-function(Mset)
{
	if(is(Mset,"RGChannelSet")) {Mset<-preprocessRaw(Mset)}
	
	qc<-getQC(Mset)
	#layout(matrix(1:2,nc=1))
	cat("  \n")
	plotQC(qc)
	cat("  \n")
	#densityPlot(Mset, sampGroups=pData(Mset)$sex)
	densityBeanPlot(Mset, sampGroups=pData(Mset)$sex)
	#controlStripPlot(RGset, controls="BISULFITE CONVERSION I")
	#controlStripPlot(RGset, controls="BISULFITE CONVERSION II")
	cat("  \n")
	return(invisible(Mset))
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

print_sex<-function(mSet, sexCol)
{
	cat("Predicted sexes:\n\n")
	gMset<-mapToGenome(mSet)
	predictedSex<-getSex(gMset, cutoff=-2)
	sample_name<-pData(mSet)$Sample_Name;
	if(is.factor(sample_name)) { sample_name<-levels(sample_name)[sample_name]; }
	obs<-pData(mSet)[[sexCol]]
	if(is.factor(obs)) { obs<-levels(obs)[obs] }
	tmp.sex<-cbind(sample_name, predicted=predictedSex$predictedSex,obs)
	print(kable(tmp.sex, caption="Predicted genders", format="markdown", align="lcc"))
	#print(pander(t(tmp.sex), caption="Predicted genders", justify="center", split.tables=80, style="markdown")) # failed
	cat("\n")
	rm(gMset)
}

mean_plot<-function(smeans, sampleMeans, sampleType, pchs=21, ...)
{    
	cellTypeNames<-factor(names(smeans))
    sampleColors <- c(rep(1, length(sampleMeans)), 1 + as.numeric(cellTypeNames))
	sampleMeans <-c(sampleMeans, smeans)
    plot(sampleMeans, pch = pchs, bg = sampleColors, ...)
    legend("bottomleft", c(sampleType, levels(cellTypeNames)),
           col = 1:(nlevels(cellTypeNames)+1), pch = 15)
}

# this function sum up the cell fractions of sub-types into one value, such as T, B, and NK cells into lymphocytes
collapse_cell_types<-function(mat, selected="all", verbose=F) # selected gives which cell type groups to be collapsed
{
	cellTypeGroup<-list(
		Lymph = c("Bcell", "CD4T", "CD8T","NK"),
		Gran  = c("Neutro", "Eosino", "Baso"),
		Tcell = c("CD4T","CD8T")
		)
	selected<-tolower(selected)
	# filter groups based on the parameter 'selected'
	if((length(selected)>1) || (selected != "all"))
	{
		selected<-sub("(.)","\\U\\1",selected, perl=T)
		cellTypeGroup<-lapply(selected, function(x) cellTypeGroup[[x]])
		names(cellTypeGroup)<-selected
	}
	keepCols<-!(colnames(mat) %in% unlist(cellTypeGroup)) # these columns won't be affected
	groups<-sapply(cellTypeGroup, function(x) any(x %in% colnames(mat)) )
	groups<-names(cellTypeGroup)[groups] # cell type groups having at least one cell type in the input data
	if(verbose) {message(sprintf("The following cell type groups are found in input data: %s \n", groups))}
	collapsedVals<-sapply(groups, function(x) {commCols<-intersect(colnames(mat),cellTypeGroup[[x]]); rowSums(mat[, commCols,drop=F],na.rm=T)} )
	results<-cbind(mat[,keepCols], collapsedVals)
	return(invisible(results))
}
