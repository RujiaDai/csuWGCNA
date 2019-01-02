hpickSoftThreshold<-function (data, dataIsExpr = TRUE, RsquaredCut = 0.85, powerVector = c(seq(1, 
    10, by = 1), seq(12, 20, by = 2)), removeFirst = FALSE, nBreaks = 10, 
    blockSize = NULL, corFnc = cor, corOptions = list(use = "p"), 
    networkType = "unsigned", moreNetworkConcepts = FALSE, verbose = 0, 
    indent = 0) 
{   networkTypes<-c('unsigned','signed','signed hybird','hybrid2','distance')
    intType = charmatch(networkType, networkTypes)
    if (is.na(intType)) 
        stop(paste("Unrecognized 'networkType'. Recognized values are", 
            paste(networkTypes, collapse = ", ")))
    nGenes = ncol(data)
    if (nGenes < 3) {
        stop("The input data data contain fewer than 3 rows (nodes).", 
            "\nThis would result in a trivial correlation network.")
    }
    if (!dataIsExpr) {
        checkSimilarity(data)
        if (any(diag(data) != 1)) 
            diag(data) = 1
    }
    if (is.null(blockSize)) {
        blockSize = blockSize(nGenes, rectangularBlocks = TRUE, 
            maxMemoryAllocation = 2^30)
        if (verbose > 0) 
            printFlush(spaste("pickSoftThreshold: will use block size ", 
                blockSize, "."))
    }
    colname1 = c("Power", "SFT.R.sq", "slope", "truncated R.sq", 
        "mean(k)", "median(k)", "max(k)")
    if (moreNetworkConcepts) {
        colname1 = c(colname1, "Density", "Centralization", "Heterogeneity")
    }
    datout = data.frame(matrix(666, nrow = length(powerVector), 
        ncol = length(colname1)))
    names(datout) = colname1
    datout[, 1] = powerVector
    spaces = indentSpaces(indent)
    if (verbose > 0) {
        cat(paste(spaces, "pickSoftThreshold: calculating connectivity for given powers..."))
        if (verbose == 1) 
            pind = initProgInd()
        else cat("\n")
    }
    corFnc = match.fun(corFnc)
    corFormals = formals(corFnc)
    if ("nThreads" %in% names(corFormals)) 
        corOptions$nThreads = 1
    datk = matrix(0, nrow = nGenes, ncol = length(powerVector))
    nThreads = WGCNAnThreads()
    nPowers = length(powerVector)
    startG = 1
    while (startG <= nGenes) {
        endG = min(startG + blockSize - 1, nGenes)
        if (verbose > 1) 
            printFlush(paste(spaces, "  ..working on genes", 
                startG, "through", endG, "of", nGenes))
        nBlockGenes = endG - startG + 1
        jobs = allocateJobs(nBlockGenes, nThreads)
        actualThreads = which(sapply(jobs, length) > 0)
        datk[c(startG:endG), ] = foreach(t = actualThreads, .combine = rbind) %dopar% 
            {
                useGenes = c(startG:endG)[jobs[[t]]]
                nGenes1 = length(useGenes)
                if (dataIsExpr) {
                  corOptions$x = data
                  corOptions$y = data[, useGenes]
                  corx = do.call(corFnc, corOptions)
                  if (intType == 1) {
                    corx = abs(corx)
                  }
                  else if (intType == 2) {
                    corx = (1 + corx)/2
                  }
                  else if (intType == 3) {
                    corx[corx < 0] = 0
                  }
                  else if (intType == 4) {
                    corx = (1 + abs(corx))/2
                  }

                  if (sum(is.na(corx)) != 0) 
                    warning(paste("Some correlations are NA in block", 
                      startG, ":", endG, "."))
                }
                else {
                  corx = data[, useGenes]
                }
                datk.local = matrix(NA, nGenes1, nPowers)
                for (j in 1:nPowers) {
                  datk.local[, j] = colSums(corx^powerVector[j], 
                    na.rm = TRUE) - 1
                }
                datk.local
            }
        startG = endG + 1
        if (verbose == 1) 
            pind = updateProgInd(endG/nGenes, pind)
    }
    if (verbose == 1) 
        printFlush("")
    for (i in c(1:length(powerVector))) {
        khelp = datk[, i]
        SFT1 = scaleFreeFitIndex(k = khelp, nBreaks = nBreaks, 
            removeFirst = removeFirst)
        datout[i, 2] = SFT1$Rsquared.SFT
        datout[i, 3] = SFT1$slope.SFT
        datout[i, 4] = SFT1$truncatedExponentialAdjRsquared
        datout[i, 5] = mean(khelp, na.rm = TRUE)
        datout[i, 6] = median(khelp, na.rm = TRUE)
        datout[i, 7] = max(khelp, na.rm = TRUE)
        if (moreNetworkConcepts) {
            Density = sum(khelp)/(nGenes * (nGenes - 1))
            datout[i, 8] = Density
            Centralization = nGenes * (max(khelp) - mean(khelp))/((nGenes - 
                1) * (nGenes - 2))
            datout[i, 9] = Centralization
            Heterogeneity = sqrt(nGenes * sum(khelp^2)/sum(khelp)^2 - 
                1)
            datout[i, 10] = Heterogeneity
        }
    }
    print(signif(data.frame(datout), 3))
    ind1 = datout[, 2] > RsquaredCut
    indcut = NA
    indcut = if (sum(ind1) > 0) 
        min(c(1:length(ind1))[ind1])
    else indcut
    powerEstimate = powerVector[indcut][[1]]
    list(powerEstimate = powerEstimate, fitIndices = data.frame(datout))
}

