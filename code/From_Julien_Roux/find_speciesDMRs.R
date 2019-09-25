## Feb 27, 2014
## Requires a lot of memory (>100g)
## Find species DMRs between pairs of species

library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))

## source("http://www.bioconductor.org/biocLite.R")
## biocLite()
## biocLite("BiocUpgrade")
## biocLite("bsseq")
## biocLite("bsseq", lib="~/R/x86_64-redhat-linux-gnu-library/2.15/") ## on the cluster
library(bsseq) # version 0.10
sessionInfo()

## Load all necessary samples (e.g., all human and chimp livers), and performs the test
## - Read granges, M, Cov and smoothed values matrices
## - Create BSseq object from them

## smoothing parameters: loading only once
load("smooth_data/combined_samples/pData.RDa")
trans <- NULL
parameters <- NULL
## load one sample to get these parameters
load(paste0("smooth_data/combined_samples/separated_samples/", pData$Condition[1], "_Smoothed.RDa"))
trans <- getBSseq(data.fit, "trans") ## the function is hard coded so it will be the same or all samples 
parameters <- getBSseq(data.fit, "parameters")
rm(data.fit) ## clean up workspace

##############################
## all pairwise comparisons ##
##############################

## for each tissue
for (tissue in unique(pData$Tissue)){
  cat("Tissue: ", tissue, "\n")

  ## for every tissue we have to subset the samples, to be able to create the BSseq object, so for the next tissue we need to reload everything
  ## infos on samples
  load("smooth_data/combined_samples/pData.RDa")
  
  ## the coordinates of the CpG sites
  load("smooth_data/combined_samples/gr.RDa")
  
  ## coverage and methylation matrices
  load("smooth_data/combined_samples/Cov.RDa")
  load("smooth_data/combined_samples/M.RDa")
  
  ## matrix of smoothed data.
  load("smooth_data/combined_samples/coef.RDa")

  ## subset the data  
  keep <- pData$Tissue == tissue
  pData <- pData[keep,]
  M <- M[,keep]
  Cov <- Cov[,keep]
  coef <- coef[,keep]
  
  allData.fit <- BSseq(gr = gr[1:10,], M = Cov[1:10,], Cov = Cov[1:10,], coef = coef[1:10,], se.coef = NULL, pData = pData, trans = trans, parameters = parameters, rmZeroCov = FALSE)
  allData.fit@rowData <- gr
  allData.fit@assays$data$coef <- coef
  rm(coef)
  allData.fit@assays$data$M <- M
  rm(M)
  allData.fit@assays$data$Cov <- Cov
  ## TO DO: remove Cov (and use getCoverage below)?

  save(allData.fit, file = paste0("DMRs/species/pairwiseSpecies_", tissue, ".RDa"))

  ## for each pair of species
  for (pair in list(c("Human", "Chimp"), c("Human", "Rhesus"), c("Chimp", "Rhesus"))){
    cat("Testing species pair: ", pair, "\n")
    keep <- pData$Species == pair[1] | pData$Species == pair[2] 
   
    ## Only keep sites with at least 2 samples with coverage of 2 in each species
    keepLoci <- which(rowSums(Cov[, pData$Species == pair[1]] >= 2) >= 2 & rowSums(Cov[, pData$Species == pair[2]] >= 2) >= 2)
    ## remove loci with NA smoothed values
    keepLoci <- keepLoci[!is.na(apply(getMeth(allData.fit[keepLoci,], type = "smooth"), 1, sum))]
    ## remove loci with mean coverage > 10X (likely to be mapped to repeated regions)
    keepLoci <- keepLoci[rowMeans(Cov[keepLoci,]) <= 10]
    ## TO DO: last 2 filterings should be made on kept columns only! This is not likely to change much but it is cleaner...
    
    cat("Tested sites: ", format(length(keepLoci), big.mark=","), "\n")
    allData.fit.subset <- allData.fit[keepLoci, keep]

    ## calculate coverage at common sites
    cat("Coverage at common sites:\n")
    print(round(colMeans(getCoverage(allData.fit.subset)), 2))

    ## compute t-statistics
    allData.tstat <- BSmooth.tstat(allData.fit.subset,
                                   group1 = rownames(pData[pData$Species == pair[1],]),
                                   group2 = rownames(pData[pData$Species == pair[2],]),
                                   estimate.var = "same",
                                   local.correct = TRUE, ## large-scale (low-frequency) mean correction. This is especially important when large-scale methylation differences between 2 conditions (e.g., cancer and normals).
                                   verbose = TRUE)
    save(allData.tstat, file = paste0("DMRs/species/", pair[1], pair[2],"_", tissue, "_tstat.RDa"))
    pdf(file = paste0("DMRs/species/", pair[1], pair[2],"_", tissue, "_histTstats.pdf"), width = 6, height = 6)
    plot(allData.tstat)
    dev.off()

    ## find DMRs
    ## with local correction: small-scale DMRs
    dmrs0 <- dmrFinder(allData.tstat, cutoff = c(-4.6, 4.6), maxGap = 300, stat = "tstat.corrected", verbose=T)
    ## without local correction: large-scale DMRs (such as in cancer).
    ## /!\ the statistic may not be symmetrical. 
    ## For cutoff, see Hansen et al. 2014
    dmrs0.uncorrected <- dmrFinder(allData.tstat, cutoff = c(-2, 2), maxGap = 300, stat = "tstat", verbose=T)

    ## TO DO: adjust cutoff? comparable between species?
    dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
    cat("DMRs found: ", nrow(dmrs), "\n")
    print(head(dmrs, n = 5))
    write.table(dmrs, file= paste0("DMRs/species/", pair[1], pair[2],"_", tissue, "_DMRs.txt"), sep="\t", quote=F, row.names=F, col.names=T)

    dmrs.uncorrected <- subset(dmrs0.uncorrected, width >= 5000)
    cat("DMRs (uncorrected t-stat) found: ", nrow(dmrs.uncorrected), "\n")
    print(head(dmrs.uncorrected, n = 5))
    write.table(dmrs.uncorrected, file= paste0("DMRs/species/", pair[1], pair[2],"_", tissue, "_DMRs.uncorrected.txt"), sep="\t", quote=F, row.names=F, col.names=T)

    ## Plotting: 2 species tested included here only
    pData$col <- rep(c(pal[1], pal[2], pal[3]), each = 4)
    pData(allData.fit) <- pData
    keepLoci <- which(rowSums(getCoverage(allData.fit)) >= 1)
    allData.fit.subset <- allData.fit[keepLoci, keep]

    n <- 100
    if (nrow(dmrs) < n){
      n <- nrow(dmrs)
    }
    pdf(file = paste0("DMRs/species/", pair[1], pair[2],"_", tissue, "_top", n, "DMRs.pdf"), width = 10, height = 5)
    plotManyRegions(allData.fit.subset, dmrs[1:n,], extend = 5000, addRegions = dmrs, addPoints=TRUE, pointsMinCov=4)
    dev.off()

    n <- 100
    if (nrow(dmrs.uncorrected) < n){
      n <- nrow(dmrs.uncorrected)
    }
    if (n > 0){
      pdf(file = paste0("DMRs/species/", pair[1], pair[2],"_", tissue, "_top", n, "DMRs.uncorrected.pdf"), width = 10, height = 5)
      plotManyRegions(allData.fit.subset, dmrs.uncorrected[1:n,], extend = 5000, addRegions = dmrs.uncorrected)
      dev.off()
    }

    ## clean up worksapce
    rm(allData.fit.subset)
    rm(allData.tstat)
  }
}
print(warnings())

## TO DO: overlap between corrected and uncorrected DMRs?
## TO DO: randomization of samples labels to estimate FDR: see Hansen 2014 Genome Res



