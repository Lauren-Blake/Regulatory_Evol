## Feb 24, 2015
## Requires a lot of memory
## Look for tissue DMRs, but on CpGs that are shared between 3 species only

library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))

## source("http://www.bioconductor.org/biocLite.R")
## biocLite()
## biocLite("BiocUpgrade")
## biocLite("bsseq")
## biocLite("bsseq", lib="~/R/x86_64-redhat-linux-gnu-library/2.15/") ## on the cluster
library(bsseq) # version 0.10
sessionInfo()

## Load smoothed data (all sampels at ~5M shared sites)
load("smooth_data/combined_samples/combinedSmoothedCommonSites.RDa")
load("smooth_data/combined_samples/pData.RDa")
allData.fit.subset$pData <- pData

## smoothing parameters: loading only once
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

## for each species
for (species in c("Human", "Rhesus", "Chimp")){
## for (species in c("Chimp")){  ## relaunch for only 1 species if there was not enough memory
  cat("Species: ", species, "\n")
  pData <- allData.fit.subset$pData
  keep <- pData$Species == species
  pData <- pData[keep, ]

  ## need to go through building BSseq object to include all parameters etc
  allData.fit <- BSseq(gr = allData.fit.subset@rowData, M = allData.fit.subset@assays$data$M[,keep], Cov = allData.fit.subset@assays$data$Cov[,keep], coef = allData.fit.subset@assays$data$coef[,keep], se.coef = NULL, pData = pData, trans = trans, parameters = parameters, rmZeroCov = FALSE)
  
  ## for each pair of tissue
  for (pair in list(c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "liver"))){
    cat("Testing tissue pair: ", pair, "\n")
    keep <- pData$Tissue == pair[1] | pData$Tissue == pair[2] 
  
    ## Only keep sites with at least 2 samples with coverage of 2 in each tissue
    keepLoci <- which(rowSums(getCoverage(allData.fit[, pData$Tissue == pair[1]]) >= 2) >= 2 & rowSums(getCoverage(allData.fit[, pData$Tissue == pair[2]]) >= 2) >= 2)
    ## remove loci with NA smoothed values
    keepLoci <- keepLoci[!is.na(apply(getMeth(allData.fit[keepLoci,], type = "smooth"), 1, sum))]
    ## remove loci with mean coverage > 10X (likely to be mapped to repeated regions)
    keepLoci <- keepLoci[rowMeans(getCoverage(allData.fit[keepLoci,])) <= 10]
    cat("Tested sites: ", format(length(keepLoci), big.mark=","), "\n")
    ## TO DO: last 2 filterings should be made on kept columns only! This is not likely to change much but it is cleaner...
    
    ## calculate coverage at common sites
    cat("Coverage at common sites:\n")
    print(round(colMeans(getCoverage(allData.fit[keepLoci, keep])), 2))

    ## compute t-statistics
    allData.tstat <- BSmooth.tstat(allData.fit[keepLoci, keep],
                                   group1 = rownames(pData[pData$Tissue == pair[1],]),
                                   group2 = rownames(pData[pData$Tissue  == pair[2],]),
                                   estimate.var = "same",
                                   local.correct = TRUE, ## large-scale (low-frequency) mean correction. This is especially important when large-scale methylation differences between 2 conditions (e.g., cancer and normals).
                                   verbose = TRUE)

    ## find DMRs
    ## with local correction: small-scale DMRs
    dmrs0 <- dmrFinder(allData.tstat, cutoff = c(-4.6, 4.6), maxGap = 300, stat = "tstat.corrected", verbose=T)
    dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
    cat("DMRs found: ", nrow(dmrs), "\n")
    print(head(dmrs, n = 5))
    write.table(dmrs, file= paste0("DMRs/tissues/", species, "_", pair[1], "_", pair[2],"_DMRs.sharedPositions.txt"), sep="\t", quote=F, row.names=F, col.names=T)

    ## clean up workspace
    rm(allData.tstat)
    rm(dmrs0)
    rm(dmrs)
  }
}
print(warnings())
