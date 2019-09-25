## Aug 21, 2014
## To have an estimate of the rate of false positives in DMR finding, this script launches several tests with samples from different tissues in the 2 groups to contrast.

library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))

## source("http://www.bioconductor.org/biocLite.R")
## biocLite()
## biocLite("BiocUpgrade")
## biocLite("bsseq")
## biocLite("bsseq", lib="~/R/x86_64-redhat-linux-gnu-library/2.15/") ## on the cluster
library(bsseq) # version 0.10
sessionInfo()

##############################
## all pairwise comparisons ##
##############################

## for each species
for (species in c("Human", "Rhesus", "Chimp")){
  cat("Species: ", species, "\n")

  ## object generated in script find_individualDMRs.R
  load(paste0("DMRs/tissues/", species, "_pairwiseTissues.RDa"))
  pData <- pData(allData.fit)
 
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
    
    ## labels permutations: see Hansen et al. 2014
    ## 4 samples vs. 4 samples: use all combinations of 2 samples of each tissue, which makes 6 combinations, but the 3 last are the reverse of the 3 first.
    ## Note: here (and in tissues DMRs eprmutations), there are less permutations than in speciesDMRs script since we use samples from the same species, we cannot use individuals randomly

    comb_1 <- combn(rownames(pData[pData$Tissue == pair[1],]), 2)
    comb_2 <- combn(rownames(pData[pData$Tissue == pair[2],]), 2)
    ## reverse the comb_2 matrix, so that in each group we have all the 4 organs represented
    comb_2 <- t(apply(t(comb_2), 2, rev))

    count <- 1
    for (i in 1:ncol(comb_1)) {
      cat("Permutation", count, ": (", paste0(c(comb_1[,i], comb_2[,i])), ") vs. (", c(rownames(pData[pData$Tissue == pair[1],])[!rownames(pData[pData$Tissue == pair[1],]) %in% comb_1[,i]], rownames(pData[pData$Tissue == pair[2],])[!rownames(pData[pData$Tissue == pair[2],]) %in% comb_2[,i]]), ")\n")
      
      ## compute t-statistics
      allData.tstat <- BSmooth.tstat(allData.fit[keepLoci, keep],
                                     group1 = c(comb_1[,i], comb_2[,i]),
                                     group2 = c(rownames(pData[pData$Tissue == pair[1],])[!rownames(pData[pData$Tissue == pair[1],]) %in% comb_1[,i]], rownames(pData[pData$Tissue == pair[2],])[!rownames(pData[pData$Tissue == pair[2],]) %in% comb_2[,i]]),
                                     estimate.var = "same",
                                     local.correct = TRUE,
                                     verbose = TRUE)
      
      ## find DMRs
      dmrs0 <- dmrFinder(allData.tstat, cutoff = c(-4.6, 4.6), maxGap = 300, stat = "tstat.corrected", verbose=T)
      dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
      cat("DMRs found: ", nrow(dmrs), "\n")
      write.table(dmrs, file= paste0("DMRs/tissues/permutations/", species, "_", pair[1], "_", pair[2],"_permut", count, "_DMRs.txt"), sep="\t", quote=F, row.names=F, col.names=T)

      ## increment permutation number
      count <- count + 1
    }
  }
}
print(warnings())
