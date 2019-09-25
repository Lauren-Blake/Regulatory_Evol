## Aug 21, 2014
## To have an estimate of the rate of false positives in DMR finding, this script launches several tests with samples from different individuals in the 2 groups to contrast.

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
## for (species in c("Human", "Chimp")){
for (species in c("Rhesus")){ ## Relaunched because there was a swapping of samples in Rhesus
  cat("Species: ", species, "\n")

  ## object generated in script find_individualDMRs.R
  load(paste0("DMRs/individuals/", species, "_pairwiseIndividuals.RDa"))
  pData <- pData(allData.fit)

  ## for each pair of individuals
  ## for (pair in list(c("1", "2"), c("1", "3"), c("1", "4"), c("2", "3"), c("2", "4"), c("3", "4"))){
  for (pair in list(c("2", "3"), c("2", "4"), c("3", "4"))){
    cat("Testing individual pair: ", pair, "\n")
    keep <- substr(pData$IndividualID, 2, 2) == pair[1] | substr(pData$IndividualID, 2, 2) == pair[2] 

    ## Only keep sites with at least 2 samples with coverage of 2 in each individual
    keepLoci <- which(rowSums(getCoverage(allData.fit[, substr(pData$IndividualID, 2, 2) == pair[1]]) >= 2) >= 2 & rowSums(getCoverage(allData.fit[, substr(pData$IndividualID, 2, 2) == pair[2]]) >= 2) >= 2)
    ## remove loci with NA smoothed values
    keepLoci <- keepLoci[!is.na(apply(getMeth(allData.fit[keepLoci,], type = "smooth"), 1, sum))]
    ## remove loci with mean coverage > 10X (likely to be mapped to repeated regions)
    keepLoci <- keepLoci[rowMeans(getCoverage(allData.fit[keepLoci,])) <= 10]
    cat("Tested sites: ", format(length(keepLoci), big.mark=","), "\n")

    ## labels permutations: see Hansen et al. 2014
    ## 4 samples vs. 4 samples: use all combinations of 2 samples of each individuals, which makes 6 combinations, but the 3 last are the reverse of teh 3 first
    ## Note: here (and in tissues DMRs permutations), there are less permutations than in speciesDMRs script since we use samples from the same species, we cannot use individuals randomly

    comb_1 <- combn(rownames(pData[substr(pData$IndividualID, 2, 2) == pair[1],]), 2)
    comb_2 <- combn(rownames(pData[substr(pData$IndividualID, 2, 2) == pair[2],]), 2)
    ## reverse the comb_2 matrix, so that in each group we have all the 4 organs represented
    comb_2 <- t(apply(t(comb_2), 2, rev))

    count <- 1
    for (i in 1:ncol(comb_1)) {
      cat("Permutation", count, ": (", paste0(c(comb_1[,i], comb_2[,i])), ") vs. (", c(rownames(pData[substr(pData$IndividualID, 2, 2) == pair[1],])[!rownames(pData[substr(pData$IndividualID, 2, 2) == pair[1],]) %in% comb_1[,i]], rownames(pData[substr(pData$IndividualID, 2, 2) == pair[2],])[!rownames(pData[substr(pData$IndividualID, 2, 2) == pair[2],]) %in% comb_2[,i]]), ")\n")
      
      ## compute t-statistics
      allData.tstat <- BSmooth.tstat(allData.fit[keepLoci, keep],
                                     group1 = c(comb_1[,i], comb_2[,i]),
                                     group2 = c(rownames(pData[substr(pData$IndividualID, 2, 2) == pair[1],])[!rownames(pData[substr(pData$IndividualID, 2, 2) == pair[1],]) %in% comb_1[,i]], rownames(pData[substr(pData$IndividualID, 2, 2) == pair[2],])[!rownames(pData[substr(pData$IndividualID, 2, 2) == pair[2],]) %in% comb_2[,i]]),
                                     estimate.var = "same",
                                     local.correct = TRUE,
                                     verbose = TRUE)
      ## save(allData.tstat, file = paste0("DMRs/individuals/", species, "_", pair[1], "_", pair[2],"_tstat.RDa"))
      
      ## find DMRs
      dmrs0 <- dmrFinder(allData.tstat, cutoff = c(-4.6, 4.6), maxGap = 300, stat = "tstat.corrected", verbose=T)
      dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
      cat("DMRs found: ", nrow(dmrs), "\n")
      write.table(dmrs, file= paste0("DMRs/individuals/permutations/", species, "_", pair[1], "_", pair[2],"_permut", count, "_DMRs.txt"), sep="\t", quote=F, row.names=F, col.names=T)

      ## increment permutation number
      count <- count + 1
    }
  }
}
print(warnings())

## TO DO: we could mix more than 2 individuals in each group? Beware of tissue (4 tissues have to be represented). E.g., H1H+H2K+H3Li+H4Lu vs. H2H+H3K+H4Li+H1Lu
