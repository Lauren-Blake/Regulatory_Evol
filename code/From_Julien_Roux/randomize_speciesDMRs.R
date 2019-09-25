## Aug 10, 2014
## Requires a lot of memory (>100g)
## To have an estimate of the rate of false positives in DMR finding, this script launches several tests with samples from both species in the 2 groups to contrast.

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

## for each tissue
for (tissue in c("heart", "kidney", "liver", "lung")){
  cat("Tissue:", tissue, "\n")

  ## object generated in script find_speciesDMRs.R
  load(paste0("DMRs/species/pairwiseSpecies_", tissue, ".RDa"))
  pData <- pData(allData.fit)
  
  ## for each pair of species: 
  for (pair in list(c("Human", "Chimp"), c("Human", "Rhesus"), c("Chimp", "Rhesus"))){
    cat("Testing species pair: ", pair, "\n")
    keep <- pData$Species == pair[1] | pData$Species == pair[2] 
   
    ## Only keep sites with at least 2 samples with coverage of 2 in each species
    keepLoci <- which(rowSums(getCoverage(allData.fit[, pData$Species == pair[1]]) >= 2) >= 2 & rowSums(getCoverage(allData.fit[, pData$Species == pair[2]]) >= 2) >= 2)
    ## remove loci with NA smoothed values
    keepLoci <- keepLoci[!is.na(apply(getMeth(allData.fit[keepLoci,], type = "smooth"), 1, sum))]
    ## remove loci with mean coverage > 10X (likely to be mapped to repeated regions)
    keepLoci <- keepLoci[rowMeans(getCoverage(allData.fit[keepLoci,])) <= 10]
    allData.fit.subset <- allData.fit[keepLoci, keep]

    ## labels permutations: see Hansen et al. 2014
    ## 4 samples vs. 4 samples: use all combinations of 2 samples of each species, which makes 6 * 6 combinations, but the 18 last are the reverse of the 18 first: 18 unique combinations
    comb_1 <-  combn(rownames(pData[pData$Species == pair[1],]), 2)
    comb_2 <-  combn(rownames(pData[pData$Species == pair[2],]), 2)

    count <- 1
    for (i in 1:ncol(comb_1)) {
      for (j in 1:ncol(comb_2)) {
        cat("Permutation", count, ": (", paste0(c(comb_1[,i], comb_2[,j])), ") vs. (", c(rownames(pData[pData$Species == pair[1],])[!rownames(pData[pData$Species == pair[1],]) %in% comb_1[,i]], rownames(pData[pData$Species == pair[2],])[!rownames(pData[pData$Species == pair[2],]) %in% comb_2[,j]]), ")\n")
       
        ## compute t-statistics
        allData.tstat <- BSmooth.tstat(allData.fit.subset,
                                       group1 = c(comb_1[,i], comb_2[,j]),
                                       group2 = c(rownames(pData[pData$Species == pair[1],])[!rownames(pData[pData$Species == pair[1],]) %in% comb_1[,i]], rownames(pData[pData$Species == pair[2],])[!rownames(pData[pData$Species == pair[2],]) %in% comb_2[,j]]),
                                       estimate.var = "same",
                                       local.correct = TRUE,
                                       verbose = TRUE)
        ## save(allData.tstat, file = paste0("DMRs/species/permutations/", pair[1], pair[2],"_", tissue, "_tstat.RDa"))

        ## find DMRs
        dmrs0 <- dmrFinder(allData.tstat, cutoff = c(-4.6, 4.6), maxGap = 300, stat = "tstat.corrected", verbose=T)
        dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
        cat("    DMRs found: ", nrow(dmrs), "\n")
        write.table(dmrs, file= paste0("DMRs/species/permutations/", pair[1], pair[2],"_", tissue, "_permut", count, "_DMRs.txt"), sep="\t", quote=F, row.names=F, col.names=T)

        ## increment permutation number
        count <- count + 1
      }
    }
  }
}
print(warnings())



