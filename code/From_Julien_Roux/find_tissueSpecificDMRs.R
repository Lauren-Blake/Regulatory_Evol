## Jul 10, 2014
## Requires a lot of memory (>100g)
## Find tissue-specific DMRs: compare 1 tissue vs. all others
## Add constraint for no difference between the other tissues

library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))

## source("http://www.bioconductor.org/biocLite.R")
## biocLite()
## biocLite("BiocUpgrade")
## biocLite("bsseq")
## biocLite("bsseq", lib="~/R/x86_64-redhat-linux-gnu-library/2.15/") ## on the cluster
library(bsseq) # version 0.10
sessionInfo()

## Load all necessary samples (e.g., all human livers and hearts), and performs the test
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

## for each species
## for (species in c("Human", "Rhesus", "Chimp")){
for (species in c("Rhesus")){ ## relaunch for only 1 species if there was not enough memory
  cat("Species: ", species, "\n")

  ## for every species we have to subset the samples, to be able to create the BSseq object, so for the next species we need to reload everything
  ## infos on samples
  load("smooth_data/combined_samples/pData.RDa")
  
  ## the coordinates of the CpG sites
  load("smooth_data/combined_samples/gr.RDa")
  
  ## coverage and methylation matrices
  load("smooth_data/combined_samples/Cov.RDa")
  load("smooth_data/combined_samples/M.RDa")
  
  ## matrix of smoothed data.
  load("smooth_data/combined_samples/coef.RDa")

  ## subset the data: keep only one species
  keep <- pData$Species == species
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
  save(allData.fit, file = paste0("DMRs/tissues/", species, "_tissueSpecific.RDa"))

  ## Only keep sites with at least 2 samples with coverage of 2 in each tissue
  keepLoci <- which(rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[1]] >= 2) >= 2 & rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[2]] >= 2) >= 2 & rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[3]] >= 2) >= 2 & rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[4]] >= 2) >= 2)
  ## remove loci with NA smoothed values
  keepLoci <- keepLoci[!is.na(apply(getMeth(allData.fit[keepLoci,], type = "smooth"), 1, sum))]
  ## remove loci with mean coverage > 10X (likely to be mapped to repeated regions)
  keepLoci <- keepLoci[rowMeans(Cov[keepLoci,]) <= 10]
  rm(Cov)
  
  cat("Tested sites: ", format(length(keepLoci), big.mark=","), "\n")
  ## allData.fit.subset <- allData.fit[keepLoci,] ## we don't create this object to save space

  ## calculate coverage at common sites
  cat("Coverage at common sites:\n")
  print(round(colMeans(getCoverage(allData.fit[keepLoci, ])), 2))
  
  ## for each tissue
  ## for (tissue in c("heart", "lung", "kidney", "liver")){
  for (tissue in c("kidney", "liver", "lung")){
    cat("Testing tissue: ", tissue, "\n")

    ## compute t-statistics between tissue of interest and other tissues
    allData.tstat <- BSmooth.tstat(allData.fit[keepLoci, ],
                                   group1 = rownames(pData[pData$Tissue == tissue,]), ## tissue of interest
                                   group2 = rownames(pData[pData$Tissue != tissue,]), ## 3 other tissues pooled
                                   estimate.var = "same",
                                   local.correct = TRUE, ## large-scale (low-frequency) mean correction. This is especially important when large-scale methylation differences between 2 conditions (e.g., cancer and normals).
                                   verbose = TRUE)
    save(allData.tstat, file = paste0("DMRs/tissues/", species, "_", tissue,"Specific_tstat.RDa"))
    pdf(file = paste0("DMRs/tissues/", species, "_", tissue ,"Specific_histTstats.pdf"), width = 6, height = 6)
    plot(allData.tstat)
    dev.off()

    ## compute other t-statistics on same sites between each pair of other tissues
    allData.tstat.others1 <- BSmooth.tstat(allData.fit[keepLoci, ],
                                           group1 = rownames(pData[pData$Tissue == unique(pData[pData$Tissue != tissue,]$Tissue)[1],]),
                                           group2 = rownames(pData[pData$Tissue == unique(pData[pData$Tissue != tissue,]$Tissue)[2],]),
                                           estimate.var = "same",
                                           local.correct = TRUE,
                                           verbose = TRUE)
    save(allData.tstat.others1, file = paste0("DMRs/tissues/", species, "_", tissue,"Specific_tstatOthers1.RDa"))

    allData.tstat.others2 <- BSmooth.tstat(allData.fit[keepLoci, ],
                                           group1 = rownames(pData[pData$Tissue == unique(pData[pData$Tissue != tissue,]$Tissue)[2],]),
                                           group2 = rownames(pData[pData$Tissue == unique(pData[pData$Tissue != tissue,]$Tissue)[3],]),
                                           estimate.var = "same",
                                           local.correct = TRUE,
                                           verbose = TRUE)
    save(allData.tstat.others2, file = paste0("DMRs/tissues/", species, "_", tissue,"Specific_tstatOthers2.RDa"))

    allData.tstat.others3 <- BSmooth.tstat(allData.fit[keepLoci, ],
                                           group1 = rownames(pData[pData$Tissue == unique(pData[pData$Tissue != tissue,]$Tissue)[1],]),
                                           group2 = rownames(pData[pData$Tissue == unique(pData[pData$Tissue != tissue,]$Tissue)[3],]),
                                           estimate.var = "same",
                                           local.correct = TRUE,
                                           verbose = TRUE)
    save(allData.tstat.others3, file = paste0("DMRs/tissues/", species, "_", tissue,"Specific_tstatOthers3.RDa"))

    ## in allData.tstat object, set tstat.corrected to NA if it is significant in any allData.tstat.others objects
    allData.tstat@stats[abs(allData.tstat.others1@stats[,"tstat.corrected"]) >= 4.6 ,"tstat.corrected"] <- NA
    allData.tstat@stats[abs(allData.tstat.others2@stats[,"tstat.corrected"]) >= 4.6 ,"tstat.corrected"] <- NA
    allData.tstat@stats[abs(allData.tstat.others3@stats[,"tstat.corrected"]) >= 4.6 ,"tstat.corrected"] <- NA

    ## TO DO? repeat without local correction? This makes more sense for pairwise comparisons

    ## find DMRs
    dmrs0 <- dmrFinder(allData.tstat, cutoff = c(-4.6, 4.6), maxGap = 300, stat = "tstat.corrected", verbose=T)

    ## TO DO: adjust cutoff? comparable between tissues?
    dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
    ## TO DO: adjust criteria? 
    cat("DMRs found: ", nrow(dmrs), "\n")
    print(head(dmrs, n = 5))
    write.table(dmrs, file= paste0("DMRs/tissues/", species, "_", tissue,"Specific_DMRs.txt"), sep="\t", quote=F, row.names=F, col.names=T)

    ## Plotting: all tissues included here
    n <- 100
    if (nrow(dmrs) < n){
      n <- nrow(dmrs)
    }
    ## To make nicer plots, we use sites observed in at least 1 sample in this species. Beware, when the species in Chimp or Rhesus, these are not only CpG sites in the human reference coordinates!
    keepLoci <- which(rowSums(getCoverage(allData.fit)) >= 1)
    cat("Number of sites used for plotting: ", format(length(keepLoci), big.mark=","), "\n")
    
    pData$col <- rep(c(pal[1], pal[2], pal[3], pal[4]), times = 4)
    pData(allData.fit) <- pData
    pdf(file = paste0("DMRs/tissues/", species, "_", tissue, "Specific_top", n, "DMRs.pdf"), width = 10, height = 5)
    plotManyRegions(allData.fit[keepLoci,], dmrs[1:n,], extend = 5000, addRegions = dmrs)
    dev.off()
    
    ## clean up worksapce
    rm(allData.tstat)
    rm(allData.tstat.others1)
    rm(allData.tstat.others2)
    rm(allData.tstat.others3)
    rm(dmrs0)
    rm(dmrs)
  }
}
print(warnings())
