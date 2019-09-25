## Aug 19, 2014
## In the initial script, there is a problem with Chimp liver, kidney and lung, and Rhesus liver and lung. Try to fix this manually in this script

library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))

## source("http://www.bioconductor.org/biocLite.R")
## biocLite()
## biocLite("BiocUpgrade")
## biocLite("bsseq")
## biocLite("bsseq", lib="~/R/x86_64-redhat-linux-gnu-library/2.15/") ## on the cluster
library(bsseq) # version 0.10
sessionInfo()

##species <- "Rhesus"
species <- "Chimp"

load("smooth_data/combined_samples/pData.RDa")
keep <- pData$Species == species
pData <- pData[keep,]

load(paste0("DMRs/tissues/", species, "_tissueSpecific.RDa"))
Cov <- getCoverage(allData.fit)

## Only keep sites with at least 2 samples with coverage of 2 in each tissue
keepLoci <- which(rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[1]] >= 2) >= 2 & rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[2]] >= 2) >= 2 & rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[3]] >= 2) >= 2 & rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[4]] >= 2) >= 2)

## remove loci with NA smoothed values
keepLoci <- keepLoci[!is.na(apply(getMeth(allData.fit[keepLoci,], type = "smooth"), 1, sum))]
## remove loci with mean coverage > 10X (likely to be mapped to repeated regions)
keepLoci <- keepLoci[rowMeans(Cov[keepLoci,]) <= 10]
rm(Cov)
cat("Tested sites: ", format(length(keepLoci), big.mark=","), "\n")

## Rhesus: 18,189,778 sites 
## Chimp: 19,006,612 sites

tissue <- "lung"
tissue <- "liver"
tissue <- "kidney"

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
## Error message for Chimp kidney: see below

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
## Error message for Chimp lung: but this does not seem to prevent the t-stat to be computed at all sites!
## Maybe this is enough to stop the script, but when launched manually it is not a problem?


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


## ## Plotting: not made on manual tests, because it can always be done later
## ## All tissues are included 
## n <- 100
## if (nrow(dmrs) < n){
##   n <- nrow(dmrs)
## }
## ## To make nicer plots, we use sites observed in at least 1 sample in this species. Beware, when the species in Chimp or Rhesus, these are not only CpG sites in the human reference coordinates!
## keepLoci <- which(rowSums(getCoverage(allData.fit)) >= 1)
## cat("Number of sites used for plotting: ", format(length(keepLoci), big.mark=","), "\n")

## pData$col <- rep(c(pal[1], pal[2], pal[3], pal[4]), times = 4)
## pData(allData.fit) <- pData
## pdf(file = paste0("DMRs/tissues/", species, "_", tissue, "Specific_top", n, "DMRs.pdf"), width = 10, height = 5)
## plotManyRegions(allData.fit[keepLoci,], dmrs[1:n,], extend = 5000, addRegions = dmrs)
## dev.off()

