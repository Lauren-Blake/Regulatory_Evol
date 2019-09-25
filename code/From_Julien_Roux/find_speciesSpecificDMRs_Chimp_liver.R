## Jul 17, 2014
## There is an error when calculating the Chimp-specific DMRs in liver.
## The error comes when calculating the Human-rhesus contrast (we don't want this contrast to be significant). This contrast was successfully computed with find_speciesDMRs.R (using more sites). We use this to fix the claclutaion of Chimp-specific DMRs in liver. 

library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))

## source("http://www.bioconductor.org/biocLite.R")
## biocLite()
## biocLite("BiocUpgrade")
## biocLite("bsseq")
## biocLite("bsseq", lib="~/R/x86_64-redhat-linux-gnu-library/2.15/") ## on the cluster
library(bsseq) # version 0.10
sessionInfo()

tissue <- "liver"
species <- "Chimp"
load("smooth_data/combined_samples/pData.RDa")
keep <- pData$Tissue == tissue
pData <- pData[keep,]

load(paste0("DMRs/species/speciesSpecific_", tissue, ".RDa"))
Cov = getCoverage(allData.fit)

## Only keep sites with at least 2 samples with coverage of 2 in each of 3 species
keepLoci <- which(rowSums(Cov[, pData$Species == unique(pData$Species)[1]] >= 2) >= 2 & rowSums(Cov[, pData$Species == unique(pData$Species)[2]] >= 2) >=2 & rowSums(Cov[, pData$Species == unique(pData$Species)[3]] >= 2) >= 2)
## remove loci with NA smoothed values
keepLoci <- keepLoci[!is.na(apply(getMeth(allData.fit[keepLoci,], type = "smooth"), 1, sum))]
## remove loci with mean coverage > 10X (likely to be mapped to repeated regions)
keepLoci <- keepLoci[rowMeans(Cov[keepLoci,]) <= 10]
rm(Cov)
  
cat("Tested sites: ", format(length(keepLoci), big.mark=","), "\n")
allData.fit.subset <- allData.fit[keepLoci,]


## compute another t-statistics on same sites between 2 others species: this is bugging!
## allData.tstat.others <- BSmooth.tstat(allData.fit.subset,
##                                       group1 = rownames(pData[pData$Species == unique(pData[pData$Species != species,]$Species)[1],]),
##                                       group2 = rownames(pData[pData$Species == unique(pData[pData$Species != species,]$Species)[2],]), 
##                                       estimate.var = "same",
##                                       local.correct = TRUE,
##                                       verbose = TRUE)

load(paste0("DMRs/species/HumanRhesus_", tissue, "_tstat.RDa"))
allData.tstat.others <- allData.tstat

## t-statistics: species of interest vs. 2 others
load(paste0("DMRs/species/", species,"Specific_", tissue, "_tstat.RDa"))

## subset allData.tstat.others to common sites
allData.tstat.others <- subsetByOverlaps(allData.tstat.others, granges(allData.tstat))
save(allData.tstat.others, file = paste0("DMRs/species/", species,"Specific_", tissue, "_tstatOthers.RDa"))

## in allData.tstat object, set tstat.corrected to NA if it is significant in allData.tstat.others
allData.tstat@stats[abs(allData.tstat.others@stats[,"tstat.corrected"]) >= 4.6 ,"tstat.corrected"] <- NA

## find DMRs
dmrs0 <- dmrFinder(allData.tstat, cutoff = c(-4.6, 4.6), maxGap = 300, stat = "tstat.corrected", verbose=T)
dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
## TO DO: adjust criteria? 
cat("DMRs found: ", nrow(dmrs), "\n")
print(head(dmrs, n = 5))
write.table(dmrs, file= paste0("DMRs/species/", species,"Specific_", tissue, "_DMRs.txt"), sep="\t", quote=F, row.names=F, col.names=T)

## Plotting: 3 species included here
## To make nicer plots, we use sites observed in at least 1 sample in this tissue. Beware, these are not only CpG sites in human!
keepLoci <- which(rowSums(getCoverage(allData.fit)) >= 1)
cat("Number of sites used for plotting: ", format(length(keepLoci), big.mark=","), "\n")

pData$col <- rep(c(pal[1], pal[2], pal[3]), each = 4)
pData(allData.fit) <- pData
pdf(file = paste0("DMRs/species/", species,"Specific_", tissue, "_top100DMRs.pdf"), width = 10, height = 5)
plotManyRegions(allData.fit[keepLoci,], dmrs[1:100,], extend = 5000, addRegions = dmrs)
dev.off()

print(warnings())

## Error chimp liver
## [BSmooth.tstat] computing stats across groups ... Error in preplot.locfit.raw(object, newdata, where, what, band) :
##  NA/NaN/Inf in foreign function call (arg 2)
##In addition: Warning messages:
##1: In lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth,  :
##  max_nr reduction problem
##2: In lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth,  :
##  max_nr reduction problem
