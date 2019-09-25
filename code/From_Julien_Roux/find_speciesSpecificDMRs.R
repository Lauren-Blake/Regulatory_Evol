## Jul 8, 2014
## Requires a lot of memory
## Find species DMRs under directional selection: compare human vs chimp+rhesus and chimp vs human+rhesus
## Add constraint for no difference between the 2 other species

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
trans <- getBSseq(data.fit, "trans") ## the function is hard coded so it will be the same for all samples 
parameters <- getBSseq(data.fit, "parameters")
rm(data.fit) ## clean up workspace

## for each tissue
for (tissue in unique(pData$Tissue)){
## for (tissue in c("liver")){
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
  save(allData.fit, file = paste0("DMRs/species/speciesSpecific_", tissue, ".RDa"))

  ## Only keep sites with at least 2 samples with coverage of 2 in each of 3 species
  keepLoci <- which(rowSums(Cov[, pData$Species == unique(pData$Species)[1]] >= 2) >= 2 & rowSums(Cov[, pData$Species == unique(pData$Species)[2]] >= 2) >=2 & rowSums(Cov[, pData$Species == unique(pData$Species)[3]] >= 2) >= 2)
  ## remove loci with NA smoothed values
  keepLoci <- keepLoci[!is.na(apply(getMeth(allData.fit[keepLoci,], type = "smooth"), 1, sum))]
  ## remove loci with mean coverage > 10X (likely to be mapped to repeated regions)
  keepLoci <- keepLoci[rowMeans(Cov[keepLoci,]) <= 10]
  rm(Cov)
  
  cat("Tested sites: ", format(length(keepLoci), big.mark=","), "\n")
  allData.fit.subset <- allData.fit[keepLoci,]

  ## calculate coverage at common sites
  cat("Coverage at common sites:\n")
  print(round(colMeans(getCoverage(allData.fit.subset)), 2))

  ## for each species
  for (species in c("Human", "Chimp")){
  ## for (species in c("Chimp")){
    cat("Testing species: ", species, "\n")

    ## compute t-statistics: species of interest vs. 2 others
    allData.tstat <- BSmooth.tstat(allData.fit.subset,
                                   group1 = rownames(pData[pData$Species == species,]),
                                   group2 = rownames(pData[pData$Species != species,]), ## 2 species mixed here
                                   estimate.var = "same",
                                   local.correct = TRUE, ## large-scale (low-frequency) mean correction. This is especially important when large-scale methylation differences between 2 conditions (e.g., cancer and normals).
                                   verbose = TRUE)
    save(allData.tstat, file = paste0("DMRs/species/", species,"Specific_", tissue, "_tstat.RDa"))
    pdf(file = paste0("DMRs/species/", species,"Specific_", tissue, "_histTstats.pdf"), width = 6, height = 6)
    plot(allData.tstat)
    dev.off()

    ## compute another t-statistics on same sites between 2 others species
    allData.tstat.others <- BSmooth.tstat(allData.fit.subset,
                                   group1 = rownames(pData[pData$Species == unique(pData[pData$Species != species,]$Species)[1],]),
                                   group2 = rownames(pData[pData$Species == unique(pData[pData$Species != species,]$Species)[2],]), 
                                   estimate.var = "same",
                                   local.correct = TRUE,
                                   verbose = TRUE)
    save(allData.tstat.others, file = paste0("DMRs/species/", species,"Specific_", tissue, "_tstatOthers.RDa"))

    ## in allData.tstat object, set tstat.corrected to NA if it is significant in allData.tstat.others
    allData.tstat@stats[abs(allData.tstat.others@stats[,"tstat.corrected"]) >= 4.6 ,"tstat.corrected"] <- NA

    ## TO DO? repeat without local correction? This makes more sense for pairwise comparisons ()!!!

    ## find DMRs
    dmrs0 <- dmrFinder(allData.tstat, cutoff = c(-4.6, 4.6), maxGap = 300, stat = "tstat.corrected", verbose=T)

    ## TO DO: adjust cutoff? comparable between species?
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
    
    ## clean up workspace
    rm(allData.tstat)
    rm(allData.tstat.others)
  }
}
print(warnings())

## Error chimp liver
## [BSmooth.tstat] computing stats across groups ... Error in preplot.locfit.raw(object, newdata, where, what, band) :
##  NA/NaN/Inf in foreign function call (arg 2)
##In addition: Warning messages:
##1: In lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth,  :
##  max_nr reduction problem
##2: In lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth,  :
##  max_nr reduction problem
