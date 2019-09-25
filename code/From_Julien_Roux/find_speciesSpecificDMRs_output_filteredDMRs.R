## Jul 22, 2014
## When we first run this script, we outputed the unfiltered DMRs
## This script just outputs the filtered DMRs in a text file

library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))

## source("http://www.bioconductor.org/biocLite.R")
## biocLite()
## biocLite("BiocUpgrade")
## biocLite("bsseq")
## biocLite("bsseq", lib="~/R/x86_64-redhat-linux-gnu-library/2.15/") ## on the cluster
library(bsseq) # version 0.10
sessionInfo()

## for each tissue
for (tissue in unique(pData$Tissue)){
  cat("Tissue: ", tissue, "\n")

  ## for every tissue we have to subset the samples, to be able to create the BSseq object, so for the next tissue we need to reload everything
  ## infos on samples
  load("smooth_data/combined_samples/pData.RDa")
  
  ## subset the data  
  keep <- pData$Tissue == tissue
  pData <- pData[keep,]

  ## for each species
  for (species in c("Human", "Chimp")){
  ## for (species in c("Chimp")){
    cat("Testing species: ", species, "\n")

##     ## compute t-statistics: species of interest vs. 2 others
##     load(paste0("DMRs/species/", species,"Specific_", tissue, "_tstat.RDa"))
##     ## compute another t-statistics on same sites between 2 others species
##     load(paste0("DMRs/species/", species,"Specific_", tissue, "_tstatOthers.RDa"))
##     ## in allData.tstat object, set tstat.corrected to NA if it is significant in allData.tstat.others
##     allData.tstat@stats[abs(allData.tstat.others@stats[,"tstat.corrected"]) >= 4.6 ,"tstat.corrected"] <- NA
##     ## find DMRs
##     dmrs0 <- dmrFinder(allData.tstat, cutoff = c(-4.6, 4.6), maxGap = 300, stat = "tstat.corrected", verbose=T)

    ## We can start by reading the output file: no need to recalculate all DMRs
    dmrs0 <- read.table(paste0("DMRs/species/", species,"Specific_", tissue, "_unfilteredDMRs.txt"), sep="\t", h=T)
    
    dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
    ## TO DO: adjust criteria? 
    cat("DMRs found: ", nrow(dmrs), "\n")
    print(head(dmrs, n = 5))
    write.table(dmrs, file= paste0("DMRs/species/", species,"Specific_", tissue, "_DMRs.txt"), sep="\t", quote=F, row.names=F, col.names=T)

  }
}
