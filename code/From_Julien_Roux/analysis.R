## Mar 10, 2014
## Diverse analyses on BS-seq data

library(RColorBrewer)
## display.brewer.all()
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))
## pal <- brewer.pal(9, "Set1")
##pal <- colorRampPalette(pal)(23)
## pie(rep(1,23), col=pal)

## source("http://www.bioconductor.org/biocLite.R")
## biocLite()
## biocLite("BiocUpgrade")
## biocLite("bsseq")
## biocLite("bsseq", lib="~/R/x86_64-redhat-linux-gnu-library/2.15/") ## on the cluster
library(bsseq) # version 0.10
sessionInfo()
setwd("~/Methylation/bsseq")

#############################################
## create matrix of samples and conditions ##
#############################################

samples <- read.table("../raw_data/SamplesDirectories.txt", sep="\t", h=T)
samples <- unique(samples[samples$Type == "BS-seq", ])
samples[,12] <- paste0("../bismark/", paste(samples$Flow.cell, samples$SampleID, sep="_"), "/")
names(samples)[12] <- "Path"
row.names(samples) <- paste(samples$Flow.cell, samples$SampleID, sep="_")
samples <- samples[,c(4,5,6,7,12)]
samples <- unique(samples)
samples[,6] <- as.factor(unlist(lapply(strsplit(as.character(samples$Condition), split=""), function(x) { return(paste0(x[1], x[2])) })))
names(samples)[6] <- "IndividualID"
samples <- samples[,c(1,2,3,4,6,5)]
samples <- samples[order(samples$SampleID),]

###########################
## Coverage: perl script ## 
###########################
coverage <- read.table(gzfile("./combine_data_perl/coverage.gz", "r"), h=T, sep="\t")
cov <- apply(coverage[,-c(1,2)], 2, mean)

##      C1H      C1K     C1Li     C1Lu      C2H      C2K     C2Li     C2Lu 
## 3.763463 3.886618 5.695652 4.539875 3.686341 2.902934 4.331456 3.430911 
##      C3H      C3K     C3Li     C3Lu      C4H      C4K     C4Li     C4Lu 
## 5.451209 4.886743 4.596521 3.070719 3.555441 3.199126 4.886564 2.603082 
##      H1H      H1K     H1Li     H1Lu      H2H      H2K     H2Li     H2Lu 
## 5.522072 1.733158 4.984468 3.808516 3.712729 3.753488 3.771592 3.600061 
##      H3H      H3K     H3Li     H3Lu      H4H      H4K     H4Li     H4Lu 
## 4.816675 5.683911 3.743829 3.350242 4.021280 3.975080 4.850055 3.696047 
##      R1H      R1K     R1Li     R1Lu      R2H      R2K     R2Li     R2Lu 
## 4.274348 3.058604 4.727085 4.052753 5.156459 4.466520 4.947335 4.509918 
##      R3H      R3K     R3Li     R3Lu      R4H      R4K     R4Li     R4Lu 
## 3.202226 4.584275 4.131124 3.385906 4.222444 3.294217 4.949179 2.811424 
## This info was added to Sup table

summary(all_coverage)
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  1.733   3.524   3.998   4.068   4.749   5.696 

pdf(file = "coverage_orthologous_sites.pdf", width = 6, height = 5)
boxplot(cov ~ c(rep("Chimp", times=16), rep("Human", times=16), rep("Rhesus", times=16)), col=pal[1:3], notch=F, pch=16, boxwex=0.7)
dev.off()

## can be taken from BSseq object too
load("smooth_data/combined_samples/combinedSmoothedCommonSites.RDa")
round(colMeans(getCoverage(allData.fit.subset)), 2)
summary(colMeans(getCoverage(allData.fit.subset)))
## same result as above
summary(getCoverage(allData.fit.subset))
summary(getCoverage(allData.fit.subset[,"H2K"]))
##  Min.   : 0.000  
##  1st Qu.: 2.000  
##  Median : 3.000  
##  Mean   : 3.753  
##  3rd Qu.: 5.000  
##  Max.   :77.000 
table(getCoverage(allData.fit.subset[,"H2K"]))
##     0      1      2      3      4      5      6      7      8      9     10  ...
## 359673 616676 800670 859777 789887 639144 463985 305425 185561 103647  54674 ...


## Histogram of coverage at all sites
pdf(file = "distributionCoverage_H2K.pdf", height = 6, width=6)
hist(getCoverage(allData.fit.subset[,"H2K"]), breaks=100, main="H2K", xlab="Coverage", xlim=c(0,80))
hist(log10(getCoverage(allData.fit.subset[,"H2K"])), breaks=20, main="H2K", xlab="log10(Coverage)", xlim=c(0,log10(80)))
dev.off()
## coverage is not continuous (1,2,3,...) so hard to plot in a density plot
## TO DO? Use all sites: here we selected sites with cov >= 1 in at least 1 sample in each species 
## TO DO? reproduce for all samples? It's likely going to be similar

## Mean coverage across all samples
pdf(file = "distributionMeanCoverage.pdf", height = 6, width=6)
plot(density(rowMeans(getCoverage(allData.fit.subset))), col=pal[1], main="", xlab="Mean coverage", lwd=1, lty=1)
abline(v=10, lty=2)
dev.off()
## log scale
pdf(file = "distributionMeanCoverageLog.pdf", height = 6, width=6)
plot(density(log10(rowMeans(getCoverage(allData.fit.subset)))), col=pal[1], main="", xlab="log10(mean coverage)", lwd=1, lty=1)
abline(v=log10(10), lty=2)
dev.off()

#########################
## Mappability of CpGs ##
#########################
load("../bsseq/smooth_data/combined_samples/Cov.RDa")
Cov = Cov[, grepl("^R", colnames(Cov), perl=T)]
summary(rowSums(Cov >= 1) == 0)
#    Mode    FALSE     TRUE     NA's 
# logical 21433025 25161372        0 
# Not possible to get some coverage on more than 21.4M in macaque

load("../bsseq/smooth_data/combined_samples/Cov.RDa")
Cov = Cov[, grepl("^H", colnames(Cov), perl=T)]
summary(rowSums(Cov >= 1) == 0)
##    Mode    FALSE     TRUE     NA's
##  logical 26186385 20408012        0

#############################################################################
## Distribution of coverage on DMR sites vs. non DMR sites: coverage bias? ##
#############################################################################

## Do it for species DMRs since it uses the 5M shared CpG
allDMRs <- vector()
for (pair in list(c("Human", "Chimp"), c("Human", "Rhesus"), c("Chimp", "Rhesus"))){
  for (tissue in c("heart", "kidney", "liver", "lung")){
    cat("Species:", pair, " / Tissue:", tissue, "\n")
    load(paste0("DMRs/species/", pair[1], pair[2], "_", tissue, "_tstat.RDa"))
    ## find DMRs
    dmrs0 <- dmrFinder(allData.tstat, cutoff = c(-4.6, 4.6), maxGap = 300, stat = "tstat.corrected", verbose=T)
    allDMRs <- sort(unique(c(allDMRs, paste(dmrs0$chr, dmrs0$start, sep=":"))))
  }
}
names <- paste(seqnames(granges(allData.fit.subset)), start(granges(allData.fit.subset)), sep=":")
summary(names %in% allDMRs)
##   Mode   FALSE    TRUE    NA's 
##logical 5045498  182971       0 

pdf(file = "distributionCoverageSpeciesDMRs.pdf", height = 6, width=6)
plot(density(log10(rowMeans(getCoverage(allData.fit.subset[names %in% allDMRs ,])))), col=pal[1], main="", xlab="log10(mean coverage)", lwd=1, lty=1)
lines(density(log10(rowMeans(getCoverage(allData.fit.subset[!names %in% allDMRs ,])))), col=pal[2], lwd=1, lty=1)
abline(v=log10(10), lty=2)
legend("topleft", legend=c("Diff. methylated", "Not diff. methylated"), lty=c(1,1), col=pal[1:2]) 
dev.off()
## There seems to be no difference: of course lower bounded because we filter a coverage >= 2X


## Do it in 1 particular contrast (Human vs. chimp in heart), and use coverage data of these particular samples
load(paste0("DMRs/species/HumanChimp_heart_tstat.RDa"))
dmrs0 <- dmrFinder(allData.tstat, cutoff = c(-4.6, 4.6), maxGap = 300, stat = "tstat.corrected", verbose=T)
DMRs <- sort(unique(paste(dmrs0$chr, dmrs0$start, sep=":")))
names <- paste(seqnames(granges(allData.fit.subset)), start(granges(allData.fit.subset)), sep=":")
length(DMRs)
info <- pData(allData.fit.subset)

## pdf(file = "distributionCoverageHumanChimpHeartDMRs.pdf", height = 6, width=6)
plot(density(log10(rowMeans(getCoverage(allData.fit.subset[names %in% DMRs, (info$Species == "Human"|info$Species == "Chimp") & info$Tissue == "heart"])))), col=pal[1], main="Human vs. Chimp in heart", xlab="log10(mean coverage)", lwd=1, lty=1)
lines(density(log10(rowMeans(getCoverage(allData.fit.subset[!(names %in% DMRs), (info$Species == "Human"|info$Species == "Chimp") & info$Tissue == "heart"])))), col=pal[2], lwd=1, lty=1)
abline(v=log10(10), lty=2)
legend("topleft", legend=c("Diff. methylated", "Not diff. methylated"), lty=c(1,1), col=pal[1:2]) 
## dev.off()

## Actually we should only plot the filtered CpG sites: they are in the Tstat object
namesFiltered <- paste(seqnames(granges(allData.tstat)), start(granges(allData.tstat)), sep=":")
length(namesFiltered) ## 15349967
pdf(file = "distributionCoverageHumanChimpHeartDMRs.pdf", height = 6, width=6)
## 7880 sites:
plot(density(log10(rowMeans(getCoverage(allData.fit.subset[names %in% DMRs, (info$Species == "Human"|info$Species == "Chimp") & info$Tissue == "heart"])))), col=pal[1], main="Human vs. Chimp in heart", xlab="log10(mean coverage)", lwd=1, lty=1)
## 4695842 sites:
lines(density(log10(rowMeans(getCoverage(allData.fit.subset[(names %in% namesFiltered & !(names %in% DMRs)), (info$Species == "Human"|info$Species == "Chimp") & info$Tissue == "heart"])))), col=pal[2], lwd=1, lty=1)
abline(v=log10(10), lty=2)
legend("topleft", legend=c("Diff. methylated", "Not diff. methylated"), lty=c(1,1), col=pal[1:2]) 

## same thing with no log scale
plot(density(rowMeans(getCoverage(allData.fit.subset[names %in% DMRs, (info$Species == "Human"|info$Species == "Chimp") & info$Tissue == "heart"]))), col=pal[1], main="Human vs. Chimp in heart", xlab="mean coverage", lwd=1, lty=1, xlim=c(0,12))
## 4695842 sites:
lines(density(rowMeans(getCoverage(allData.fit.subset[(names %in% namesFiltered & !(names %in% DMRs)), (info$Species == "Human"|info$Species == "Chimp") & info$Tissue == "heart"]))), col=pal[2], lwd=1, lty=1)
legend("topright", legend=c("Diff. methylated", "Not diff. methylated"), lty=c(1,1), col=pal[1:2]) 

## make smoother plots
plot(density(rowMeans(getCoverage(allData.fit.subset[names %in% DMRs, (info$Species == "Human"|info$Species == "Chimp") & info$Tissue == "heart"])), adjust=2), col=pal[1], main="Human vs. Chimp in heart", xlab="mean coverage", lwd=1, lty=1, xlim=c(0,12))
## 4695842 sites:
lines(density(rowMeans(getCoverage(allData.fit.subset[(names %in% namesFiltered & !(names %in% DMRs)), (info$Species == "Human"|info$Species == "Chimp") & info$Tissue == "heart"])), adjust=2), col=pal[2], lwd=1, lty=1)
legend("topright", legend=c("Diff. methylated", "Not diff. methylated"), lty=c(1,1), col=pal[1:2]) 
dev.off()

## stat test to check significance of difference:
ks.test(rowMeans(getCoverage(allData.fit.subset[names %in% DMRs, (info$Species == "Human"|info$Species == "Chimp") & info$Tissue == "heart"])), rowMeans(getCoverage(allData.fit.subset[(names %in% namesFiltered & !(names %in% DMRs)), (info$Species == "Human"|info$Species == "Chimp") & info$Tissue == "heart"])))
## D = 0.051, p-value < 2.2e-16
## alternative hypothesis: two-sided


## This is taking into account all detected CpGs. Maybe the difference attenuates when we filter the DMRs?
dmrs1 <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1) ## filtered
DMRs.filtered <- sort(unique(paste(dmrs1$chr, dmrs1$start, sep=":")))
pdf(file = "distributionCoverageHumanChimpHeartFilteredDMRs.pdf", height = 6, width=6)
## 5320 sites
plot(density(rowMeans(getCoverage(allData.fit.subset[names %in% DMRs.filtered, (info$Species == "Human"|info$Species == "Chimp") & info$Tissue == "heart"])), adjust=2), col=pal[1], main="Human vs. Chimp in heart", xlab="mean coverage", lwd=1, lty=1, xlim=c(0,12))
## 4695842 sites:
lines(density(rowMeans(getCoverage(allData.fit.subset[(names %in% namesFiltered & !(names %in% DMRs)), (info$Species == "Human"|info$Species == "Chimp") & info$Tissue == "heart"])), adjust=2), col=pal[2], lwd=1, lty=1)
legend("topright", legend=c("Within DMRs", "Not diff. methylated"), lty=c(1,1), col=pal[1:2]) 
dev.off()


## Let's try a tissue contrast:
## Here we only use the shared positions between species, but we don't have to.
## We could use the object already saved for Human, but it quite big and we need a session with a lot of memory
## load("DMRs/tissues/Human_pairwiseTissues.RDa")
## Anyway the results are going to be similar
load(paste0("DMRs/tissues/Human_liver_lung_tstat.RDa"))
dmrs0 <- dmrFinder(allData.tstat, cutoff = c(-4.6, 4.6), maxGap = 300, stat = "tstat.corrected", verbose=T)
DMRs <- sort(unique(paste(dmrs0$chr, dmrs0$start, sep=":")))
length(DMRs) ##16753
dmrs1 <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1) ## filtered
DMRs.filtered <- sort(unique(paste(dmrs1$chr, dmrs1$start, sep=":")))
length(DMRs.filtered) ##13010
info <- pData(allData.fit.subset)
## plot the filtered (based on coverage mainly) CpG sites: they are in the Tstat object
namesFiltered <- paste(seqnames(granges(allData.tstat)), start(granges(allData.tstat)), sep=":")
length(namesFiltered) ## 21406393
pdf(file = "distributionCoverageHumanLiverLungFilteredDMRs.pdf", height = 6, width=6)
## 3691 sites
plot(density(rowMeans(getCoverage(allData.fit.subset[names %in% DMRs.filtered, info$Species == "Human" & (info$Tissue == "heart" | info$Tissue == "lung")])), adjust=2), col=pal[1], main="Human liver vs. lung", xlab="mean coverage", lwd=1, lty=1, xlim=c(0,12))
## 4591686 sites:
lines(density(rowMeans(getCoverage(allData.fit.subset[(names %in% namesFiltered & !(names %in% DMRs)), info$Species == "Human" & (info$Tissue == "heart" | info$Tissue == "lung")])), adjust=2), col=pal[2], lwd=1, lty=1)
legend("topright", legend=c("Within DMRs", "Not diff. methylated"), lty=c(1,1), col=pal[1:2]) 
dev.off()

#########################################################################################
## Export the coordinates of the CpG sites used for DMR finding (to calculate density) ##
#########################################################################################
## Tissue DMRs
for (species in c("Human", "Rhesus", "Chimp")){
  tissuePairs <- list(c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "liver"), c("heart", "Specific"), c("liver", "Specific"), c("kidney", "Specific"), c("lung", "Specific"))
  for (pair in tissuePairs){
    cat("Testing tissue pair: ", pair, "in", species, "\n")
    ## Reading DMR file
    if (pair[2] == "Specific"){
      load(paste0("DMRs/tissues/", species, "_", pair[1], pair[2],"_tstat.RDa"))
    }
    else {
      load(paste0("DMRs/tissues/", species, "_", pair[1], "_", pair[2],"_tstat.RDa"))
    }
    gr <- granges(allData.tstat)
    ## Granges object is 0-based (start and end position): convert to Bed file (start is 0-based, end is 1-based)    
    df <- data.frame(seqnames=seqnames(gr), starts=start(gr), ends=end(gr)+1)
    ## disable scientific notations
    options("scipen"=100)
    
    if (pair[2] == "Specific"){
      write.table(df, file=paste0("DMRs/tissues/", species, "_", pair[1], pair[2],"_tstat_allSites.bed"), quote=F, sep="\t", row.names=F, col.names=F)

    } else {
      write.table(df, file=paste0("DMRs/tissues/", species, "_", pair[1], "_", pair[2],"_tstat_allSites.bed"), quote=F, sep="\t", row.names=F, col.names=F)
    }
  }
}

## Command line to get number of CpG sites per DMR: 
## coverageBed -a DMRs/tissues/Human_heart_lung_tstat_allSites.bed -b ../annotation/DMRs/tissues/Human_heart_lung_DMRs.bed -counts | less
## coverageBed -a DMRs/tissues/Human_heart_lung_tstat_allSites.bed -b ../annotation/DMRs/tissues/randomized/Human_heart_lung_randomDMRs_73.bed -counts > DMRs/tissues/number_CpG_sites/Human_heart_lung_randomDMRs_73_number_CpG_sites.bed
## Launched with launch_all_calculate_DMRs_CpG_density.pl

## Reading the results and plotting them:
library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))
for (species in c("Human", "Chimp", "Rhesus")){
  for (pair in list(c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "liver"), c("heart", "Specific"), c("liver", "Specific"), c("kidney", "Specific"), c("lung", "Specific"))){
    cat("Testing tissue pair: ", pair, "in", species, "\n")
    if (pair[2] != "Specific"){
      baseName <- paste0(species, "_", pair[1], "_", pair[2])
    } else {
      baseName <- paste0(species, "_", pair[1], "Specific")
    }
    fileName <- paste0("./DMRs/tissues/number_CpG_sites/", baseName, "_DMRs_number_CpG_sites.bed")
    realFeatures <- read.table(fileName, h=F, sep="\t", check.names = F)
    realFeatures[,5] <-  ( realFeatures[,3]- realFeatures[,2] ) / realFeatures[,4]
    nDMRs <- length(realFeatures[,1])
    ymax <- max(density(realFeatures[,5])$y)
    xdens <- density(realFeatures[,5])$x
    xmax <- max(xdens[is.finite(xdens)])
      
    ## control for length only
    randomFeatures <- list()
    for (i in 1:100){
      print(i)
      fileName <- paste0("./DMRs/tissues/number_CpG_sites/", baseName, "_randomDMRs_", i, "_number_CpG_sites.bed")
      randomFeatures[[i]] <- read.table(fileName, h=F, sep="\t", check.names = F)
      randomFeatures[[i]][,5] <-  ( randomFeatures[[i]][,3]- randomFeatures[[i]][,2] ) / randomFeatures[[i]][,4]
      if (max(density(randomFeatures[[i]][,5])$y) > ymax){
        ymax <- max(density(randomFeatures[[i]][,5])$y)
      }
      xdens <- density(randomFeatures[[i]][,5])$x
      if (max(xdens[is.finite(xdens)]) > xmax){
        xmax <- max(xdens[is.finite(xdens)])
      }
    }

    plotTitle <- vector()
    if (pair[2] == "Specific"){
      plotTitle <- paste0(pair[1], "-specific in ", species, ": ", nDMRs, " DMRs")
    } else {
      plotTitle <- paste0(pair[1], " vs. ", pair[2], " in ", species, ": ", nDMRs, " DMRs")
    }

    pdf(file = paste0("./DMRs/tissues/number_CpG_sites/", baseName, "_DMRs_number_CpG_sites.pdf"), width = 6, height = 5)
    ## xlimit <- c(0, xmax) ## sometimes resulting scale is hard to read
    xlimit <- c(0, 1000)
    ylimit <- c(0, ymax)
    ## first plot the random DMRs density plots
    plot(density(randomFeatures[[1]][,5]), xlim=xlimit, ylim=ylimit, main=plotTitle, xlab=expression(DMR ~ CpG ~ density^{-1}), col="gray", lwd=0.5)
    for (i in 2:100){
      lines(density(randomFeatures[[i]][,5]), col="gray", lwd=0.5)
    }
    ## draw real DMRs density
    lines(density(realFeatures[,5]), col=pal[1], lwd=2)
    dev.off()
  }
}

## species DMRs
for (pair in list(c("Human", "Chimp"), c("Human", "Rhesus"), c("Chimp", "Rhesus"), c("Human", "Specific"), c("Chimp", "Specific"))){
  for (tissue in c("heart", "kidney", "liver", "lung")){
    cat("Species:", pair, " / Tissue:", tissue, "\n")

    ## Reading DMR file
    load(paste0("DMRs/species/", pair[1], pair[2], "_", tissue, "_tstat.RDa"))
    gr <- granges(allData.tstat)
    ## Granges object is 0-based (start and end position): convert to Bed file (start is 0-based, end is 1-based)    
    df <- data.frame(seqnames=seqnames(gr), starts=start(gr), ends=end(gr)+1)
    ## disable scientific notations
    options("scipen"=100)
    
    write.table(df, file=paste0("DMRs/species/", pair[1], pair[2], "_", tissue, "_tstat_allSites.bed"), quote=F, sep="\t", row.names=F, col.names=F)
  }
}

## Command line to get number of CpG sites per DMR: 
# coverageBed -a DMRs/species/HumanChimp_lung_tstat_allSites.bed -b ../annotation/DMRs/species/HumanChimp_lung_DMRs.bed -counts | less
## Launched with launch_all_calculate_DMRs_CpG_density.pl

## Reading the results and plotting them:
library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))

for (pair in list(c("Human", "Chimp"), c("Human", "Rhesus"), c("Chimp", "Rhesus"), c("Human", "Specific"), c("Chimp", "Specific"))){
  for (tissue in c("heart", "kidney", "liver", "lung")){
    cat("Species:", pair, " / Tissue:", tissue, "\n")

    baseName <- paste0(pair[1], pair[2], "_", tissue)

    fileName <- paste0("./DMRs/species/number_CpG_sites/", baseName, "_DMRs_number_CpG_sites.bed")
    realFeatures <- read.table(fileName, h=F, sep="\t", check.names = F)
    realFeatures[,5] <-  ( realFeatures[,3]- realFeatures[,2] ) / realFeatures[,4]
    nDMRs <- length(realFeatures[,1])
    ymax <- max(density(realFeatures[,5])$y)
    xdens <- density(realFeatures[,5])$x
    xmax <- max(xdens[is.finite(xdens)])
      
    ## control for length only
    randomFeatures <- list()
    for (i in 1:100){
      print(i)
      fileName <- paste0("./DMRs/species/number_CpG_sites/", baseName, "_randomDMRs_", i, "_number_CpG_sites.bed")
      randomFeatures[[i]] <- read.table(fileName, h=F, sep="\t", check.names = F)
      randomFeatures[[i]][,5] <-  ( randomFeatures[[i]][,3]- randomFeatures[[i]][,2] ) / randomFeatures[[i]][,4]
      if (max(density(randomFeatures[[i]][,5])$y) > ymax){
        ymax <- max(density(randomFeatures[[i]][,5])$y)
      }
      xdens <- density(randomFeatures[[i]][,5])$x
      if (max(xdens[is.finite(xdens)]) > xmax){
        xmax <- max(xdens[is.finite(xdens)])
      }
    }
    plotTitle <- vector()
    if (pair[2] == "Specific"){
      plotTitle <- paste0(pair[1], "-specific in ", tissue, ": ", nDMRs, " DMRs")
    } else {
      plotTitle <- paste0(pair[1], " vs. ", pair[2], " in ", tissue, ": ", nDMRs, " DMRs")
    }

    pdf(file = paste0("./DMRs/species/number_CpG_sites/", baseName, "_DMRs_number_CpG_sites.pdf"), width = 6, height = 5)
    ## xlimit <- c(0, xmax) ## sometimes resulting scale is hard to read
    xlimit <- c(0, 1000)
    ylimit <- c(0, ymax)
    ## first plot the random DMRs density plots
    plot(density(randomFeatures[[1]][,5]), xlim=xlimit, ylim=ylimit, main=plotTitle, xlab=expression(DMR ~ CpG ~ density^{-1}), col="gray", lwd=0.5)
    for (i in 2:100){
      lines(density(randomFeatures[[i]][,5]), col="gray", lwd=0.5)
    }
    ## draw real DMRs density
    lines(density(realFeatures[,5]), col=pal[1], lwd=2)
    dev.off()
  }
}

##################################################################################
## Plotting the density of sequence conservation scores for DMRs and randomDMRs ##
##################################################################################

## Reading the results and plotting them:
library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))
for (species in c("Human", "Chimp", "Rhesus")){
  for (pair in list(c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "liver"), c("heart", "Specific"), c("liver", "Specific"), c("kidney", "Specific"), c("lung", "Specific"))){
    cat("Testing tissue pair: ", pair, "in", species, "\n")
    if (pair[2] != "Specific"){
      baseName <- paste0(species, "_", pair[1], "_", pair[2])
    } else {
      baseName <- paste0(species, "_", pair[1], "Specific")
    }
    fileName <- paste0("../annotation/DMRs/tissues/", baseName, "_DMRs_meanSeqConservation.bed")
    realFeatures <- read.table(fileName, h=F, sep="\t", check.names = F)
    nDMRs <- length(realFeatures[,1])
    ymax <- max(density(realFeatures[,5])$y)
    xdens <- density(realFeatures[,5])$x
    xmin <- min(xdens[is.finite(xdens)])
    xmax <- max(xdens[is.finite(xdens)])
      
    ## control for length only
    randomFeatures <- list()
    for (i in 1:100){
      print(i)
      fileName <- paste0("../annotation/DMRs/tissues/randomized/", baseName, "_randomDMRs_", i, "_meanSeqConservation.bed")
      randomFeatures[[i]] <- read.table(fileName, h=F, sep="\t", check.names = F)
      if (max(density(randomFeatures[[i]][,5])$y) > ymax){
        ymax <- max(density(randomFeatures[[i]][,5])$y)
      }
      xdens <- density(randomFeatures[[i]][,5])$x
      if (min(xdens[is.finite(xdens)]) < xmin){
        xmin <- min(xdens[is.finite(xdens)])
      }
      if (max(xdens[is.finite(xdens)]) > xmax){
        xmax <- max(xdens[is.finite(xdens)])
      }
    }

    ## control for length and CpG density
    randomFeaturesControl <- list()
    for (i in 1:100){
      print(i)
      fileName <- paste0("../annotation/DMRs/tissues/randomized_control_CpG_density/", baseName, "_randomDMRs_", i, "_meanSeqConservation.bed")
      randomFeaturesControl[[i]] <- read.table(fileName, h=F, sep="\t", check.names = F)
      if (max(density(randomFeaturesControl[[i]][,5])$y) > ymax){
        ymax <- max(density(randomFeaturesControl[[i]][,5])$y)
      }
      xdens <- density(randomFeaturesControl[[i]][,5])$x
      if (min(xdens[is.finite(xdens)]) < xmin){
        xmin <- min(xdens[is.finite(xdens)])
      }
      if (max(xdens[is.finite(xdens)]) > xmax){
        xmax <- max(xdens[is.finite(xdens)])
      }
    }

    plotTitle <- vector()
    if (pair[2] == "Specific"){
      plotTitle <- paste0(pair[1], "-specific in ", species, ": ", nDMRs, " DMRs")
    } else {
      plotTitle <- paste0(pair[1], " vs. ", pair[2], " in ", species, ": ", nDMRs, " DMRs")
    }

    pdf(file = paste0("../annotation/DMRs/tissues/", baseName, "_DMRs_meanSeqConservation.pdf"), width = 6, height = 5)
    ## xlimit <- c(xmin, xmax)
    xlimit <- c(-5, 5) ## not necessary to render the full range
    ylimit <- c(0, ymax)
    ## first plot the random DMRs density plots
    plot(density(randomFeatures[[1]][,5]), xlim=xlimit, ylim=ylimit, main=plotTitle, xlab="mean RS score", col=gray(0.8), lwd=0.5)
    for (i in 2:100){
      lines(density(randomFeatures[[i]][,5]), col=gray(0.8), lwd=0.5)
    }
    for (i in 2:100){
      lines(density(randomFeaturesControl[[i]][,5]), col=gray(0.4), lwd=0.5)
    }
    ## draw real DMRs density
    lines(density(realFeatures[,5]), col=pal[1], lwd=2)
    dev.off()
  }
}

## species DMRs
for (pair in list(c("Human", "Chimp"), c("Human", "Rhesus"), c("Chimp", "Rhesus"), c("Human", "Specific"), c("Chimp", "Specific"))){
  for (tissue in c("heart", "kidney", "liver", "lung")){
    cat("Species:", pair, " / Tissue:", tissue, "\n")

    baseName <- paste0(pair[1], pair[2], "_", tissue)
    fileName <- paste0("../annotation/DMRs/species/", baseName, "_DMRs_meanSeqConservation.bed")

    realFeatures <- read.table(fileName, h=F, sep="\t", check.names = F)
    nDMRs <- length(realFeatures[,1])
    ymax <- max(density(realFeatures[,5])$y)
    xdens <- density(realFeatures[,5])$x
    xmin <- min(xdens[is.finite(xdens)])
    xmax <- max(xdens[is.finite(xdens)])
      
    ## control for length only
    randomFeatures <- list()
    for (i in 1:100){
      print(i)
      fileName <- paste0("../annotation/DMRs/species/randomized/", baseName, "_randomDMRs_", i, "_meanSeqConservation.bed")
      randomFeatures[[i]] <- read.table(fileName, h=F, sep="\t", check.names = F)
      if (max(density(randomFeatures[[i]][,5])$y) > ymax){
        ymax <- max(density(randomFeatures[[i]][,5])$y)
      }
      xdens <- density(randomFeatures[[i]][,5])$x
      if (min(xdens[is.finite(xdens)]) < xmin){
        xmin <- min(xdens[is.finite(xdens)])
      }
      if (max(xdens[is.finite(xdens)]) > xmax){
        xmax <- max(xdens[is.finite(xdens)])
      }
    }

    ## control for length and CpG density
    randomFeaturesControl <- list()
    for (i in 1:100){
      print(i)
      fileName <- paste0("../annotation/DMRs/species/randomized_control_CpG_density/", baseName, "_randomDMRs_", i, "_meanSeqConservation.bed")
      randomFeaturesControl[[i]] <- read.table(fileName, h=F, sep="\t", check.names = F)
      if (max(density(randomFeaturesControl[[i]][,5])$y) > ymax){
        ymax <- max(density(randomFeaturesControl[[i]][,5])$y)
      }
      xdens <- density(randomFeaturesControl[[i]][,5])$x
      if (min(xdens[is.finite(xdens)]) < xmin){
        xmin <- min(xdens[is.finite(xdens)])
      }
      if (max(xdens[is.finite(xdens)]) > xmax){
        xmax <- max(xdens[is.finite(xdens)])
      }
    }

    plotTitle <- vector()
    if (pair[2] == "Specific"){
      plotTitle <- paste0(pair[1], "-specific in ", tissue, ": ", nDMRs, " DMRs")
    } else {
      plotTitle <- paste0(pair[1], " vs. ", pair[2], " in ", tissue, ": ", nDMRs, " DMRs")
    }

    pdf(file = paste0("../annotation/DMRs/species/", baseName, "_DMRs_meanSeqConservation.pdf"), width = 6, height = 5)
    ## xlimit <- c(xmin, xmax)
    xlimit <- c(-5, 5) ## not necessary to render the full range
    ylimit <- c(0, ymax)
    ## first plot the random DMRs density plots
    plot(density(randomFeatures[[1]][,5]), xlim=xlimit, ylim=ylimit, main=plotTitle, xlab="mean RS score", col=gray(0.8), lwd=0.5)
    for (i in 2:100){
      lines(density(randomFeatures[[i]][,5]), col=gray(0.8), lwd=0.5)
    }
    for (i in 2:100){
      lines(density(randomFeaturesControl[[i]][,5]), col=gray(0.4), lwd=0.5)
    }
    ## draw real DMRs density
    lines(density(realFeatures[,5]), col=pal[1], lwd=2)
    dev.off()
  }
}

##########################################################
## Plotting the distance to TSS for DMRs and randomDMRs ##
##########################################################
## Only plotting for human DMRs here (since it is human annotation)

## Reading the results and plotting them:
library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))
for (species in c("Human")){
  for (pair in list(c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "liver"), c("heart", "Specific"), c("liver", "Specific"), c("kidney", "Specific"), c("lung", "Specific"))){
    cat("Testing tissue pair: ", pair, "in", species, "\n")
    if (pair[2] != "Specific"){
      baseName <- paste0(species, "_", pair[1], "_", pair[2])
    } else {
      baseName <- paste0(species, "_", pair[1], "Specific")
    }
    fileName <- paste0("../annotation/DMRs/tissues/", baseName, "_DMRs_distanceToTSS.bed")
    realFeatures <- read.table(fileName, h=F, sep="\t", check.names = F)
    nDMRs <- length(realFeatures[,1])
      
    ## control for length only
    randomFeatures <- list()
    for (i in 1:100){
      print(i)
      fileName <- paste0("../annotation/DMRs/tissues/randomized/", baseName, "_randomDMRs_", i, "_distanceToTSS.bed")
      randomFeatures[[i]] <- read.table(fileName, h=F, sep="\t", check.names = F)
    }

    ## control for length and CpG density
    randomFeaturesControl <- list()
    for (i in 1:100){
      print(i)
      fileName <- paste0("../annotation/DMRs/tissues/randomized_control_CpG_density/", baseName, "_randomDMRs_", i, "_distanceToTSS.bed")
      randomFeaturesControl[[i]] <- read.table(fileName, h=F, sep="\t", check.names = F)
    }

    plotTitle <- vector()
    if (pair[2] == "Specific"){
      plotTitle <- paste0(pair[1], "-specific in ", species, ": ", nDMRs, " DMRs")
    } else {
      plotTitle <- paste0(pair[1], " vs. ", pair[2], " in ", species, ": ", nDMRs, " DMRs")
    }

    pdf(file = paste0("../annotation/DMRs/tissues/", baseName, "_DMRs_distanceToTSS.pdf"), width = 6, height = 5)
    xlimit <- c(0, 7) ## better rendering in log scale
    ylimit <- c(0, 0.8) ## better rendering in log scale
    ## first plot the random DMRs density plots
    plot(density(log10(randomFeatures[[1]][,10]+1)), xlim=xlimit, ylim=ylimit, main=plotTitle, xlab="log10(distance to TSS + 1)", col=gray(0.8), lwd=0.5)
    for (i in 2:100){
      lines(density(log10(randomFeatures[[i]][,10]+1)), col=gray(0.8), lwd=0.5)
    }
    for (i in 2:100){
      lines(density(log10(randomFeaturesControl[[i]][,10]+1)), col=gray(0.4), lwd=0.5)
    }
    ## draw real DMRs density
    lines(density(log10(realFeatures[,10]+1)), col=pal[1], lwd=2)
    dev.off()
  }
}
## Deconvolute 3 distributions (gaussians) with mclust
species <- "Human"
pair <- c("heart", "liver")
baseName <- paste0(species, "_", pair[1], "_", pair[2])
fileName <- paste0("../annotation/DMRs/tissues/", baseName, "_DMRs_distanceToTSS.bed")
realFeatures <- read.table(fileName, h=F, sep="\t", check.names = F)
nDMRs <- length(realFeatures[,1])

## install.packages("mclust") ## doesn't work
## In ~/R/x86_64-redhat-linux-gnu-library/2.15/:
## wget http://cran.r-project.org/src/contrib/Archive/mclust/mclust_4.2.tar.gz
## R CMD INSTALL mclust_4.2.tar.gz 
library(mclust)
mod = densityMclust(log10(realFeatures[,10] + 1), G=3)
summary(mod, parameters = TRUE)
pdf(file = paste0("../annotation/DMRs/tissues/", baseName, "_DMRs_distanceToTSS_deconvolution.pdf"), width = 6, height = 5)
plot(mod, what = "density", data = log10(realFeatures[,10] + 1), breaks = 100)
## Remove 0s:
summary(realFeatures[,10]==0) ## 0.093
mod = densityMclust(log10(realFeatures[realFeatures[,10]!=0, 10]), G=2)
summary(mod, parameters = TRUE)
plot(mod, what = "density", data = log10(realFeatures[realFeatures[,10]!=0, 10] + 1), breaks = 100)
plot(density(log10(realFeatures[realFeatures[,10]!=0, 10] + 1), adjust=1), ylim=c(0, 0.6), lwd=2, xlab="log10(distance to TSS + 1)", main="", xlim=c(0, 7))
lines(density(c(rep(-2, times=sum(mod$classification == 2)), rnorm(sum(mod$classification == 1), mean = mod$parameters$mean[1], sd = sqrt(mod$parameters$variance$sigmasq)[1])), adjust=1), col=pal[4], lwd=2)
lines(density(c(rep(-2, times=sum(mod$classification == 1)), rnorm(sum(mod$classification == 2), mean = mod$parameters$mean[2], sd = sqrt(mod$parameters$variance$sigmasq)[2])), adjust=1), col=pal[8], lwd=2)
eval.points <- seq(from = mod$range[1], to = mod$range[2], length = 1000)
d <- predict.densityMclust(mod, eval.points)
lines(eval.points, d, col=pal[2], lwd=2, lty=2)
## just identify the 2 modes
plot(density(log10(realFeatures[realFeatures[,10]!=0, 10] + 1), adjust=1), ylim=c(0, 0.6), lwd=2, xlab="log10(distance to TSS + 1)", main="", xlim=c(0, 7))
abline(v=3)
abline(v=4.1)
dev.off()
## Not very convincing: fit is poor. The distribution is probably not a sum of 2 simple gaussians


## species DMRs
## Only plotting for DMRs involving human here (since it is human annotation)
for (pair in list(c("Human", "Chimp"), c("Human", "Rhesus"), c("Human", "Specific"))){
  for (tissue in c("heart", "kidney", "liver", "lung")){
    cat("Species:", pair, " / Tissue:", tissue, "\n")

    baseName <- paste0(pair[1], pair[2], "_", tissue)
    fileName <- paste0("../annotation/DMRs/species/", baseName, "_DMRs_distanceToTSS.bed")
    realFeatures <- read.table(fileName, h=F, sep="\t", check.names = F)
    nDMRs <- length(realFeatures[,1])
      
    ## control for length only
    randomFeatures <- list()
    for (i in 1:100){
      print(i)
      fileName <- paste0("../annotation/DMRs/species/randomized/", baseName, "_randomDMRs_", i, "_distanceToTSS.bed")
      randomFeatures[[i]] <- read.table(fileName, h=F, sep="\t", check.names = F)
    }

    ## control for length and CpG density
    randomFeaturesControl <- list()
    for (i in 1:100){
      print(i)
      fileName <- paste0("../annotation/DMRs/species/randomized_control_CpG_density/", baseName, "_randomDMRs_", i, "_distanceToTSS.bed")
      randomFeaturesControl[[i]] <- read.table(fileName, h=F, sep="\t", check.names = F)
    }

    plotTitle <- vector()
    if (pair[2] == "Specific"){
      plotTitle <- paste0(pair[1], "-specific in ", tissue, ": ", nDMRs, " DMRs")
    } else {
      plotTitle <- paste0(pair[1], " vs. ", pair[2], " in ", tissue, ": ", nDMRs, " DMRs")
    }

    pdf(file = paste0("../annotation/DMRs/species/", baseName, "_DMRs_distanceToTSS.pdf"), width = 6, height = 5)
    xlimit <- c(0, 7) ## better rendering in log scale
    ylimit <- c(0, 0.8) ## better rendering in log scale
    ## first plot the random DMRs density plots
    plot(density(log10(randomFeatures[[1]][,10]+1)), xlim=xlimit, ylim=ylimit, main=plotTitle, xlab="log10(distance to TSS + 1)", col=gray(0.8), lwd=0.5)
    for (i in 2:100){
      lines(density(log10(randomFeatures[[i]][,10]+1)), col=gray(0.8), lwd=0.5)
    }
    for (i in 2:100){
      lines(density(log10(randomFeaturesControl[[i]][,10]+1)), col=gray(0.4), lwd=0.5)
    }
    ## draw real DMRs density
    lines(density(log10(realFeatures[,10]+1)), col=pal[1], lwd=2)
    dev.off()
  }
}

#################################################################
## Sup table: number of sites with >2X coverage in each sample ##
#################################################################
## For Chimp and Rhesus, use the number of sites maped to human genome (could be non-CpGs actually)
## We read the matrix of coverage combined for all samples:
load("smooth_data/combined_samples/Cov.RDa")
dim(Cov)
## 46594397       48
## number of sites:
numSites <- colSums(Cov >= 1)
numSites <- cbind(numSites, colSums(Cov >= 2))
numSites <- cbind(numSites, colSums(Cov >= 4))
colnames(numSites) <- c("1X", "2X", "4X")
write.table(numSites, file="number_sites_used.txt", sep="\t", quote=F, row.names=T, col.names=F)
## Added this to SupTables.xls (2nd tab: BS-seq - Conditions stats)

########################################
## histogram of methylation estimates ##
########################################
## Use data at ~5M common sites
load("./smooth_data/combined_samples/combinedSmoothedCommonSites.RDa")
info <- pData(allData.fit.subset)

pdf(file = "distributionMethylationSmoothed.pdf", height = 6, width=6)
for (species in unique(info$Species)){
  i <- 0
  for (tissue in unique(info$Tissue)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      if (i == 0){
        plot(density(getMeth(allData.fit.subset[,sample])), col=pal[as.integer(as.factor(info[sample,]$Tissue))], main=species, xlab="Methylation percentage", lwd=1, lty=j)
        i <- 1
      }
      else {
        lines(density(getMeth(allData.fit.subset[,sample])), col=pal[as.integer(as.factor(info[sample,]$Tissue))], lwd=1, lty=j)
      }
    }
  }
  legend("top", legend=c(as.character(unique(info$Tissue)), as.character(unique(info[info$Species == species,]$IndividualID))), lty=c(1,1,1,1,1,2,3,4), col=c(pal[1:4], rep("black",4))) 
}
dev.off()
## Same plot, but gathering same tissue in different species on different plots
pdf(file = "distributionMethylationSmoothed_by_tissues.pdf", height = 6, width=24)
par(mfrow=c(1,4), ps = 14, cex = 1, cex.main = 1) ## Used to keep labels at Font 14pt
for (tissue in unique(info$Tissue)){
  i <- 0
  for (species in unique(info$Species)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      if (i == 0){
        plot(density(getMeth(allData.fit.subset[,sample])), col=pal[as.integer(as.factor(info[sample,]$Species))+9], main=tissue, xlab="Methylation percentage", lwd=2, lty=j, ylim=c(0,4.5))
        i <- 1
      }
      else {
        lines(density(getMeth(allData.fit.subset[,sample])), col=pal[as.integer(as.factor(info[sample,]$Species))+9], lwd=2, lty=j)
      }
    }
  }
  if (tissue==unique(info$Tissue)[1]){
    legend("top", legend=c(as.character(unique(info$Species)), as.character(unique(info[info$Tissue == tissue,]$IndividualID))), lty=c(1,1,1,rep(c(1,2,3,4), times=3)), lwd=2, col=c(pal[10:12], rep("black",12)), bty = "n") 
  }
}
dev.off()

## density on raw estimates: depends a lot on coverage
pdf(file = "distributionMethylationRaw.pdf", height = 6, width=6)
for (species in unique(info$Species)){
  i <- 0
  for (tissue in unique(info$Tissue)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      if (i == 0){
        plot(density(na.omit(getMeth(allData.fit.subset[,sample], type="raw"))), col=pal[as.integer(as.factor(info[sample,]$Tissue))], main=species, xlab="Methylation percentage", lwd=1, lty=j)
        i <- 1
      }
      else {
        lines(density(na.omit(getMeth(allData.fit.subset[,sample], type="raw"))), col=pal[as.integer(as.factor(info[sample,]$Tissue))], lwd=1, lty=j)
      }
    }
  }
  legend("top", legend=unique(info$Tissue), lty=c(1,1), col=pal[1:4]) 
}
dev.off()

## density on raw estimates: using high coverage sites only
pdf(file = "distributionMethylationRaw_highCoverage.pdf", height = 6, width=6)
for (species in unique(info$Species)){
  i <- 0
  for (tissue in unique(info$Tissue)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      Cov <- getCoverage(allData.fit.subset[,sample])
      keepLoci <- which(Cov >= 5 & Cov <= 10)
      if (i == 0){
        plot(density(na.omit(getMeth(allData.fit.subset[keepLoci,sample], type="raw"))), col=pal[as.integer(as.factor(info[sample,]$Tissue))], main=species, xlab="Methylation percentage", lwd=1, ylim=c(0,11), lty=j)
        i <- 1
      }
      else {
        lines(density(na.omit(getMeth(allData.fit.subset[keepLoci,sample], type="raw"))), col=pal[as.integer(as.factor(info[sample,]$Tissue))], lwd=1, lty=j)
      }
    }
  }
  legend("top", legend=unique(info$Tissue), lty=c(1,1), col=pal[1:4]) 
}
dev.off()
## by tissue
pdf(file = "distributionMethylationRaw_highCoverage_by_tissues.pdf", height = 6, width=24)
par(mfrow=c(1,4), ps = 14, cex = 1, cex.main = 1) ## Used to keep labels at Font 14pt
for (tissue in unique(info$Tissue)){
  i <- 0
  for (species in unique(info$Species)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      Cov <- getCoverage(allData.fit.subset[,sample])
      keepLoci <- which(Cov >= 5 & Cov <= 10)
      if (i == 0){
        plot(density(na.omit(getMeth(allData.fit.subset[keepLoci,sample], type="raw"))), col=pal[as.integer(as.factor(info[sample,]$Species))+9], main=tissue, xlab="Methylation percentage", lwd=2, ylim=c(0,11), lty=j)
        i <- 1
      }
      else {
        lines(density(na.omit(getMeth(allData.fit.subset[keepLoci,sample], type="raw"))), col=pal[as.integer(as.factor(info[sample,]$Species))+9], lwd=2, lty=j)
      }
    }
  }
  if (tissue==unique(info$Tissue)[1]){
    legend("top", legend=c(as.character(unique(info$Species)), as.character(unique(info[info$Tissue == tissue,]$IndividualID))), lty=c(1,1,1,rep(c(1,2,3,4), times=3)), lwd=2, col=c(pal[10:12], rep("black",12)), bty = "n") 
  }
}
dev.off()

## density on sites overlapping orthologous exons
dim(allData.fit.subset) ## 5228469      48
head(granges(allData.fit.subset))

## first load orthologous exons hg19 coordinates
orthoExons <- read.table("../orthoExon/orthoExons_hg19.bed", h=F, sep="\t")
orthoExons <- GRanges(orthoExons[,1], IRanges(start=orthoExons[,2], end=orthoExons[,3]))
length(orthoExons) ## 187889
## overalp
length(findOverlaps(granges(allData.fit.subset), orthoExons)) ## 681371 sites
## keep subset of sites
allData.fit.subset.orthoExons <- subsetByOverlaps(allData.fit.subset, orthoExons)
## 660533 sites. Why the differences?
pdf(file = "distributionMethylationSmoothed_overlapOrthoExons.pdf", height = 6, width=6)
for (species in unique(info$Species)){
  i <- 0
  for (tissue in unique(info$Tissue)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset.orthoExons[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      if (i == 0){
        plot(density(getMeth(allData.fit.subset.orthoExons[,sample])), col=pal[as.integer(as.factor(info[sample,]$Tissue))], main=species, xlab="Methylation percentage", lwd=1, lty=j)
        i <- 1
      }
      else {
        lines(density(getMeth(allData.fit.subset.orthoExons[,sample])), col=pal[as.integer(as.factor(info[sample,]$Tissue))], lwd=1, lty=j)
      }
    }
  }
  legend("top", legend=c(as.character(unique(info$Tissue)), as.character(unique(info[info$Species == species,]$IndividualID))), lty=c(1,1,1,1,1,2,3,4), col=c(pal[1:4], rep("black",4))) 
}
dev.off()
## by tissue
pdf(file = "distributionMethylationSmoothed_overlapOrthoExons_by_tissues.pdf", height = 6, width=24)
par(mfrow=c(1,4), ps = 14, cex = 1, cex.main = 1) ## Used to keep labels at Font 14pt
for (tissue in unique(info$Tissue)){
  i <- 0
  for (species in unique(info$Species)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset.orthoExons[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      if (i == 0){
        plot(density(getMeth(allData.fit.subset.orthoExons[,sample])), col=pal[as.integer(as.factor(info[sample,]$Species))+9], main=tissue, xlab="Methylation percentage", lwd=2, lty=j, ylim=c(0,4.5))
        i <- 1
      }
      else {
        lines(density(getMeth(allData.fit.subset.orthoExons[,sample])), col=pal[as.integer(as.factor(info[sample,]$Species))+9], lwd=2, lty=j)
      }
    }
  }
  if (tissue==unique(info$Tissue)[1]){
    legend("top", legend=c(as.character(unique(info$Species)), as.character(unique(info[info$Tissue == tissue,]$IndividualID))), lty=c(1,1,1,rep(c(1,2,3,4), times=3)), lwd=2, col=c(pal[10:12], rep("black",12)), bty = "n") 
  }
}
dev.off()

## Filter exonic positions with C>T SNPs
sites <- read.table(gzfile("../annotation/orthologous_promoters/overlap_promoters_sharedPositions.gz"), h=T, sep="\t")[,1:2]
## load SNP data for each individual. We will filter positions in all individuals
## sites are 0-based
SNPs <- read.table("../genotypes/C2T/C2T_all.bed", h=F, sep="\t") ## 412345 sites
names(SNPs) <- c("chr", "start", "end")
summary(paste(sites$chromosome, sites$position, sep=":") %in% paste(SNPs$chr, SNPs$start, sep=":"))
## filter SNPs and intersect with orthoExons

## allData.fit.subset.orthoExons.noSNP <- subsetByOverlaps(allData.fit.subset[!(paste(sites$chromosome, sites$position, sep=":") %in% paste(SNPs$chr, SNPs$start, sep=":")), ], orthoExons)
## First, load plot_by_tissue function definition below
plot_by_tissue(sites[!(paste(sites$chromosome, sites$position, sep=":") %in% paste(SNPs$chr, SNPs$start, sep=":")), ], allData.fit.subset.orthoExons, info, "distributionMethylationSmoothed_overlapOrthoExonsNoSNP_by_tissues.pdf", c(0,4.5))
## 653364 sites used


## density using only sites with high coverage: Cov >= 5 & Cov <= 10
pdf(file = "distributionMethylationSmoothed_highCoverage.pdf", height = 6, width=6)
for (species in unique(info$Species)){
  i <- 0
  for (tissue in unique(info$Tissue)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      Cov <- getCoverage(allData.fit.subset[,sample])
      ## keepLoci <- which(Cov > quantile(Cov, probs = 0.9)) ## includes many very high coverage suspicious regions
      keepLoci <- which(Cov >= 5 & Cov <= 10)
      if (i == 0){
        plot(density(getMeth(allData.fit.subset[keepLoci,sample])), col=pal[as.integer(as.factor(info[sample,]$Tissue))], main=species, xlab="Methylation percentage", lwd=1, lty=j, ylim=c(0,4))
        i <- 1
      }
      else {
        lines(density(getMeth(allData.fit.subset[keepLoci,sample])), col=pal[as.integer(as.factor(info[sample,]$Tissue))], lwd=1, lty=j)
      }
    }
  }
  legend("top", legend=c(as.character(unique(info$Tissue)), as.character(unique(info[info$Species == species,]$IndividualID))), lty=c(1,1,1,1,1,2,3,4), col=c(pal[1:4], rep("black",4))) 
}
dev.off()

## density on sites with medium coverage (not too low, not too high)
pdf(file = "distributionMethylationSmoothed_mediumCoverage.pdf", height = 6, width=6)
for (species in unique(info$Species)){
  i <- 0
  for (tissue in unique(info$Tissue)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      Cov <- getCoverage(allData.fit.subset[,sample])
      keepLoci <- which(Cov >= 3 & Cov <= 7)
      if (i == 0){
        plot(density(getMeth(allData.fit.subset[keepLoci,sample])), col=pal[as.integer(as.factor(info[sample,]$Tissue))], main=species, xlab="Methylation percentage", lwd=1, lty=j, ylim=c(0,4))
        i <- 1
      }
      else {
        lines(density(getMeth(allData.fit.subset[keepLoci,sample])), col=pal[as.integer(as.factor(info[sample,]$Tissue))], lwd=1, lty=j)
      }
    }
  }
  legend("top", legend=c(as.character(unique(info$Tissue)), as.character(unique(info[info$Species == species,]$IndividualID))), lty=c(1,1,1,1,1,2,3,4), col=c(pal[1:4], rep("black",4))) 
}
dev.off()

## Yoav's suggestion for coverage filtering:
Cov <- getCoverage(allData.fit.subset)
colnames(Cov) <- colnames(allData.fit.subset)
## At least 2 out of 4 individuals per species/tissue with 10X > coverage >= 2x
keepLoci <- which(
                  ## Chimp
                  rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[1] & pData$Species == unique(pData$Species)[1]] >= 2) >= 2 & rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[2] & pData$Species == unique(pData$Species)[1]] >= 2) >= 2 & rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[3] & pData$Species == unique(pData$Species)[1]] >= 2) >= 2 & rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[4] & pData$Species == unique(pData$Species)[1]] >= 2) >= 2 &
                  ## Human
                  rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[1] & pData$Species == unique(pData$Species)[2]] >= 2) >= 2 & rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[2] & pData$Species == unique(pData$Species)[2]] >= 2) >= 2 & rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[3] & pData$Species == unique(pData$Species)[2]] >= 2) >= 2 & rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[4] & pData$Species == unique(pData$Species)[2]] >= 2) >= 2 &
                  ## Macaque
                  rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[1] & pData$Species == unique(pData$Species)[3]] >= 2) >= 2 & rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[2] & pData$Species == unique(pData$Species)[3]] >= 2) >= 2 & rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[3] & pData$Species == unique(pData$Species)[3]] >= 2) >= 2 & rowSums(Cov[, pData$Tissue == unique(pData$Tissue)[4] & pData$Species == unique(pData$Species)[3]] >= 2) >= 2 & 
## + Average coverage between 2x and 10x across all 4 individuals per species/tissue
                  rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[1] & pData$Species == unique(pData$Species)[1]]) >= 2 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[2] & pData$Species == unique(pData$Species)[1]]) >= 2 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[3] & pData$Species == unique(pData$Species)[1]]) >= 2 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[4] & pData$Species == unique(pData$Species)[1]]) >= 2 &
                  rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[1] & pData$Species == unique(pData$Species)[1]]) < 10 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[2] & pData$Species == unique(pData$Species)[1]]) < 10 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[3] & pData$Species == unique(pData$Species)[1]]) < 10 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[4] & pData$Species == unique(pData$Species)[1]]) < 10 &
                  rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[1] & pData$Species == unique(pData$Species)[2]]) >= 2 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[2] & pData$Species == unique(pData$Species)[2]]) >= 2 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[3] & pData$Species == unique(pData$Species)[2]]) >= 2 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[4] & pData$Species == unique(pData$Species)[2]]) >= 2 &
                  rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[1] & pData$Species == unique(pData$Species)[1]]) < 10 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[2] & pData$Species == unique(pData$Species)[1]]) < 10 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[3] & pData$Species == unique(pData$Species)[1]]) < 10 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[4] & pData$Species == unique(pData$Species)[1]]) < 10 &
                  rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[1] & pData$Species == unique(pData$Species)[3]]) >= 2 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[2] & pData$Species == unique(pData$Species)[3]]) >= 2 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[3] & pData$Species == unique(pData$Species)[3]]) >= 2 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[4] & pData$Species == unique(pData$Species)[3]]) >= 2 &
                  rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[1] & pData$Species == unique(pData$Species)[3]]) < 10 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[2] & pData$Species == unique(pData$Species)[3]]) < 10 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[3] & pData$Species == unique(pData$Species)[3]]) < 10 & rowMeans(Cov[, pData$Tissue == unique(pData$Tissue)[4] & pData$Species == unique(pData$Species)[3]]) <10)

length(keepLoci) ## 2928421 sites: not very stringent
allData.fit.subset.filtered <- allData.fit.subset[keepLoci,]

pdf(file = "distributionMethylationSmoothed_evenCoverage.pdf", height = 6, width=6)
for (species in unique(info$Species)){
  i <- 0
  for (tissue in unique(info$Tissue)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      if (i == 0){
        plot(density(getMeth(allData.fit.subset.filtered[,sample])), col=pal[as.integer(as.factor(info[sample,]$Tissue))], main=species, xlab="Methylation percentage", lwd=1, lty=j, ylim=c(0,4))
        i <- 1
      }
      else {
        lines(density(getMeth(allData.fit.subset.filtered[,sample])), col=pal[as.integer(as.factor(info[sample,]$Tissue))], lwd=1, lty=j)
      }
    }
  }
  legend("top", legend=c(as.character(unique(info$Tissue)), as.character(unique(info[info$Species == species,]$IndividualID))), lty=c(1,1,1,1,1,2,3,4), col=c(pal[1:4], rep("black",4))) 
}
dev.off()
## by tissue
pdf(file = "distributionMethylationSmoothed_evenCoverage_by_tissues.pdf", height = 6, width=24)
par(mfrow=c(1,4), ps = 14, cex = 1, cex.main = 1) ## Used to keep labels at Font 14pt
for (tissue in unique(info$Tissue)){
  i <- 0
  for (species in unique(info$Species)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset.filtered[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      if (i == 0){
        plot(density(getMeth(allData.fit.subset.filtered[,sample])), col=pal[as.integer(as.factor(info[sample,]$Species))+9], main=tissue, xlab="Methylation percentage", lwd=2, lty=j, ylim=c(0,4.5))
        i <- 1
      }
      else {
        lines(density(getMeth(allData.fit.subset.filtered[,sample])), col=pal[as.integer(as.factor(info[sample,]$Species))+9], lwd=2, lty=j)
      }
    }
  }
  if (tissue==unique(info$Tissue)[1]){
    legend("top", legend=c(as.character(unique(info$Species)), as.character(unique(info[info$Tissue == tissue,]$IndividualID))), lty=c(1,1,1,rep(c(1,2,3,4), times=3)), lwd=2, col=c(pal[10:12], rep("black",12)), bty = "n") 
  }
}
dev.off()

## We can be even more stringent! 10X>Coverage>=2X in all samples?
keepLoci <- which(rowSums(Cov >= 2) == 48)
## 136160 sites. Looking at coverage it seems thta we need an upper bound!
keepLoci <- which(rowSums(Cov >= 2) == 48 & rowSums(Cov < 10) == 48)
## 3363 sites only: requiring all sites to >=2X selects mostly repeats / CNVs
allData.fit.subset.filtered <- allData.fit.subset[keepLoci,]
pdf(file = "distributionMethylationSmoothed_evenCoverage2.pdf", height = 6, width=6)
for (species in unique(info$Species)){
  i <- 0
  for (tissue in unique(info$Tissue)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      if (i == 0){
        plot(density(getMeth(allData.fit.subset.filtered[,sample])), col=pal[as.integer(as.factor(info[sample,]$Tissue))], main=species, xlab="Methylation percentage", lwd=1, lty=j, ylim=c(0,4))
        i <- 1
      }
      else {
        lines(density(getMeth(allData.fit.subset.filtered[,sample])), col=pal[as.integer(as.factor(info[sample,]$Tissue))], lwd=1, lty=j)
      }
    }
  }
  legend("top", legend=c(as.character(unique(info$Tissue)), as.character(unique(info[info$Species == species,]$IndividualID))), lty=c(1,1,1,1,1,2,3,4), col=c(pal[1:4], rep("black",4))) 
}
dev.off()
## by tissue
pdf(file = "distributionMethylationSmoothed_evenCoverage2_by_tissues.pdf", height = 6, width=24)
par(mfrow=c(1,4), ps = 14, cex = 1, cex.main = 1) ## Used to keep labels at Font 14pt
for (tissue in unique(info$Tissue)){
  i <- 0
  for (species in unique(info$Species)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset.filtered[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      if (i == 0){
        plot(density(getMeth(allData.fit.subset.filtered[,sample])), col=pal[as.integer(as.factor(info[sample,]$Species))+9], main=tissue, xlab="Methylation percentage", lwd=2, lty=j, ylim=c(0,4.5))
        i <- 1
      }
      else {
        lines(density(getMeth(allData.fit.subset.filtered[,sample])), col=pal[as.integer(as.factor(info[sample,]$Species))+9], lwd=2, lty=j)
      }
    }
  }
  if (tissue==unique(info$Tissue)[1]){
    legend("top", legend=c(as.character(unique(info$Species)), as.character(unique(info[info$Tissue == tissue,]$IndividualID))), lty=c(1,1,1,rep(c(1,2,3,4), times=3)), lwd=2, col=c(pal[10:12], rep("black",12)), bty = "n") 
  }
}
dev.off()
## Still with this reduced number of sites, the trend is still clear!


## Density plots using sites overlaping transposon-free regions
## http://www.biomedcentral.com/content/supplementary/1471-2164-8-470-s2.txt
## Converted to hg19 with https://genome.ucsc.edu/cgi-bin/hgLiftOver
TFfree <- read.table("../annotation/Transposon_free_regions/transposon_free_regions_hg19.bed", h=F, sep="\t")
TFfree <- GRanges(TFfree[,1], IRanges(start=TFfree[,2], end=TFfree[,3]))
length(TFfree) ## 9199
## overalp
length(findOverlaps(granges(allData.fit.subset), TFfree)) ## 471892 sites
## keep subset of sites
allData.fit.subset.TFfree <- subsetByOverlaps(allData.fit.subset, TFfree)
## 471873 sites. Why the differences?
pdf(file = "distributionMethylationSmoothed_overlapTransposonFree.pdf", height = 6, width=6)
for (species in unique(info$Species)){
  i <- 0
  for (tissue in unique(info$Tissue)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset.TFfree[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      if (i == 0){
        plot(density(getMeth(allData.fit.subset.TFfree[,sample])), col=pal[as.integer(as.factor(info[sample,]$Tissue))], main=species, xlab="Methylation percentage", lwd=1, lty=j, ylim=c(0,6))
        i <- 1
      }
      else {
        lines(density(getMeth(allData.fit.subset.TFfree[,sample])), col=pal[as.integer(as.factor(info[sample,]$Tissue))], lwd=1, lty=j)
      }
    }
  }
  legend("top", legend=c(as.character(unique(info$Tissue)), as.character(unique(info[info$Species == species,]$IndividualID))), lty=c(1,1,1,1,1,2,3,4), col=c(pal[1:4], rep("black",4))) 
}
dev.off()
## by tissue
pdf(file = "distributionMethylationSmoothed_overlapTransposonFree_by_tissues.pdf", height = 6, width=24)
par(mfrow=c(1,4), ps = 14, cex = 1, cex.main = 1) ## Used to keep labels at Font 14pt
for (tissue in unique(info$Tissue)){
  i <- 0
  for (species in unique(info$Species)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset.TFfree[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      if (i == 0){
        plot(density(getMeth(allData.fit.subset.TFfree[,sample])), col=pal[as.integer(as.factor(info[sample,]$Species))+9], main=tissue, xlab="Methylation percentage", lwd=2, lty=j, ylim=c(0,6))
        i <- 1
      }
      else {
        lines(density(getMeth(allData.fit.subset.TFfree[,sample])), col=pal[as.integer(as.factor(info[sample,]$Species))+9], lwd=2, lty=j)
      }
    }
  }
  if (tissue==unique(info$Tissue)[1]){
    legend("top", legend=c(as.character(unique(info$Species)), as.character(unique(info[info$Tissue == tissue,]$IndividualID))), lty=c(1,1,1,rep(c(1,2,3,4), times=3)), lwd=2, col=c(pal[10:12], rep("black",12)), bty = "n") 
  }
}
dev.off()
## Problem: we have no guarantee that Transposon-free regions are conserved in chimp and macaque. A better apporach in to lookat regions where repeats are perfectly conserevd (shared) across the three species:

## Conserved repeats (conserved regions only)
repeats <- read.table(gzfile("../annotation/orthologous_repeats/overlap_repeats_sharedPositions.gz"), h=T, sep="\t")
repeats <- repeats[!is.na(repeats["repeats_conserved_panTro3.rheMac2"]), ]
dim(repeats) ## 759431 sites
plot_by_tissue(repeats[,1:2], allData.fit.subset, info, "distributionMethylationSmoothed_overlapConservedRepeats_by_tissues.pdf", c(0,6))

repeats_unique <- aggregate(repeats[,1:2], list(repeats[,"repeats_conserved_panTro3.rheMac2"]), function(x){ return( sample(x, 1) )})
dim(repeats_unique) ## 
plot_by_tissue(repeats_unique[,2:3], allData.fit.subset, info, "distributionMethylationSmoothed_overlapConservedRepeatsUnique_by_tissues.pdf", c(0,6))


## Density plots using sites overlaping TSS +/2kb regions
promoter <- read.table(gzfile("../annotation/15_features_sharedPositions.gz"), h=T, sep="\t")[,c(1,2,15)]
promoter <- GRanges(promoter[promoter$promoter==1, 1], IRanges(start=promoter[promoter$promoter==1, 2], end=promoter[promoter$promoter==1, 2]+1))
length(promoter) ## 1719690
## keep subset of sites
allData.fit.subset.promoter <- subsetByOverlaps(allData.fit.subset, promoter)
pdf(file = "distributionMethylationSmoothed_overlapPromoter.pdf", height = 6, width=6)
for (species in unique(info$Species)){
  i <- 0
  for (tissue in unique(info$Tissue)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset.promoter[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      if (i == 0){
        plot(density(getMeth(allData.fit.subset.promoter[,sample])), col=pal[as.integer(as.factor(info[sample,]$Tissue))], main=species, xlab="Methylation percentage", lwd=1, lty=j, ylim=c(0,8))
        i <- 1
      }
      else {
        lines(density(getMeth(allData.fit.subset.promoter[,sample])), col=pal[as.integer(as.factor(info[sample,]$Tissue))], lwd=1, lty=j)
      }
    }
  }
  legend("top", legend=c(as.character(unique(info$Tissue)), as.character(unique(info[info$Species == species,]$IndividualID))), lty=c(1,1,1,1,1,2,3,4), col=c(pal[1:4], rep("black",4))) 
}
dev.off()
## by tissue
pdf(file = "distributionMethylationSmoothed_overlapPromoter_by_tissues.pdf", height = 6, width=24)
par(mfrow=c(1,4), ps = 14, cex = 1, cex.main = 1) ## Used to keep labels at Font 14pt
for (tissue in unique(info$Tissue)){
  i <- 0
  for (species in unique(info$Species)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset.promoter[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      if (i == 0){
        plot(density(getMeth(allData.fit.subset.promoter[,sample])), col=pal[as.integer(as.factor(info[sample,]$Species))+9], main=tissue, xlab="Methylation percentage", lwd=2, lty=j, ylim=c(0,8))
        i <- 1
      }
      else {
        lines(density(getMeth(allData.fit.subset.promoter[,sample])), col=pal[as.integer(as.factor(info[sample,]$Species))+9], lwd=2, lty=j)
      }
    }
  }
  if (tissue==unique(info$Tissue)[1]){
    legend("top", legend=c(as.character(unique(info$Species)), as.character(unique(info[info$Tissue == tissue,]$IndividualID))), lty=c(1,1,1,rep(c(1,2,3,4), times=3)), lwd=2, col=c(pal[10:12], rep("black",12)), bty = "n") 
  }
}
dev.off()


## Density plots using sites overlaping proximal promoters +/- 250 nt from TSS
promoter <- read.table(gzfile("../annotation/15_features_sharedPositions.gz"), h=T, sep="\t")[,c(1,2,17)]
promoter <- GRanges(promoter[promoter$promoter==1, 1], IRanges(start=promoter[promoter$promoter==1, 2], end=promoter[promoter$promoter==1, 2]+1))
length(promoter) ## 668118
## keep subset of sites
allData.fit.subset.promoter <- subsetByOverlaps(allData.fit.subset, promoter)
pdf(file = "distributionMethylationSmoothed_overlapProximalPromoter.pdf", height = 6, width=6)
for (species in unique(info$Species)){
  i <- 0
  for (tissue in unique(info$Tissue)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset.promoter[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      if (i == 0){
        plot(density(getMeth(allData.fit.subset.promoter[,sample])), col=pal[as.integer(as.factor(info[sample,]$Tissue))], main=species, xlab="Methylation percentage", lwd=1, lty=j, ylim=c(0,11))
        i <- 1
      }
      else {
        lines(density(getMeth(allData.fit.subset.promoter[,sample])), col=pal[as.integer(as.factor(info[sample,]$Tissue))], lwd=1, lty=j)
      }
    }
  }
  legend("top", legend=c(as.character(unique(info$Tissue)), as.character(unique(info[info$Species == species,]$IndividualID))), lty=c(1,1,1,1,1,2,3,4), col=c(pal[1:4], rep("black",4))) 
}
dev.off()
## by tissue
pdf(file = "distributionMethylationSmoothed_overlapProximalPromoter_by_tissues.pdf", height = 6, width=24)
par(mfrow=c(1,4), ps = 14, cex = 1, cex.main = 1) ## Used to keep labels at Font 14pt
for (tissue in unique(info$Tissue)){
  i <- 0
  for (species in unique(info$Species)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset.promoter[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      if (i == 0){
        plot(density(getMeth(allData.fit.subset.promoter[,sample])), col=pal[as.integer(as.factor(info[sample,]$Species))+9], main=tissue, xlab="Methylation percentage", lwd=2, lty=j, ylim=c(0,11))
        i <- 1
      }
      else {
        lines(density(getMeth(allData.fit.subset.promoter[,sample])), col=pal[as.integer(as.factor(info[sample,]$Species))+9], lwd=2, lty=j)
      }
    }
  }
  if (tissue==unique(info$Tissue)[1]){
    legend("top", legend=c(as.character(unique(info$Species)), as.character(unique(info[info$Tissue == tissue,]$IndividualID))), lty=c(1,1,1,rep(c(1,2,3,4), times=3)), lwd=2, col=c(pal[10:12], rep("black",12)), bty = "n") 
  }
}
dev.off()

## Conserved promoters (conserved regions only)
promoter <- read.table(gzfile("../annotation/orthologous_promoters/overlap_promoters_sharedPositions.gz"), h=T, sep="\t")
promoter <- GRanges(promoter[!is.na(promoter["promoter_conserved_panTro3.rheMac2"]), 1], IRanges(start=promoter[!is.na(promoter["promoter_conserved_panTro3.rheMac2"]), 2], end=promoter[!is.na(promoter["promoter_conserved_panTro3.rheMac2"]), 2]+1))
length(promoter) ## 185270
## keep subset of sites
allData.fit.subset.promoter <- subsetByOverlaps(allData.fit.subset, promoter)
pdf(file = "distributionMethylationSmoothed_overlapConservedPromoter_by_tissues.pdf", height = 6, width=24)
par(mfrow=c(1,4), ps = 14, cex = 1, cex.main = 1) ## Used to keep labels at Font 14pt
for (tissue in unique(info$Tissue)){
  i <- 0
  for (species in unique(info$Species)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset.promoter[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      if (i == 0){
        plot(density(getMeth(allData.fit.subset.promoter[,sample])), col=pal[as.integer(as.factor(info[sample,]$Species))+9], main=tissue, xlab="Methylation percentage", lwd=2, lty=j, ylim=c(0,8))
        i <- 1
      }
      else {
        lines(density(getMeth(allData.fit.subset.promoter[,sample])), col=pal[as.integer(as.factor(info[sample,]$Species))+9], lwd=2, lty=j)
      }
    }
  }
  if (tissue==unique(info$Tissue)[1]){
    legend("top", legend=c(as.character(unique(info$Species)), as.character(unique(info[info$Tissue == tissue,]$IndividualID))), lty=c(1,1,1,rep(c(1,2,3,4), times=3)), lwd=2, col=c(pal[10:12], rep("black",12)), bty = "n") 
  }
}
dev.off()

## Conserved promoters (conserved regions only) + taking one CpG per promoter randomly
promoter <- read.table(gzfile("../annotation/orthologous_promoters/overlap_promoters_sharedPositions.gz"), h=T, sep="\t")
promoter <- promoter[!is.na(promoter["promoter_conserved_panTro3.rheMac2"]), ]
promoter_unique <- aggregate(promoter[,1:2], list(promoter[,"promoter_conserved_panTro3.rheMac2"]), function(x){ return( sample(x, 1) )})
dim(promoter_unique) ## 5257
## Function to make plotting simpler (less text)
plot_by_tissue(promoter_unique[,2:3], allData.fit.subset, info, "distributionMethylationSmoothed_overlapConservedPromoterUnique_by_tissues.pdf", c(0,3))

plot_by_tissue <- function(sites, data, info, pdfFile, ylimit=NA){
  subset <- GRanges(sites[, 1], IRanges(start=sites[, 2], end=sites[, 2]+1))
  ## keep subset of sites
  data <- subsetByOverlaps(data, subset)
  cat(length(data), "sites used\n")
  
  ## Plotting
  if(is.na(ylimit[1])){
    ylimit <- c(0, max(density(getMeth(data[,"H1H"]))$y) * 1.2)
  }
  pdf(file = pdfFile, height = 6, width=24)
  par(mfrow=c(1,4), ps = 14, cex = 1, cex.main = 1) ## Used to keep labels at Font 14pt
  for (tissue in unique(info$Tissue)){
    i <- 0
    for (species in unique(info$Species)){
      j <- 0
      for (sample in sampleNames(data[, info$Species == species & info$Tissue == tissue])){
        cat(sample, "\n")
        j <- j+1
        if (i == 0){
          plot(density(getMeth(data[,sample])), col=pal[as.integer(as.factor(info[sample,]$Species))+9], main=tissue, xlab="Methylation percentage", lwd=2, lty=j, ylim=ylimit)
          i <- 1
        }
        else {
          lines(density(getMeth(data[,sample])), col=pal[as.integer(as.factor(info[sample,]$Species))+9], lwd=2, lty=j)
        }
      }
    }
    if (tissue==unique(info$Tissue)[1]){
      legend("top", legend=c(as.character(unique(info$Species)), as.character(unique(info[info$Tissue == tissue,]$IndividualID))), lty=c(1,1,1,rep(c(1,2,3,4), times=3)), lwd=2, col=c(pal[10:12], rep("black",12)), bty = "n") 
    }
  }
  dev.off()
}

## TO DO? promoters conserved only in human/chimp?

## Conserved CGI (conserved regions only)
CGI <- read.table(gzfile("../annotation/orthologous_CGI/overlap_CGI_sharedPositions.gz"), h=T, sep="\t")
CGI <- CGI[!is.na(CGI["CGI_conserved_panTro3.rheMac2"]), ]
dim(CGI) ## 589693 sites
plot_by_tissue(CGI[,1:2], allData.fit.subset, info, "distributionMethylationSmoothed_overlapConservedCGI_by_tissues.pdf", c(0,50))

CGI_unique <- aggregate(CGI[,1:2], list(CGI[,"CGI_conserved_panTro3.rheMac2"]), function(x){ return( sample(x, 1) )})
dim(CGI_unique) ## 12996
plot_by_tissue(CGI_unique[,2:3], allData.fit.subset, info, "distributionMethylationSmoothed_overlapConservedCGIUnique_by_tissues.pdf", c(0,25))

## Conserved promoters and conserved CGI
promoter <- read.table(gzfile("../annotation/orthologous_promoters/overlap_promoters_sharedPositions.gz"), h=T, sep="\t")
CGI <- read.table(gzfile("../annotation/orthologous_CGI/overlap_CGI_sharedPositions.gz"), h=T, sep="\t")
promoter <- promoter[!is.na(promoter["promoter_conserved_panTro3.rheMac2"]) & !is.na(CGI["CGI_conserved_panTro3.rheMac2"]), ] ## 89492 sites
plot_by_tissue(promoter[,1:2], allData.fit.subset, info, "distributionMethylationSmoothed_overlapConservedPromoterConservedCGI_by_tissues.pdf", NA)
## Taking a unique sites per conserved promoter
promoter_unique <- aggregate(promoter[,1:2], list(promoter[,"promoter_conserved_panTro3.rheMac2"]), function(x){ return( sample(x, 1) )})
dim(promoter_unique) ## 1830
## Function to make plotting simpler (less text)
plot_by_tissue(promoter[,1:2], allData.fit.subset, info, "distributionMethylationSmoothed_overlapConservedPromoterConservedCGIUnique_by_tissues.pdf", NA)


## Filter most variable sites
allSds <- apply(getMeth(allData.fit.subset), 1, sd) ## long!
pdf(file = "temp.pdf", height = 6, width=6)
hist(allSds, breaks=100)
hist(log10(allSds), breaks=100)
dev.off()

sites <- read.table(gzfile("../annotation/orthologous_promoters/overlap_promoters_sharedPositions.gz"), h=T, sep="\t")[,1:2]
plot_by_tissue(sites[allSds < 0.02, ], allData.fit.subset, info, "distributionMethylationSmoothed_lowSd_by_tissues.pdf", NA) ## 542741 sites
## Only unmethylated sites. We still see a very slight difference between species in liver and kidney
plot_by_tissue(sites[allSds > 0.03 & allSds < 0.1, ], allData.fit.subset, info, "distributionMethylationSmoothed_mediumSd_by_tissues.pdf", c(0,5.5)) ## 3542368 sites


## Density plot taking average per 1kb window
sites <- read.table(gzfile("../annotation/sharedPositions.bed.gz"), h=T, sep="\t")[,1:2]
sites <- GRanges(sites[, 1], IRanges(start=sites[, 2], end=sites[, 2]))

## keep autosomal chromosomes only
seqLength <- read.table("../annotation/hg19.genome", h=T)[c(1:7, 9:20, 22:24),2]
names(seqLength) <- read.table("../annotation/hg19.genome", h=T)[c(1:7, 9:20, 22:24),1]
## create 1kb tiles 
tiles <- unlist(tileGenome(seqLength, tilewidth=1000))
## overlap with sites
overlaps <- findOverlaps(sites, tiles)
sites <- sites[overlaps@queryHits]
## subset data
allData.fit.subset.sites <- subsetByOverlaps(allData.fit.subset, sites)

pdf(file = "distributionMethylationSmoothed_average1kbWindows.pdf", height = 6, width=24)
par(mfrow=c(1,4), ps = 14, cex = 1, cex.main = 1) ## Used to keep labels at Font 14pt
for (tissue in unique(info$Tissue)){
  i <- 0
  for (species in unique(info$Species)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset.sites[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      if (i == 0){
        plot(density(aggregate(getMeth(allData.fit.subset.sites[,sample]), list(overlaps@subjectHits), mean)[,2]), col=pal[as.integer(as.factor(info[sample,]$Species))+9], main=tissue, xlab="Methylation percentage", lwd=2, lty=j, ylim=c(0,6))
        i <- 1
      }
      else {
        lines(density(aggregate(getMeth(allData.fit.subset.sites[,sample]), list(overlaps@subjectHits), mean)[,2]), col=pal[as.integer(as.factor(info[sample,]$Species))+9], lwd=2, lty=j)
      }
    }
  }
  if (tissue==unique(info$Tissue)[1]){
    legend("top", legend=c(as.character(unique(info$Species)), as.character(unique(info[info$Tissue == tissue,]$IndividualID))), lty=c(1,1,1,rep(c(1,2,3,4), times=3)), lwd=2, col=c(pal[10:12], rep("black",12)), bty = "n") 
  }
}
dev.off()
## TO DO: change window size?

## Density plot taking 1 random CpG site per 1kb window

## Different ways to proceed :
## randomSites <- aggregate(as.data.frame(sites)[,1:2], list(overlaps@subjectHits), function(x){ return( sample(x, 1) )})
## randomSites <- tapply(overlaps@queryHits, list(overlaps@subjectHits), function(x){ return( sample(x, 1) )})
## randomSites <- aggregate(overlaps@queryHits, list(overlaps@subjectHits), function(x){ return( sample(x, 1) )})
## These seem to give strange results. For example first site is chosen in 2 different categories!?! No time to investigate... But this works fine:
randomSites <- do.call( rbind, lapply( split(overlaps@queryHits, overlaps@subjectHits), function(x) x[sample(1:length(x), 1)] ))
length(randomSites) ## 1186867
## subset data
allData.fit.subset.sites <- allData.fit.subset[randomSites,]

pdf(file = "distributionMethylationSmoothed_randomSiteIn1kbWindows.pdf", height = 6, width=24)
par(mfrow=c(1,4), ps = 14, cex = 1, cex.main = 1) ## Used to keep labels at Font 14pt
for (tissue in unique(info$Tissue)){
  i <- 0
  for (species in unique(info$Species)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset.sites[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      if (i == 0){
        plot(density(getMeth(allData.fit.subset.sites[,sample])), col=pal[as.integer(as.factor(info[sample,]$Species))+9], main=tissue, xlab="Methylation percentage", lwd=2, lty=j, ylim=c(0,8))
        i <- 1
      }
      else {
        lines(density(getMeth(allData.fit.subset.sites[,sample])), col=pal[as.integer(as.factor(info[sample,]$Species))+9], lwd=2, lty=j)
      }
    }
  }
  if (tissue==unique(info$Tissue)[1]){
    legend("top", legend=c(as.character(unique(info$Species)), as.character(unique(info[info$Tissue == tissue,]$IndividualID))), lty=c(1,1,1,rep(c(1,2,3,4), times=3)), lwd=2, col=c(pal[10:12], rep("black",12)), bty = "n") 
  }
}
dev.off()

## Density plot taking average per 1kb window, but using only selected windows #########################
## - where CpG density didn't change
## - where CpG density stayed above a critical threshold
## - ...

## First load data summed on 1kb tiles
density1kb <- read.table("../CpG_maps_to_hg19/density_allCpGs_1kb.txt.gz", h=T, sep="\t")
## 3,012,613 tiles

## Only use windows with a least 1 human, 1 chimp and 1 rhesus CpG detected
## At least we can be sure there is some orthology
density1kb <- density1kb[density1kb$numHuman > 0 & density1kb$numChimp > 0 & density1kb$numRhesus > 0, ]
## 2,345,232 tiles left

## Density human: distribution
pdf(file = "density_allCpGs_1kb.pdf", height = 6, width=24)
par(mfrow=c(1,4), ps = 14, cex = 1, cex.main = 1) ## Used to keep labels at Font 14pt
hist(density1kb$numHuman, breaks=100, main="")
hist(log10(density1kb$numHuman), breaks=100, main="")
## peak at n=5
hist(log10(density1kb$numChimp), breaks=100, main="")
hist(log10(density1kb$numRhesus), breaks=100, main="")
## Distributions are quite similar for chimp and rhesus!

## New line:
## Number of CpGs compared to rhesus (lost in human/chimp or gained in rhesus)
plot(density(log10(density1kb$numHumanRefRhesus)), col=pal[10], lwd=2, xlab="", main="")
lines(density(log10(density1kb$numChimpRefRhesus)), col=pal[11], lwd=2)
lines(density(log10(density1kb$numRhesus)), col=pal[12], lwd=2)
## Virtually identical between human and chimp: good!

## Evolution of absolute number of CpGs compared to rhesus
plot(density(density1kb$numRhesus - density1kb$numHumanRefRhesus), col=pal[10], lwd=2, xlab="", main="")
lines(density(density1kb$numRhesus - density1kb$numChimpRefRhesus), col=pal[11], lwd=2, xlab="", main="")
lines(density(density1kb$numRhesus), col=pal[12], lwd=2)

## Evolution of relative number of CpGs compared to rhesus
plot(density((density1kb$numRhesus - density1kb$numHumanRefRhesus)/density1kb$numRhesus), col=pal[10], lwd=2, xlab="", main="")
lines(density((density1kb$numRhesus - density1kb$numChimpRefRhesus)/density1kb$numRhesus), col=pal[11], lwd=2, xlab="", main="")
## 0: no change compared to rhesus
## 1: all positions changed
## We cannot go to negative values since we start from CpG positions in macaque that have orthology
## -> Slight differences between human and chimp, but it seems minor. Apparently we got rid of the orthology problem

## Same thing but basic coverage filtering: at least 5 macaque CpGs
plot(density((density1kb$numRhesus[density1kb$numRhesus >= 5] - density1kb$numHumanRefRhesus[density1kb$numRhesus >= 5])/density1kb$numRhesus[density1kb$numRhesus >= 5]), col=pal[10], lwd=2, xlab="", main="")
lines(density((density1kb$numRhesus[density1kb$numRhesus >= 5] - density1kb$numChimpRefRhesus[density1kb$numRhesus >= 5])/density1kb$numRhesus[density1kb$numRhesus >= 5]), col=pal[11], lwd=2, xlab="", main="")
## Removes a large part of 0 and 1. Human and chimp still very similar

## How many tiles left when filtering based on these criteria?
summary(density1kb$numRhesus >= 5)
## 1,767,417 tiles
summary((density1kb$numRhesus[density1kb$numRhesus >= 5] - density1kb$numHumanRefRhesus[density1kb$numRhesus >= 5])/density1kb$numRhesus[density1kb$numRhesus >= 5] < 0.2)
## 38,690 tiles only
summary((density1kb$numRhesus[density1kb$numRhesus >= 5] - density1kb$numHumanRefRhesus[density1kb$numRhesus >= 5])/density1kb$numRhesus[density1kb$numRhesus >= 5] < 0.2 & (density1kb$numRhesus[density1kb$numRhesus >= 5] - density1kb$numChimpRefRhesus[density1kb$numRhesus >= 5])/density1kb$numRhesus[density1kb$numRhesus >= 5] < 0.2)
## 28,961

## Load these tiles in a granges object
tiles <- density1kb[density1kb$numRhesus >= 5, ]
tiles <- tiles[(density1kb$numRhesus[density1kb$numRhesus >= 5] - density1kb$numHumanRefRhesus[density1kb$numRhesus >= 5])/density1kb$numRhesus[density1kb$numRhesus >= 5] < 0.2 & (density1kb$numRhesus[density1kb$numRhesus >= 5] - density1kb$numChimpRefRhesus[density1kb$numRhesus >= 5])/density1kb$numRhesus[density1kb$numRhesus >= 5] < 0.2, ]
write.table(tiles, "../CpG_maps_to_hg19/1kb_tiles_stable_density.txt", sep="\t", quote=F, row.names=F, col.names=T)
tiles <- GRanges(tiles[, 1], IRanges(start=tiles[, 2], end=tiles[, 2]))

## New line:
## Evolution of relative number of CpGs compared to rhesus
plot(density(density1kb$numHumanRefRhesus/density1kb$numRhesus), col=pal[10], lwd=2, xlab="", main="")
lines(density(density1kb$numChimpRefRhesus/density1kb$numRhesus), col=pal[11], lwd=2, xlab="", main="")
## 1: no change compared to rhesus
## 0: all positions changed
## We cannot go above 1 since we start from CpG positions in macaque that have orthology

## Same thing but basic coverage filtering: at least 10 macaque CpGs
plot(density(density1kb$numHumanRefRhesus[density1kb$numRhesus >= 10]/density1kb$numRhesus[density1kb$numRhesus >= 10]), col=pal[10], lwd=2, xlab="", main="")
lines(density(density1kb$numChimpRefRhesus[density1kb$numRhesus >= 10]/density1kb$numRhesus[density1kb$numRhesus >= 10]), col=pal[11], lwd=2, xlab="", main="")
## divergence between human and chimp, with rhesus as reference
plot(density(density1kb$numChimpRefRhesus[density1kb$numRhesus >= 10] - density1kb$numHumanRefRhesus[density1kb$numRhesus >= 10]), col=pal[10], lwd=2, xlab="", main="")
## divergence between human and chimp, with rhesus as reference
plot(density(density1kb$numChimpRefRhesus[density1kb$numRhesus >= 10]/density1kb$numRhesus[density1kb$numRhesus >= 10] - density1kb$numHumanRefRhesus[density1kb$numRhesus >= 10]/density1kb$numRhesus[density1kb$numRhesus >= 10]), col=pal[10], lwd=2, xlab="", main="")

## restrict divergence between human and chimp to 0: what is divergence to rhesus looking like?
plot(density(density1kb$numHumanRefRhesus[density1kb$numRhesus >= 10 & density1kb$numHumanRefRhesus == density1kb$numChimpRefRhesus]/density1kb$numRhesus[density1kb$numRhesus >= 10 & density1kb$numHumanRefRhesus == density1kb$numChimpRefRhesus]), col=pal[10], lwd=2, xlab="", main="")
lines(density(density1kb$numChimpRefRhesus[density1kb$numRhesus >= 10 & density1kb$numHumanRefRhesus == density1kb$numChimpRefRhesus]/density1kb$numRhesus[density1kb$numRhesus >= 10 & density1kb$numHumanRefRhesus == density1kb$numChimpRefRhesus]), col=pal[11], lwd=2, xlab="", main="")
## pretty much the same

## what are distributions looking like?
plot(density(density1kb$numHuman[density1kb$numRhesus >= 10 & density1kb$numHumanRefRhesus == density1kb$numChimpRefRhesus]), col=pal[10], lwd=2, xlab="", main="")
lines(density(density1kb$numChimp[density1kb$numRhesus >= 10 & density1kb$numHumanRefRhesus == density1kb$numChimpRefRhesus]), col=pal[11], lwd=2, xlab="", main="")
lines(density(density1kb$numRhesus[density1kb$numRhesus >= 10 & density1kb$numHumanRefRhesus == density1kb$numChimpRefRhesus]), col=pal[12], lwd=2, xlab="", main="")
dev.off()



##        study distribtuion for gradient of changes in CpG density
##        restrict to CpGs thta are shared betweeb human and chimp: mor elikely to be gaine din rhesus




## Then overlap selected tiles with shared CpG sites
## allSites <- read.table(gzfile("../annotation/sharedPositions.bed.gz"), h=T, sep="\t")[,1:2]
## allSites <- GRanges(sites[, 1], IRanges(start=sites[, 2], end=sites[, 2]))
overlaps <- findOverlaps(allSites, tiles)
sites <- allSites[overlaps@queryHits]
## subset data
allData.fit.subset.sites <- subsetByOverlaps(allData.fit.subset, sites)
pdf(file = "distributionMethylationSmoothed_average1kbWindows_considerDensity.pdf", height = 6, width=24)
par(mfrow=c(1,4), ps = 14, cex = 1, cex.main = 1) ## Used to keep labels at Font 14pt
## Record mean methylation at tile positions to plot heatmap below
meanMeth <- matrix(nrow=length(overlaps@subjectHits), ncol=48)
colnames(meanMeth) <- unique(info$Condition)
for (tissue in unique(info$Tissue)){
  i <- 0
  for (species in unique(info$Species)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset.sites[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      if (i == 0){
        plot(density(aggregate(getMeth(allData.fit.subset.sites[,sample]), list(overlaps@subjectHits), mean)[,2]), col=pal[as.integer(as.factor(info[sample,]$Species))+9], main=tissue, xlab="Methylation percentage", lwd=2, lty=j, ylim=c(0,8))
        i <- 1
      }
      else {
        lines(density(aggregate(getMeth(allData.fit.subset.sites[,sample]), list(overlaps@subjectHits), mean)[,2]), col=pal[as.integer(as.factor(info[sample,]$Species))+9], lwd=2, lty=j)
      }
      meanMeth[,sample] <- aggregate(getMeth(allData.fit.subset.sites[,sample]), list(overlaps@subjectHits), mean)[,2]
    }
  }
  if (tissue==unique(info$Tissue)[1]){
    legend("top", legend=c(as.character(unique(info$Species)), as.character(unique(info[info$Tissue == tissue,]$IndividualID))), lty=c(1,1,1,rep(c(1,2,3,4), times=3)), lwd=2, col=c(pal[10:12], rep("black",12)), bty = "n") 
  }
}
dev.off()
## Few sites overlap these tiles (664), so the profiles is noisier than for previous plots, but there seem to be no or very few differences across species: good!
## Indeed the changes in CpG content drove the differences in methylation levels distribution across species
## The distribution is shifted to the left, suggesting that many of these regions were conserved because of low methylation...
## TO DO: we might want to exclude tiles overlapping CGI, promoters, etc since these are highy conserved parts of the genome because of function... 

## Plot heatmap mean methylation on tiles where CpG content stable: should be clustering more by tissue
library(gplots)
colors <- colorRampPalette(c(brewer.pal(9, "Blues")[1],brewer.pal(9, "Blues")[9]))(100)
pdf(file = "heatmap_1kb_tiles_stable_density.pdf", width = 12, height = 8)
h <- heatmap.2( meanMeth, scale="none", col = colors, margins = c(12, 12), trace='none', denscol="white", labCol=paste(info$Species, info$Tissue, sep=" "), labRow=NA, ColSideColors=pal[as.integer(as.factor(info$Species))+9], cexCol = 1.5, dendrogram='col', key=T)
## heatmap base on autocorrelation matrix
cors <- cor(meanMeth, method="spearman", use="pairwise.complete.obs") 
h <- heatmap.2(cors , scale="none", col = colors, margins = c(12, 12), trace='none', denscol="white", labCol=paste(info$Species, info$Tissue, sep=" "), labRow=info$Condition, ColSideColors=pal[as.integer(as.factor(info$Species))+9], RowSideColors=pal[as.integer(as.factor(info$Tissue))], cexCol = 1.5)
dev.off()
## Disappointing! At these sites, the methylation is so conserved / stable that there is very low tissue structure...
## Macaque samples still cluster out :(


## Absolute number of changes low (and number rhesus >= 5)
tiles <- density1kb[density1kb$numRhesus >= 5, ]
tiles <- tiles[(density1kb$numRhesus[density1kb$numRhesus >= 5] - density1kb$numHumanRefRhesus[density1kb$numRhesus >= 5]) < 4 & (density1kb$numRhesus[density1kb$numRhesus >= 5] - density1kb$numChimpRefRhesus[density1kb$numRhesus >= 5]) < 4, ]
## 143,321
## ... (not shown) ..
## -> not really meaningful...


## CpG density stayed above a critical threshold
tiles <- density1kb[density1kb$numChimp >= 20 & density1kb$numChimp >= 20 & density1kb$numRhesus >= 20, ]
tiles <- GRanges(tiles[, 1], IRanges(start=tiles[, 2], end=tiles[, 2]))
## 126,068 tiles

## Then overlap selected tiles with shared CpG sites
overlaps <- findOverlaps(allSites, tiles)
sites <- allSites[overlaps@queryHits]
## subset data
allData.fit.subset.sites <- subsetByOverlaps(allData.fit.subset, sites)
pdf(file = "distributionMethylationSmoothed_average1kbWindows_considerDensity2.pdf", height = 6, width=24)
par(mfrow=c(1,4), ps = 14, cex = 1, cex.main = 1) ## Used to keep labels at Font 14pt
for (tissue in unique(info$Tissue)){
  i <- 0
  for (species in unique(info$Species)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset.sites[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      if (i == 0){
        plot(density(aggregate(getMeth(allData.fit.subset.sites[,sample]), list(overlaps@subjectHits), mean)[,2]), col=pal[as.integer(as.factor(info[sample,]$Species))+9], main=tissue, xlab="Methylation percentage", lwd=2, lty=j, ylim=c(0,4))
        i <- 1
      }
      else {
        lines(density(aggregate(getMeth(allData.fit.subset.sites[,sample]), list(overlaps@subjectHits), mean)[,2]), col=pal[as.integer(as.factor(info[sample,]$Species))+9], lwd=2, lty=j)
      }
    }
  }
  if (tissue==unique(info$Tissue)[1]){
    legend("top", legend=c(as.character(unique(info$Species)), as.character(unique(info[info$Tissue == tissue,]$IndividualID))), lty=c(1,1,1,rep(c(1,2,3,4), times=3)), lwd=2, col=c(pal[10:12], rep("black",12)), bty = "n") 
  }
}
dev.off()
## high density tiles show reduced differences. At one point, it doesn't matter if CpG disappear, there are enough of them to maintain the methylation levels stable
## TO DO: repeat for different cutoffs of density (10, 20, 30, etc)


## Use direct density estimates numHuman, numChimp and numRhesus (more data). Large changes will be seen for regions with missing orthology, so these will be eliminated in practice...
## Also, if a CpG is replaced by another one within the tile, the density will not change

## Evolution of relative number of CpGs compared to human
## Only tiles with at least 5 CpGs in human are considered 
plot(density((density1kb$numHuman[density1kb$numHuman >= 5] - density1kb$numRhesus[density1kb$numHuman >= 5])/density1kb$numHuman[density1kb$numHuman >= 5]), col=pal[12], lwd=2, xlab="", main="", ylim=c(0, 3), xlim=c(-2, 1))
lines(density((density1kb$numHuman[density1kb$numHuman >= 5] - density1kb$numChimp[density1kb$numHuman >= 5])/density1kb$numHuman[density1kb$numHuman >= 5]), col=pal[11], lwd=2, xlab="", main="")
## 0: no change compared to rhesus
## 1: all human positions no found in chimp/macaque
## <0: more positions found in chimp/macaque
## The distribution is peaking at 0: good! A lot more signal than when we start from rhesus CpGs and took care of the orthology of CpGs...
densChimp <- (density1kb$numHuman[density1kb$numHuman >= 5] - density1kb$numChimp[density1kb$numHuman >= 5])/density1kb$numHuman[density1kb$numHuman >= 5]
densRhesus <- (density1kb$numHuman[density1kb$numHuman >= 5] - density1kb$numRhesus[density1kb$numHuman >= 5])/density1kb$numHuman[density1kb$numHuman >= 5]

## How many tiles left when filtering based on these criteria?
summary(densChimp > -0.1 & densChimp < 0.1 & densRhesus > -0.1 & densRhesus < 0.1) 
## 125,085 tiles

## Load tiles in a granges object
tiles <- density1kb[density1kb$numRhesus >= 5, ]
tiles <- tiles[densChimp > -0.1 & densChimp < 0.1 & densRhesus > -0.1 & densRhesus < 0.1, ]
tiles <- GRanges(tiles[, 1], IRanges(start=tiles[, 2], end=tiles[, 2]))
## Then overlap selected tiles with shared CpG sites
overlaps <- findOverlaps(allSites, tiles)
sites <- allSites[overlaps@queryHits]
## subset data
allData.fit.subset.sites <- subsetByOverlaps(allData.fit.subset, sites)

pdf(file = "distributionMethylationSmoothed_average1kbWindows_considerDensity3.pdf", height = 6, width=24)
par(mfrow=c(1,4), ps = 14, cex = 1, cex.main = 1) ## Used to keep labels at Font 14pt
for (tissue in unique(info$Tissue)){
  i <- 0
  for (species in unique(info$Species)){
    j <- 0
    for (sample in sampleNames(allData.fit.subset.sites[, info$Species == species & info$Tissue == tissue])){
      cat(sample, "\n")
      j <- j+1
      if (i == 0){
        plot(density(aggregate(getMeth(allData.fit.subset.sites[,sample]), list(overlaps@subjectHits), mean)[,2]), col=pal[as.integer(as.factor(info[sample,]$Species))+9], main=tissue, xlab="Methylation percentage", lwd=2, lty=j, ylim=c(0,4))
        i <- 1
      }
      else {
        lines(density(aggregate(getMeth(allData.fit.subset.sites[,sample]), list(overlaps@subjectHits), mean)[,2]), col=pal[as.integer(as.factor(info[sample,]$Species))+9], lwd=2, lty=j)
      }
    }
  }
  if (tissue==unique(info$Tissue)[1]){
    legend("top", legend=c(as.character(unique(info$Species)), as.character(unique(info[info$Tissue == tissue,]$IndividualID))), lty=c(1,1,1,rep(c(1,2,3,4), times=3)), lwd=2, col=c(pal[10:12], rep("black",12)), bty = "n") 
  }
}
dev.off()
## Less good if density is kept but not orthology (one lost CpG can be replaced by another one), at least with this window size...
## CCL: only a strict control, using regions were very few of the identified ortgologous sites evolved, allows to attentuate species differences... But this is a minority of regions, so across the genome this expalins major differences across species


######################################################################
## What is the proportion of UMR (<10%), FMR (>50%), LMR (or cytosines showing intermediate methylation levels, see Stadler et al. 2011)
## apply(getMeth(allData.fit.subset), 2, function(x){ cat(sum(x<=0.1)/length(x), "\t", sum(x>0.1 & x<0.5)/length(x), "\t", sum(x>=0.5)/length(x), "\n") })
load("./smooth_data/combined_samples/combinedSmoothedCommonSites.RDa")
info <- pData(allData.fit.subset)

props <- cbind(sampleNames(allData.fit.subset), t(rbind(apply(getMeth(allData.fit.subset), 2, function(x){ return(c(sum(x<=0.1)/length(x), sum(x>0.1 & x<0.5)/length(x), sum(x>=0.5)/length(x))) }))))
colnames(props) <- c("Sample", "Low methylation", "Intermediate methylation", "High methylation")
write.table(props, file = "distributionMethylationClasses.txt", sep="\t", quote=F, row.names=F, col.names=T)
# This info was added to Sup table

## mean methylation level of all samples
M <- getMeth(allData.fit.subset)
methylation <- colMeans(M)
## This info was added to Sup table

## boxplot of mean percent methylation per tissue, and per species
pdf(file = paste0("percent_methylation_tissue.pdf"), width = 6, height = 6)
boxplot(methylation ~ as.factor(info$Tissue), col=pal[1:4], notch=F, pch=16, boxwex=0.7)
dev.off()
pdf(file = paste0("percent_methylation_species.pdf"), width = 6, height = 6)
boxplot(methylation ~ as.factor(info$Species), col=pal[1:3], notch=F, pch=16, boxwex=0.7)
dev.off()
pdf(file = paste0("percent_methylation_species_tissues.pdf"), width = 10, height = 6)
boxplot(methylation ~ as.factor(paste(info$Species, info$Tissue)), col=pal[1:4], notch=F, pch=16, boxwex=0.7, ylab="Mean methylation", names=rep(c("heart", "kidney", "liver", "lung"), 3), at =c(1,2,3,4, 6,7,8,9, 11,12,13,14))
mtext(unique(info$Species)[1], side = 1, line = 3, at=2.5)
mtext(unique(info$Species)[2], side = 1, line = 3, at=7.5)
mtext(unique(info$Species)[3], side = 1, line = 3, at=12.5)
dev.off()
## Overall the same trend in each species
## This is consistent with density plots (see above)

######################################################################
## DMR finding with bsseq: Test case between human and chimp hearts ##
######################################################################

## ## make up the bsseq object (too big to be constructed with BSseq() function directly)
## keep <- grepl("^[HC]\\dH", colnames(Cov), perl=T)
## pData <- pData[keep,]
## M <- M[,keep]
## Cov <- Cov[,keep]
## coef <- coef[,keep]
## allData.fit <- BSseq(gr = gr[1:10,], M = Cov[1:10,], Cov = Cov[1:10,], coef = coef[1:10,], se.coef = NULL, pData = pData, trans = trans, parameters = parameters, rmZeroCov = FALSE)
## allData.fit@rowData <- gr
## allData.fit@assays$data$coef <- coef
## rm(coef)
## allData.fit@assays$data$M <- M
## rm(M)
## allData.fit@assays$data$Cov <- Cov
## ## subset object to common sites
## keepLoci <- which(rowSums(Cov) >= 1) ## 31,997,653 loci
## Cov <- Cov[keepLoci,]
## allData.fit <- allData.fit[keepLoci,]

## ## Only keep sites with at least 2 samples with coverage of 2 in each species
## keepLoci <- which(rowSums(Cov[, grepl("^H", colnames(Cov), perl=T)] >= 2) >= 2 & rowSums(Cov[, grepl("^C", colnames(Cov), perl=T)] >= 2) >= 2)
## ## remove loci with NA smoothed values
## keepLoci <- keepLoci[!is.na(apply(getMeth(allData.fit[keepLoci,], type = "smooth"), 1, sum))]
## keepLoci <- keepLoci[rowMeans(Cov[keepLoci,]) <= 10]
## length(keepLoci) ## 15,355,691
## allData.fit.subset <- allData.fit[keepLoci,]

## ## calculate coverage at common sites
## round(colMeans(getCoverage(allData.fit.subset)), 2)

## ## compute t-statistics
## allData.tstat <- BSmooth.tstat(allData.fit.subset,
##                                group1 = c("H1H", "H2H", "H3H", "H4H"),
##                                group2 = c("C1H", "C2H", "C3H", "C4H"),
##                                estimate.var = "same",
##                                local.correct = TRUE, ## large-scale (low-frequency) mean correction. This is especially important when large-scale methylation differences between 2 conditions (e.g., cancer and normals).
##                                verbose = TRUE)
## allData.tstat
## pdf(file = "DMRs/species/HumanChimp_heart_histTstats.pdf", width = 6, height = 6)
## plot(allData.tstat)
## dev.off()

## ## find DMRs
## dmrs0 <- dmrFinder(allData.tstat, cutoff = c(-4.6, 4.6))
## dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
## nrow(dmrs) ## 15337
## head(dmrs, n = 3)

## ## Plotting
## pData$col <- rep(c(pal[1], pal[2]), each = 4)
## pData(allData.fit) <- pData
## ## pdf(file = "DMRs/species/HumanChimp_heart_DMR1.pdf", width = 10, height = 5)
## ## plotRegion(allData.fit, dmrs[1,], extend = 10000, addRegions = dmrs)
## ## dev.off()
## pdf(file = "DMRs/species/HumanChimp_heart_top100DMRs.pdf", width = 10, height = 5)
## plotManyRegions(allData.fit, dmrs[1:100,], extend = 5000, addRegions = dmrs)
## dev.off()


## Quicker way: load the data (RData files) already generated and recreate the plots (for example)
## Human specific DMRs in heart:
species <- "Human"
tissue <- "heart"
load("DMRs/species/HumanSpecific_heart.RDa")
pData <- pData(allData.fit)
Cov <- getCoverage(allData.fit)
hist(log10(rowSums(Cov)), breaks=100)
abline(v=log10(120))
summary(rowMeans(Cov[keepLoci,]) <= 10)

keepLoci <- which(rowSums(Cov[, pData$Species == unique(pData$Species)[1]] >= 2) >= 2 & rowSums(Cov[, pData$Species == unique(pData$Species)[2]] >= 2) >=2 & rowSums(Cov[, pData$Species == unique(pData$Species)[3]] >= 2) >= 2)
keepLoci <- keepLoci[!is.na(apply(getMeth(allData.fit[keepLoci,], type = "smooth"), 1, sum))]
keepLoci <- keepLoci[rowMeans(Cov[keepLoci,]) <= 10]
allData.fit.subset <- allData.fit[keepLoci,]

load("DMRs/species/HumanSpecific_heart_tstat.RDa")
allData.tstat.others <- BSmooth.tstat(allData.fit.subset,
                                      group1 = rownames(pData[pData$Species == unique(pData[pData$Species != species,]$Species)[1],]),
                                      group2 = rownames(pData[pData$Species == unique(pData[pData$Species != species,]$Species)[2],]), 
                                      estimate.var = "same",
                                      local.correct = TRUE,
                                      verbose = TRUE)
allData.tstat@stats[abs(allData.tstat.others@stats[,"tstat.corrected"]) >= 4.6 ,"tstat.corrected"] <- NA
dmrs0 <- dmrFinder(allData.tstat, cutoff = c(-4.6, 4.6), maxGap = 300, stat = "tstat.corrected", verbose=T)
dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
summary(dmrs$direction == "hyper")

## Plotting
keepLoci <- which(rowSums(getCoverage(allData.fit)) >= 1)
pData$col <- rep(c(pal[1], pal[2], pal[3]), each = 4)
pData(allData.fit) <- pData
## pdf(file = paste0("DMRs/species/", species,"Specific_", tissue, "_top100DMRs.pdf"), width = 10, height = 5)
plotManyRegions(allData.fit[keepLoci, ], dmrs[1:100,], extend = 5000, addRegions = dmrs)
dev.off()

########################
## Recreate DMR plots ##
########################
## For example if we want to change some plotting options

pair <- c("Human", "Chimp")
tissue <- "liver"
## load the Tstat object already generated
load("DMRs/species/HumanChimp_liver_tstat.RDa")
dmrs0 <- dmrFinder(allData.tstat, cutoff = c(-4.6, 4.6), maxGap = 300, stat = "tstat.corrected", verbose=T)
dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)

## The BSseq object has to be reloaded
load("smooth_data/combined_samples/pData.RDa")
load("DMRs/species/pairwiseSpecies_liver.RDa")
keep <- pData$Species == pair[1] | pData$Species == pair[2] 

# Plotting: 2 species tested included here only
pData$col <- rep(c(pal[11], pal[10]), each = 4)
pData(allData.fit) <- pData
keepLoci <- which(rowSums(getCoverage(allData.fit)) >= 1)
allData.fit.subset <- allData.fit[keepLoci, keep]

pdf(file = paste0("DMRs/species/", pair[1], pair[2],"_", tissue, "_top", n, "DMRs_replotted.pdf"), width = 10, height = 5)
plotManyRegions(allData.fit.subset, dmrs[1:n,], extend = 5000, addRegions = dmrs, addPoints=TRUE, pointsMinCov=4)
dev.off()


#######################################
## tissue-specificity of speciesDMRs ##
#######################################
## see find_tissue_specificity_of_speciesDMRs.pl
spec <- read.table("tissue_specificity_of_speciesDMRs.txt", h=T, sep="\t")

## proportion of unique DMRs
summary(spec[,4]/spec[,3])
cbind(spec[,1:2], spec[,4]/spec[,3])

## proportion of shared DMRs
summary(spec[,5]/spec[,3])
cbind(spec[,1:2], spec[,5]/spec[,3])

################################
## conservation of tissueDMRs ##
################################
## First and easier criterion: a tDMR is conserved if a tDMR is found in the orthologous region in another species. In out case all species are mapped to the human coordiantes so it's very easy: just intercept the coordinates. See find_conservation_of_tissueDMRs.pl and find_conservation_of_tissueDMRs_allDMRs.pl 
cons <- read.table("conservation_of_tissueDMRs.txt", h=T, sep="\t")

## Problem: sometimes, in the orthologous regions, the CpG sites are not conserved. We need a reference list with the possible orthology status for each DMR for all DMR types.
## For each DMR type, take all orthologus regions. Load t-stats file in other species (only sites with decent coverage are there). If >=3 sites distant less than 300bp there, this region is considered valid for orthology

## How to best (and quickly) do the overlap between DMRs and t-stat objects? Granges function subsetByOverlaps is slow and it needs to be done DMR by DMR. Instead let's use findOverlaps which gives the many to many correspondance between tstats and all DMRs
library("GenomicRanges")
 
tissuePairs <- list(c("heart", "liver"), c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "Specific"), c("liver", "Specific"), c("kidney", "Specific"), c("lung", "Specific"))
allSpecies <- c("Human", "Rhesus", "Chimp")
for (species in allSpecies){
  for (pair in tissuePairs){
    if (pair[2] == "Specific"){
      baseName <- paste0("./DMRs/tissues/", species, "_", pair[1], pair[2])
    } else {
      baseName <- paste0("./DMRs/tissues/", species, "_", pair[1], "_", pair[2])
    }
    ## Reading DMRs objects
    dmrs <- read.table(paste0(baseName, "_DMRs.txt"), sep="\t", h=T)

    ## Preparing final orthology data frame to export
    allOrth <- data.frame(cbind(rep(NA, length(dmrs[,1])), rep(NA, length(dmrs[,1]))))
    names(allOrth) <- allSpecies[-grep(species, allSpecies)]
    
    ## Comparing the DMRs in species under consideration to tstats of the 2 other species
    for (species2 in allSpecies[-grep(species, allSpecies)]){
      cat("Testing orthology for", length(dmrs[,1]), species, pair, "DMRs in", species2, "\n")
      ## Reading tstat objects
      if (pair[2] == "Specific"){
        load(paste0("./DMRs/tissues/", species2, "_", pair[1], pair[2],"_tstat.RDa"))
      } else {
        load(paste0("./DMRs/tissues/", species2, "_", pair[1], "_", pair[2],"_tstat.RDa"))
      }
      tstat <- granges(allData.tstat)
      cat("\t Tstat object loaded\n")
      ## Converting DMRs to Granges
      dmrsGr <- GRanges(dmrs$chr, IRanges(dmrs$start, dmrs$end))

      ## Overlap tstat object and all DMRs. Resulting object has query hits (row in tstat) and subject hits (row in DMRs) correspondance. 
      ##   queryHits subjectHits 
      ## 1       386       18834 
      ## 2       387       18834 
      ## 3       530        5056 
      overlap <- findOverlaps(tstat, dmrsGr)
      overlap <- as.data.frame(overlap)
      cat("\t Overlap object created\n")   

      ## Replacing query hits with coordinates of CpG sites
      overlap$queryHits <- start(tstat[overlap$queryHits])
      ## create vector for orthology status
      orth <- rep(0, times=length(dmrs[,1]))
      names(orth) <- 1:length(dmrs[,1])
      ## Create a list with keys = all DMRs (subject hit), and values = orthology or not
      orth2 <- tapply(overlap$queryHits, list(overlap$subjectHits),
                     function(x){
                       ## We should have at least 3 CpG, distant less than 300 bp
                       ## Logical vector indicating if intervals lengths between consective elements is less than 300bp:
                       allDiffs <- diff(x) - 1 <= 300;
                       ## if there is a run of at least 2 TRUE elements, return 1, ortherwise return 0
                       if (length(allDiffs) < 2){ return(0) } else {
                         return(ifelse(any(rle(allDiffs)$lengths[rle(allDiffs)$values] >= 2), 1, 0))
                       }
                     })
      ## update the orthology status vector (which, contrary to orth2, includes DMRs with no orthologous sites, set to 0)
      orth[names(orth2[orth2 == 1])] <- 1
      cat("\t", sum(orth == 1), "out of", length(orth), "tDMRs have orthology to", species2, "\n")
      allOrth[[species2]] <- orth
    }
    ## Export orthology results into file
    write.table(cbind(dmrs[,1:3], allOrth), file=paste0(baseName,"_DMRs_orthology.txt"), sep="\t", quote=F, row.names=F, col.names=T)
  }
}
## tapply optimized the calculation as lot (from >1 h per DMR type to a few seconds)
## Check with Human heart-liver:
##   - 18205 out of 22932 with orthology to rhesus (79.4%)
##   - 20605 to chimp (89.9%)
##   - 1231 DMRs with no orthology to chimp or rhesus
##   - 17109 DMRs with orthology human/chimp/rhesus
##   - 3496 DMRs with orthology human/chimp 
##   - 1096 DMRs with orthology human/rhesus 

###################
## Now we can calculate the proportion of species-specific DMRs, relative to the list of all tDMRs with possible orthology
## Loop across all contrasts to reproduce table similar to conservation_of_tissueDMRs.txt
make_summary_table_cons <- function(filterFileExtension=NA, filterColumn=NA){
  counts <- as.data.frame(matrix(ncol=13, nrow=30)) ## 10 contrasts  * 3 species
  names(counts) <- c("Contrast", "Species", "Number of DMRs", "Number of DMRs with no orthology", "Number of DMRs with human-chimp orthology", "Number of DMRs with human-rhesus orthology", "Number of DMRs with chimp-rhesus orthology", "Number of DMRs with human-chimp-rhesus orthology", "Number of species-specific DMRs", "Number of human-chimp DMRs", "Number of human-rhesus DMRs", "Number of chimp-rhesus DMRs", "Number of human-chimp-rhesus DMRs")
  tissuePairs <- list(c("heart", "liver"), c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "Specific"), c("liver", "Specific"), c("kidney", "Specific"), c("lung", "Specific"))
  allSpecies <- c("Human", "Chimp", "Rhesus")
  i <- 1
  for (species in allSpecies){
    for (pair in tissuePairs){
      cat(species, ":", pair, "\n")
      if (pair[2] == "Specific"){
        baseName <- paste0("./DMRs/tissues/", species, "_", pair[1], pair[2])
        counts[i, 1] <- paste0(pair[1], pair[2])
      } else {
        baseName <- paste0("./DMRs/tissues/", species, "_", pair[1], "_", pair[2])
        counts[i, 1] <- paste0(pair[1], "_", pair[2])
      }
      counts[i, 2] <- species

      ## Conservation
      cons <- read.table(paste0(baseName,"_DMRs_conservation.txt"), h=T)
      ## Sort DMRs
      cons <- cons[order(cons$chr, cons$start, cons$end),]
      ## Orthology
      orth <- read.table(paste0(baseName,"_DMRs_orthology.txt"), h=T)
      ## Sort DMRs
      orth <- orth[order(orth$chr, orth$start, orth$end),]

      ## vector with the 2 other species
      allSpecies2 <- allSpecies[-grep(species, allSpecies)]
      
      if(!is.na(filterFileExtension) & !is.na(filterColumn[1])){
        ## Filter list of DMRs based on some feature annotation
        ## Several filter columns allowed (e.g., promoter & CGI)
        features <- read.table(paste0("../annotation/", baseName, filterFileExtension), h=T, sep="\t", check.names = F)
        features <- features[order(features$chr, features$start, features$end),]
        ## Problem: in some files, we have 0/1 encoding, in some others we have NA/name of feature: replace 0s by NAs
        filter <- rep(TRUE, times=length(orth[,1]))
        for (j in 1:length(filterColumn)){
          if (!is.null(features[[filterColumn[j]]])){
            features[[filterColumn[j]]][features[[filterColumn[j]]] == 0] <- NA
            filter[is.na(features[[filterColumn[j]]])] <- FALSE
          } else {
            cat(paste0("  Warning: there seem to be no ", filterColumn[j], " column in the file used for filtering\n"))
          }
        }
        ## Keeps only overlapping DMRs
        cons <- cons[filter,]
        orth <- orth[filter,]
      }

      ## Total number of DMRs
      counts[i, 3] <- length(cons[, 1])
      ## Number of DMRs with no orthology to any other species
      counts[i, 4] <- sum(apply(orth[,4:5], 1, sum) == 0)
      ## Number of DMRs with human-chimp orthology
      if (species == "Human" | species == "Chimp"){
        counts[i, 5] <- sum(orth[[allSpecies2[1]]] == 1)
      } else if (species == "Rhesus"){
        counts[i, 5] <- NA
      }
      ## Number of DMRs with human-rhesus orthology
      if (species == "Human"){
        counts[i, 6] <- sum(orth[[allSpecies2[2]]] == 1)
      } else if (species == "Chimp"){
        counts[i, 6] <- NA
      } else if (species == "Rhesus"){
        counts[i, 6] <- sum(orth[[allSpecies2[1]]] == 1)
      }
      ## Number of DMRs with chimp-rhesus orthology
      if (species == "Rhesus" | species == "Chimp"){
        counts[i, 7] <- sum(orth[[allSpecies2[2]]] == 1)
      } else if (species == "Human"){
        counts[i, 7] <- NA
      }
      ## Number of DMRs with human-chimp-rhesus orthology
      counts[i, 8] <- sum(apply(orth[,4:5], 1, sum) == 2)

      ## Conservation encoding:
      ## for human:  unique (0), human-chimp (1),  human-rhesus (2) and human-chimp-rhesus (4)
      ## for chimp:  unique (0), human-chimp (1),  chimp-rhesus (3) and human-chimp-rhesus (4)
      ## for rhesus: unique (0), human-rhesus (2), chimp-rhesus (3) and human-chimp-rhesus (4)
      ## This numbers are filtered for orthology status. Simplest is to consider only DMRs with HCR orthology status. Otherwise we get into complex considerations sicne the reference is changing. Additionally conservation can be found even without strict ortholgoy criteria met
      
      ## Number of species-specific DMRs
      counts[i, 9] <- sum(cons$conservation == 0 & apply(orth[,4:5], 1, sum) == 2)
      ## Number of human-chimp DMRs
      if (species == "Human" | species == "Chimp"){
        counts[i, 10] <- sum(cons$conservation == 1 & apply(orth[,4:5], 1, sum) == 2)
      } else if (species == "Rhesus"){
        counts[i, 10] <- NA
      } 
      ## Number of human-rhesus DMRs
      if (species == "Human" | species == "Rhesus"){
        counts[i, 11] <- sum(cons$conservation == 2 & apply(orth[,4:5], 1, sum) == 2)
      } else if (species == "Chimp"){
        counts[i, 11] <- NA
      }
      ## Number of chimp-rhesus DMRs
      if (species == "Human"){
        counts[i, 12] <- NA
      } else if ((species == "Chimp") | (species == "Rhesus")){
        counts[i, 12] <- sum(cons$conservation == 3 & apply(orth[,4:5], 1, sum) == 2)
      }
      ## Number of human-chimp-rhesus DMRs 
      counts[i, 13] <- sum(cons$conservation == 4 & apply(orth[,4:5], 1, sum) == 2)

      ## increment row
      i <- i + 1
    }
  }
  return(counts)
}
# Run it on full set of DMRs with HCR orthology
counts <- make_summary_table_cons(NA, NA)
write.table(counts, file="conservation_of_tissueDMRs_with_HCR_orthology.txt", sep="\t", quote=F, row.names=F, col.names=T)

## Plotting: make plots automatically with all contrasts and 3 species in same plot
make_conservation_plots <- function(counts, suffix=NA){
  ## Transform to relative counts
  counts[, 9:13] <- counts[, 9:13] / counts[,8]

  ## Open graphics
  if (is.na(suffix)){
    pdf(file = "conservation_of_tissueDMRs.pdf", width = 18, height = 5)
  } else {
    pdf(file = paste0("conservation_of_tissueDMRs_", suffix, ".pdf"), width = 18, height = 5)
  }
  par(mfrow=c(1,3),
      mar=c(6, 5, 0.2, 2), ## let some room for x labels
      ps = 14, cex = 1, cex.main = 1) ## keep labels at 14pt size
  allSpecies <- c("Chimp", "Human", "Rhesus")
  for (species in allSpecies){
    ## select species
    counts_to_plot <- counts[counts$Species == species, ]
    row.names(counts_to_plot) <- counts_to_plot$Contrast
    number_DMRs_with_orthology <- counts[counts$Species == species, 8]

    ## select columns of proportions and order then from conserved to species-specific
    counts_to_plot <- as.data.frame(t(counts_to_plot)[13:9,])
    ## Remove NA rows
    counts_to_plot <- counts_to_plot[!is.na(counts_to_plot[,1]),]
    ## Sort pairwise contrasts / tissue-specific contrasts from less species-specific to more species-specific
    number_DMRs_with_orthology <- c(number_DMRs_with_orthology[1:6][order(counts_to_plot[4,1:6])], number_DMRs_with_orthology[7:10][order(counts_to_plot[4,7:10])]) ## first sort number of DMRs
    counts_to_plot <- cbind(counts_to_plot[,1:6][, order(counts_to_plot[4,1:6])], counts_to_plot[,7:10][, order(counts_to_plot[4,7:10])])
    ## Plot
    if (species == "Human" | species == "Chimp"){
      mp <- barplot(as.matrix(counts_to_plot), col=seq[c(2,2,4,6)], density=c(NA, 55, NA, NA), angle=c(NA, 45, NA, NA), ylab="", xlab="", xaxt="n", space=c(rep(0.2,times=6), 0.5, rep(0.2, times=3)), ylim = c(0, 1.4), axes=F)
    } else if (species == "Rhesus"){
      mp <- barplot(as.matrix(counts_to_plot), col=seq[c(2,2,2,6)], density=c(NA, 55, 55, NA), angle=c(NA, 45, 135, NA), ylab="", xlab="", xaxt="n", space=c(rep(0.2,times=6), 0.5, rep(0.2, times=3)), ylim = c(0, 1.4), axes=F)
    }
    axis(2, at=seq(0, 1, 0.2))
    mtext(side = 2, paste0("Proportion of ", tolower(species), " tDMRs"), line = 3, at=0.5)
    ## Format x labels
    labels <- colnames(counts_to_plot)
    simpleCap <- function(x) {
      s <- strsplit(x, " ")[[1]]
      paste(toupper(substring(s, 1,1)), substring(s, 2),
            sep="", collapse=" ")
    }
    labels <- unlist(lapply(labels, simpleCap))
    labels <- gsub("_", " vs. ", labels, perl=TRUE)
    labels <- gsub("Specific", "-specific", labels, perl=TRUE)
    text(mp, -0.05, srt = 45, adj = 1, labels = labels, xpd = TRUE, cex=1)
    ## Add legend on top
    if (species == "Human"){
      legend("topleft", rev(c("Conserved human to macaque (+chimp)", "Conserved human to macaque (-chimp)", "Conserved human to chimp", "Human-specific")), fill=rev(seq[c(2,2,4,6)]), density=rev(c(NA, 55, NA, NA)), angle=rev(c(NA, 45, NA, NA)), bty="n")
    } else if (species == "Chimp"){
      legend("topleft", rev(c("Conserved chimp to macaque (+human)", "Conserved chimp to macaque (-human)", "Conserved chimp to human", "Chimp-specific")), fill=rev(seq[c(2,2,4,6)]), density=rev(c(NA, 55, NA, NA)), angle=rev(c(NA, 45, NA, NA)), bty="n")
    } else if (species == "Rhesus"){
      legend("topleft", rev(c("Conserved macaque to human/chimp", "Conserved macaque to chimp (-human)", "Conserved macaque to human (-chimp)", "Macaque-specific")), fill=rev(seq[c(2,2,2,6)]), density=rev(c(NA, 55, 55, NA)), angle=rev(c(NA, 45, 135, NA)), bty="n")
    }

    ## Add number of DMRs considered for each contrast
    text(mp, 0.98, srt = 90, adj = 1, labels = number_DMRs_with_orthology, xpd = TRUE, cex=1)
  }
  dev.off()
}

## Plotting:
library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))
seq <- rev(brewer.pal(9, "Blues")) ## sequential colors
make_conservation_plots(counts)

## Conservation for subset of DMRs overlapping specific features
## Promoter:
counts <- make_summary_table_cons("_DMRs_15_features_0.2.gz", "promoter")
make_conservation_plots(counts, "promoter")
## Slightly higher than over all tDMRs

## Proximal promoter:
counts <- make_summary_table_cons("_DMRs_15_features_0.2.gz", "promoter_proximal")
make_conservation_plots(counts, "promoter_proximal")

## Conserved promoters
counts <- make_summary_table_cons("_DMRs_conserved_features.gz", "promoter_conserved_panTro3+rheMac2")
make_conservation_plots(counts, "promoter_conserved")

## Conserved CGI
counts <- make_summary_table_cons("_DMRs_conserved_features.gz", "CGI_conserved_panTro3+rheMac2")
make_conservation_plots(counts, "CGI_conserved")

## Conserved repeats
counts <- make_summary_table_cons("_DMRs_conserved_features.gz", "repeats_conserved_panTro3+rheMac2")
make_conservation_plots(counts, "repeats_conserved")

## TO DO: tissueDMRs next to conserved TSSs of 1-to-1 orthologous genes that show DE
## TO DO: tissueDMRs next to conserved TSSs of 1-to-1 orthologous genes that show conserved DE

##########################################################################
## ## Older code:
## ## absolute counts
## ## focus on one particular contrast: heart vs. liver
## counts <- cons[cons$Species == "Human" & cons$Contrast == "heart_liver", -c(1,2)]
## seq <- rev(brewer.pal(9, "Blues")) ## sequential colors
## pdf(file = "conservation_of_tissueDMRs_human_heart_liver.pdf", width = 6, height = 5)
## barplot(rep(counts$Number.of.DMRs, 3), col=c(seq[6], NA, NA), names.arg=c("Human", "Human to chimp", "Human to macaque"))
## barplot(c(0, counts$Number.of.human.chimp.DMRs+ counts$Number.of.human.chimp.rhesus.DMRs, 0), col=c(NA, seq[4], NA), add=T, axes=F)
## barplot(c(0, 0, counts$Number.of.human.chimp.rhesus.DMRs+counts$Number.of.human.rhesus.DMRs), col=c(NA, NA, seq[2]), add=T, axes=F)
## dev.off()
## ## Used in PPT presentation (barplots below good for paper, but too intensive for oral presentation)

## ## proportions / Human
## counts <- t(as.matrix(cons[cons$Species == "Human", c(8,6,5,4)]/cons[cons$Species == "Human", 3]))
## colnames(counts) <- cons[cons$Species == "Human", 1]
## counts <- counts[, c(5:10, 1:4)]
## ## colSums(counts) ## check, should be 1
## ## Sort contrasts (% species-specific)
## counts <- cbind(counts[,1:6][, order(counts[4,1:6])], counts[,7:10][, order(counts[4,7:10])])
  
## ## plotting
## pdf(file = "conservation_of_tissueDMRs_human.pdf", width = 6, height = 5)
## seq <- rev(brewer.pal(9, "Blues")) ## sequential colors
## par(mar=c(6, 5, 0.2, 2)) ## let some room for x labels
## mp <- barplot(counts, col=seq[c(2,2,4,6)], density=c(NA, 55, NA, NA), angle=c(NA, 45, NA, NA), ylab="", xlab="", xaxt="n", space=c(rep(0.2,times=6), 0.5, rep(0.2, times=3)), ylim = c(0, 1.4), axes=F)
## axis(2, at=seq(0, 1, 0.2))
## mtext(side = 2, "Proportion of human DMRs", line = 3, at=0.5)
## labels <- colnames(counts)
## simpleCap <- function(x) {
##   s <- strsplit(x, " ")[[1]]
##   paste(toupper(substring(s, 1,1)), substring(s, 2),
##       sep="", collapse=" ")
## }
## labels <- unlist(lapply(labels, simpleCap))
## labels <- gsub("_", " vs. ", labels, perl=TRUE)
## labels <- gsub("Specific", "-specific", labels, perl=TRUE)
## text(mp, -0.05, srt = 45, adj = 1, labels = labels, xpd = TRUE, cex=1)
## legend("topleft", rev(c("Conserved human to macaque (+chimp)", "Conserved human to macaque (-chimp)", "Conserved human to chimp", "Human-specific")), fill=rev(seq[c(2,2,4,6)]), density=rev(c(NA, 55, NA, NA)), angle=rev(c(NA, 45, NA, NA)), bty="n")
## dev.off()

## ## Less stringent conservation criterion ###################  
## ## Within orthologous region in chimp or macaque, load t-stat and look for at least 1 CpG is diff methylated
## ## First try on 1 contrast: human heart-liver
## dmrs <- read.table("./DMRs/tissues/Human_heart_liver_DMRs.txt", sep="\t", h=T, stringsAsFactors=F, as.is=T)
## ## library("GenomicRanges")
## ## dmrsGr <- GRanges(dmrs$chr, IRanges(dmrs$start, dmrs$end))

## ## Load chimp heart-liver tstats
## load("DMRs/tissues/Chimp_heart_liver_tstat.RDa")
## allData.tstat ## 20,944,276 CpG sites translated in human coordinates. Not all of them are corresponding to human CpGs!
## names.tstat <- data.frame(chr=seqnames(granges(allData.tstat)), stringsAsFactors=FALSE)
## names.tstat$start <- start(granges(allData.tstat))
## names.tstat$position <- paste(names.tstat$chr, names.tstat$start, sep=":")
## dim(names.tstat) 
## names.tstat$chr <- as.character(names.tstat$chr)

## ## Foreach DMR, retrieve tstat of all CpGs in chimp object
## ## consChimp <- rep(0, length(dmrs[,1]))
## ## for (i in 1:length(dmrs[,1])){
## ##   print(i)
## ##   allData.tstat.subset <- subsetByOverlaps(allData.tstat, dmrsGr[i,])
## ##   if (sum(abs(allData.tstat.subset@stats[,"tstat.corrected"]) >= 4.6) > 0){
## ##     consChimp[i] <- 1
## ##   }
## ## }
## ## Very long!!!

## ## Faster + takes care of direction, but still takes ~12h!:
## allStats <- allData.tstat@stats[, "tstat.corrected"]
## names(allStats) <-  names.tstat$start ## helpful to know CpG site position. problem if several DMRs with same start position  
## consChimp <- rep(0, length(dmrs[,1]))
## ## Relaxed t-stat cutoff: (-3;3). Given distribution of t-stats, this seems reasonable (similar to 2-cutoff approaches for DE)
## ## open ~/clusterhome/Methylation/bsseq/DMRs/tissues/Human_heart_liver_histTstats.pdf 
## consChimpRelaxed <- rep(0, length(dmrs[,1]))
## ## Look around DMR: are there differentially methylated CpGs in same direction around (1kb upstream and downstream)?
## consChimpExtended <- rep(0, length(dmrs[,1]))
## consChimpExtendedRelaxed <- rep(0, length(dmrs[,1]))
## for (i in 1:length(dmrs[,1])){
##   print(i)
##   ## subset to right chromosome
##   allStats.subset <- allStats[names.tstat$chr == dmrs$chr[i] & names.tstat$start >= dmrs$start[i] & names.tstat$start <= dmrs$end[i]]
##   allStats.subset.extended <- allStats[names.tstat$chr == dmrs$chr[i] & names.tstat$start >= dmrs$start[i]-1000 & names.tstat$start <= dmrs$end[i]+1000]  
##   if (dmrs$direction[i] == "hyper" && (sum(allStats.subset >= 4.6) > 0)){
##     consChimp[i] <- 1
##   } else if (dmrs$direction[i] == "hypo" && (sum(allStats.subset <= -4.6) > 0)){
##     consChimp[i] <- 1
##   } 
##   if (dmrs$direction[i] == "hyper" && (sum(allStats.subset >= 3) > 0)){
##     consChimpRelaxed[i] <- 1
##   } else if (dmrs$direction[i] == "hypo" && (sum(allStats.subset <= -3) > 0)){
##     consChimpRelaxed[i] <- 1
##   } 
##   if (dmrs$direction[i] == "hyper" && (sum(allStats.subset.extended >= 4.6) > 0)){
##     consChimpExtended[i] <- 1
##   } else if (dmrs$direction[i] == "hypo" && (sum(allStats.subset.extended <= -4.6) > 0)){
##     consChimpExtended[i] <- 1
##   } 
##   if (dmrs$direction[i] == "hyper" && (sum(allStats.subset.extended >= 3) > 0)){
##     consChimpExtendedRelaxed[i] <- 1
##   } else if (dmrs$direction[i] == "hypo" && (sum(allStats.subset.extended <= -3) > 0)){
##     consChimpExtendedRelaxed[i] <- 1
##   } 
## }
## summary(consChimp == 1)
## summary(consChimpExtended == 1)

## ## We expect something higher than the 48% (what we find with just overlap of DMRs):
## ## We get 11587/22932 = 50.5%
## ## For extended: we get something bigger, as expected, but not so different: tends to indicate that there is no sliding/replacement of tDMRs by nearby DMRs 
## ## 12538/22932 = 54.7%

## ## [Similar chunk of code removed for rhesus]
## ## We expect something higher than the 37% (what we find with just overlap of DMRs)
## ## Only a bit more: 8771/22932 = 38.2%
## ## Extended: we get something bigger, but not a lot: a prioiri there is no systematic sliding/replacement of tDMRs to nearby DMRs 
## ## 10166/22932 = 44.3%
## ## Export conservation results into file
## write.table(cbind(dmrs[,1:3],consChimp, consChimpExtended, consRhesus, consRhesusExtended), file="conservation_of_tissueDMRs_human_heart_liver_relaxed.txt", sep="\t", quote=F, row.names=F, col.names=T)

## Less stringent conservation criterion for all contrasts ######################  
## Takes ~15min per contrast to compare with 1 other species (Granges function subsetByOverlaps is slow and it needs to be done DMR by DMR. Instead let's use findOverlaps which gives the many to many correspondance between tstats and all DMRs)
library("GenomicRanges")
 
tissuePairs <- list(c("heart", "liver"), c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "Specific"), c("liver", "Specific"), c("kidney", "Specific"), c("lung", "Specific"))
allSpecies <- c("Human", "Rhesus", "Chimp")
for (species in allSpecies){
  for (pair in tissuePairs){
    if (pair[2] == "Specific"){
      baseName <- paste0("./DMRs/tissues/", species, "_", pair[1], pair[2])
    }
    else {
      baseName <- paste0("./DMRs/tissues/", species, "_", pair[1], "_", pair[2])
    }
    ## Reading DMRs objects
    dmrs <- read.table(paste0(baseName, "_DMRs.txt"), sep="\t", h=T)
    
    ## Comparing the DMRs in species under consideration to tstats of the 2 other species
    for (species2 in allSpecies[-grep(species, allSpecies)]){
      cat("Testing conservation of", species, pair, "DMRs in", species2, "\n")
      ## Reading tstat objects
      if (pair[2] == "Specific"){
        load(paste0("./DMRs/tissues/", species2, "_", pair[1], pair[2],"_tstat.RDa"))
      }
      else {
        load(paste0("./DMRs/tissues/", species2, "_", pair[1], "_", pair[2],"_tstat.RDa"))
      }
      cat("\t Tstat object loaded\n")

      ## Converting DMRs to Granges
      dmrsGr <- GRanges(dmrs$chr, IRanges(dmrs$start, dmrs$end))
      dmrsGrExtended <- GRanges(dmrs$chr, IRanges(dmrs$start-1000, dmrs$end+1000))
      
      ## add t-stat as metadata to granges(allData.tstat)
      tstat <- granges(allData.tstat)
      tstat$tstat <- allData.tstat@stats[, "tstat.corrected"]

      ## Overlap tstat object and all DMRs. Resulting object has query hits (row in tstat) and subject hits (row in DMRs) correspondance. 
      ##   queryHits subjectHits 
      ## 1       386       18834 
      ## 2       387       18834 
      ## 3       530        5056 
      overlap <- findOverlaps(tstat, dmrsGr)
      ## For extended DMRs, they can overlap, bu this is not a problem since many to many hits are kept (option select="all" by default)
      overlapExtended <- findOverlaps(tstat, dmrsGrExtended)
      cat("\t Overlap objects created\n")
      
      cons <- rep(0, length(dmrs[,1]))
      ## Relaxed t-stat cutoff: (-3;3). Given distribution of t-stats, this seems reasonable (similar to 2-cutoff approaches for DE)
      ## open ~/clusterhome/Methylation/bsseq/DMRs/tissues/Human_heart_liver_histTstats.pdf 
      consRelaxed <- rep(0, length(dmrs[,1]))
      ## Look around DMR: are there differentially methylated CpGs in same direction around (1kb upstream and downstream)?
      consExtended <- rep(0, length(dmrs[,1]))
      consExtendedRelaxed <- rep(0, length(dmrs[,1]))
      for (i in 1:length(dmrs[,1])){
        if (i %% 1000 == 0) { cat("\t", i, "/", length(dmrs[,1]), "\n") }

        ## for each DMR (subject hit), get tstats corresponding to all query hits (overlapped positions) 
        allTstats <- tstat$tstat[overlap@queryHits[overlap@subjectHits == i]]
        if (dmrs$direction[i] == "hyper" && (sum(allTstats >= 4.6) > 0)){
          cons[i] <- 1
        } else if (dmrs$direction[i] == "hypo" && (sum(allTstats <= -4.6) > 0)){
          cons[i] <- 1
        } 
        if (dmrs$direction[i] == "hyper" && (sum(allTstats >= 3) > 0)){
          consRelaxed[i] <- 1
        } else if (dmrs$direction[i] == "hypo" && (sum(allTstats <= -3) > 0)){
          consRelaxed[i] <- 1
        } 

        ## Same thing for extended regions
        allTstatsExtended <- tstat$tstat[overlapExtended@queryHits[overlapExtended@subjectHits == i]]
        if (dmrs$direction[i] == "hyper" && (sum(allTstatsExtended >= 4.6) > 0)){
          consExtended[i] <- 1
        } else if (dmrs$direction[i] == "hypo" && (sum(allTstatsExtended <= -4.6) > 0)){
          consExtended[i] <- 1
        } 
        if (dmrs$direction[i] == "hyper" && (sum(allTstatsExtended >= 3) > 0)){
          consExtendedRelaxed[i] <- 1
        } else if (dmrs$direction[i] == "hypo" && (sum(allTstatsExtended <= -3) > 0)){
          consExtendedRelaxed[i] <- 1
        }   
      }
      ##summary(cons == 1)
      ##summary(consRelaxed == 1)
      ##summary(consExtended == 1)
      ##summary(consExtendedRelaxed == 1)

      ## Export conservation results into file
      write.table(cbind(dmrs[,1:3],cons, consRelaxed, consExtended, consExtendedRelaxed), file=paste0(baseName,"_DMRs_conservation_to_", species2, "_relaxed+extended.txt"), sep="\t", quote=F, row.names=F, col.names=T)
    }
  }
}
## TO DO: could probably be optimized even further by using tapply instead of looping though DMRs (see above for non relaxed conservation)

## Now we need to recalculate the orthology criterion: only 1 sites in orthologous region in t-stat object
## - get it for orthologous region
## - get it for extended region: +/- 1 kb from region
for (species in allSpecies){
  for (pair in tissuePairs){
    if (pair[2] == "Specific"){
      baseName <- paste0("./DMRs/tissues/", species, "_", pair[1], pair[2])
    } else {
      baseName <- paste0("./DMRs/tissues/", species, "_", pair[1], "_", pair[2])
    }
    ## Reading DMRs objects
    dmrs <- read.table(paste0(baseName, "_DMRs.txt"), sep="\t", h=T)

    allSpecies2 <- allSpecies[-grep(species, allSpecies)]
    ## Preparing final orthology data frame to export
    allOrth <- data.frame(a=rep(NA, times=length(dmrs[,1])), b=rep(NA, times=length(dmrs[,1])), c=rep(NA, times=length(dmrs[,1])), d=rep(NA, times=length(dmrs[,1])))
    names(allOrth) <- c(paste0(allSpecies2, "Relaxed"), paste0(allSpecies2, "ExtendedRelaxed"))
   
    ## Comparing the DMRs in species under consideration to tstats of the 2 other species
    for (species2 in allSpecies2){
      cat("Testing orthology for", length(dmrs[,1]), species, pair, "DMRs in", species2, "\n")
      ## Reading tstat objects
      if (pair[2] == "Specific"){
        load(paste0("./DMRs/tissues/", species2, "_", pair[1], pair[2],"_tstat.RDa"))
      } else {
        load(paste0("./DMRs/tissues/", species2, "_", pair[1], "_", pair[2],"_tstat.RDa"))
      }
      tstat <- granges(allData.tstat)
      cat("\t Tstat object loaded\n")
      ## Converting DMRs to Granges
      dmrsGr <- GRanges(dmrs$chr, IRanges(dmrs$start, dmrs$end))
      dmrsGrExtended <- GRanges(dmrs$chr, IRanges(dmrs$start-1000, dmrs$end+1000))

      ## Overlap tstat object and all DMRs. Resulting object has query hits (row in tstat) and subject hits (row in DMRs) correspondance. 
      overlap <- findOverlaps(tstat, dmrsGr)
      overlap <- as.data.frame(overlap)
      overlapExtended <- findOverlaps(tstat, dmrsGrExtended)
      overlapExtended <- as.data.frame(overlapExtended)
      cat("\t Overlap objects created\n")   

      ## create vector for orthology status
      orth <- rep(0, times=length(dmrs[,1]))
      names(orth) <- 1:length(dmrs[,1])
      orthExtended <- rep(0, times=length(dmrs[,1]))
      names(orthExtended) <- 1:length(dmrs[,1])

      ## Orthology only necessitates 1 site with t-stat to overlap = at least 1 subjectHit in overlap object 
      orth[sort(unique(overlap$subjectHits))] <- 1
      orthExtended[sort(unique(overlapExtended$subjectHits))] <- 1

      cat("\t", sum(orth == 1), "out of", length(orth), "tDMRs have relaxed orthology to", species2, "\n")
      allOrth[[paste0(species2, "Relaxed")]] <- orth
      cat("\t", sum(orthExtended == 1), "out of", length(orthExtended), "tDMRs have relaxed+extended orthology to", species2, "\n")
      allOrth[[paste0(species2, "ExtendedRelaxed")]] <- orthExtended
    }
    ## Export orthology results into file
    write.table(cbind(dmrs[,1:3], allOrth), file=paste0(baseName,"_DMRs_orthology_relaxed+extended.txt"), sep="\t", quote=F, row.names=F, col.names=T)
  }
}

## Loop across all contrasts to reproduce table similar to conservation_of_tissueDMRs.txt
make_summary_table_cons_relaxed <- function(filterFileExtension=NA, filterColumn=NA){
  counts <- as.data.frame(matrix(ncol=14, nrow=120)) ## 10 contrasts  * 3 species * 4 criteria (cons, consRelaxed, consExtended, consExtendedRelaxed)
  names(counts) <- c("Contrast", "Species", "Criterion", "Number of DMRs", "Number of DMRs with no orthology", "Number of DMRs with human-chimp orthology", "Number of DMRs with human-rhesus orthology", "Number of DMRs with chimp-rhesus orthology", "Number of DMRs with human-chimp-rhesus orthology", "Number of species-specific DMRs", "Number of human-chimp DMRs", "Number of human-rhesus DMRs", "Number of chimp-rhesus DMRs", "Number of human-chimp-rhesus DMRs")
  tissuePairs <- list(c("heart", "liver"), c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "Specific"), c("liver", "Specific"), c("kidney", "Specific"), c("lung", "Specific"))
  allSpecies <- c("Human", "Chimp", "Rhesus")
  i <- 1
  for (species in allSpecies){
    for (pair in tissuePairs){
      cat(species, ":", pair, "\n")
      if (pair[2] == "Specific"){
        baseName <- paste0("./DMRs/tissues/", species, "_", pair[1], pair[2])
        counts[(4*i - 3):(4*i), 1] <- paste0(pair[1], pair[2])
      } else {
        baseName <- paste0("./DMRs/tissues/", species, "_", pair[1], "_", pair[2])
        counts[(4*i - 3):(4*i), 1] <- paste0(pair[1], "_", pair[2])
      }
      counts[(4*i - 3):(4*i), 2] <- species

      ## vector with the 2 other species
      allSpecies2 <- allSpecies[-grep(species, allSpecies)]

      ## Conservation    
      cons <- list()
      cons[[1]] <- read.table(paste0(baseName,"_DMRs_conservation_to_", allSpecies2[1], "_relaxed+extended.txt"), h=T, sep="\t")
      cons[[2]] <- read.table(paste0(baseName,"_DMRs_conservation_to_", allSpecies2[2], "_relaxed+extended.txt"), h=T, sep="\t")
      names(cons) <- allSpecies2
      ## Sort DMRs
      cons[[1]] <- cons[[1]][order(cons[[1]]$chr, cons[[1]]$start, cons[[1]]$end),]
      cons[[2]] <- cons[[2]][order(cons[[2]]$chr, cons[[2]]$start, cons[[2]]$end),]
      
      ## Orthology
      orth <- read.table(paste0(baseName,"_DMRs_orthology_relaxed+extended.txt"), h=T)
      ## Sort DMRs
      orth <- orth[order(orth$chr, orth$start, orth$end),]

      if(!is.na(filterFileExtension) & !is.na(filterColumn[1])){
        ## Filter list of DMRs based on some feature annotation
        ## Several filter columns allowed (e.g., promoter & CGI)
        features <- read.table(paste0("../annotation/", baseName, filterFileExtension), h=T, sep="\t", check.names = F)
        features <- features[order(features$chr, features$start, features$end),]
        ## Problem: in some files, we have 0/1 encoding, in some others we have NA/name of feature: replace 0s by NAs
        filter <- rep(TRUE, times=length(orth[,1]))
        for (j in 1:length(filterColumn)){
          if (!is.null(features[[filterColumn[j]]])){
            features[[filterColumn[j]]][features[[filterColumn[j]]] == 0] <- NA
            filter[is.na(features[[filterColumn[j]]])] <- FALSE
          } else {
            cat(paste0("  Warning: there seem to be no ", filterColumn[j], " column in the file used for filtering\n"))
          }
        }
        ## Keeps only overlapping DMRs
        cons[[1]] <- cons[[1]][filter,]
        cons[[2]] <- cons[[2]][filter,]
        orth <- orth[filter,]
      }

      counts[4*i - 3, 3] <- "cons"
      counts[4*i - 2, 3] <- "consRelaxed"
      counts[4*i - 1, 3] <- "consExtended"
      counts[4*i    , 3] <- "consExtendedRelaxed"

      ## Total number of DMRs
      counts[(4*i - 3):(4*i), 4] <- length(cons[[allSpecies2[1]]][, 1])
      ## Number of DMRs with no orthology to any other species
      counts[(4*i - 3):(4*i - 2), 5] <- sum(apply(orth[,4:5], 1, sum) == 0)
      counts[(4*i - 1):(4*i)    , 5] <- sum(apply(orth[,6:7], 1, sum) == 0)
      ## Number of DMRs with human-chimp orthology
      if (species == "Human" | species == "Chimp"){
        counts[(4*i - 3):(4*i - 2), 6] <- sum(orth[[paste0(allSpecies2[2], "Relaxed")]] == 1)
        counts[(4*i - 1):(4*i)    , 6] <- sum(orth[[paste0(allSpecies2[2], "ExtendedRelaxed")]] == 1)
      } else if (species == "Rhesus"){
        counts[(4*i - 3):(4*i), 6] <- NA
      }
      ## Number of DMRs with human-rhesus orthology
      if (species == "Human"){
        counts[(4*i - 3):(4*i - 2), 7] <- sum(orth[[paste0(allSpecies2[2], "Relaxed")]] == 1)
        counts[(4*i - 1):(4*i)    , 7] <- sum(orth[[paste0(allSpecies2[2], "ExtendedRelaxed")]] == 1)
      } else if (species == "Chimp"){
        counts[(4*i - 3):(4*i), 7] <- NA
      } else if (species == "Rhesus"){
        counts[(4*i - 3):(4*i - 2), 7] <- sum(orth[[paste0(allSpecies2[1], "Relaxed")]] == 1)
        counts[(4*i - 1):(4*i)    , 7] <- sum(orth[[paste0(allSpecies2[1], "ExtendedRelaxed")]] == 1)
      }
      ## Number of DMRs with chimp-rhesus orthology
      if (species == "Rhesus" | species == "Chimp"){
        counts[(4*i - 3):(4*i -2 ), 8] <- sum(orth[[paste0(allSpecies2[2], "Relaxed")]] == 1)
        counts[(4*i - 1):(4*i)    , 8] <- sum(orth[[paste0(allSpecies2[2], "ExtendedRelaxed")]] == 1)
      } else if (species == "Human"){
        counts[(4*i - 3):(4*i), 8] <- NA
      }
      ## Number of DMRs with human-chimp-rhesus orthology
      counts[(4*i - 3):(4*i - 2), 9] <- sum(apply(orth[,4:5], 1, sum) == 2)
      counts[(4*i - 1):(4*i)    , 9] <- sum(apply(orth[,6:7], 1, sum) == 2)

      ## Conservation encoding:
      ## for human:  unique (0), human-chimp (1),  human-rhesus (2) and human-chimp-rhesus (4)
      ## for chimp:  unique (0), human-chimp (1),  chimp-rhesus (3) and human-chimp-rhesus (4)
      ## for rhesus: unique (0), human-rhesus (2), chimp-rhesus (3) and human-chimp-rhesus (4)
      ## This numbers are filtered for orthology status. Simplest is to consider only DMRs with HCR orthology status. Otherwise we get into complex considerations sicne the reference is changing. Additionally conservation can be found even without strict orthology criteria met
      
      ## Number of species-specific DMRs
      counts[4*i - 3, 10] <- sum(cons[[allSpecies2[1]]]$cons == 0 & cons[[allSpecies2[2]]]$cons == 0 & apply(orth[,4:5], 1, sum) == 2)
      counts[4*i - 2, 10] <- sum(cons[[allSpecies2[1]]]$consRelaxed == 0 & cons[[allSpecies2[2]]]$consRelaxed == 0 & apply(orth[,4:5], 1, sum) == 2)
      counts[4*i - 1, 10] <- sum(cons[[allSpecies2[1]]]$consExtended == 0 & cons[[allSpecies2[2]]]$consExtended == 0 & apply(orth[,6:7], 1, sum) == 2)
      counts[4*i,     10] <- sum(cons[[allSpecies2[1]]]$consExtendedRelaxed == 0 & cons[[allSpecies2[2]]]$consExtendedRelaxed == 0 & apply(orth[,6:7], 1, sum) == 2)
      ## Number of human-chimp DMRs
      if (species == "Human" | species == "Chimp"){
        counts[4*i - 3, 11] <- sum(cons[[allSpecies2[1]]]$cons == 1 & apply(orth[,4:5], 1, sum) == 2)
        counts[4*i - 2, 11] <- sum(cons[[allSpecies2[1]]]$consRelaxed == 1 & apply(orth[,4:5], 1, sum) == 2)
        counts[4*i - 1, 11] <- sum(cons[[allSpecies2[1]]]$consExtended == 1 & apply(orth[,6:7], 1, sum) == 2)
        counts[4*i    , 11] <- sum(cons[[allSpecies2[1]]]$consExtendedRelaxed == 1 & apply(orth[,6:7], 1, sum) == 2)
      } else if (species == "Rhesus"){
        counts[(4*i - 3):(4*i), 11] <- NA
      } 
      ## Number of human-rhesus DMRs
      if (species == "Human"){
        counts[4*i - 3, 12] <- sum(cons[[allSpecies2[2]]]$cons == 1 & apply(orth[,4:5], 1, sum) == 2)
        counts[4*i - 2, 12] <- sum(cons[[allSpecies2[2]]]$consRelaxed == 1 & apply(orth[,4:5], 1, sum) == 2)
        counts[4*i - 1, 12] <- sum(cons[[allSpecies2[2]]]$consExtended == 1 & apply(orth[,6:7], 1, sum) == 2)
        counts[4*i    , 12] <- sum(cons[[allSpecies2[2]]]$consExtendedRelaxed == 1 & apply(orth[,6:7], 1, sum) == 2)
      } else if (species == "Rhesus"){
        counts[4*i - 3, 12] <- sum(cons[[allSpecies2[1]]]$cons == 1 & apply(orth[,4:5], 1, sum) == 2)
        counts[4*i - 2, 12] <- sum(cons[[allSpecies2[1]]]$consRelaxed == 1 & apply(orth[,4:5], 1, sum) == 2)
        counts[4*i - 1, 12] <- sum(cons[[allSpecies2[1]]]$consExtended == 1 & apply(orth[,6:7], 1, sum) == 2)
        counts[4*i    , 12] <- sum(cons[[allSpecies2[1]]]$consExtendedRelaxed == 1 & apply(orth[,6:7], 1, sum) == 2)
      } else if (species == "Chimp"){
        counts[(4*i - 3):(4*i), 12] <- NA
      }
      ## Number of chimp-rhesus DMRs
      if (species == "Human"){
        counts[(4*i - 3):(4*i), 13] <- NA
      } else if (species == "Chimp" | species == "Rhesus"){
        counts[4*i - 3, 13] <- sum(cons[[allSpecies2[2]]]$cons == 1 & apply(orth[,4:5], 1, sum) == 2)
        counts[4*i - 2, 13] <- sum(cons[[allSpecies2[2]]]$consRelaxed == 1 & apply(orth[,4:5], 1, sum) == 2)
        counts[4*i - 1, 13] <- sum(cons[[allSpecies2[2]]]$consExtended == 1 & apply(orth[,6:7], 1, sum) == 2)
        counts[4*i    , 13] <- sum(cons[[allSpecies2[2]]]$consExtendedRelaxed == 1 & apply(orth[,6:7], 1, sum) == 2)
      } 
      ## Number of human-chimp-rhesus DMRs 
      counts[(4*i - 3):(4*i), 14] <- sum(cons$conservation == 4 & apply(orth[,4:5], 1, sum) == 2)

      counts[4*i - 3, 14] <- sum(cons[[allSpecies2[1]]]$cons == 1 & cons[[allSpecies2[2]]]$cons == 1 & apply(orth[,4:5], 1, sum) == 2)
      counts[4*i - 2, 14] <- sum(cons[[allSpecies2[1]]]$consRelaxed == 1 & cons[[allSpecies2[2]]]$consRelaxed == 1 & apply(orth[,4:5], 1, sum) == 2)
      counts[4*i - 1, 14] <- sum(cons[[allSpecies2[1]]]$consExtended == 1 & cons[[allSpecies2[2]]]$consExtended == 1 & apply(orth[,6:7], 1, sum) == 2)
      counts[4*i,     14] <- sum(cons[[allSpecies2[1]]]$consExtendedRelaxed == 1 & cons[[allSpecies2[2]]]$consExtendedRelaxed == 1 & apply(orth[,6:7], 1, sum) == 2)
      
      ## increment row
      i <- i + 1
    }
  }
  ## We want the rows to sum up to 1: subtract HCR to HC, HR and CR
  counts[,11] <- counts[,11] - counts[,14]
  counts[,12] <- counts[,12] - counts[,14]
  counts[,13] <- counts[,13] - counts[,14]
  
  return(counts)
}

counts <- make_summary_table_cons_relaxed(NA, NA)
write.table(counts, file="conservation_of_tissueDMRs_with_HCR_orthology_relaxed.txt", sep="\t", quote=F, row.names=F, col.names=T)
## counts <- read.table("conservation_of_tissueDMRs_with_HCR_orthology_relaxed.txt", sep="\t", h=T, check.names=F)

## Plots
## We need to input a matrix with 13 column: remove "Criterion" column
make_conservation_plots(counts[counts$Criterion == "cons", -3], "1site")
## Very similar to strict conservation (above), but actually slightly lower for pairwise contrasts! Slightly higher for tissue-specific contrasts 
make_conservation_plots(counts[counts$Criterion == "consRelaxed", -3], "1site_relaxed")
## A lot higher!
make_conservation_plots(counts[counts$Criterion == "consExtended", -3], "1site_extended")
## A bit higher than "1site"
make_conservation_plots(counts[counts$Criterion == "consExtendedRelaxed", -3], "1site_relaxed+extended")
## A bit higher than relaxed

## It is kind of surprising that Actually we get slightly lower conservation for "1site" compared to strict conservation (overlap of tDMRs): this is because the total set of DMRs with orthology is higher.
## To compare conservation on the total set of DMRs with strict HCR orthology (3 CpG sites, etc), we restrict set of DMRs using orthology file:
counts <- make_summary_table_cons_relaxed("_DMRs_orthology.txt", c("Human", "Chimp", "Rhesus"))
## Warnings are issued, because there is always orthology to 2 species, not 3 as specified in the columns
make_conservation_plots(counts[counts$Criterion == "cons", -3], "1site_with_strict_HCR_orthology")
## Indeed, the conservation is higher with this relaxed criterion (1 site) 

## repeat plots, for subset of DMRs overlapping specific features
counts <- make_summary_table_cons_relaxed("_DMRs_15_features_0.2.gz", "promoter")
make_conservation_plots(counts[counts$Criterion == "cons", -3], "1site_promoter")

counts <- make_summary_table_cons_relaxed("_DMRs_15_features_0.2.gz", "promoter_proximal")
make_conservation_plots(counts[counts$Criterion == "cons", -3], "1site_promoter_proximal")

## Conserved promoters
counts <- make_summary_table_cons_relaxed("_DMRs_conserved_features.gz", "promoter_conserved_panTro3+rheMac2")
make_conservation_plots(counts[counts$Criterion == "cons", -3], "1site_promoter_conserved")

## Conserved CGI
counts <- make_summary_table_cons_relaxed("_DMRs_conserved_features.gz", "CGI_conserved_panTro3+rheMac2")
make_conservation_plots(counts[counts$Criterion == "cons", -3], "CGI_conserved")

## Conserved repeats
counts <- make_summary_table_cons_relaxed("_DMRs_conserved_features.gz", "repeats_conserved_panTro3+rheMac2")
make_conservation_plots(counts[counts$Criterion == "cons", -3], "repeats_conserved")

## TO DO: tissueDMRs next to conserved TSSs of 1-to-1 orthologous genes that show DE
## TO DO: tissueDMRs next to conserved TSSs of 1-to-1 orthologous genes that show conserved DE


###################################################################
## TO DO: what is below may be a bit finer: look at distance to TSS
##        repeat same as below but for all contrasts + relaxed conservation criteria
## 
## ## Look at conservation of human tDMRs that are close to genes ##
## ## human heart vs. liver (22932 DMRs)
## ## The code below is similar to code used in ../DE_vs_DMRs/analysis.R
## closest <- read.table("../DE_vs_DMRs/Human_heart_liver_DMRs_closest_promoter_2kb.gz", h=F, sep="\t")
## ## merge with DMR table
## dmrs <- read.table("./DMRs/tissues/Human_heart_liver_DMRs.txt", sep="\t", h=T)
## dmrs <- cbind(dmrs, closest[, 4:7])
## names(dmrs)[17:20] <- c("geneID","geneStart","geneEnd","distance")
## ## Merge with features (to be sure DMR overlaps 20% of promoter)
## features <- read.table("../annotation/DMRs/tissues/Human_heart_liver_DMRs_15_features_0.2.gz", sep="\t", h=T, check.names=F)
## ## not same ordering as dmrs
## row.names(features) <- paste(features[,1], features[,2], sep=":")
## features <- features[paste(dmrs[,1], dmrs[,2], sep=":"),]
## summary(paste(features[,1], features[,2], sep=":") == paste(dmrs[,1], dmrs[,2], sep=":"))## now it's good
## dmrs <- cbind(dmrs, features[, 4:18])
## ## Conservation info
## conservation <- read.table("../bsseq/DMRs/tissues/Human_heart_liver_DMRs_conservation.txt", sep="\t", h=T)
## dim(conservation)
## summary(paste(conservation[,1], conservation[,2], sep=":") == paste(dmrs[,1], dmrs[,2], sep=":"))  
## row.names(conservation) <- paste(conservation[,1], conservation[,2], sep=":")
## dmrs <- cbind(dmrs, conservation[, 4])
## names(dmrs)[36] <- "conservation"
## summary(as.factor(dmrs$conservation))
## ##   0    1    2    4 
## ## 9175 5345 2732 5680
## ## 0: H only
## ## 1: HC
## ## 2: HR
## ## 4: HCR
## ## DMRs conserved in chimp
## sum(dmrs$conservation == 1 | dmrs$conservation == 4)/length(dmrs[,1]) # 48%
## ## DMRs conserved in rhesus
## sum(dmrs$conservation == 2 | dmrs$conservation == 4)/length(dmrs[,1]) # 37%

## ## Relaxed conservation info (see above)
## consChimp <- read.table("conservation_of_tissueDMRs_human_heart_liver_relaxed.txt", h=T, sep="\t")[,4]
## consChimpExtended <- read.table("conservation_of_tissueDMRs_human_heart_liver_relaxed.txt", h=T, sep="\t")[,5]
## consRhesus <- read.table("conservation_of_tissueDMRs_human_heart_liver_relaxed.txt", h=T, sep="\t")[,6]
## consRhesusExtended <- read.table("conservation_of_tissueDMRs_human_heart_liver_relaxed.txt", h=T, sep="\t")[,7]
## dmrs <- cbind(dmrs, consChimp, consChimpExtended, consRhesus, consRhesusExtended)
## ## DMRs conserved in chimp
## sum(dmrs$consChimp)/length(dmrs[,1]) # 50.5%
## ## DMRs conserved in rhesus
## sum(dmrs$consRhesus)/length(dmrs[,1]) # 38%

## ## check consistency of both conservation measures
## summary((dmrs$conservation == 1 | dmrs$conservation == 4) &  dmrs$consChimp == 1)
## summary((dmrs$conservation == 2 | dmrs$conservation == 4) &  dmrs$consRhesus == 1)
## summary(dmrs$conservation == 0 & dmrs$consChimp == 1)
## ## 375 newly discovered conserved in chimp DMRs only
## summary(dmrs$conservation == 0 & dmrs$consRhesus == 1 & dmrs$consChimp == 1)
## ## 34 newly discovered conserved DMRs only
## summary(dmrs$conservation == 0 & dmrs$consRhesusExtended == 1 & dmrs$consChimpExtended == 1)
## ## 284 newly discovered conserved DMRs including neighboring CpGs

## ## Overlap with DE genes, with and without considering direction (i.e., hyper methylation associted with lower expression)
## DE <- read.table("../limma/gene_lists/genes/DMRs_comparison/HumanHeartvsLiver.txt", sep="\t", h=T)
## row.names(DE) <- DE[,1]
## ## Conservation of DE
## cons <- read.table("../limma/gene_lists/genes/DMRs_comparison/HumanHeartvsLiver_conservation.txt", sep="\t", row.names=1, h=T, check.names=F)
## summary(row.names(cons) == row.names(DE))
## DE <- cbind(DE, cons[, 1:5])
## ## Merge tables
## mergedDMRsExpr <- merge(dmrs, DE, by.x="geneID", by.y="genes")
## dim(mergedDMRsExpr) ## 14143 DMRs next to genes with expression data


## ## Look at conservation of human tDMRs that are close to genes (overlap promoter) ################
## ## proximal promoter: same as promoter
## subset <- dmrs$promoter == 1 ## 9709
## ## DMRs conserved in chimp
## sum(dmrs$conservation[subset] == 1 | dmrs$conservation[subset] == 4)/length(dmrs[subset,1]) # 51%
## ## DMRs conserved in rhesus
## sum(dmrs$conservation[subset] == 2 | dmrs$conservation[subset] == 4)/length(dmrs[subset,1]) # 40%
## ## Slightly better!
## seq <- rev(brewer.pal(9, "Blues")) ## sequential colors
## pdf(file = "conservation_of_tissueDMRs_human_heart_liver_promoter.pdf", width = 6, height = 5)
## barplot(rep(sum(subset), 3), col=c(seq[6], NA, NA), names.arg=c("Human", "Human to chimp", "Human to macaque"))
## barplot(c(0, sum(dmrs$conservation[subset] == 1 | dmrs$conservation[subset] == 4), 0), col=c(NA, seq[4], NA), add=T, axes=F)
## barplot(c(0, 0, sum(dmrs$conservation[subset] == 2 | dmrs$conservation[subset] == 4)), col=c(NA, NA, seq[2]), add=T, axes=F)
## dev.off()

## ## Look at conservation of human tDMRs that are close to genes with expression
## ## These genes are 1-to-1 orthologous, so in practice  we're restricting to tDMRs next to 1-to-1 orthologus genes
## subset <- mergedDMRsExpr$promoter == 1 ## 6991
## sum(mergedDMRsExpr$conservation[subset] == 1 | mergedDMRsExpr$conservation[subset] == 4)/length(mergedDMRsExpr[subset,1]) # 54%
## sum(mergedDMRsExpr$conservation[subset] == 2 | mergedDMRsExpr$conservation[subset] == 4)/length(mergedDMRsExpr[subset,1]) # 43%
## ## Slightly better next to genes which have expression
## seq <- rev(brewer.pal(9, "Blues")) ## sequential colors
## pdf(file = "conservation_of_tissueDMRs_human_heart_liver_promoter_expression.pdf", width = 6, height = 5)
## barplot(rep(sum(subset), 3), col=c(seq[6], NA, NA), names.arg=c("Human", "Human to chimp", "Human to macaque"))
## barplot(c(0, sum(mergedDMRsExpr$conservation[subset] == 1 | mergedDMRsExpr$conservation[subset] == 4), 0), col=c(NA, seq[4], NA), add=T, axes=F)
## barplot(c(0, 0, sum(mergedDMRsExpr$conservation[subset] == 2 | mergedDMRsExpr$conservation[subset] == 4)), col=c(NA, NA, seq[2]), add=T, axes=F)
## dev.off()


## ## Look at conservation of human tDMRs that are close to DE genes
## subset <- mergedDMRsExpr$promoter == 1 & mergedDMRsExpr$adj.P.Val < 0.01 ## 3662
## sum(mergedDMRsExpr$conservation[subset] == 1 | mergedDMRsExpr$conservation[subset] == 4)/length(mergedDMRsExpr[subset,1]) # 58%
## sum(mergedDMRsExpr$conservation[subset] == 2 | mergedDMRsExpr$conservation[subset] == 4)/length(mergedDMRsExpr[subset,1]) # 47%
## ## Slightly better next to genes which are DE

## ## Look at conservation of human tDMRs that are close to DE genes (check direction)
## subset <- mergedDMRsExpr$promoter == 1 & mergedDMRsExpr$adj.P.Val < 0.01 & ifelse(mergedDMRsExpr$logFC > 0, "hypo", "hyper") == mergedDMRsExpr$direction ## 2626
## sum(mergedDMRsExpr$conservation[subset] == 1 | mergedDMRsExpr$conservation[subset] == 4)/length(mergedDMRsExpr[subset,1]) # 61%
## sum(mergedDMRsExpr$conservation[subset] == 2 | mergedDMRsExpr$conservation[subset] == 4)/length(mergedDMRsExpr[subset,1]) # 50%
## ## Even better
## pdf(file = "conservation_of_tissueDMRs_human_heart_liver_promoter_DE.pdf", width = 6, height = 5)
## barplot(rep(sum(subset), 3), col=c(seq[6], NA, NA), names.arg=c("Human", "Human to chimp", "Human to macaque"))
## barplot(c(0, sum(mergedDMRsExpr$conservation[subset] == 1 | mergedDMRsExpr$conservation[subset] == 4), 0), col=c(NA, seq[4], NA), add=T, axes=F)
## barplot(c(0, 0, sum(mergedDMRsExpr$conservation[subset] == 2 | mergedDMRsExpr$conservation[subset] == 4)), col=c(NA, NA, seq[2]), add=T, axes=F)
## dev.off()

## ## Look at conservation of human tDMRs that are close to conserved DE genes
## subset <- mergedDMRsExpr$promoter == 1 & mergedDMRsExpr$adj.P.Val < 0.01 & ifelse(mergedDMRsExpr$logFC > 0, "hypo", "hyper") == mergedDMRsExpr$direction & mergedDMRsExpr$"ChimpDE_10%" == 1 & mergedDMRsExpr$"RhesusDE_10%" == 1 ## 2181

## sum(mergedDMRsExpr$conservation[subset] == 1 | mergedDMRsExpr$conservation[subset] == 4)/length(mergedDMRsExpr[subset,1]) # 62.5%
## sum(mergedDMRsExpr$conservation[subset] == 2 | mergedDMRsExpr$conservation[subset] == 4)/length(mergedDMRsExpr[subset,1]) # 53%
## pdf(file = "conservation_of_tissueDMRs_human_heart_liver_promoter_DE_conserved.pdf", width = 6, height = 5)
## barplot(rep(sum(subset), 3), col=c(seq[6], NA, NA), names.arg=c("Human", "Human to chimp", "Human to macaque"))
## barplot(c(0, sum(mergedDMRsExpr$conservation[subset] == 1 | mergedDMRsExpr$conservation[subset] == 4), 0), col=c(NA, seq[4], NA), add=T, axes=F)
## barplot(c(0, 0, sum(mergedDMRsExpr$conservation[subset] == 2 | mergedDMRsExpr$conservation[subset] == 4)), col=c(NA, NA, seq[2]), add=T, axes=F)
## dev.off()

## ## Same thing but using relaxed conservation
## subset <- mergedDMRsExpr$promoter == 1 & mergedDMRsExpr$adj.P.Val < 0.01 & ifelse(mergedDMRsExpr$logFC > 0, "hypo", "hyper") == mergedDMRsExpr$direction & mergedDMRsExpr$"ChimpDE_10%" == 1 & mergedDMRsExpr$"RhesusDE_10%" == 1 ## 2181
## sum(mergedDMRsExpr$consChimp[subset] == 1)/length(mergedDMRsExpr[subset,1]) # 65.9%
## sum(mergedDMRsExpr$consRhesus[subset] == 1)/length(mergedDMRsExpr[subset,1]) # 55.6%
## sum(mergedDMRsExpr$consChimpExtended[subset] == 1)/length(mergedDMRsExpr[subset,1]) # 71%
## sum(mergedDMRsExpr$consRhesusExtended[subset] == 1)/length(mergedDMRsExpr[subset,1]) # 62.7%

## ## TO DO: calculate % conservation next to genes/DE genes/conserved DE genes for relaxed and relaxed+extended
## ## pdf(file = "conservation_of_tissueDMRs_human_heart_liver_promoter_relaxed.pdf", width = 6, height = 5)
## ## pdf(file = "conservation_of_tissueDMRs_human_heart_liver_promoter_DE_relaxed.pdf", width = 6, height = 5)
## ## pdf(file = "conservation_of_tissueDMRs_human_heart_liver_promoter_DE_conserved_relaxed.pdf", width = 6, height = 5)
## ## barplot(rep(counts$Number.of.DMRs, 3), col=c(seq[6], NA, NA), names.arg=c("Human", "Human to chimp", "Human to macaque"))
## ## barplot(c(0, ..., 0), col=c(NA, seq[4], NA), add=T, axes=F)
## ## barplot(c(0, 0, ...), col=c(NA, NA, seq[2]), add=T, axes=F)
## ## dev.off()

## ## pdf(file = "conservation_of_tissueDMRs_human_heart_liver_promoter_relaxed_extended.pdf", width = 6, height = 5)
## ## pdf(file = "conservation_of_tissueDMRs_human_heart_liver_promoter_DE_relaxed_extended.pdf", width = 6, height = 5)
## ## pdf(file = "conservation_of_tissueDMRs_human_heart_liver_promoter_DE_conserved_relaxed_extended.pdf", width = 6, height = 5)
## ## barplot(rep(counts$Number.of.DMRs, 3), col=c(seq[6], NA, NA), names.arg=c("Human", "Human to chimp", "Human to macaque"))
## ## barplot(c(0, ..., 0), col=c(NA, seq[4], NA), add=T, axes=F)
## ## barplot(c(0, 0, ...), col=c(NA, NA, seq[2]), add=T, axes=F)
## ## dev.off()


## ... ###################  
## TO DO: tie to orthologous genes: the closest gene in chimp should be orthologous to closest gene in human
##        this would also help to test if there is "replacement" of tDMRs associated to same gene


###############################################
## DMRs: length, effect size and CpG density ##
###############################################

## length
plot_speciesDMRs_characteristics("~/clusterhome/Methylation/bsseq/length_speciesDMRs", "width")
## species-specific are shorter: because intersection of 2 pairwise contrasts
## For final figures in paper
plot_speciesDMRs_characteristics2("~/Methylation/bsseq/length_speciesDMRs", "width", legend=T)


## CpG density
plot_speciesDMRs_characteristics("~/clusterhome/Methylation/bsseq/CpG_density_speciesDMRs", "invdensity")
plot_speciesDMRs_characteristics2("~/Methylation/bsseq/CpG_density_speciesDMRs", "invdensity", legend=F)

## Mean diff in % methylation
plot_speciesDMRs_characteristics("~/clusterhome/Methylation/bsseq/mean_methylation_difference_speciesDMRs", "meanDiff")
plot_speciesDMRs_characteristics2("~/Methylation/bsseq/mean_methylation_difference_speciesDMRs", "meanDiff", legend=F)
## liver has higher effect size
## assymetry over/under methylation. Especially for species-specific changes. Makes sense since they are oriented!
## species-specific: higher effect size: methodological?
## human-chimp lower effect size: increase of difference in methylation with evolutionary age?
## TO DO? take absolute value of mean diff in % meth? No, we would loose information!


## Tissue DMRs
allSpecies <- c("Chimp", "Human", "Rhesus")
tissuePairs <- list(c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "liver"), c("heart", "Specific"), c("liver", "Specific"), c("kidney", "Specific"), c("lung", "Specific"))
for (species in allSpecies){
  for (pair in tissuePairs){
    ## cat("Testing tissue pair: ", pair, "in", species, "\n")
    ## Reading DMR file
    if (pair[2] == "Specific"){
      dmrs <- read.table(paste0("./DMRs/tissues/", species, "_", pair[1], pair[2], "_DMRs.txt"), sep="\t", h=T)
    }
    else {
      dmrs <- read.table(paste0("./DMRs/tissues/", species, "_", pair[1], "_", pair[2], "_DMRs.txt"), sep="\t", h=T)   
    }
    ## cat(pair, "in", species, "median length: ", median(dmrs$width), "\n")
    ## cat(pair, "in", species, "density: ", median(dmrs$invdensity), "\n")
    cat(pair, "in", species, "absolute effect size: ", median(abs(dmrs$meanDiff)), "\n")
  }
}
  
## length
plot_tissueDMRs_characteristics("~/clusterhome/Methylation/bsseq/length_tissueDMRs", "width")
## Tissue-specific slightly shorter than pairwise comparisons. Especially, lung-specific DMRs are shorter
## Kidney-specific DMRs in Chimp: very short (strange)
## Patterns very similar between species
plot_tissueDMRs_characteristics2("~/Methylation/bsseq/length_tissueDMRs", "width", legend=T, 2)


## CpG density
plot_tissueDMRs_characteristics("~/clusterhome/Methylation/bsseq/CpG_density_tissueDMRs", "invdensity")
plot_tissueDMRs_characteristics2("~/Methylation/bsseq/CpG_density_tissueDMRs", "invdensity", legend=F, 1.5)
## overall similar between species and between contrasts... liver contrasts are a bit strange in Human

## Mean diff in % methylation
plot_tissueDMRs_characteristics("~/clusterhome/Methylation/bsseq/mean_methylation_difference_tissueDMRs", "meanDiff")
plot_tissueDMRs_characteristics2("~/Methylation/bsseq/mean_methylation_tissueDMRs", "meanDiff", legend=F, 1.7)
## Overall similar between species.
## tissue-specific contrasts very enriched for under-methylation (some more than others)
## pairwise contrasts not always symetrical: depending on characteristics of tissue compared


## ## Individual DMRs
## ## length
## plot_individualDMRs_characteristics("~/clusterhome/Methylation/bsseq/length_individualDMRs", "width")
## ## quite short DMRs!
## ## strange things happening with R2/3
## ## CpG density
## plot_individualDMRs_characteristics("~/clusterhome/Methylation/bsseq/CpG_density_individualDMRs", "invdensity")
## ## Mean diff in % methylation
## plot_individualDMRs_characteristics("~/clusterhome/Methylation/bsseq/mean_methylation_difference_individualDMRs", "meanDiff")
## ## quite symetrical, except contrasts involving C3


###########################################
## overlap of DMRs with genomic features ##
###########################################
## See folder ~/Methylation/annotation/
## We have for each set of DMRs
##   - the number of DMRs that overlap different genomic features
##   - 100 sets of randomized DMRs of same length (on same chromosome)  and their overlap different genomic features

## Idea: for each type of DMR, and for each feature, draw a boxplot of number of overlaps from random sets of DMRs. Compare to the "real" values for tested DMR. This gives an idea of significance.

## first try on 1 particular DMR: Human vs. Chimp in liver ##########################
fileName <- "../annotation/DMRs/species/HumanChimp_liver_DMRs_15_features.gz"
realFeatures <- read.table(fileName, h=T, sep="\t", check.names = F)
nDMRs <- length(realFeatures[,1])
summary(realFeatures[,-c(1:3)] == 1)
##    exon         first_exon         CDS            3_UTR        
##  FALSE:8093      FALSE:8757      FALSE:8824      FALSE:9298     
##  TRUE :1703      TRUE :1039      TRUE :972       TRUE :498      

##    5_UTR           intron        intergenic      conserved      
##  FALSE:9190      FALSE:4914      FALSE:6108      FALSE:6854     
##  TRUE :606       TRUE :4882      TRUE :3688      TRUE :2942     

##  CpG_island      CpG_island_shore CpG_island_shelf   repeat       
##  FALSE:8751      FALSE:5390       FALSE:8931       FALSE:3989     
##  TRUE :1045      TRUE :4406       TRUE :865        TRUE :5807     

##   promoter       promoter_coding promoter_proximal
##  FALSE:4973      FALSE:5799      FALSE:8465       
##  TRUE :4823      TRUE :3997      TRUE :1331       

realFeatures <- apply(realFeatures[,-c(1:3)], 2, sum)

randomFeatures <- matrix(nrow=100, ncol=15)
colnames(randomFeatures) = names(realFeatures)
for (i in 1:100){
  print(i)
  randomFeatures[i,] <- apply(read.table(paste0("../annotation/DMRs/species/randomized/HumanChimp_liver_randomDMRs_", i, "_15_features.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
}

## Reordering features based on their "natural" order 
## randomFeatures <- randomFeatures[, c(9, 10, 11, 13, 14, 15, 5, 2, 1, 3, 6, 4, 7, 12, 8)]
## realFeatures <- realFeatures[c(9, 10, 11, 13, 14, 15, 5, 2, 1, 3, 6, 4, 7, 12, 8)]

## Reordering features based on fold enrichment / depletion
ordered <- order(realFeatures / apply(randomFeatures, 2, median), decreasing=T)
randomFeatures <- randomFeatures[, ordered]
realFeatures <- realFeatures[ordered]

## Plotting boxplots:
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))

pdf(file = "test.pdf", width = 6, height = 5)
par(mar=c(7, 5, 4, 2) + 0.1, mgp=c(4, 1, 0)) 
## let 7% room above for fold-change indication
ylimit <- c(0, max(realFeatures, randomFeatures)*1.08)
## draw empty plot to be able to add backrgoudn grey lines below boxes
plot(1:15, type="n", axes=F, xlab="", ylab="", ylim=ylimit, yaxs="i", xlim=c(1,15))
## add vertical gray dashed lines to help reading
abline(v=c(1:15), col="gray90", lty=3)
## add boxplot
mp <- boxplot(randomFeatures, pch=20, cex=0.8, notch=F, ylab="Number of overlapping DMRs", xlab="", xaxt="n", ylim=ylimit, yaxs="i", whisklty=1, medlty=0, col="dimgray", border="dimgray", add=T, las=2)
## text 6% below x-axis / substitute _ for spaces in colnames
axis(1, labels=NA, at=c(1:15))
featNames <- names(realFeatures)
featNames <- gsub("_UTR", "'_UTR", featNames, perl=TRUE)
featNames <- gsub("coding", "(coding)", featNames, perl=TRUE)
featNames <- gsub("proximal", "(proximal)", featNames, perl=TRUE)
featNames <- gsub("_", " ", featNames, perl=TRUE)
text(c(1:15), par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.06, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)
## plot the observed overlap for DMRs (diamond)
## red for enrichment, blue for depletion compared to median of randomized data)
points(realFeatures, col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), pch=18, cex=1.3)
text(c(1:15), par("usr")[4]*0.97, labels=round(realFeatures / apply(randomFeatures, 2, median),1), col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), cex=0.8)
## header: number of DMRs and name of contrast
## species contrast string treatment
titleToPrint <- gsub("_DMRs_15_features.gz", "", fileName, perl=TRUE)
titleToPrint <- gsub("\\.\\.\\/annotation\\/DMRs\\/.+\\/", "", titleToPrint, perl=TRUE)
titleToPrint <- gsub("(?!^)(?=[[:upper:]])", " vs. ", titleToPrint, perl=TRUE) ## see http://stackoverflow.com/questions/7988959/splitting-string-based-on-letters-case
titleToPrint <- gsub("_", " in ", titleToPrint, perl=TRUE)
title(paste0(titleToPrint, ": ", nDMRs, " DMRs"))
dev.off()


## TO DO?
## - log scale for y-axis?


## overlap at east 1bp and random DMRs of same length (but not same CpG density) ##########################
## Species DMRs
for (tissue in c("heart", "kidney", "liver", "lung")){
  for (pair in list(c("Human", "Chimp"), c("Human", "Rhesus"), c("Chimp", "Rhesus"), c("Human", "Specific"), c("Chimp", "Specific"))){
    cat("Testing species pair: ", pair, "in", tissue, "\n")
    fileName <- paste0("../annotation/DMRs/species/", pair[1], pair[2], "_", tissue, "_DMRs_15_features.gz")
    realFeatures <- read.table(fileName, h=T, sep="\t", check.names = F)
    nDMRs <- length(realFeatures[,1])
    realFeatures <- apply(realFeatures[,-c(1:3)], 2, sum)

    randomFeatures <- matrix(nrow=100, ncol=15)
    colnames(randomFeatures) = names(realFeatures)
    for (i in 1:100){
      randomFeatures[i,] <- apply(read.table(paste0("../annotation/DMRs/species/randomized/", pair[1], pair[2], "_", tissue, "_randomDMRs_", i, "_15_features.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
    }

    ## Reordering features based on fold enrichment / depletion
    ordered <- order(realFeatures / apply(randomFeatures, 2, median), decreasing=T)
    randomFeatures <- randomFeatures[, ordered]
    realFeatures <- realFeatures[ordered]

    ## Plotting boxplots:
    pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))

    pdf(file = paste0("../annotation/DMRs/species/", pair[1], pair[2], "_", tissue, "_DMRs_boxplot_15_features.pdf"), width = 6, height = 5)
    par(mar=c(7, 5, 4, 2) + 0.1, mgp=c(4, 1, 0)) 
    ## let 7% room above for fold-change indication
    ylimit <- c(0, max(realFeatures, randomFeatures)*1.08)
    ## draw empty plot to be able to add backrgoudn grey lines below boxes
    plot(1:15, type="n", axes=F, xlab="", ylab="", ylim=ylimit, yaxs="i", xlim=c(1,15))
    ## add vertical gray dashed lines to help reading
    abline(v=c(1:15), col="gray90", lty=3)
    ## add boxplot
    mp <- boxplot(randomFeatures, pch=20, cex=0.8, notch=F, ylab="Number of overlapping DMRs", xlab="", xaxt="n", ylim=ylimit, yaxs="i", whisklty=1, medlty=0, col="dimgray", border="dimgray", add=T, las=2)
    ## text 6% below x-axis / substitute _ for spaces in colnames
    axis(1, labels=NA, at=c(1:15))
    featNames <- names(realFeatures)
    featNames <- gsub("_UTR", "'_UTR", featNames, perl=TRUE)
    featNames <- gsub("coding", "(coding)", featNames, perl=TRUE)
    featNames <- gsub("proximal", "(proximal)", featNames, perl=TRUE)
    featNames <- gsub("_", " ", featNames, perl=TRUE)
    text(c(1:15), par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.06, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)
    ## plot the observed overlap for DMRs (diamond)
    ## red for enrichment, blue for depletion compared to median of randomized data)
    points(realFeatures, col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), pch=18, cex=1.3)
    text(c(1:15), par("usr")[4]*0.97, labels=round(realFeatures / apply(randomFeatures, 2, median),1), col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), cex=0.8)
    ## header: number of DMRs and name of contrast
    ## species contrast string treatment
    if (pair[2] == "Specific"){
      title(paste0(pair[1], "-specific in ", tissue, ": ", nDMRs, " DMRs"))
    }
    else {
      title(paste0(pair[1], " vs. ", pair[2], " in ", tissue, ": ", nDMRs, " DMRs"))
    }
    dev.off()
  }
}

## Tissue DMRs
for (species in c("Human", "Chimp", "Rhesus")){
  for (pair in list(c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "liver"), c("heart", "Specific"), c("liver", "Specific"), c("kidney", "Specific"), c("lung", "Specific"))){

    cat("Testing tissue pair: ", pair, "in", species, "\n")
    if (pair[2] != "Specific"){
      baseName <- paste0(species, "_", pair[1], "_", pair[2])
    }
    else {
      baseName <- paste0(species, "_", pair[1], "Specific")
    }      
    realFeatures <- read.table(paste0("../annotation/DMRs/tissues/", baseName, "_DMRs_15_features.gz"), h=T, sep="\t", check.names = F)
    nDMRs <- length(realFeatures[,1])
    realFeatures <- apply(realFeatures[,-c(1:3)], 2, sum)

    randomFeatures <- matrix(nrow=100, ncol=15)
    colnames(randomFeatures) = names(realFeatures)
    for (i in 1:100){
      randomFeatures[i,] <- apply(read.table(paste0("../annotation/DMRs/tissues/randomized/", baseName, "_randomDMRs_", i, "_15_features.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
    }

    ## Reordering features based on fold enrichment / depletion
    ordered <- order(realFeatures / apply(randomFeatures, 2, median), decreasing=T)
    randomFeatures <- randomFeatures[, ordered]
    realFeatures <- realFeatures[ordered]

    ## Plotting boxplots:
    pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))

    pdf(file = paste0("../annotation/DMRs/tissues/", baseName, "_DMRs_boxplot_15_features.pdf"), width = 6, height = 5)
    par(mar=c(7, 5, 4, 2) + 0.1, mgp=c(4, 1, 0)) 
    ## let 7% room above for fold-change indication
    ylimit <- c(0, max(realFeatures, randomFeatures)*1.08)
    ## draw empty plot to be able to add backrgoudn grey lines below boxes
    plot(1:15, type="n", axes=F, xlab="", ylab="", ylim=ylimit, yaxs="i", xlim=c(1,15))
    ## add vertical gray dashed lines to help reading
    abline(v=c(1:15), col="gray90", lty=3)
    ## add boxplot
    mp <- boxplot(randomFeatures, pch=20, cex=0.8, notch=F, ylab="Number of overlapping DMRs", xlab="", xaxt="n", ylim=ylimit, yaxs="i", whisklty=1, medlty=0, col="dimgray", border="dimgray", add=T, las=2)
    ## text 6% below x-axis / substitute _ for spaces in colnames
    axis(1, labels=NA, at=c(1:15))
    featNames <- names(realFeatures)
    featNames <- gsub("_UTR", "'_UTR", featNames, perl=TRUE)
    featNames <- gsub("coding", "(coding)", featNames, perl=TRUE)
    featNames <- gsub("proximal", "(proximal)", featNames, perl=TRUE)
    featNames <- gsub("_", " ", featNames, perl=TRUE)
    text(c(1:15), par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.06, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)
    ## plot the observed overlap for DMRs (diamond)
    ## red for enrichment, blue for depletion compared to median of randomized data)
    points(realFeatures, col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), pch=18, cex=1.3)
    text(c(1:15), par("usr")[4]*0.97, labels=round(realFeatures / apply(randomFeatures, 2, median),1), col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), cex=0.8)
    ## header: number of DMRs and name of contrast
    ## tissue contrast string treatment
    if (pair[2] == "Specific"){
      title(paste0(pair[1], "-specific in ", species, ": ", nDMRs, " DMRs"))
    }
    else {
      title(paste0(pair[1], " vs. ", pair[2], " in ", species, ": ", nDMRs, " DMRs"))
    }
    dev.off()
  }
}

## Individual DMRs
for (species in c("Human", "Chimp", "Rhesus")){
  for (pair in list(c("1", "2"), c("1", "3"), c("1", "4"), c("2", "3"), c("2", "4"), c("3", "4"))){
    cat("Testing individual pair: ", pair, "in", species, "\n")
    baseName <- paste0(species, "_", pair[1], "_", pair[2])

    realFeatures <- read.table(paste0("../annotation/DMRs/individuals/", baseName, "_DMRs_15_features.gz"), h=T, sep="\t", check.names = F)
    nDMRs <- length(realFeatures[,1])
    realFeatures <- apply(realFeatures[,-c(1:3)], 2, sum)

    randomFeatures <- matrix(nrow=100, ncol=15)
    colnames(randomFeatures) = names(realFeatures)
    for (i in 1:100){
      randomFeatures[i,] <- apply(read.table(paste0("../annotation/DMRs/individuals/randomized/", baseName, "_randomDMRs_", i, "_15_features.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
    }

    ## Reordering features based on fold enrichment / depletion
    ordered <- order(realFeatures / apply(randomFeatures, 2, median), decreasing=T)
    randomFeatures <- randomFeatures[, ordered]
    realFeatures <- realFeatures[ordered]

    ## Plotting boxplots:
    pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))

    pdf(file = paste0("../annotation/DMRs/individuals/", baseName, "_DMRs_boxplot_15_features.pdf"), width = 6, height = 5)
    par(mar=c(7, 5, 4, 2) + 0.1, mgp=c(4, 1, 0)) 
    ## let 7% room above for fold-change indication
    ylimit <- c(0, max(realFeatures, randomFeatures)*1.08)
    ## draw empty plot to be able to add backrgoudn grey lines below boxes
    plot(1:15, type="n", axes=F, xlab="", ylab="", ylim=ylimit, yaxs="i", xlim=c(1,15))
    ## add vertical gray dashed lines to help reading
    abline(v=c(1:15), col="gray90", lty=3)
    ## add boxplot
    mp <- boxplot(randomFeatures, pch=20, cex=0.8, notch=F, ylab="Number of overlapping DMRs", xlab="", xaxt="n", ylim=ylimit, yaxs="i", whisklty=1, medlty=0, col="dimgray", border="dimgray", add=T, las=2)
    ## text 6% below x-axis / substitute _ for spaces in colnames
    axis(1, labels=NA, at=c(1:15))
    featNames <- names(realFeatures)
    featNames <- gsub("_UTR", "'_UTR", featNames, perl=TRUE)
    featNames <- gsub("coding", "(coding)", featNames, perl=TRUE)
    featNames <- gsub("proximal", "(proximal)", featNames, perl=TRUE)
    featNames <- gsub("_", " ", featNames, perl=TRUE)
    text(c(1:15), par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.06, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)
    ## plot the observed overlap for DMRs (diamond)
    ## red for enrichment, blue for depletion compared to median of randomized data)
    points(realFeatures, col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), pch=18, cex=1.3)
    text(c(1:15), par("usr")[4]*0.97, labels=round(realFeatures / apply(randomFeatures, 2, median),1), col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), cex=0.8)
    ## header: number of DMRs and name of contrast
    ## individual contrast string treatment
    title(paste0("Individual ", pair[1], " vs. ", pair[2], " in ", species, ": ", nDMRs, " DMRs"))
    dev.off()
  }
}


## overlap of at least 20% of DMR length / random DMRs of same length and same CpG density #####################

library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))


## Species DMRs
for (tissue in c("heart", "kidney", "liver", "lung")){
  for (pair in list(c("Human", "Chimp"), c("Human", "Rhesus"), c("Chimp", "Rhesus"), c("Human", "Specific"), c("Chimp", "Specific"))){
    cat("Testing species pair: ", pair, "in", tissue, "\n")
    fileName <- paste0("../annotation/DMRs/species/", pair[1], pair[2], "_", tissue, "_DMRs_15_features_0.2.gz")
    realFeatures <- read.table(fileName, h=T, sep="\t", check.names = F)
    nDMRs <- length(realFeatures[,1])
    realFeatures <- apply(realFeatures[,-c(1:3)], 2, sum)

    ## control for length only
    randomFeatures <- matrix(nrow=100, ncol=15)
    colnames(randomFeatures) = names(realFeatures)
    for (i in 1:100){
      randomFeatures[i,] <- apply(read.table(paste0("../annotation/DMRs/species/randomized/", pair[1], pair[2], "_", tissue, "_randomDMRs_", i, "_15_features_0.2.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
    }

    ## control for length and CpG density
    randomFeaturesControlDensity <- matrix(nrow=100, ncol=15)
    colnames(randomFeaturesControlDensity) = names(realFeatures)
    for (i in 1:100){
      randomFeaturesControlDensity[i,] <- apply(read.table(paste0("../annotation/DMRs/species/randomized_control_CpG_density/", pair[1], pair[2], "_", tissue, "_randomDMRs_", i, "_15_features_0.2.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
    }

    ## Common plot for control for length and comtrol for length+density (2 pages):
    pdf(file = paste0("../annotation/DMRs/species/", pair[1], pair[2], "_", tissue, "_DMRs_boxplot_15_features_0.2.pdf"), width = 6, height = 5)

    ## Reordering features based on fold enrichment compared to randomDMRs of same length
    ordered <- order(realFeatures / apply(randomFeatures, 2, median), decreasing=T)
    randomFeaturesControlDensity <- randomFeaturesControlDensity[, ordered]
    randomFeatures <- randomFeatures[, ordered]
    realFeatures <- realFeatures[ordered]

    ## first boxplot
    par(mar=c(7, 5, 4, 2) + 0.1, mgp=c(4, 1, 0)) 
    ## let 7% room above for fold-change indication
    ylimit <- c(0, max(realFeatures, randomFeatures)*1.08)
    ## draw empty plot to be able to add background grey lines below boxes
    plot(1:15, type="n", axes=F, xlab="", ylab="", ylim=ylimit, yaxs="i", xlim=c(1,15))
    ## add vertical gray dashed lines to help reading
    abline(v=c(1:15), col="gray90", lty=3)
    ## add boxplot
    mp <- boxplot(randomFeatures, pch=20, cex=0.8, notch=F, ylab="Number of overlapping DMRs", xlab="", xaxt="n", ylim=ylimit, yaxs="i", whisklty=1, medlty=0, col="dimgray", border="dimgray", add=T, las=2)
    ## text 6% below x-axis / substitute _ for spaces in colnames
    axis(1, labels=NA, at=c(1:15))
    featNames <- names(realFeatures)
    featNames <- gsub("_UTR", "'_UTR", featNames, perl=TRUE)
    featNames <- gsub("coding", "(coding)", featNames, perl=TRUE)
    featNames <- gsub("proximal", "(proximal)", featNames, perl=TRUE)
    featNames <- gsub("_", " ", featNames, perl=TRUE)
    text(c(1:15), par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.06, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)
    ## plot the observed overlap for DMRs (diamond)
    ## red for enrichment, blue for depletion compared to median of randomized data)
    points(realFeatures, col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), pch=18, cex=1.3)
    text(c(1:15), par("usr")[4]*0.97, labels=round(realFeatures / apply(randomFeatures, 2, median),1), col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), cex=0.8)
    ## header: number of DMRs and name of contrast
    ## species contrast string treatment
    plotTitle <- vector()
    if (pair[2] == "Specific"){
      plotTitle <- paste0(pair[1], "-specific in ", tissue, ": ", nDMRs, " DMRs")
    } else {
      plotTitle <- paste0(pair[1], " vs. ", pair[2], " in ", tissue, ": ", nDMRs, " DMRs")
    }
    title(plotTitle)

    
    ## 2nd boxplot with control for CpG density
    ## Reorder based on fold enrichment
    ordered <- order(realFeatures / apply(randomFeaturesControlDensity, 2, median), decreasing=T)
    randomFeaturesControlDensity <- randomFeaturesControlDensity[, ordered]
    randomFeatures <- randomFeatures[, ordered]
    realFeatures <- realFeatures[ordered]
    featNames <- featNames[ordered]
    
    par(mar=c(7, 5, 4, 2) + 0.1, mgp=c(4, 1, 0)) 
    ## let 7% room above for fold-change indication
    ylimit <- c(0, max(realFeatures, randomFeaturesControlDensity)*1.08)
    ## draw empty plot to be able to add backgroumd grey lines below boxes
    plot(1:15, type="n", axes=F, xlab="", ylab="", ylim=ylimit, yaxs="i", xlim=c(1,15))
    ## add vertical gray dashed lines to help reading
    abline(v=c(1:15), col="gray90", lty=3)
    ## add boxplot
    mp <- boxplot(randomFeaturesControlDensity, pch=20, cex=0.8, notch=F, ylab="Number of overlapping DMRs", xlab="", xaxt="n", ylim=ylimit, yaxs="i", whisklty=1, medlty=0, col="dimgray", border="dimgray", add=T, las=2)
    ## text 6% below x-axis / substitute _ for spaces in colnames
    axis(1, labels=NA, at=c(1:15))
    text(c(1:15), par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.06, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)
    ## plot the observed overlap for DMRs (diamond)
    ## red for enrichment, blue for depletion compared to median of randomized data)
    points(realFeatures, col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), pch=18, cex=1.3)
    text(c(1:15), par("usr")[4]*0.97, labels=round(realFeatures / apply(randomFeaturesControlDensity, 2, median),1), col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), cex=0.8)
    ## header: number of DMRs and name of contrast
    title(plotTitle)
    dev.off()
  }
}

## Tissue DMRs
for (species in c("Human", "Chimp", "Rhesus")){
  for (pair in list(c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "liver"), c("heart", "Specific"), c("liver", "Specific"), c("kidney", "Specific"), c("lung", "Specific"))){

    cat("Testing tissue pair: ", pair, "in", species, "\n")
    if (pair[2] != "Specific"){
      baseName <- paste0(species, "_", pair[1], "_", pair[2])
    }
    else {
      baseName <- paste0(species, "_", pair[1], "Specific")
    }
    fileName <- paste0("../annotation/DMRs/tissues/", baseName, "_DMRs_15_features_0.2.gz")
    realFeatures <- read.table(fileName, h=T, sep="\t", check.names = F)
    nDMRs <- length(realFeatures[,1])
    realFeatures <- apply(realFeatures[,-c(1:3)], 2, sum)

    ## control for length only
    randomFeatures <- matrix(nrow=100, ncol=15)
    colnames(randomFeatures) = names(realFeatures)
    for (i in 1:100){
      randomFeatures[i,] <- apply(read.table(paste0("../annotation/DMRs/tissues/randomized/", baseName, "_randomDMRs_", i, "_15_features_0.2.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
    }

    ## control for length and CpG density
    randomFeaturesControlDensity <- matrix(nrow=100, ncol=15)
    colnames(randomFeaturesControlDensity) = names(realFeatures)
    for (i in 1:100){
      randomFeaturesControlDensity[i,] <- apply(read.table(paste0("../annotation/DMRs/tissues/randomized_control_CpG_density/", baseName, "_randomDMRs_", i, "_15_features_0.2.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
    }

    ## Common plot for control for length and comtrol for length+density (2 pages):
    pdf(file = paste0("../annotation/DMRs/tissues/", baseName, "_DMRs_boxplot_15_features_0.2.pdf"), width = 6, height = 5)

    ## Reordering features based on fold enrichment compared to randomDMRs of same length
    ordered <- order(realFeatures / apply(randomFeatures, 2, median), decreasing=T)
    randomFeaturesControlDensity <- randomFeaturesControlDensity[, ordered]
    randomFeatures <- randomFeatures[, ordered]
    realFeatures <- realFeatures[ordered]

    ## first boxplot
    par(mar=c(7, 5, 4, 2) + 0.1, mgp=c(4, 1, 0)) 
    ## let 7% room above for fold-change indication
    ylimit <- c(0, max(realFeatures, randomFeatures)*1.08)
    ## draw empty plot to be able to add background grey lines below boxes
    plot(1:15, type="n", axes=F, xlab="", ylab="", ylim=ylimit, yaxs="i", xlim=c(1,15))
    ## add vertical gray dashed lines to help reading
    abline(v=c(1:15), col="gray90", lty=3)
    ## add boxplot
    mp <- boxplot(randomFeatures, pch=20, cex=0.8, notch=F, ylab="Number of overlapping DMRs", xlab="", xaxt="n", ylim=ylimit, yaxs="i", whisklty=1, medlty=0, col="dimgray", border="dimgray", add=T, las=2)
    ## text 6% below x-axis / substitute _ for spaces in colnames
    axis(1, labels=NA, at=c(1:15))
    featNames <- names(realFeatures)
    featNames <- gsub("_UTR", "'_UTR", featNames, perl=TRUE)
    featNames <- gsub("coding", "(coding)", featNames, perl=TRUE)
    featNames <- gsub("proximal", "(proximal)", featNames, perl=TRUE)
    featNames <- gsub("_", " ", featNames, perl=TRUE)
    text(c(1:15), par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.06, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)
    ## plot the observed overlap for DMRs (diamond)
    ## red for enrichment, blue for depletion compared to median of randomized data)
    points(realFeatures, col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), pch=18, cex=1.3)
    text(c(1:15), par("usr")[4]*0.97, labels=round(realFeatures / apply(randomFeatures, 2, median),1), col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), cex=0.8)
    ## header: number of DMRs and name of contrast
    ## species contrast string treatment
    plotTitle <- vector()
    if (pair[2] == "Specific"){
      plotTitle <- paste0(pair[1], "-specific in ", species, ": ", nDMRs, " DMRs")
    } else {
      plotTitle <- paste0(pair[1], " vs. ", pair[2], " in ", species, ": ", nDMRs, " DMRs")
    }
    title(plotTitle)

    
    ## 2nd boxplot with control for CpG density
    ## Reorder based on fold enrichment
    ordered <- order(realFeatures / apply(randomFeaturesControlDensity, 2, median), decreasing=T)
    randomFeaturesControlDensity <- randomFeaturesControlDensity[, ordered]
    randomFeatures <- randomFeatures[, ordered]
    realFeatures <- realFeatures[ordered]
    featNames <- featNames[ordered]
    
    par(mar=c(7, 5, 4, 2) + 0.1, mgp=c(4, 1, 0)) 
    ## let 7% room above for fold-change indication
    ylimit <- c(0, max(realFeatures, randomFeaturesControlDensity)*1.08)
    ## draw empty plot to be able to add backgroumd grey lines below boxes
    plot(1:15, type="n", axes=F, xlab="", ylab="", ylim=ylimit, yaxs="i", xlim=c(1,15))
    ## add vertical gray dashed lines to help reading
    abline(v=c(1:15), col="gray90", lty=3)
    ## add boxplot
    mp <- boxplot(randomFeaturesControlDensity, pch=20, cex=0.8, notch=F, ylab="Number of overlapping DMRs", xlab="", xaxt="n", ylim=ylimit, yaxs="i", whisklty=1, medlty=0, col="dimgray", border="dimgray", add=T, las=2)
    ## text 6% below x-axis / substitute _ for spaces in colnames
    axis(1, labels=NA, at=c(1:15))
    text(c(1:15), par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.06, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)
    ## plot the observed overlap for DMRs (diamond)
    ## red for enrichment, blue for depletion compared to median of randomized data)
    points(realFeatures, col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), pch=18, cex=1.3)
    text(c(1:15), par("usr")[4]*0.97, labels=round(realFeatures / apply(randomFeaturesControlDensity, 2, median),1), col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), cex=0.8)
    ## header: number of DMRs and name of contrast
    title(plotTitle)
    dev.off()
  }
}

## individual DMRs
for (species in c("Human", "Chimp", "Rhesus")){
  for (pair in list(c("1", "2"), c("1", "3"), c("1", "4"), c("2", "3"), c("2", "4"), c("3", "4"))){
    cat("Testing individual pair: ", pair, "in", species, "\n")
    fileName <- paste0("../annotation/DMRs/individuals/", species, "_", pair[1], "_", pair[2], "_DMRs_15_features_0.2.gz")
    realFeatures <- read.table(fileName, h=T, sep="\t", check.names = F)
    nDMRs <- length(realFeatures[,1])
    realFeatures <- apply(realFeatures[,-c(1:3)], 2, sum)

    ## control for length only
    randomFeatures <- matrix(nrow=100, ncol=15)
    colnames(randomFeatures) = names(realFeatures)
    for (i in 1:100){
      randomFeatures[i,] <- apply(read.table(paste0("../annotation/DMRs/individuals/randomized/", species, "_", pair[1], "_", pair[2], "_randomDMRs_", i, "_15_features_0.2.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
    }

    ## control for length and CpG density
    randomFeaturesControlDensity <- matrix(nrow=100, ncol=15)
    colnames(randomFeaturesControlDensity) = names(realFeatures)
    for (i in 1:100){
      randomFeaturesControlDensity[i,] <- apply(read.table(paste0("../annotation/DMRs/individuals/randomized_control_CpG_density/", species, "_", pair[1], "_", pair[2], "_randomDMRs_", i, "_15_features_0.2.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
    }

    ## Common plot for control for length and comtrol for length+density (2 pages):
    pdf(file = paste0("../annotation/DMRs/individuals/", species, "_", pair[1], "_", pair[2], "_DMRs_boxplot_15_features_0.2.pdf"), width = 6, height = 5)

    ## Reordering features based on fold enrichment compared to randomDMRs of same length
    ordered <- order(realFeatures / apply(randomFeatures, 2, median), decreasing=T)
    randomFeaturesControlDensity <- randomFeaturesControlDensity[, ordered]
    randomFeatures <- randomFeatures[, ordered]
    realFeatures <- realFeatures[ordered]

    ## first boxplot
    par(mar=c(7, 5, 4, 2) + 0.1, mgp=c(4, 1, 0)) 
    ## let 7% room above for fold-change indication
    ylimit <- c(0, max(realFeatures, randomFeatures)*1.08)
    ## draw empty plot to be able to add background grey lines below boxes
    plot(1:15, type="n", axes=F, xlab="", ylab="", ylim=ylimit, yaxs="i", xlim=c(1,15))
    ## add vertical gray dashed lines to help reading
    abline(v=c(1:15), col="gray90", lty=3)
    ## add boxplot
    mp <- boxplot(randomFeatures, pch=20, cex=0.8, notch=F, ylab="Number of overlapping DMRs", xlab="", xaxt="n", ylim=ylimit, yaxs="i", whisklty=1, medlty=0, col="dimgray", border="dimgray", add=T, las=2)
    ## text 6% below x-axis / substitute _ for spaces in colnames
    axis(1, labels=NA, at=c(1:15))
    featNames <- names(realFeatures)
    featNames <- gsub("_UTR", "'_UTR", featNames, perl=TRUE)
    featNames <- gsub("coding", "(coding)", featNames, perl=TRUE)
    featNames <- gsub("proximal", "(proximal)", featNames, perl=TRUE)
    featNames <- gsub("_", " ", featNames, perl=TRUE)
    text(c(1:15), par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.06, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)
    ## plot the observed overlap for DMRs (diamond)
    ## red for enrichment, blue for depletion compared to median of randomized data)
    points(realFeatures, col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), pch=18, cex=1.3)
    text(c(1:15), par("usr")[4]*0.97, labels=round(realFeatures / apply(randomFeatures, 2, median),1), col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), cex=0.8)
    ## header: number of DMRs and name of contrast
    ## species contrast string treatment
    plotTitle <- paste0(pair[1], " vs. ", pair[2], " in ", species, ": ", nDMRs, " DMRs")
    title(plotTitle)

    ## 2nd boxplot with control for CpG density
    ## Reorder based on fold enrichment
    ordered <- order(realFeatures / apply(randomFeaturesControlDensity, 2, median), decreasing=T)
    randomFeaturesControlDensity <- randomFeaturesControlDensity[, ordered]
    randomFeatures <- randomFeatures[, ordered]
    realFeatures <- realFeatures[ordered]
    featNames <- featNames[ordered]
    
    par(mar=c(7, 5, 4, 2) + 0.1, mgp=c(4, 1, 0)) 
    ## let 7% room above for fold-change indication
    ylimit <- c(0, max(realFeatures, randomFeaturesControlDensity)*1.08)
    ## draw empty plot to be able to add backgroumd grey lines below boxes
    plot(1:15, type="n", axes=F, xlab="", ylab="", ylim=ylimit, yaxs="i", xlim=c(1,15))
    ## add vertical gray dashed lines to help reading
    abline(v=c(1:15), col="gray90", lty=3)
    ## add boxplot
    mp <- boxplot(randomFeaturesControlDensity, pch=20, cex=0.8, notch=F, ylab="Number of overlapping DMRs", xlab="", xaxt="n", ylim=ylimit, yaxs="i", whisklty=1, medlty=0, col="dimgray", border="dimgray", add=T, las=2)
    ## text 6% below x-axis / substitute _ for spaces in colnames
    axis(1, labels=NA, at=c(1:15))
    text(c(1:15), par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.06, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)
    ## plot the observed overlap for DMRs (diamond)
    ## red for enrichment, blue for depletion compared to median of randomized data)
    points(realFeatures, col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), pch=18, cex=1.3)
    text(c(1:15), par("usr")[4]*0.97, labels=round(realFeatures / apply(randomFeaturesControlDensity, 2, median),1), col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), cex=0.8)
    ## header: number of DMRs and name of contrast
    title(plotTitle)
    dev.off()
  }
}

## overlap of at least 20% of DMR length / random DMRs of same length and same CpG density #####################
## Hyper and hypo methylated DMRs separated. Only tissue-specific and species-specific comparisons considered

library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))

## Species DMRs
for (tissue in c("heart", "kidney", "liver", "lung")){
  for (pair in list(c("Human", "Specific"), c("Chimp", "Specific"))){

    ## Common plot for control for length and comtrol for length+density and 2 directions (4 pages):
    pdf(file = paste0("../annotation/DMRs/species/", pair[1], pair[2], "_", tissue, "_polarized_DMRs_boxplot_15_features_0.2.pdf"), width = 6, height = 5)


    for (direction in c("hyper", "hypo")){
      cat("Testing ", pair, direction, "methylated DMRs in", tissue, "\n")
      fileName <- paste0("../annotation/DMRs/species/", pair[1], pair[2], "_", tissue, "_", direction, "_DMRs_15_features_0.2.gz")
      realFeatures <- read.table(fileName, h=T, sep="\t", check.names = F)
      nDMRs <- length(realFeatures[,1])
      realFeatures <- apply(realFeatures[,-c(1:3)], 2, sum)

      ## control for length only
      randomFeatures <- matrix(nrow=100, ncol=15)
      colnames(randomFeatures) = names(realFeatures)
      for (i in 1:100){
        randomFeatures[i,] <- apply(read.table(paste0("../annotation/DMRs/species/randomized/", pair[1], pair[2], "_", tissue, "_", direction, "_randomDMRs_", i, "_15_features_0.2.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
      }

      ## control for length and CpG density
      randomFeaturesControlDensity <- matrix(nrow=100, ncol=15)
      colnames(randomFeaturesControlDensity) = names(realFeatures)
      for (i in 1:100){
        randomFeaturesControlDensity[i,] <- apply(read.table(paste0("../annotation/DMRs/species/randomized_control_CpG_density/", pair[1], pair[2], "_", tissue, "_", direction, "_randomDMRs_", i, "_15_features_0.2.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
      }

      ## Reordering features based on fold enrichment compared to randomDMRs of same length
      ordered <- order(realFeatures / apply(randomFeatures, 2, median), decreasing=T)
      randomFeaturesControlDensity <- randomFeaturesControlDensity[, ordered]
      randomFeatures <- randomFeatures[, ordered]
      realFeatures <- realFeatures[ordered]

      ## first boxplot
      par(mar=c(7, 5, 4, 2) + 0.1, mgp=c(4, 1, 0)) 
      ## let 7% room above for fold-change indication
      ylimit <- c(0, max(realFeatures, randomFeatures)*1.08)
      ## draw empty plot to be able to add background grey lines below boxes
      plot(1:15, type="n", axes=F, xlab="", ylab="", ylim=ylimit, yaxs="i", xlim=c(1,15))
      ## add vertical gray dashed lines to help reading
      abline(v=c(1:15), col="gray90", lty=3)
      ## add boxplot
      mp <- boxplot(randomFeatures, pch=20, cex=0.8, notch=F, ylab="Number of overlapping DMRs", xlab="", xaxt="n", ylim=ylimit, yaxs="i", whisklty=1, medlty=0, col="dimgray", border="dimgray", add=T, las=2)
      ## text 6% below x-axis / substitute _ for spaces in colnames
      axis(1, labels=NA, at=c(1:15))
      featNames <- names(realFeatures)
      featNames <- gsub("_UTR", "'_UTR", featNames, perl=TRUE)
      featNames <- gsub("coding", "(coding)", featNames, perl=TRUE)
      featNames <- gsub("proximal", "(proximal)", featNames, perl=TRUE)
      featNames <- gsub("_", " ", featNames, perl=TRUE)
      text(c(1:15), par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.06, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)
      ## plot the observed overlap for DMRs (diamond)
      ## red for enrichment, blue for depletion compared to median of randomized data)
      points(realFeatures, col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), pch=18, cex=1.3)
      text(c(1:15), par("usr")[4]*0.97, labels=round(realFeatures / apply(randomFeatures, 2, median),1), col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), cex=0.8)
      ## header: number of DMRs and name of contrast
      ## species contrast string treatment
      plotTitle <- vector()
      plotTitle <- paste0(pair[1], "-specific ", direction, "-methylated DMRs in ", tissue, ": ", nDMRs, " DMRs")
      title(plotTitle)

      
      ## 2nd boxplot with control for CpG density
      ## Reorder based on fold enrichment
      ordered <- order(realFeatures / apply(randomFeaturesControlDensity, 2, median), decreasing=T)
      randomFeaturesControlDensity <- randomFeaturesControlDensity[, ordered]
      randomFeatures <- randomFeatures[, ordered]
      realFeatures <- realFeatures[ordered]
      featNames <- featNames[ordered]
      
      par(mar=c(7, 5, 4, 2) + 0.1, mgp=c(4, 1, 0)) 
      ## let 7% room above for fold-change indication
      ylimit <- c(0, max(realFeatures, randomFeaturesControlDensity)*1.08)
      ## draw empty plot to be able to add backgroumd grey lines below boxes
      plot(1:15, type="n", axes=F, xlab="", ylab="", ylim=ylimit, yaxs="i", xlim=c(1,15))
      ## add vertical gray dashed lines to help reading
      abline(v=c(1:15), col="gray90", lty=3)
      ## add boxplot
      mp <- boxplot(randomFeaturesControlDensity, pch=20, cex=0.8, notch=F, ylab="Number of overlapping DMRs", xlab="", xaxt="n", ylim=ylimit, yaxs="i", whisklty=1, medlty=0, col="dimgray", border="dimgray", add=T, las=2)
      ## text 6% below x-axis / substitute _ for spaces in colnames
      axis(1, labels=NA, at=c(1:15))
      text(c(1:15), par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.06, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)
      ## plot the observed overlap for DMRs (diamond)
      ## red for enrichment, blue for depletion compared to median of randomized data)
      points(realFeatures, col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), pch=18, cex=1.3)
      text(c(1:15), par("usr")[4]*0.97, labels=round(realFeatures / apply(randomFeaturesControlDensity, 2, median),1), col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), cex=0.8)
      ## header: number of DMRs and name of contrast
      title(plotTitle)
    }
    dev.off()
  }
}

## Tissue DMRs
for (species in c("Human", "Chimp", "Rhesus")){
  for (pair in list(c("heart", "Specific"), c("liver", "Specific"), c("kidney", "Specific"), c("lung", "Specific"))){
    ## Common plot for control for length and control for length+density and 2 directions (4 pages):
    pdf(file = paste0("../annotation/DMRs/tissues/", species, "_", pair[1], pair[2], "_polarized_DMRs_boxplot_15_features_0.2.pdf"), width = 6, height = 5)

    for (direction in c("hyper", "hypo")){
      cat("Testing", pair, direction, "methylated DMRs in", species, "\n")
      baseName <- paste0(species, "_", pair[1], "Specific")
      fileName <- paste0("../annotation/DMRs/tissues/", baseName, "_", direction, "_DMRs_15_features_0.2.gz")
      realFeatures <- read.table(fileName, h=T, sep="\t", check.names = F)
      nDMRs <- length(realFeatures[,1])
      realFeatures <- apply(realFeatures[,-c(1:3)], 2, sum)

      ## control for length only
      randomFeatures <- matrix(nrow=100, ncol=15)
      colnames(randomFeatures) = names(realFeatures)
      for (i in 1:100){
        randomFeatures[i,] <- apply(read.table(paste0("../annotation/DMRs/tissues/randomized/", baseName, "_", direction, "_randomDMRs_", i, "_15_features_0.2.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
      }

      ## control for length and CpG density
      randomFeaturesControlDensity <- matrix(nrow=100, ncol=15)
      colnames(randomFeaturesControlDensity) = names(realFeatures)
      for (i in 1:100){
        randomFeaturesControlDensity[i,] <- apply(read.table(paste0("../annotation/DMRs/tissues/randomized_control_CpG_density/", baseName, "_", direction, "_randomDMRs_", i, "_15_features_0.2.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
      }

      ## Reordering features based on fold enrichment compared to randomDMRs of same length
      ordered <- order(realFeatures / apply(randomFeatures, 2, median), decreasing=T)
      randomFeaturesControlDensity <- randomFeaturesControlDensity[, ordered]
      randomFeatures <- randomFeatures[, ordered]
      realFeatures <- realFeatures[ordered]

      ## first boxplot
      par(mar=c(7, 5, 4, 2) + 0.1, mgp=c(4, 1, 0)) 
      ## let 7% room above for fold-change indication
      ylimit <- c(0, max(realFeatures, randomFeatures)*1.08)
      ## draw empty plot to be able to add background grey lines below boxes
      plot(1:15, type="n", axes=F, xlab="", ylab="", ylim=ylimit, yaxs="i", xlim=c(1,15))
      ## add vertical gray dashed lines to help reading
      abline(v=c(1:15), col="gray90", lty=3)
      ## add boxplot
      mp <- boxplot(randomFeatures, pch=20, cex=0.8, notch=F, ylab="Number of overlapping DMRs", xlab="", xaxt="n", ylim=ylimit, yaxs="i", whisklty=1, medlty=0, col="dimgray", border="dimgray", add=T, las=2)
      ## text 6% below x-axis / substitute _ for spaces in colnames
      axis(1, labels=NA, at=c(1:15))
      featNames <- names(realFeatures)
      featNames <- gsub("_UTR", "'_UTR", featNames, perl=TRUE)
      featNames <- gsub("coding", "(coding)", featNames, perl=TRUE)
      featNames <- gsub("proximal", "(proximal)", featNames, perl=TRUE)
      featNames <- gsub("_", " ", featNames, perl=TRUE)
      text(c(1:15), par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.06, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)
      ## plot the observed overlap for DMRs (diamond)
      ## red for enrichment, blue for depletion compared to median of randomized data)
      points(realFeatures, col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), pch=18, cex=1.3)
      text(c(1:15), par("usr")[4]*0.97, labels=round(realFeatures / apply(randomFeatures, 2, median),1), col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), cex=0.8)
      ## header: number of DMRs and name of contrast
      ## species contrast string treatment
      plotTitle <- paste0(pair[1], "-specific ", direction, "-methylated DMRs in ", species, ": ", nDMRs, " DMRs")
      title(plotTitle)

      
      ## 2nd boxplot with control for CpG density
      ## Reorder based on fold enrichment
      ordered <- order(realFeatures / apply(randomFeaturesControlDensity, 2, median), decreasing=T)
      randomFeaturesControlDensity <- randomFeaturesControlDensity[, ordered]
      randomFeatures <- randomFeatures[, ordered]
      realFeatures <- realFeatures[ordered]
      featNames <- featNames[ordered]
      
      par(mar=c(7, 5, 4, 2) + 0.1, mgp=c(4, 1, 0)) 
      ## let 7% room above for fold-change indication
      ylimit <- c(0, max(realFeatures, randomFeaturesControlDensity)*1.08)
      ## draw empty plot to be able to add backgroumd grey lines below boxes
      plot(1:15, type="n", axes=F, xlab="", ylab="", ylim=ylimit, yaxs="i", xlim=c(1,15))
      ## add vertical gray dashed lines to help reading
      abline(v=c(1:15), col="gray90", lty=3)
      ## add boxplot
      mp <- boxplot(randomFeaturesControlDensity, pch=20, cex=0.8, notch=F, ylab="Number of overlapping DMRs", xlab="", xaxt="n", ylim=ylimit, yaxs="i", whisklty=1, medlty=0, col="dimgray", border="dimgray", add=T, las=2)
      ## text 6% below x-axis / substitute _ for spaces in colnames
      axis(1, labels=NA, at=c(1:15))
      text(c(1:15), par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.06, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)
      ## plot the observed overlap for DMRs (diamond)
      ## red for enrichment, blue for depletion compared to median of randomized data)
      points(realFeatures, col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), pch=18, cex=1.3)
      text(c(1:15), par("usr")[4]*0.97, labels=round(realFeatures / apply(randomFeaturesControlDensity, 2, median),1), col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), cex=0.8)
      ## header: number of DMRs and name of contrast
      title(plotTitle)
    }
    dev.off()
  }
}


## Only plot barplot of % overlap of each feature (no comparison with randomDMRs) ###########################

## Species DMRs
for (tissue in c("heart", "kidney", "liver", "lung")){
  for (pair in list(c("Human", "Chimp"), c("Human", "Rhesus"), c("Chimp", "Rhesus"), c("Human", "Specific"), c("Chimp", "Specific"))){
    cat("Testing species pair: ", pair, "in", tissue, "\n")
    fileName <- paste0("../annotation/DMRs/species/", pair[1], pair[2], "_", tissue, "_DMRs_15_features_0.2.gz")
    realFeatures <- read.table(fileName, h=T, sep="\t", check.names = F)
    nDMRs <- length(realFeatures[,1])
    realFeatures <- apply(realFeatures[,-c(1:3)], 2, sum)

    ## Barplot
    pdf(file = paste0("../annotation/DMRs/species/", pair[1], pair[2], "_", tissue, "_DMRs_barplot_15_features_0.2.pdf"), width = 6, height = 5)

    ## Reordering features based on % overlap
    ordered <- order(realFeatures / nDMRs, decreasing=T)
    realFeatures <- realFeatures[ordered]

    par(mar=c(7, 5, 4, 2) + 0.1) 
    mp <- barplot(realFeatures/nDMRs, ylab="Percentage of DMRs", xlab="", col="darkgrey", xaxt="n")

    featNames <- names(realFeatures)
    featNames <- gsub("_UTR", "'_UTR", featNames, perl=TRUE)
    featNames <- gsub("coding", "(coding)", featNames, perl=TRUE)
    featNames <- gsub("proximal", "(proximal)", featNames, perl=TRUE)
    featNames <- gsub("_", " ", featNames, perl=TRUE)
    text(mp[,1], par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.03, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)

    ## header: number of DMRs and name of contrast
    ## species contrast string treatment
    plotTitle <- vector()
    if (pair[2] == "Specific"){
      plotTitle <- paste0(pair[1], "-specific in ", tissue, ": ", nDMRs, " DMRs")
    } else {
      plotTitle <- paste0(pair[1], " vs. ", pair[2], " in ", tissue, ": ", nDMRs, " DMRs")
    }
    title(plotTitle)
    dev.off()
  }
}


## Tissue DMRs
for (species in c("Human", "Chimp", "Rhesus")){
  for (pair in list(c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "liver"), c("heart", "Specific"), c("liver", "Specific"), c("kidney", "Specific"), c("lung", "Specific"))){

    cat("Testing tissue pair: ", pair, "in", species, "\n")
    if (pair[2] != "Specific"){
      baseName <- paste0(species, "_", pair[1], "_", pair[2])
    }
    else {
      baseName <- paste0(species, "_", pair[1], "Specific")
    }
    fileName <- paste0("../annotation/DMRs/tissues/", baseName, "_DMRs_15_features_0.2.gz")

    realFeatures <- read.table(fileName, h=T, sep="\t", check.names = F)
    nDMRs <- length(realFeatures[,1])
    realFeatures <- apply(realFeatures[,-c(1:3)], 2, sum)

    ## Barplot
    pdf(file = paste0("../annotation/DMRs/tissues/", baseName, "_DMRs_barplot_15_features_0.2.pdf"), width = 6, height = 5)

    ## Reordering features based on % overlap
    ordered <- order(realFeatures / nDMRs, decreasing=T)
    realFeatures <- realFeatures[ordered]

    par(mar=c(7, 5, 4, 2) + 0.1) 
    mp <- barplot(realFeatures/nDMRs, ylab="Percentage of DMRs", xlab="", col="darkgrey", xaxt="n")

    featNames <- names(realFeatures)
    featNames <- gsub("_UTR", "'_UTR", featNames, perl=TRUE)
    featNames <- gsub("coding", "(coding)", featNames, perl=TRUE)
    featNames <- gsub("proximal", "(proximal)", featNames, perl=TRUE)
    featNames <- gsub("_", " ", featNames, perl=TRUE)
    text(mp[,1], par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.03, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)

    ## header: number of DMRs and name of contrast
    ## species contrast string treatment
    plotTitle <- vector()
    if (pair[2] == "Specific"){
      plotTitle <- paste0(pair[1], "-specific in ", species, ": ", nDMRs, " DMRs")
    } else {
      plotTitle <- paste0(pair[1], " vs. ", pair[2], " in ", species, ": ", nDMRs, " DMRs")
    }
    title(plotTitle)
    dev.off()
  }
}


## Individual DMRs
for (species in c("Human", "Chimp", "Rhesus")){
  for (pair in list(c("1", "2"), c("1", "3"), c("1", "4"), c("2", "3"), c("2", "4"), c("3", "4"))){
    cat("Testing individual pair: ", pair, "in", species, "\n")
    fileName <- paste0("../annotation/DMRs/individuals/", species, "_", pair[1], "_", pair[2], "_DMRs_15_features_0.2.gz")

    realFeatures <- read.table(fileName, h=T, sep="\t", check.names = F)
    nDMRs <- length(realFeatures[,1])
    realFeatures <- apply(realFeatures[,-c(1:3)], 2, sum)

    ## Barplot
    pdf(file = paste0("../annotation/DMRs/individuals/", species, "_", pair[1], "_", pair[2], "_DMRs_barplot_15_features_0.2.pdf"), width = 6, height = 5)

    ## Reordering features based on % overlap
    ordered <- order(realFeatures / nDMRs, decreasing=T)
    realFeatures <- realFeatures[ordered]

    par(mar=c(7, 5, 4, 2) + 0.1) 
    mp <- barplot(realFeatures/nDMRs, ylab="Percentage of DMRs", xlab="", col="darkgrey", xaxt="n")

    featNames <- names(realFeatures)
    featNames <- gsub("_UTR", "'_UTR", featNames, perl=TRUE)
    featNames <- gsub("coding", "(coding)", featNames, perl=TRUE)
    featNames <- gsub("proximal", "(proximal)", featNames, perl=TRUE)
    featNames <- gsub("_", " ", featNames, perl=TRUE)
    text(mp[,1], par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.03, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)

    ## header: number of DMRs and name of contrast
    ## species contrast string treatment
    plotTitle <- paste0(pair[1], " vs. ", pair[2], " in ", species, ": ", nDMRs, " DMRs")
    title(plotTitle)
    dev.off()
  }
}


## Repeat analysis: different classes and families #######################################
## overlap of at least 1bp / random DMRs of same length / same length and same CpG density 
library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))

## Species DMRs
for (tissue in c("heart", "kidney", "liver", "lung")){
  for (pair in list(c("Human", "Chimp"), c("Human", "Rhesus"), c("Chimp", "Rhesus"), c("Human", "Specific"), c("Chimp", "Specific"))){
    cat("Testing species pair: ", pair, "in", tissue, "\n")
    fileName <- paste0("../annotation/DMRs/species/", pair[1], pair[2], "_", tissue, "_DMRs_repeats.gz")
    realFeatures <- read.table(fileName, h=T, sep="\t", check.names = F)
    nDMRs <- length(realFeatures[,1])
    realFeatures <- apply(realFeatures[,-c(1:3)], 2, sum)

    ## control for length only
    randomFeatures <- matrix(nrow=100, ncol=13)
    colnames(randomFeatures) = names(realFeatures)
    for (i in 1:100){
      ## print(i)
      randomFeatures[i,] <- apply(read.table(paste0("../annotation/DMRs/species/randomized/", pair[1], pair[2], "_", tissue, "_randomDMRs_", i, "_repeats.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
    }

    ## control for length and CpG density
    randomFeaturesControlDensity <- matrix(nrow=100, ncol=13)
    colnames(randomFeaturesControlDensity) = names(realFeatures)
    for (i in 1:100){
      randomFeaturesControlDensity[i,] <- apply(read.table(paste0("../annotation/DMRs/species/randomized_control_CpG_density/", pair[1], pair[2], "_", tissue, "_randomDMRs_", i, "_repeats.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
    }

    ## Common plot for control for length and comtrol for length+density (2 pages):
    pdf(file = paste0("../annotation/DMRs/species/", pair[1], pair[2], "_", tissue, "_DMRs_boxplot_repeats.pdf"), width = 6, height = 5)

    ## Reordering features based on fold enrichment compared to randomDMRs of same length
    ordered <- order(realFeatures / apply(randomFeatures, 2, median), decreasing=T)
    randomFeaturesControlDensity <- randomFeaturesControlDensity[, ordered]
    randomFeatures <- randomFeatures[, ordered]
    realFeatures <- realFeatures[ordered]

    ## first boxplot
    par(mar=c(7, 5, 4, 2) + 0.1, mgp=c(4, 1, 0)) 
    ## let 7% room above for fold-change indication
    ylimit <- c(0, max(realFeatures, randomFeatures)*1.08)
    ## draw empty plot to be able to add background grey lines below boxes
    plot(1:13, type="n", axes=F, xlab="", ylab="", ylim=ylimit, yaxs="i", xlim=c(1,13))
    ## add vertical gray dashed lines to help reading
    abline(v=c(1:13), col="gray90", lty=3)
    ## add boxplot
    mp <- boxplot(randomFeatures, pch=20, cex=0.8, notch=F, ylab="Number of overlapping DMRs", xlab="", xaxt="n", ylim=ylimit, yaxs="i", whisklty=1, medlty=0, col="dimgray", border="dimgray", add=T, las=2)
    ## text 6% below x-axis / substitute _ for spaces in colnames
    axis(1, labels=NA, at=c(1:13))
    featNames <- names(realFeatures)
    featNames <- gsub("_", " ", featNames, perl=TRUE)
    text(c(1:13), par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.06, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)
    ## plot the observed overlap for DMRs (diamond)
    ## red for enrichment, blue for depletion compared to median of randomized data)
    points(realFeatures, col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), pch=18, cex=1.3)
    text(c(1:13), par("usr")[4]*0.97, labels=round(realFeatures / apply(randomFeatures, 2, median), 2), col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), cex=0.8)
    ## header: number of DMRs and name of contrast
    ## species contrast string treatment
    plotTitle <- vector()
    if (pair[2] == "Specific"){
      plotTitle <- paste0(pair[1], "-specific in ", tissue, ": ", nDMRs, " DMRs")
    } else {
      plotTitle <- paste0(pair[1], " vs. ", pair[2], " in ", tissue, ": ", nDMRs, " DMRs")
    }
    title(plotTitle)

    
    ## 2nd boxplot with control for CpG density
    ## Reorder based on fold enrichment
    ordered <- order(realFeatures / apply(randomFeaturesControlDensity, 2, median), decreasing=T)
    randomFeaturesControlDensity <- randomFeaturesControlDensity[, ordered]
    randomFeatures <- randomFeatures[, ordered]
    realFeatures <- realFeatures[ordered]
    featNames <- featNames[ordered]
    
    par(mar=c(7, 5, 4, 2) + 0.1, mgp=c(4, 1, 0)) 
    ## let 7% room above for fold-change indication
    ylimit <- c(0, max(realFeatures, randomFeaturesControlDensity)*1.08)
    ## draw empty plot to be able to add backgroumd grey lines below boxes
    plot(1:13, type="n", axes=F, xlab="", ylab="", ylim=ylimit, yaxs="i", xlim=c(1,13))
    ## add vertical gray dashed lines to help reading
    abline(v=c(1:13), col="gray90", lty=3)
    ## add boxplot
    mp <- boxplot(randomFeaturesControlDensity, pch=20, cex=0.8, notch=F, ylab="Number of overlapping DMRs", xlab="", xaxt="n", ylim=ylimit, yaxs="i", whisklty=1, medlty=0, col="dimgray", border="dimgray", add=T, las=2)
    ## text 6% below x-axis / substitute _ for spaces in colnames
    axis(1, labels=NA, at=c(1:13))
    text(c(1:13), par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.06, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)
    ## plot the observed overlap for DMRs (diamond)
    ## red for enrichment, blue for depletion compared to median of randomized data)
    points(realFeatures, col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), pch=18, cex=1.3)
    text(c(1:13), par("usr")[4]*0.97, labels=round(realFeatures / apply(randomFeaturesControlDensity, 2, median), 2), col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), cex=0.8)
    ## header: number of DMRs and name of contrast
    title(plotTitle)

    dev.off()
  }
}

## Tissue DMRs
for (species in c("Human", "Chimp", "Rhesus")){
  for (pair in list(c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "liver"), c("heart", "Specific"), c("liver", "Specific"), c("kidney", "Specific"), c("lung", "Specific"))){

    cat("Testing tissue pair: ", pair, "in", species, "\n")
    if (pair[2] != "Specific"){
      baseName <- paste0(species, "_", pair[1], "_", pair[2])
    }
    else {
      baseName <- paste0(species, "_", pair[1], "Specific")
    }
    fileName <- paste0("../annotation/DMRs/tissues/", baseName, "_DMRs_repeats.gz")
    realFeatures <- read.table(fileName, h=T, sep="\t", check.names = F)
    nDMRs <- length(realFeatures[,1])
    realFeatures <- apply(realFeatures[,-c(1:3)], 2, sum)

    ## control for length only
    randomFeatures <- matrix(nrow=100, ncol=13)
    colnames(randomFeatures) = names(realFeatures)
    for (i in 1:100){
      randomFeatures[i,] <- apply(read.table(paste0("../annotation/DMRs/tissues/randomized/", baseName, "_randomDMRs_", i, "_repeats.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
    }

    ## control for length and CpG density
    randomFeaturesControlDensity <- matrix(nrow=100, ncol=13)
    colnames(randomFeaturesControlDensity) = names(realFeatures)
    for (i in 1:100){
      randomFeaturesControlDensity[i,] <- apply(read.table(paste0("../annotation/DMRs/tissues/randomized_control_CpG_density/", baseName, "_randomDMRs_", i, "_repeats.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
    }

    ## Common plot for control for length and comtrol for length+density (2 pages):
    pdf(file = paste0("../annotation/DMRs/tissues/", baseName, "_DMRs_boxplot_repeats.pdf"), width = 6, height = 5)

    ## Reordering features based on fold enrichment compared to randomDMRs of same length
    ordered <- order(realFeatures / apply(randomFeatures, 2, median), decreasing=T)
    randomFeaturesControlDensity <- randomFeaturesControlDensity[, ordered]
    randomFeatures <- randomFeatures[, ordered]
    realFeatures <- realFeatures[ordered]

    ## first boxplot
    par(mar=c(7, 5, 4, 2) + 0.1, mgp=c(4, 1, 0)) 
    ## let 7% room above for fold-change indication
    ylimit <- c(0, max(realFeatures, randomFeatures)*1.08)
    ## draw empty plot to be able to add background grey lines below boxes
    plot(1:13, type="n", axes=F, xlab="", ylab="", ylim=ylimit, yaxs="i", xlim=c(1,13))
    ## add vertical gray dashed lines to help reading
    abline(v=c(1:13), col="gray90", lty=3)
    ## add boxplot
    mp <- boxplot(randomFeatures, pch=20, cex=0.8, notch=F, ylab="Number of overlapping DMRs", xlab="", xaxt="n", ylim=ylimit, yaxs="i", whisklty=1, medlty=0, col="dimgray", border="dimgray", add=T, las=2)
    ## text 6% below x-axis / substitute _ for spaces in colnames
    axis(1, labels=NA, at=c(1:13))
    featNames <- names(realFeatures)
    featNames <- gsub("_", " ", featNames, perl=TRUE)
    text(c(1:13), par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.06, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)
    ## plot the observed overlap for DMRs (diamond)
    ## red for enrichment, blue for depletion compared to median of randomized data)
    points(realFeatures, col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), pch=18, cex=1.3)
    text(c(1:13), par("usr")[4]*0.97, labels=round(realFeatures / apply(randomFeatures, 2, median), 2), col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), cex=0.8)
    ## header: number of DMRs and name of contrast
    ## species contrast string treatment
    plotTitle <- vector()
    if (pair[2] == "Specific"){
      plotTitle <- paste0(pair[1], "-specific in ", species, ": ", nDMRs, " DMRs")
    } else {
      plotTitle <- paste0(pair[1], " vs. ", pair[2], " in ", species, ": ", nDMRs, " DMRs")
    }
    title(plotTitle)

    
    ## 2nd boxplot with control for CpG density
    ## Reorder based on fold enrichment
    ordered <- order(realFeatures / apply(randomFeaturesControlDensity, 2, median), decreasing=T)
    randomFeaturesControlDensity <- randomFeaturesControlDensity[, ordered]
    randomFeatures <- randomFeatures[, ordered]
    realFeatures <- realFeatures[ordered]
    featNames <- featNames[ordered]
    
    par(mar=c(7, 5, 4, 2) + 0.1, mgp=c(4, 1, 0)) 
    ## let 7% room above for fold-change indication
    ylimit <- c(0, max(realFeatures, randomFeaturesControlDensity)*1.08)
    ## draw empty plot to be able to add backgroumd grey lines below boxes
    plot(1:13, type="n", axes=F, xlab="", ylab="", ylim=ylimit, yaxs="i", xlim=c(1,13))
    ## add vertical gray dashed lines to help reading
    abline(v=c(1:13), col="gray90", lty=3)
    ## add boxplot
    mp <- boxplot(randomFeaturesControlDensity, pch=20, cex=0.8, notch=F, ylab="Number of overlapping DMRs", xlab="", xaxt="n", ylim=ylimit, yaxs="i", whisklty=1, medlty=0, col="dimgray", border="dimgray", add=T, las=2)
    ## text 6% below x-axis / substitute _ for spaces in colnames
    axis(1, labels=NA, at=c(1:13))
    text(c(1:13), par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.06, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)
    ## plot the observed overlap for DMRs (diamond)
    ## red for enrichment, blue for depletion compared to median of randomized data)
    points(realFeatures, col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), pch=18, cex=1.3)
    text(c(1:13), par("usr")[4]*0.97, labels=round(realFeatures / apply(randomFeaturesControlDensity, 2, median), 2), col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), cex=0.8)
    ## header: number of DMRs and name of contrast
    title(plotTitle)
    dev.off()
  }
}
## TO DO? do not plot Human-specific repeats in tissue plots


## Contrast-specific features: DE exons, DS exons, FANTOM enhancers, and other DMRs ###################
library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))

plotFeatures <- function(realFeatures, randomFeatures, type){
  par(mar=c(9, 5, 4, 2) + 0.1, mgp=c(4, 1, 0)) 
  ## let 7% room above for fold-change indication
  ylimit <- c(0, max(realFeatures, randomFeatures)*1.08)
  ## draw empty plot to be able to add background grey lines below boxes
  plot(1:length(realFeatures), type="n", axes=F, xlab="", ylab="", ylim=ylimit, yaxs="i", xlim=c(0,length(realFeatures)+1))
  ## add vertical gray dashed lines to help reading
  abline(v=c(1:length(realFeatures)), col="gray90", lty=3)
  ## add boxplot
  mp <- boxplot(randomFeatures, pch=20, cex=0.8, notch=F, ylab="Number of overlapping DMRs", xlab="", xaxt="n", ylim=ylimit, yaxs="i", whisklty=1, medlty=0, col="dimgray", border="dimgray", add=T, las=2, boxwex=0.8)
  ## text 6% below x-axis / substitute _ for spaces in colnames
  axis(1, labels=NA, at=c(1:length(realFeatures)))
  featNames <- names(realFeatures)
  featNames <- gsub("_", " ", featNames, perl=TRUE)
  text(c(1:length(realFeatures)), par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.06, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)
  ## plot the observed overlap for DMRs (diamond)
  ## red for enrichment, blue for depletion compared to median of randomized data)
  points(realFeatures, col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), pch=18, cex=1.3)
  text(c(1:length(realFeatures)), par("usr")[4]*0.97, labels=round(realFeatures / apply(randomFeatures, 2, median), 2), col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), cex=0.8)
  ## header: number of DMRs and name of contrast
  plotTitle <- vector()
  if (type == 'tissue'){
    plotTitle <- paste0(pair[1], " vs. ", pair[2], " in ", species, ": ", nDMRs, " DMRs")
  }
  if (type == 'species'){
    plotTitle <- paste0(pair[1], " vs. ", pair[2], " in ", tissue, ": ", nDMRs, " DMRs")
  }
  title(plotTitle)
  ## TO DO: to be cleaner, pair and tissue should be passed as argument
}

## Species DMRs
for (tissue in c("heart", "kidney", "liver", "lung")){
  for (pair in list(c("Human", "Chimp"), c("Human", "Rhesus"), c("Chimp", "Rhesus"))){
    cat("Testing species pair: ", pair, "in", tissue, "\n")
    fileName <- paste0("../annotation/DMRs/species/", pair[1], pair[2], "_", tissue, "_DMRs_contrast_specific_features.gz")
    realFeatures <- read.table(fileName, h=T, sep="\t", check.names = F)
    nDMRs <- length(realFeatures[,1])
    realFeatures <- apply(realFeatures[,-c(1:3)], 2, sum)

    ## control for length only
    randomFeatures <- matrix(nrow=100, ncol=19)
    colnames(randomFeatures) = names(realFeatures)
    for (i in 1:100){
      print(i)
      randomFeatures[i,] <- apply(read.table(paste0("../annotation/DMRs/species/randomized/", pair[1], pair[2], "_", tissue, "_randomDMRs_", i, "_contrast_specific_features.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
    }

    ## control for length and CpG density
    randomFeaturesControlDensity <- matrix(nrow=100, ncol=19)
    colnames(randomFeaturesControlDensity) = names(realFeatures)
    for (i in 1:100){      randomFeaturesControlDensity[i,] <- apply(read.table(paste0("../annotation/DMRs/species/randomized_control_CpG_density/", pair[1], pair[2], "_", tissue, "_randomDMRs_", i, "_contrast_specific_features.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
    }

    ## Common plot for control for length and comtrol for length+density (2 pages):
    pdf(file = paste0("../annotation/DMRs/species/", pair[1], pair[2], "_", tissue, "_DMRs_boxplot_contrast_specific_features.pdf"), width = 6, height = 5)

    ## Three separate boxplots: exons / enhancers / DMRs
    plotFeatures(realFeatures[1:2], randomFeatures[,1:2], 'species')
    ## TO DO: add all exons to compare
    ordered <- order(realFeatures[3:11] / apply(randomFeatures[, 3:11], 2, median), decreasing=T)
    plotFeatures(realFeatures[3:11][ordered], randomFeatures[,3:11][, ordered], 'species')
    ordered <- order(realFeatures[12:19] / apply(randomFeatures[, 12:19], 2, median), decreasing=T)
    plotFeatures(realFeatures[12:19][ordered], randomFeatures[,12:19][, ordered], 'species')
    
    ## 2nd boxplot with control for CpG density
    plotFeatures(realFeatures[1:2], randomFeaturesControlDensity[,1:2], 'species')
    ordered <- order(realFeatures[3:11] / apply(randomFeaturesControlDensity[, 3:11], 2, median), decreasing=T)
    plotFeatures(realFeatures[3:11][ordered], randomFeaturesControlDensity[,3:11][, ordered], 'species')
    ordered <- order(realFeatures[12:19] / apply(randomFeaturesControlDensity[, 12:19], 2, median), decreasing=T)
    plotFeatures(realFeatures[12:19][ordered], randomFeaturesControlDensity[,12:19][, ordered], 'species')
 
    dev.off()
  }
}

## Tissue DMRs
for (species in c("Human", "Chimp", "Rhesus")){
  for (pair in list(c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "liver"))){

    cat("Testing tissue pair: ", pair, "in", species, "\n")
    baseName <- paste0(species, "_", pair[1], "_", pair[2])
    fileName <- paste0("../annotation/DMRs/tissues/", baseName, "_DMRs_contrast_specific_features.gz")
    realFeatures <- read.table(fileName, h=T, sep="\t", check.names = F)
    nDMRs <- length(realFeatures[,1])
    realFeatures <- apply(realFeatures[,-c(1:3)], 2, sum)

    ## control for length only
    randomFeatures <- matrix(nrow=100, ncol=length(realFeatures))
    colnames(randomFeatures) = names(realFeatures)
    for (i in 1:100){
      randomFeatures[i,] <- apply(read.table(paste0("../annotation/DMRs/tissues/randomized/", baseName, "_randomDMRs_", i, "_contrast_specific_features.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
    }

    ## control for length and CpG density
    randomFeaturesControlDensity <- matrix(nrow=100, ncol=length(realFeatures))
    colnames(randomFeaturesControlDensity) = names(realFeatures)
    for (i in 1:100){
      randomFeaturesControlDensity[i,] <- apply(read.table(paste0("../annotation/DMRs/tissues/randomized_control_CpG_density/", baseName, "_randomDMRs_", i, "_contrast_specific_features.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
    }

    ## Common plot for control for length and comtrol for length+density (2 pages):
    pdf(file = paste0("../annotation/DMRs/tissues/", baseName, "_DMRs_boxplot_contrast_specific_features.pdf"), width = 6, height = 5)

    ## Three separate boxplots: exons / enhancers / DMRs
    plotFeatures(realFeatures[1:2], randomFeatures[,1:2], 'tissue')
    ordered <- order(realFeatures[3:11] / apply(randomFeatures[, 3:11], 2, median), decreasing=T)
    plotFeatures(realFeatures[3:11][ordered], randomFeatures[,3:11][, ordered], 'tissue')
    ordered <- order(realFeatures[12:length(realFeatures)] / apply(randomFeatures[, 12:length(realFeatures)], 2, median), decreasing=T)
    plotFeatures(realFeatures[12:length(realFeatures)][ordered], randomFeatures[,12:length(realFeatures)][, ordered], 'tissue')
    
    ## 2nd boxplot with control for CpG density
    plotFeatures(realFeatures[1:2], randomFeaturesControlDensity[,1:2], 'tissue')
    ordered <- order(realFeatures[3:11] / apply(randomFeaturesControlDensity[, 3:11], 2, median), decreasing=T)
    plotFeatures(realFeatures[3:11][ordered], randomFeaturesControlDensity[,3:11][, ordered], 'tissue')
    ordered <- order(realFeatures[12:length(realFeatures)] / apply(randomFeaturesControlDensity[, 12:length(realFeatures)], 2, median), decreasing=T)
    plotFeatures(realFeatures[12:length(realFeatures)][ordered], randomFeaturesControlDensity[,12:length(realFeatures)][, ordered], 'tissue')
    dev.off()
  }
}


## Promoters, enhancers (many sources), chromatin states, etc
library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))


plotFeatures <- function(realFeatures, randomFeatures, type, plotTitle){
  par(mar=c(9, 5, 4, 2) + 0.1, mgp=c(4, 1, 0)) 
  ## let 7% room above for fold-change indication
  ylimit <- c(0, max(realFeatures, randomFeatures)*1.08)
  ## draw empty plot to be able to add background grey lines below boxes
  plot(1:length(realFeatures), type="n", axes=F, xlab="", ylab="", ylim=ylimit, yaxs="i", xlim=c(0,length(realFeatures)+1))
  ## add vertical gray dashed lines to help reading
  abline(v=c(1:length(realFeatures)), col="gray90", lty=3)
  ## add boxplot
  mp <- boxplot(randomFeatures, pch=20, cex=0.8, notch=F, ylab="Number of overlapping DMRs", xlab="", xaxt="n", ylim=ylimit, yaxs="i", whisklty=1, medlty=0, col="dimgray", border="dimgray", add=T, las=2, boxwex=0.8)
  ## text 6% below x-axis / substitute _ for spaces in colnames
  axis(1, labels=NA, at=c(1:length(realFeatures)))
  featNames <- names(realFeatures)
  featNames <- gsub("_", " ", featNames, perl=TRUE)
  text(c(1:length(realFeatures)), par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.06, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)
  ## plot the observed overlap for DMRs (diamond)
  ## red for enrichment, blue for depletion compared to median of randomized data)
  points(realFeatures, col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), pch=18, cex=1.3)
  text(c(1:length(realFeatures)), par("usr")[4]*0.97, labels=round(realFeatures / apply(randomFeatures, 2, median), 1), col=ifelse(realFeatures >= mp$stats[3,], pal[1], pal[2]), cex=0.8)
  ## header: number of DMRs and name of contrast
  title(plotTitle)
  ## TO DO: to be cleaner, pair and tissue should be passed as argument
}


## Species DMRs
for (tissue in c("heart", "kidney", "liver", "lung")){
  for (pair in list(c("Human", "Chimp"), c("Human", "Rhesus"), c("Chimp", "Rhesus"), c("Human", "Specific"), c("Chimp", "Specific"))){
    cat("Testing species pair: ", pair, "in", tissue, "\n")
    fileName <- paste0("../annotation/DMRs/species/", pair[1], pair[2], "_", tissue, "_DMRs_promoters_enhancers_features.gz")
    realFeatures <- read.table(fileName, h=T, sep="\t", check.names = F)
    nDMRs <- length(realFeatures[,1])
    realFeatures <- apply(realFeatures[,-c(1:3)], 2, sum)

    ## control for length only
    randomFeatures <- matrix(nrow=100, ncol=length(realFeatures))
    colnames(randomFeatures) = names(realFeatures)
    for (i in 1:100){
      #print(i)
      randomFeatures[i,] <- apply(read.table(paste0("../annotation/DMRs/species/randomized/", pair[1], pair[2], "_", tissue, "_randomDMRs_", i, "_promoters_enhancers_features.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
    }

    ## control for length and CpG density
    randomFeaturesControlDensity <- matrix(nrow=100, ncol=length(realFeatures))
    colnames(randomFeaturesControlDensity) = names(realFeatures)
    for (i in 1:100){
      #print(i)
      randomFeaturesControlDensity[i,] <- apply(read.table(paste0("../annotation/DMRs/species/randomized_control_CpG_density/", pair[1], pair[2], "_", tissue, "_randomDMRs_", i, "_promoters_enhancers_features.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
    }

    ## Common plot for control for length and comtrol for length+density (2 pages):
    pdf(file = paste0("../annotation/DMRs/species/", pair[1], pair[2], "_", tissue, "_DMRs_boxplot_promoters_enhancers_features.pdf"), width = 6, height = 5)

    plotTitle <- vector()
    if (pair[2] == "Specific"){
      plotTitle <- paste0(pair[1], "-specific in ", tissue, ": ", nDMRs, " DMRs")
    } else {
      plotTitle <- paste0(pair[1], " vs. ", pair[2], " in ", tissue, ": ", nDMRs, " DMRs")
    }

    for(columns in list(c(1:9), c(10:14), c(15,16), c(17:19), c(20:22), c(23:37), c(38:52), c(53:67), c(68:82), c(83:97), c(98:112), c(113:127), c(128:142))){
      ordered <- order(realFeatures[columns] / apply(randomFeatures[, columns], 2, median), decreasing=T)
      plotFeatures(realFeatures[columns][ordered], randomFeatures[,columns][, ordered], 'species', plotTitle)
    }
    ## 2nd boxplot with control for CpG density
    for(columns in list(c(1:9), c(10:14), c(15,16), c(17:19), c(20:22), c(23:37), c(38:52), c(53:67), c(68:82), c(83:97), c(98:112), c(113:127), c(128:142))){
      ordered <- order(realFeatures[columns] / apply(randomFeaturesControlDensity[, columns], 2, median), decreasing=T)
      plotFeatures(realFeatures[columns][ordered], randomFeaturesControlDensity[,columns][, ordered], 'species', plotTitle)
    }
    dev.off()
  }
}

## Tissue DMRs
for (species in c("Human", "Chimp", "Rhesus")){
  for (pair in list(c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "liver"), c("heart", "Specific"), c("liver", "Specific"), c("kidney", "Specific"), c("lung", "Specific"))){
    cat("Testing tissue pair: ", pair, "in", species, "\n")
    if (pair[2] != "Specific"){
      baseName <- paste0(species, "_", pair[1], "_", pair[2])
    } else {
      baseName <- paste0(species, "_", pair[1], "Specific")
    }
    fileName <- paste0("../annotation/DMRs/tissues/", baseName, "_DMRs_promoters_enhancers_features.gz")
    realFeatures <- read.table(fileName, h=T, sep="\t", check.names = F)
    nDMRs <- length(realFeatures[,1])
    realFeatures <- apply(realFeatures[,-c(1:3)], 2, sum)

    ## control for length only
    randomFeatures <- matrix(nrow=100, ncol=length(realFeatures))
    colnames(randomFeatures) = names(realFeatures)
    for (i in 1:100){
      #print(i)
      randomFeatures[i,] <- apply(read.table(paste0("../annotation/DMRs/tissues/randomized/", baseName, "_randomDMRs_", i, "_promoters_enhancers_features.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
    }

    ## control for length and CpG density
    randomFeaturesControlDensity <- matrix(nrow=100, ncol=length(realFeatures))
    colnames(randomFeaturesControlDensity) = names(realFeatures)
    for (i in 1:100){
      #print(i)
      randomFeaturesControlDensity[i,] <- apply(read.table(paste0("../annotation/DMRs/tissues/randomized_control_CpG_density/", baseName, "_randomDMRs_", i, "_promoters_enhancers_features.gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
    }

    ## Common plot for control for length and comtrol for length+density (2 pages):
    pdf(file = paste0("../annotation/DMRs/tissues/", baseName, "_DMRs_boxplot_promoters_enhancers_features.pdf"), width = 6, height = 5)

    plotTitle <- vector()
    if (pair[2] == "Specific"){
      plotTitle <- paste0(pair[1], "-specific in ", species, ": ", nDMRs, " DMRs")
    } else {
      plotTitle <- paste0(pair[1], " vs. ", pair[2], " in ", species, ": ", nDMRs, " DMRs")
    }

    for( columns in list(c(1:9), c(10:14), c(15,16), c(17:19), c(20:22), c(23:37), c(38:52), c(53:67), c(68:82), c(83:97), c(98:112), c(113:127), c(128:142))){
      ordered <- order(realFeatures[columns] / apply(randomFeatures[, columns], 2, median), decreasing=T)
      plotFeatures(realFeatures[columns][ordered], randomFeatures[,columns][, ordered], 'tissue', plotTitle)
    }
    ## 2nd boxplot with control for CpG density
    for( columns in list(c(1:9), c(10:14), c(15,16), c(17:19), c(20:22), c(23:37), c(38:52), c(53:67), c(68:82), c(83:97), c(98:112), c(113:127), c(128:142))){
      ordered <- order(realFeatures[columns] / apply(randomFeaturesControlDensity[, columns], 2, median), decreasing=T)
      plotFeatures(realFeatures[columns][ordered], randomFeaturesControlDensity[,columns][, ordered], 'tissue', plotTitle)
    }
    dev.off()
  }
}

###########################################################################
## Boxplots of genomic features overlap, but for some sDMRs or tDMRs only #
###########################################################################
## Try to put random and random+controlCpGdensity DMRs side by side

library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))
plotFeaturesSideBySide <- function(realFeatures, randomFeatures, randomFeaturesControlCpG, plotFile){

  ## We want plot window to be proportional to number of features (so that different lots could be juxtaposed using Illustrator)
  ## http://www.r-bloggers.com/setting-graph-margins-in-r-using-the-par-function-and-lots-of-cow-milk/
  ## 1 margin line = 0.2 inches
  ## Let's aim for 9 features -> width = 6 inches (whole plot, looks good)
  ## Window width in inches = 6 - 5*0.2 - 2*0.2 = 4.6 inches (see below par() command for margins size)
  ## pdf width = (9 features, so 8 intervals + 0.5 left + 0.5 right)*0.46 + 1.4
  ## /!\ xaxs = "i" option should be added to plot functions to stick to the xlimits
  pdf(file = plotFile, width = (length(realFeatures)*0.46 + 1.4), height = 5)

  ## Larger margins below plot, ylabel at 3.5 lines from y-axis (mgp), y-axis numbers at 0.7 line from y-axis
  par(mar=c(10, 5, 2, 2) + 0.1, mgp=c(3.5, 0.7, 0)) 
  ## let 12% room above for fold-change indication
  ylimit <- c(0, max(realFeatures, randomFeatures, randomFeaturesControlCpG)*1.12)
  ## draw empty plot to be able to add background grey lines below boxes
  plot(1:length(realFeatures), type="n", axes=F, xlab="", ylab="Number of overlapping DMRs", ylim=ylimit, yaxs="i", xaxs = "i", xlim=c(0.5,length(realFeatures)+0.5))
  ## add vertical gray dashed lines to make reading easier
  abline(v=c(1:length(realFeatures)), col="gray90", lty=3)

  ## add boxplots of random DMRs
  mp1 <- boxplot(randomFeatures, pch=20, cex=0.5, notch=F, ylab="", xlab="", xaxt="n", ylim=ylimit, xaxs = "i", yaxs="i", whisklty=1, medlty=0, col=gray(0.8), border=gray(0.8), add=T, las=2, boxwex=0.3, at=1:length(realFeatures)-0.15)
  ## add boxplots of randomFeaturesControlCpG DMRs
  mp2 <- boxplot(randomFeaturesControlCpG, pch=20, cex=0.5, notch=F, ylab="", xlab="", xaxt="n", ylim=ylimit, xaxs = "i", yaxs="i", whisklty=1, medlty=0, col=gray(0.4), border=gray(0.4), add=T, las=2, boxwex=0.3, at=1:length(realFeatures)+0.15)

  ## plot the observed overlap for DMRs (diamond)
  points(realFeatures, col=ifelse(realFeatures >= mp2$stats[3,], pal[1], pal[2]), pch=18, cex=1.3)
  ## red for enrichment, blue for depletion compared to median of randomized data + control CpG)
 
  ## Enrichment compared to random features: on enrichment next to the other in different colors
  ## 2 significant digits
  text(c(1:length(realFeatures)), par("usr")[4]*0.96, labels=signif(realFeatures / apply(randomFeatures, 2, median), 2), col=gray(0.8), cex=0.8, adj=c(1.1, 0.9))
  text(c(1:length(realFeatures)), par("usr")[4]*0.96, labels="\\", col="black", cex=0.8, adj=c(0.5, 0.5)) ## separating slash
  text(c(1:length(realFeatures)), par("usr")[4]*0.96, labels=signif(realFeatures / apply(randomFeaturesControlCpG, 2, median), 2), col=gray(0.4), cex=0.8, adj=c(-0.1, 0.1))
  
  ## feature labels 6% below x-axis / substitute _ for spaces in colnames
  axis(1, labels=NA, at=c(1:length(realFeatures)))
  featNames <- names(realFeatures)
  featNames <- gsub("_", " ", featNames, perl=TRUE)
  text(c(1:length(realFeatures)), par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.06, srt = 45, adj = 1, labels = featNames, xpd = TRUE, cex=1)

  dev.off()
}
## ## Test: 
## columns <- c(1:9)
## ordered <- order(realFeatures[columns] / apply(randomFeatures[, columns], 2, median), decreasing=T)
## plotFeaturesSideBySide(realFeatures[columns][ordered], randomFeatures[,columns][, ordered], randomFeaturesControlDensity[,columns][, ordered], "../annotation/test.pdf")
 
## Human heart vs. liver tDMRs
plotFeaturesHumanHeartVsLiver <- function(features, columns){
  species <-  "Human"
  pair <- c("heart", "liver")
  baseName <- paste0(species, "_", pair[1], "_", pair[2])
  fileName <- paste0("../annotation/DMRs/tissues/", baseName, "_DMRs_", features, ".gz")
  realFeatures <- read.table(fileName, h=T, sep="\t", check.names = F)
  nDMRs <- length(realFeatures[,1])
  realFeatures <- apply(realFeatures[,-c(1:3)], 2, sum)

  ## randomDMRs: control for length only
  print("Reading random DMRs...")
  randomFeatures <- matrix(nrow=100, ncol=length(realFeatures))
  colnames(randomFeatures) = names(realFeatures)
  for (i in 1:100){
    print(i)
    randomFeatures[i,] <- apply(read.table(paste0("../annotation/DMRs/tissues/randomized/", baseName, "_randomDMRs_", i, "_", features, ".gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
  }
  ## randomDMRs: control for length and CpG density
  print("Reading random DMRs + control for CpG density...")
  randomFeaturesControlDensity <- matrix(nrow=100, ncol=length(realFeatures))
  colnames(randomFeaturesControlDensity) = names(realFeatures)
  for (i in 1:100){
    print(i)
    randomFeaturesControlDensity[i,] <- apply(read.table(paste0("../annotation/DMRs/tissues/randomized_control_CpG_density/", baseName, "_randomDMRs_", i, "_", features, ".gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
  }

  ## Plot boxplot of features overlap for selected columns
  print("Plotting...")
  ordered <- order(realFeatures[columns] / apply(randomFeaturesControlDensity[, columns], 2, median), decreasing=T)
  columnsForFileName <- ""
  if (length(columns) == columns[length(columns)] - columns[1] + 1){
    columnsForFileName <- paste0("columns_", columns[1], "-", columns[length(columns)])
  } else {
    columnsForFileName <- paste0("columns_", paste(columns, collapse="_"))
  }
  plotFeaturesSideBySide(realFeatures[columns][ordered], randomFeatures[,columns][, ordered], randomFeaturesControlDensity[,columns][, ordered], paste0("../annotation/DMRs/tissues/", baseName, "_DMRs_boxplot_", features, "_side_by_side_", columnsForFileName, ".pdf"))
}

plotFeaturesHumanHeartVsLiver("15_features_0.2", c(1:15))
plotFeaturesHumanHeartVsLiver("15_features_0.2", c(1:7, 9:15)) ## removed conserved elements
plotFeaturesHumanHeartVsLiver("15_features_0.2", c(1:11, 13:15)) ## removed repeats

plotFeaturesHumanHeartVsLiver("repeats", c(1:13))

for( columns in list(c(1:2), c(12:17))){
  plotFeaturesHumanHeartVsLiver("contrast_specific_features", columns)
}

for( columns in list(c(1:9), c(10:14), c(15,16), c(17:19), c(20:22), c(23:37), c(38:52), c(53:67), c(68:82), c(83:97), c(98:112), c(113:127), c(128:142))){
  plotFeaturesHumanHeartVsLiver("promoters_enhancers_features", columns)
}
plotFeaturesHumanHeartVsLiver("promoters_enhancers_features", c(1:14)) ## FANTOM enhancers and promoters together
plotFeaturesHumanHeartVsLiver("promoters_enhancers_features", c(1,3:5,9:14)) ## FANTOM tissue-specific enhancers and promoters together
plotFeaturesHumanHeartVsLiver("promoters_enhancers_features", c(17:22)) ## ENCODE enhancers and promoters together
plotFeaturesHumanHeartVsLiver("promoters_enhancers_features", c(1:22)) ## FANTOM, VISTA, ENCODE enhancers and promoters together



## repeat for Human vs. Chimp sDMRs in heart
plotFeaturesHumanVsChimpHeart <- function(features, columns){
  tissue <-  "heart"
  pair <- c("Human", "Chimp")
  baseName <- paste0(pair[1], pair[2], "_", tissue)
  fileName <- paste0("../annotation/DMRs/species/", baseName, "_DMRs_", features, ".gz")
  realFeatures <- read.table(fileName, h=T, sep="\t", check.names = F)
  nDMRs <- length(realFeatures[,1])
  realFeatures <- apply(realFeatures[,-c(1:3)], 2, sum)

  ## randomDMRs: control for length only
  print("Reading random DMRs...")
  randomFeatures <- matrix(nrow=100, ncol=length(realFeatures))
  colnames(randomFeatures) = names(realFeatures)
  for (i in 1:100){
    print(i)
    randomFeatures[i,] <- apply(read.table(paste0("../annotation/DMRs/species/randomized/", baseName, "_randomDMRs_", i, "_", features, ".gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
  }
  ## randomDMRs: control for length and CpG density
  print("Reading random DMRs + control for CpG density...")
  randomFeaturesControlDensity <- matrix(nrow=100, ncol=length(realFeatures))
  colnames(randomFeaturesControlDensity) = names(realFeatures)
  for (i in 1:100){
    print(i)
    randomFeaturesControlDensity[i,] <- apply(read.table(paste0("../annotation/DMRs/species/randomized_control_CpG_density/", baseName, "_randomDMRs_", i, "_", features, ".gz"), h=T, sep="\t", check.names = F)[,-c(1:3)], 2, sum)
  }

  ## Plot boxplot of features overlap for selected columns
  print("Plotting...")
  ordered <- order(realFeatures[columns] / apply(randomFeaturesControlDensity[, columns], 2, median), decreasing=T)
  columnsForFileName <- ""
  if (length(columns) == columns[length(columns)] - columns[1] + 1){
    columnsForFileName <- paste0("columns_", columns[1], "-", columns[length(columns)])
  } else {
    columnsForFileName <- paste0("columns_", paste(columns, collapse="_"))
  }
  plotFeaturesSideBySide(realFeatures[columns][ordered], randomFeatures[,columns][, ordered], randomFeaturesControlDensity[,columns][, ordered], paste0("../annotation/DMRs/species/", baseName, "_DMRs_boxplot_", features, "_side_by_side_", columnsForFileName, ".pdf"))
}

plotFeaturesHumanVsChimpHeart("15_features_0.2", c(1:15))
plotFeaturesHumanVsChimpHeart("15_features_0.2", c(1:7, 9:15)) ## removed conserved elements
plotFeaturesHumanVsChimpHeart("15_features_0.2", c(1:11, 13:15)) ## removed repeats

plotFeaturesHumanVsChimpHeart("repeats", c(1:13))

for( columns in list(c(1:2), c(12:19))){
  plotFeaturesHumanVsChimpHeart("contrast_specific_features", columns)
}

for( columns in list(c(1:9), c(10:14), c(15,16), c(17:19), c(20:22), c(23:37), c(38:52), c(53:67), c(68:82), c(83:97), c(98:112), c(113:127), c(128:142))){
  plotFeaturesHumanVsChimpHeart("promoters_enhancers_features", columns)
}
plotFeaturesHumanVsChimpHeart("promoters_enhancers_features", c(1:14)) ## FANTOM enhancers and promoters together
plotFeaturesHumanVsChimpHeart("promoters_enhancers_features", c(1,3:5,9:14)) ## FANTOM tissue-specific enhancers and promoters together
plotFeaturesHumanVsChimpHeart("promoters_enhancers_features", c(17:22)) ## ENCODE enhancers and promoters together
plotFeaturesHumanVsChimpHeart("promoters_enhancers_features", c(1:22)) ## FANTOM, VISTA, ENCODE enhancers and promoters together

## TO DO: another contrast involving more distant species? Human vs. Rhesus in heart?
##        Maybe it would be good to show another tissue? Human vs. Rhesus in liver?


#################################################################
## Heatmap of mean methylation in all samples at DMR locations ##
#################################################################
## Similarly to Irizarry et al. 2009, Figure 6 and 7
library(bsseq)
load("smooth_data/combined_samples/pData.RDa")
samples <- pData

## which contrast to study
contrast <- "tissues/Human_heart_liver"
## contrast <- "tissues/Human_lung_kidney"
## contrast <- "tissues/Rhesus_lung_kidney"
## contrast <- "species/HumanChimp_liver"
## contrast <- "species/HumanChimp_lung"

## Load a DMR file and record coordinates, for example: 
dmrs <- read.table(paste0("./DMRs/", contrast, "_DMRs.txt"), sep="\t", h=T)
## create GRanges object
library("GenomicRanges")
dmrsGr <- GRanges(dmrs$chr, IRanges(dmrs$start, dmrs$end))

## Get mean methylation at DMR positions for all samples (long!)
meanMeth <- matrix(nrow=length(dmrs[,1]), ncol=48)
colnames(meanMeth) <- unique(samples$Condition)
for (sample in unique(samples$Condition)){
  cat(sample, "\n")
  ## load sample and calculate the mean methylation at each DMR location
  load(paste0("./smooth_data/smoothed_samples/", sample, "_Smoothed.RDa"))

  meanMeth[,sample] <- getMeth(data.fit, regions = dmrsGr, type = "smooth", what = "perRegion")
}

## plot heatmap
library(gplots)
## colors <- rev(colorRampPalette(brewer.pal(10,"RdBu"))(100)) ## We want red to be higher methylation and blue lower
colors <- colorRampPalette(c(brewer.pal(9, "Blues")[1],brewer.pal(9, "Blues")[9]))(100) ## Actually it is better to use color scale using white to dark blue

## pdf(file = paste0("DMRs/", contrast, "_DMRs_heatmap.pdf"), width = 12, height = 8)
pdf(file = "test.pdf", width = 12, height = 8)
## testing on 100 DMRs
h <- heatmap.2( meanMeth[1:100,], scale="none", col = colors, margins = c(12, 12), trace='none', denscol="white", labCol=paste(samples$Species, samples$Tissue, sep=" "), labRow=NA, ColSideColors=pal[as.integer(as.factor(samples$Species))+9], cexCol = 1.5, dendrogram='col', key=T)
dev.off()
## problem: we would like the color range to span 0 to 1 fully

## Create breaks manually
palette.breaks <- seq(0, 1, 0.01)
colors <- colorRampPalette(c(brewer.pal(9, "Blues")[1],brewer.pal(9, "Blues")[9]))(100) ## Actually it is better to use color scale using white to dark blue
pdf(file = "test.pdf", width = 12, height = 8)
## testing on 100 DMRs
h <- heatmap.2( meanMeth[1:100,], scale="none", col = colors, breaks = palette.breaks, margins = c(12, 12), trace='none', denscol="white", labCol=paste(samples$Species, samples$Tissue, sep=" "), labRow=NA, ColSideColors=pal[as.integer(as.factor(samples$Species))+9], cexCol = 1.5, dendrogram='col', key=T)
dev.off()
## TO DO: modify function plot_dmrs to accomodate the solution

## Heatmap using only DMRs conserved to rhesus
cons <- read.table(paste0("DMRs/", contrast, "_DMRs_conservation.txt"), h=T)
## Loading mean Meth object
load(paste0("DMRs/", contrast, "_DMRs_mean_Methylation.RDa"))
## conserved DMRs
head(meanMeth[cons$conservation == 4,]) 
pdf(file = "test.pdf", width = 12, height = 8)
h <- heatmap.2( meanMeth[cons$conservation == 4,], scale="none", col = colors, breaks = palette.breaks, margins = c(12, 12), trace='none', denscol="white", labCol=paste(samples$Species, samples$Tissue, sep=" "), labRow=NA, ColSideColors=pal[as.integer(as.factor(samples$Species))+9], cexCol = 1.5, dendrogram='col', key=T)
dev.off()
## Cluster by tissue! Makes sense ;)


## Repeat for all DMRs types and contrasts ###############################################################
allSpecies <- c("Human", "Chimp", "Rhesus")
tissuePairs <- list(c("heart", "liver"), c("heart", "Specific"), c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("liver", "Specific"), c("kidney", "Specific"), c("lung", "Specific"))
for (species in allSpecies){
  for (pair in tissuePairs){
    cat("Testing tissue pair: ", pair, "in", species, "\n")
    if (pair[2] == "Specific"){
      contrast <- paste0("tissues/", species, "_", pair[1], pair[2])
    } else {
      contrast <- paste0("tissues/", species, "_", pair[1], "_", pair[2])
    }
    plot_dmrs(contrast)
  }
}

speciesPairs <- list(c("Human", "Chimp"), c("Human", "Rhesus"), c("Chimp", "Rhesus"), c("Human", "Specific"), c("Chimp", "Specific"))
allTissues <- c("heart", "kidney", "liver", "lung")
  
for (tissue in allTissues){
  for (pair in speciesPairs){
    cat("Testing species pair: ", pair, "in", tissue, "\n")
    contrast <- paste0("species/", pair[1], pair[2], "_", tissue)
    plot_dmrs(contrast)
  }
}

## allSpecies <- c("Human", "Chimp", "Rhesus")
## individualPairs <- list(c("1", "2"), c("1", "3"), c("1", "4"), c("2", "3"), c("2", "4"), c("3", "4"))
## for (species in allSpecies){
##   for (pair in individualPairs){
##     cat("Testing individual pair: ", pair, "in", species, "\n")
##     contrast <- paste0("tissues/", species, "_", pair[1], "_", pair[2])
##     plot_dmrs(contrast)
##   }
## }


## Repeat for all DMRs types and contrasts ###############################################################
## Heatmap on correlation matrix
## Color scale: white to dark blue
source("functions.R")
library(gplots)
colors <- colorRampPalette(c(brewer.pal(9, "Blues")[1],brewer.pal(9, "Blues")[9]))(100)

allSpecies <- c("Human", "Chimp", "Rhesus")
tissuePairs <- list(c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "liver"), c("heart", "Specific"), c("liver", "Specific"), c("kidney", "Specific"), c("lung", "Specific"))
for (species in allSpecies){
  for (pair in tissuePairs){
    cat("Testing tissue pair: ", pair, "in", species, "\n")
    if (pair[2] == "Specific"){
      contrast <- paste0("tissues/", species, "_", pair[1], pair[2])
    } else {
      contrast <- paste0("tissues/", species, "_", pair[1], "_", pair[2])
    }
    plot_dmrs_correlation(contrast)
  }
}

speciesPairs <- list(c("Human", "Chimp"), c("Human", "Rhesus"), c("Chimp", "Rhesus"), c("Human", "Specific"), c("Chimp", "Specific"))
allTissues <- c("heart", "kidney", "liver", "lung")
 
for (tissue in allTissues){
  for (pair in speciesPairs){
    cat("Testing species pair: ", pair, "in", tissue, "\n")
    contrast <- paste0("species/", pair[1], pair[2], "_", tissue)
    plot_dmrs_correlation(contrast)
  }
}

## Repeat for conserved tDMRs ###############################################################
## Heatmap on correlation matrix
## Color scale: white to dark blue
source("functions.R")
library(gplots)
colors <- colorRampPalette(c(brewer.pal(9, "Blues")[1],brewer.pal(9, "Blues")[9]))(100)

allSpecies <- c("Human", "Chimp", "Rhesus")
tissuePairs <- list(c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "liver"), c("heart", "Specific"), c("liver", "Specific"), c("kidney", "Specific"), c("lung", "Specific"))
for (species in allSpecies){
  for (pair in tissuePairs){
    cat("Testing tissue pair: ", pair, "in", species, "\n")
    if (pair[2] == "Specific"){
      contrast <- paste0("tissues/", species, "_", pair[1], pair[2])
    } else {
      contrast <- paste0("tissues/", species, "_", pair[1], "_", pair[2])
    }
    ## Conserved human-chimp-macaque
    plot_dmrs_correlation_conserved(contrast, 4)
  }
}


## As in Irizarry 2009, pool together all tissue-specific DMRs in 1 species and plot heatmap #################
## Do we see perfect tissue clustering as they do?

## Use Rhesus DMRs at first
tissuePairs <- list(c("heart", "Specific"), c("kidney", "Specific"), c("lung", "Specific"), c("liver", "Specific"))
i <- 0
for (pair in tissuePairs){
  cat("Testing tissue pair: ", pair, "in Human\n")
  contrast <- paste0("tissues/Rhesus_", pair[1], pair[2])
  load(paste0("DMRs/", contrast, "_DMRs_mean_Methylation.RDa"))
  if ( i==0 ){
    allDMRs <- meanMeth
  } else {
    allDMRs <- rbind(allDMRs, meanMeth)
  }
  i <- i+1
}
dim(allDMRs) ## 38950

## plot heatmap: too large. Sample 5,000 rows at random
library(gplots)
colors <- colorRampPalette(c(brewer.pal(9, "Blues")[1],brewer.pal(9, "Blues")[9]))(100)
palette.breaks <- seq(0, 1, 0.01)

pdf(file = paste0("heatmap_all_Rhesus_tissue-specific_tDMRs.pdf"), width = 12, height = 8)
h <- heatmap.2( allDMRs[  sample(1:length(allDMRs[,1]), 5000),], scale="none", col = colors, breaks = palette.breaks, margins = c(12, 12), trace='none', denscol="white", labCol=paste(samples$Species, samples$Tissue, sep=" "), labRow=NA, ColSideColors=pal[as.integer(as.factor(samples$Species))+9], cexCol = 1.5, dendrogram='col', key=T)

## calculate correlation matrix
cors <- cor(allDMRs, method="spearman", use="pairwise.complete.obs") 

h <- heatmap.2( cors, scale="none", col = colors, margins = c(12, 12), trace='none', denscol="white", labCol=paste(samples$Species, samples$Tissue, sep=" "), labRow=samples$Condition, ColSideColors=pal[as.integer(as.factor(samples$Species))+9], RowSideColors=pal[as.integer(as.factor(samples$Tissue))], cexCol = 1.5)
dev.off()
## Rhesus clusters out
## On correlation matrix heatmap, liver samples cluster out. Because more liverSpecific DMRs probably: 

## Use 1000 of each tissue-specific DMR type only:
tissuePairs <- list(c("heart", "Specific"), c("kidney", "Specific"), c("lung", "Specific"), c("liver", "Specific"))
i <- 0
for (pair in tissuePairs){
  cat("Testing tissue pair: ", pair, "in Rhesus\n")
  contrast <- paste0("tissues/Rhesus_", pair[1], pair[2])
  load(paste0("DMRs/", contrast, "_DMRs_mean_Methylation.RDa"))
  if ( i==0 ){
    allDMRs <- meanMeth[sample(1:length(meanMeth[,1]), 1000), ]
  } else {
    allDMRs <- rbind(allDMRs, meanMeth[sample(1:length(meanMeth[,1]), 1000), ])
  }
  i <- i+1
}
dim(allDMRs) ## 4000
## In theory tissue-specific DMRs are not supposed to overlap, so we're not dealing with this here

## plot heatmaps: 
pdf(file = paste0("heatmap_all_Rhesus_tissue-specific_tDMRs_1000_each.pdf"), width = 12, height = 8)
h <- heatmap.2( allDMRs, scale="none", col = colors, breaks = palette.breaks, margins = c(12, 12), trace='none', denscol="white", labCol=paste(samples$Species, samples$Tissue, sep=" "), labRow=NA, ColSideColors=pal[as.integer(as.factor(samples$Species))+9], cexCol = 1.5, dendrogram='col', key=T)

## calculate correlation matrix
cors <- cor(allDMRs, method="spearman", use="pairwise.complete.obs") 

h <- heatmap.2( cors, scale="none", col = colors, margins = c(12, 12), trace='none', denscol="white", labCol=paste(samples$Species, samples$Tissue, sep=" "), labRow=samples$Condition, ColSideColors=pal[as.integer(as.factor(samples$Species))+9], RowSideColors=pal[as.integer(as.factor(samples$Tissue))], cexCol = 1.5)
dev.off()
## Rhesus clusters out: in contradiction with Irizarry et al.
## Used in sup Figure

## TO DO: repeat as they do: center-nromalzie data by species

## TO DO: repeat for conserved DMRs

## Same thing in human and in chimp
i <- 0
for (pair in tissuePairs){
  ## contrast <- paste0("tissues/Chimp_", pair[1], pair[2])
  contrast <- paste0("tissues/Human_", pair[1], pair[2])
  load(paste0("DMRs/", contrast, "_DMRs_mean_Methylation.RDa"))
  if ( i==0 ){
    allDMRs <- meanMeth[sample(1:length(meanMeth[,1]), 1000), ]
  } else {
    allDMRs <- rbind(allDMRs, meanMeth[sample(1:length(meanMeth[,1]), 1000), ])
  }
  i <- i+1
}
## pdf(file = paste0("heatmap_all_Chimp_tissue-specific_tDMRs_1000_each.pdf"), width = 12, height = 8)
pdf(file = paste0("heatmap_all_Human_tissue-specific_tDMRs_1000_each.pdf"), width = 12, height = 8)
h <- heatmap.2( allDMRs, scale="none", col = colors, breaks = palette.breaks, margins = c(12, 12), trace='none', denscol="white", labCol=paste(samples$Species, samples$Tissue, sep=" "), labRow=NA, ColSideColors=pal[as.integer(as.factor(samples$Species))+9], cexCol = 1.5, dendrogram='col', key=T)
cors <- cor(allDMRs, method="spearman", use="pairwise.complete.obs") 
h <- heatmap.2( cors, scale="none", col = colors, margins = c(12, 12), trace='none', denscol="white", labCol=paste(samples$Species, samples$Tissue, sep=" "), labRow=samples$Condition, ColSideColors=pal[as.integer(as.factor(samples$Species))+9], RowSideColors=pal[as.integer(as.factor(samples$Tissue))], cexCol = 1.5)
dev.off()

## Pool together all pairwise tissue DMRs in 1 species and plot heatmap ###############################
## Do we expect perfect tissue clustering? for human-chimp at least?

## Use human DMRs at first
tissuePairs <- list(c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "liver"))
allDMRs <- GRanges()
i <- 0
for (pair in tissuePairs){
  cat("Testing tissue pair: ", pair, "in Human\n")

  dmrs <- read.table(paste0("./DMRs/tissues/Human_", pair[1], "_", pair[2], "_DMRs.txt"), sep="\t", h=T)
  ## create GRanges object
  dmrsGr <- GRanges(dmrs$chr, IRanges(dmrs$start, dmrs$end))
  allDMRs <- union(allDMRs, dmrsGr)
}
length(allDMRs) ## 61629

## Get mean methylation at DMR positions for all samples
meanMeth <- matrix(nrow=length(allDMRs), ncol=48)
colnames(meanMeth) <- unique(samples$Condition)
for (sample in unique(samples$Condition)){
  cat(sample, "\n")
  ## load sample and calculate the mean methylation at each DMR location
  load(paste0("./smooth_data/smoothed_samples/", sample, "_Smoothed.RDa"))
  meanMeth[,sample] <- getMeth(data.fit, regions = allDMRs, type = "smooth", what = "perRegion")
}
save(meanMeth, file="mean_methylation_all_Human_tDMRs.RDa")
## plot heatmap: too large. Sample 10,000 rows at random
pdf(file = paste0("heatmap_all_Human_tDMRs.pdf"), width = 12, height = 8)
h <- heatmap.2( meanMeth[  sample(1:length(meanMeth[,1]), 10000),], scale="none", col = colors, margins = c(12, 12), trace='none', denscol="white", labCol=paste(samples$Species, samples$Tissue, sep=" "), labRow=NA, ColSideColors=pal[as.integer(as.factor(samples$Species))], cexCol = 1.5, dendrogram='col', key=FALSE)
dev.off()
## Here we used all human tDMRs, even those with no orthology in Chimp or Macaque
## TO DO: retsrict to tDMRs with full orthology

## TO DO: heatmap using all CpGs falling within any tDMR

## Use only the top 5000 top DMRs for each tissue contrast ##################################
tissuePairs <- list(c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "liver"))
allDMRs <- GRanges()
i <- 0
for (pair in tissuePairs){
  cat("Testing tissue pair: ", pair, "in Human\n")

  dmrs <- read.table(paste0("./DMRs/tissues/Human_", pair[1], "_", pair[2], "_DMRs.txt"), sep="\t", h=T)
  ## create GRanges object
  dmrsGr <- GRanges(dmrs$chr, IRanges(dmrs$start, dmrs$end))
  allDMRs <- union(allDMRs, dmrsGr[1:5000])
}
length(allDMRs)

## Get mean methylation at DMR positions for all samples
meanMeth <- matrix(nrow=length(allDMRs), ncol=48)
colnames(meanMeth) <- unique(samples$Condition)
for (sample in unique(samples$Condition)){
  cat(sample, "\n")
  ## load sample and calculate the mean methylation at each DMR location
  load(paste0("./smooth_data/smoothed_samples/", sample, "_Smoothed.RDa"))
  meanMeth[,sample] <- getMeth(data.fit, regions = allDMRs, type = "smooth", what = "perRegion")
}
save(meanMeth, file="mean_methylation_top5000_Human_tDMRs.RDa")
pdf(file = paste0("heatmap_top5000_Human_tDMRs.pdf"), width = 12, height = 8)
h <- heatmap.2( meanMeth, scale="none", col = colors, margins = c(12, 12), trace='none', denscol="white", labCol=paste(samples$Species, samples$Tissue, sep=" "), labRow=NA, ColSideColors=pal[as.integer(as.factor(samples$Species))], cexCol = 1.5, dendrogram='col', key=FALSE)
dev.off()


## Use only 10000 random DMRs for each tissue contrast ############################################
tissuePairs <- list(c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "liver"))
allDMRs <- GRanges()
i <- 0
for (pair in tissuePairs){
  cat("Testing tissue pair: ", pair, "in Human\n")

  dmrs <- read.table(paste0("./DMRs/tissues/Human_", pair[1], "_", pair[2], "_DMRs.txt"), sep="\t", h=T)
  ## create GRanges object
  dmrsGr <- GRanges(dmrs$chr, IRanges(dmrs$start, dmrs$end))
  allDMRs <- union(allDMRs, dmrsGr[ sample(1:length(dmrsGr), 5000)])
}
length(allDMRs) ## 25199
## Random DMRs overalp less between tissue contrasts than the top DMRs (see above)

## Get mean methylation at DMR positions for all samples
meanMeth <- matrix(nrow=length(allDMRs), ncol=48)
colnames(meanMeth) <- unique(samples$Condition)
for (sample in unique(samples$Condition)){
  cat(sample, "\n")
  ## load sample and calculate the mean methylation at each DMR location
  load(paste0("./smooth_data/smoothed_samples/", sample, "_Smoothed.RDa"))
  meanMeth[,sample] <- getMeth(data.fit, regions = allDMRs, type = "smooth", what = "perRegion")
}
save(meanMeth, file="mean_methylation_random5000_Human_tDMRs.RDa")
pdf(file = paste0("heatmap_random5000_Human_tDMRs.pdf"), width = 12, height = 8)
h <- heatmap.2( meanMeth, scale="none", col = colors, margins = c(12, 12), trace='none', denscol="white", labCol=paste(samples$Species, samples$Tissue, sep=" "), labRow=NA, ColSideColors=pal[as.integer(as.factor(samples$Species))], cexCol = 1.5, dendrogram='col', key=FALSE)
dev.off()
## Similar than using all DMRs


## limit to tDMRs located > 2kb from annotated TSS ##########################################
contrast <- "tissues/Human_heart_liver"
dmrs <- read.table(paste0("./DMRs/", contrast, "_DMRs.txt"), sep="\t", h=T)

features <- read.table(paste0("../annotation/DMRs/", contrast, "_DMRs_15_features_0.2.gz"), h=T, sep="\t", check.names = F)
## not same ordering
row.names(features) <- paste(features[,1], features[,2], sep=":")
features <- features[paste(dmrs[,1], dmrs[,2], sep=":"),]
summary(paste(features[,1], features[,2], sep=":") == paste(dmrs[,1], dmrs[,2], sep=":"))## now it's good

dmrs <- dmrs[features$intergenic == 1 & features$promoter == 0,]
dmrsGr <- GRanges(dmrs$chr, IRanges(dmrs$start, dmrs$end))
dim(dmrs) #4868/22932

meanMeth <- matrix(nrow=length(dmrs[,1]), ncol=48)
colnames(meanMeth) <- unique(samples$Condition)
for (sample in unique(samples$Condition)){
  cat(sample, "\n")
  ## load sample and calculate the mean methylation at each DMR location
  load(paste0("./smooth_data/smoothed_samples/", sample, "_Smoothed.RDa"))
  meanMeth[,sample] <- getMeth(data.fit, regions = dmrsGr, type = "smooth", what = "perRegion")
}

## plot heatmap
pdf(file = paste0("DMRs/", contrast, "_DMRs_2kb_upstream_heatmap.pdf"), width = 12, height = 8)
h <- heatmap.2( meanMeth, scale="none", col = colors, margins = c(12, 12), trace='none', denscol="white", labCol=paste(samples$Species, samples$Tissue, sep=" "), labRow=NA, ColSideColors=pal[as.integer(as.factor(samples$Species))], cexCol = 1.5, dendrogram='col', key=FALSE)
dev.off()
## Not very different than full set
## TO DO: 2kb threshold really meaningful? Maybe it should be larger? We should also remove conserved regions
## TO DO: subset full meanMeth object to conserved only DMRs 2kb upstream (see above)


######################################
## Correlation with Athma's results ##
######################################

athma <- read.table("../Athma/table_S1.txt", h=T, sep="\t", check.names = F)
## Data for 9911 autosomal probes

colnames(athma)
## What are the different columns meaning? What is H2, H1 and H0? Reg1? Reg2?
# In the simplest model (H0), the region’s methylation value is assumed to be constant across all three tissues, while in the second alternative (H2) the methylation value is allowed to differ between all three tissues. The first alternative (H1) models the situation where the methylation level at the site of interest is constant in the two non-target tissues but differs in the target tissue

## Good correlation with methylation levels? This would indicate that smoothing is working ##########
library(bsseq)
load("./smooth_data/combined_samples/combinedSmoothedCommonSites.RDa")
names <- paste(seqnames(granges(allData.fit.subset)), start(granges(allData.fit.subset)), sep=":")

## maping of Athma's data to our data
athma$Chr <- paste0("chr", athma$Chr)
## summary(paste(athma$Chr, athma$CGloc, sep=":") %in% names)
## ##    Mode   FALSE    TRUE    NA's 
## ## logical    9879      32       0 
## summary(paste(athma$Chr, athma$CGloc-1, sep=":") %in% names)
## ##    Mode   FALSE    TRUE    NA's
## ## logical    9893      18       0
## ## This is not a problem of 0-based or 1-based coordinates. The coordinates are actually from hg18!

## Remapping to hg19:
## ## Loading annotation for the illumina array from BioC
## source("http://bioconductor.org/biocLite.R")
## biocLite("IlluminaHumanMethylation27k.db", suppressUpdates=FALSE, lib="~/R/x86_64-redhat-linux-gnu-library/2.15/")
## ## install.packages("IlluminaHumanMethylation27k.db", repos="~/Methylation/Athma/IlluminaHumanMethylation27k.db/", lib="~/R/x86_64-redhat-linux-gnu-library/2.15/")
## ## installation didn't work :( see https://support.bioconductor.org/p/34220/
## ## But I don't want to reinstall R
## ## -> installation on local Mac was OK
## library(IlluminaHumanMethylation27k.db)
## ls("package:IlluminaHumanMethylation27k.db")
## ## Annotation here is on hg18 :(

## Use liftover tool to remap (Athma's coordinates are 1-based):
## Export to file:
## head(paste(paste(athma$Chr, athma$CGloc, sep=":"), athma$CGloc, sep="-"))
## write.table(paste(paste(athma$Chr, athma$CGloc, sep=":"), athma$CGloc, sep="-"), file = "../Athma/methylation_probes.txt", sep="\t", quote=F, row.names=F, col.names=F)
## remap file using http://www.genome.ucsc.edu/cgi-bin/hgLiftOver
## remove unmapped lines from initial table:
## chr17:59817773-59817773
## chr7:142178732-142178732
## chr17:59817613-59817613
summary(paste(athma$Chr, athma$CGloc, sep=":") %in% c("chr17:59817773", "chr7:142178732", "chr17:59817613"))
athma <- athma[!(paste(athma$Chr, athma$CGloc, sep=":") %in% c("chr17:59817773", "chr7:142178732", "chr17:59817613")),]
## Add new coordinates to table:
athma$remapped <- read.table("../Athma/methylation_probes_mapped.txt", h=F, sep="-", check.names = F)[,1]

## summary(athma$remapped %in% names)
##    Mode   FALSE    TRUE    NA's 
## logical    9904       4       0 
## Our coordinates are 0-based!

athma$remapped_0_based <- paste(unlist(lapply(strsplit(as.character(athma$remapped), split=":"), function(x) { return(x[1]) })), as.numeric(unlist(lapply(strsplit(as.character(athma$remapped), split=":"), function(x) { return(x[2]) }))) - 1, sep=":")
summary(athma$remapped_0_based %in% names)
##  Mode   FALSE    TRUE    NA's 
## logical    5063    4845       0 
## Better! (although not perfect, but here we used on 5M sites conserved HCR. Maybe Athma has sites conserved only in human/chimp)

## TO DO: maybe limited by 5M sites hsared also by macaque. This is not useful, we could take sites shared between huamdn and chimp only

athma_all <- athma ## store the full data

########### Correlation of methylation signal ##########
## Keep only these positions
athma <- athma[athma$remapped_0_based %in% names , ]
## sort by chromosome and position
athma <- athma[order(gsub('chr', '', unlist(lapply(strsplit(as.character(athma$remapped), split=":"), function(x) { return(x[1]) })), perl=TRUE), as.numeric(unlist(lapply(strsplit(as.character(athma$remapped), split=":"), function(x) { return(x[2]) })))),]
## sorting of chromosomes is not by numeric value by by characters (similar as our data)
unique(unlist(lapply(strsplit(as.character(athma$remapped), split=":"), function(x) { return(x[1]) })))
unique(seqnames(granges(allData.fit.subset)))

## check correspondance
summary(names[names %in% athma$remapped_0_based] == athma$remapped_0_based)
## filter out data
tab_corresp_athma <- getMeth(allData.fit.subset[names %in% athma$remapped_0_based, ])
tab_corresp_athma_raw <- getMeth(allData.fit.subset[names %in% athma$remapped_0_based, ], type="raw")
colnames(tab_corresp_athma) <- colnames(allData.fit.subset)
colnames(tab_corresp_athma_raw) <- colnames(allData.fit.subset)

## check smoothed and not smoothed data
## chimp heart
pdf(file = "test.pdf", width = 6, height = 5)
plot(athma$ChiHrt.mean, apply(tab_corresp_athma[,c("C1H", "C2H", "C3H", "C4H")], 1, mean), pch=16, col=rgb(0,0,0,0.1), ylab="BS-seq data", xlab="Microarray data")
plot(athma$ChiHrt.mean, apply(tab_corresp_athma_raw[,c("C1H", "C2H", "C3H", "C4H")], 1, mean), pch=16, col=rgb(0,0,0,0.1), ylab="BS-seq data", xlab="Microarray data")
dev.off()
cor.test(athma$ChiHrt.mean, apply(tab_corresp_athma[,c("C1H", "C2H", "C3H", "C4H")], 1, mean), method="spearman")
## rho = 0.7602853
cor.test(athma$ChiHrt.mean, apply(tab_corresp_athma_raw[,c("C1H", "C2H", "C3H", "C4H")], 1, mean), method="spearman")
## rho = 0.7548931 
## Quite good too!

## add density plots in margins of figure
## Adapted from function found on http://sas-and-r.blogspot.ch/2012/09/example-103-enhanced-scatterplot-with.html
pdf(file = "comparisonAthmaMicroarray_mean_ChimpHeart.pdf", width = 6, height = 5)
scatterhist(athma$ChiHrt.mean, apply(tab_corresp_athma[,c("C1H", "C2H", "C3H", "C4H")], 1, mean), ylab="BS-seq data", xlab="Microarray data", xsize=1, cleanup=F)
clip(0,1,0,1)
abline(a=0, b=1, col="darkgrey", lty=2)
lines(lowess(athma$ChiHrt.mean, apply(tab_corresp_athma[,c("C1H", "C2H", "C3H", "C4H")], 1, mean)), col="darkgrey")

scatterhist(athma$ChiHrt.mean, apply(tab_corresp_athma_raw[,c("C1H", "C2H", "C3H", "C4H")], 1, mean), ylab="BS-seq data", xlab="Microarray data", xsize=1, cleanup=F)
clip(0,1,0,1)
abline(a=0, b=1, col="darkgrey", lty=2)
## lowess was too noisy here

dev.off()


## Chimp liver
cor.test(athma$ChiLiv.mean, apply(tab_corresp_athma[,c("C1Li", "C2Li", "C3Li", "C4Li")], 1, mean), method="spearman")
## rho = 0.7738892 
cor.test(athma$ChiLiv.mean, apply(tab_corresp_athma_raw[,c("C1Li", "C2Li", "C3Li", "C4Li")], 1, mean), method="spearman")
## rho = 0.7722345
pdf(file = "comparisonAthmaMicroarray_mean_ChimpLiver.pdf", width = 6, height = 5)
scatterhist(athma$ChiLiv.mean, apply(tab_corresp_athma[,c("C1Li", "C2Li", "C3Li", "C4Li")], 1, mean), ylab="BS-seq data", xlab="Microarray data", xsize=1, cleanup=F)
clip(0,1,0,1)
abline(a=0, b=1, col="darkgrey", lty=2)
lines(lowess(athma$ChiLiv.mean, apply(tab_corresp_athma[,c("C1Li", "C2Li", "C3Li", "C4Li")], 1, mean)), col="darkgrey")
scatterhist(athma$ChiLiv.mean, apply(tab_corresp_athma_raw[,c("C1Li", "C2Li", "C3Li", "C4Li")], 1, mean), ylab="BS-seq data", xlab="Microarray data", xsize=1, cleanup=F)
clip(0,1,0,1)
abline(a=0, b=1, col="darkgrey", lty=2)
dev.off()


## Chimp kidney
cor.test(athma$ChiKid.mean, apply(tab_corresp_athma[,c("C1K", "C2K", "C3K", "C4K")], 1, mean), method="spearman")
## rho = 0.7828201
cor.test(athma$ChiKid.mean, apply(tab_corresp_athma_raw[,c("C1K", "C2K", "C3K", "C4K")], 1, mean), method="spearman")
## rho = 0.7831497
pdf(file = "comparisonAthmaMicroarray_mean_ChimpKidney.pdf", width = 6, height = 5)
scatterhist(athma$ChiKid.mean, apply(tab_corresp_athma[,c("C1K", "C2K", "C3K", "C4K")], 1, mean), ylab="BS-seq data", xlab="Microarray data", xsize=1, cleanup=F)
clip(0,1,0,1)
abline(a=0, b=1, col="darkgrey", lty=2)
lines(lowess(athma$ChiKid.mean, apply(tab_corresp_athma[,c("C1K", "C2K", "C3K", "C4K")], 1, mean)), col="darkgrey")
scatterhist(athma$ChiKid.mean, apply(tab_corresp_athma_raw[,c("C1K", "C2K", "C3K", "C4K")], 1, mean), ylab="BS-seq data", xlab="Microarray data", xsize=1, cleanup=F)
clip(0,1,0,1)
abline(a=0, b=1, col="darkgrey", lty=2)
dev.off()


## Human heart
cor.test(athma$HumHrt.mean, apply(tab_corresp_athma[,c("H1H", "H2H", "H3H", "H4H")], 1, mean), method="spearman")
## rho = 0.7565347 
cor.test(athma$HumHrt.mean, apply(tab_corresp_athma_raw[,c("H1H", "H2H", "H3H", "H4H")], 1, mean), method="spearman")
## rho = 0.7733894
pdf(file = "comparisonAthmaMicroarray_mean_HumanHeart.pdf", width = 6, height = 5)
scatterhist(athma$HumHrt.mean[athma$HumHrt.mean > 0], apply(tab_corresp_athma[athma$HumHrt.mean > 0 ,c("H1H", "H2H", "H3H", "H4H")], 1, mean), ylab="BS-seq data", xlab="Microarray data", xsize=1, cleanup=F)
clip(0,1,0,1)
abline(a=0, b=1, col="darkgrey", lty=2)
lines(lowess(athma$HumHrt.mean, apply(tab_corresp_athma[,c("H1H", "H2H", "H3H", "H4H")], 1, mean)), col="darkgrey")
scatterhist(athma$HumHrt.mean[athma$HumHrt.mean > 0], apply(tab_corresp_athma_raw[athma$HumHrt.mean > 0, c("H1H", "H2H", "H3H", "H4H")], 1, mean), ylab="BS-seq data", xlab="Microarray data", xsize=1, cleanup=F)
clip(0,1,0,1)
abline(a=0, b=1, col="darkgrey", lty=2)
dev.off()
## Strange: there are some negativemethylation values! Thi screws up the histogram plots

## Human liver
cor.test(athma$HumLiv.mean[athma$HumLiv.mean > 0], apply(tab_corresp_athma[athma$HumLiv.mean > 0, c("H1Li", "H2Li", "H3Li", "H4Li")], 1, mean), method="spearman")
## rho = 0.8069868
cor.test(athma$HumLiv.mean[athma$HumLiv.mean > 0], apply(tab_corresp_athma_raw[athma$HumLiv.mean > 0, c("H1Li", "H2Li", "H3Li", "H4Li")], 1, mean), method="spearman")
## rho = 0.7884884
pdf(file = "comparisonAthmaMicroarray_mean_HumanLiver.pdf", width = 6, height = 5)
scatterhist(athma$HumLiv.mean[athma$HumLiv.mean > 0], apply(tab_corresp_athma[athma$HumLiv.mean > 0, c("H1Li", "H2Li", "H3Li", "H4Li")], 1, mean), ylab="BS-seq data", xlab="Microarray data", xsize=1, cleanup=F)
clip(0,1,0,1)
abline(a=0, b=1, col="darkgrey", lty=2)
lines(lowess(athma$HumLiv.mean, apply(tab_corresp_athma[,c("H1Li", "H2Li", "H3Li", "H4Li")], 1, mean)), col="darkgrey")
scatterhist(athma$HumLiv.mean[athma$HumLiv.mean > 0], apply(tab_corresp_athma_raw[athma$HumLiv.mean > 0, c("H1Li", "H2Li", "H3Li", "H4Li")], 1, mean), ylab="BS-seq data", xlab="Microarray data", xsize=1, cleanup=F)
clip(0,1,0,1)
abline(a=0, b=1, col="darkgrey", lty=2)
dev.off()


## Human kidney
cor.test(athma$HumKid.mean, apply(tab_corresp_athma[,c("H1K", "H2K", "H3K", "H4K")], 1, mean), method="spearman")
## rho = 0.7807188
cor.test(athma$HumKid.mean, apply(tab_corresp_athma_raw[,c("H1K", "H2K", "H3K", "H4K")], 1, mean), method="spearman")
## rho =0.790441
pdf(file = "comparisonAthmaMicroarray_mean_HumanKidney.pdf", width = 6, height = 5)
scatterhist(athma$HumKid.mean, apply(tab_corresp_athma[,c("H1K", "H2K", "H3K", "H4K")], 1, mean), ylab="BS-seq data", xlab="Microarray data", xsize=1, cleanup=F)
clip(0,1,0,1)
abline(a=0, b=1, col="darkgrey", lty=2)
lines(lowess(athma$HumKid.mean, apply(tab_corresp_athma[,c("H1K", "H2K", "H3K", "H4K")], 1, mean)), col="darkgrey")
scatterhist(athma$HumKid.mean, apply(tab_corresp_athma_raw[,c("H1K", "H2K", "H3K", "H4K")], 1, mean), ylab="BS-seq data", xlab="Microarray data", xsize=1, cleanup=F)
clip(0,1,0,1)
abline(a=0, b=1, col="darkgrey", lty=2)
dev.off()


########### Correlation of tissue DMRs locations #############
## We can use p-value H1 vs. H0 (tissue-specific)
## We could also filter H2 vs. H1 low p-values (not tissue-specific).
## Unfortunately we have no indication on which tissue differ from the others
athma <- athma_all
summary(athma$HumLivTDMR_H1H0_pvalue < 0.001)
##   Mode   FALSE    TRUE    NA's 
## logical    7926     1982       0 
summary(athma$HumLivTDMR_H1H0_pvalue < 0.001 & athma$HumLivTDMR_H2H1_pvalue < 0.05)
##   Mode   FALSE    TRUE    NA's 
## logical    8986     922       0 

load("./DMRs/tissues/Human_liverSpecific_tstat.RDa") ## we only consider the score of the tests for tissue of interest vs. all other tissues. Comparisons between other tissues (which had to be NS) are ignored here
names.tstat <- paste(seqnames(granges(allData.tstat)), start(granges(allData.tstat)), sep=":")
length(names.tstat)

summary(athma$remapped_0_based %in% names.tstat) ## all but 1 is true
athma <- athma[athma$remapped_0_based %in% names.tstat, ]

## sort by chromosome and position
athma <- athma[order(gsub('chr', '', unlist(lapply(strsplit(as.character(athma$remapped), split=":"), function(x) { return(x[1]) })), perl=TRUE), as.numeric(unlist(lapply(strsplit(as.character(athma$remapped), split=":"), function(x) { return(x[2]) })))),]
## sorting of chromosomes is not by numeric value by by characters (similar as our data)

## there should be no need to sort data: check correspondance
summary(names.tstat[names.tstat %in% athma$remapped_0_based] == athma$remapped_0_based) ## OK

## Append the t-stat for these positions
athma$liver.specific.tstat <- allData.tstat@stats[names.tstat %in% athma$remapped_0_based, "tstat.corrected"]


## test the agreement of t-stats and microarray p-values
## overall the agreement is low
cor.test(abs(athma$liver.specific.tstat), athma$HumLivTDMR_H1H0_pvalue, method="spearman")
## rho = -0.25
cor.test(abs(athma$liver.specific.tstat[athma$HumLivTDMR_H1H0_pvalue < 0.01]), athma$HumLivTDMR_H1H0_pvalue[athma$HumLivTDMR_H1H0_pvalue < 0.01], method="spearman")
## rho = -0.33

## plots:
pdf(file = "comparisonAthmaMicroarray_DMRs_liverSpecific.pdf", width = 6, height = 5)
## Absolute t-stat vs. p-value
plot(abs(athma$liver.specific.tstat), athma$HumLivTDMR_H1H0_pvalue, pch=16, col=rgb(0,0,0,0.1), xlab="BS-seq tstat", ylab="Microarray p-value", log="y")
## hard to read: too many points. Also because of log scale, all what we see is significant for microarray

## not on log scale
plot(athma$liver.specific.tstat, athma$HumLivTDMR_H1H0_pvalue, pch=16, col=rgb(0,0,0,0.1), xlab="BS-seq tstat", ylab="Microarray p-value", log="")
plot(abs(athma$liver.specific.tstat), athma$HumLivTDMR_H1H0_pvalue, pch=16, col=rgb(0,0,0,0.1), xlab="BS-seq tstat", ylab="Microarray p-value", log="")
## Here we see an enrichment of low p-values for high t-stat: good

## Remove H2H1 significant
plot(athma$liver.specific.tstat[athma$HumLivTDMR_H2H1_pvalue > 0.05], athma$HumLivTDMR_H1H0_pvalue[athma$HumLivTDMR_H2H1_pvalue > 0.05], pch=16, col=rgb(0,0,0,0.1), xlab="BS-seq tstat", ylab="Microarray p-value", log="")
## not changing much

##  Remove H2H1 significant + only H1H0 significant at FDR 1% (not same criteria in Athma's paper)
plot(athma$liver.specific.tstat[athma$HumLivTDMR_H2H1_pvalue > 0.05 & p.adjust(athma$HumLivTDMR_H1H0_pvalue, method="BH") < 0.01], athma$HumLivTDMR_H1H0_pvalue[athma$HumLivTDMR_H2H1_pvalue > 0.05 & p.adjust(athma$HumLivTDMR_H1H0_pvalue, method="BH") < 0.01], pch=16, col=rgb(0,0,0,0.1), xlab="BS-seq tstat", ylab="Microarray p-value", log="")
## y-axis in log scale
plot(athma$liver.specific.tstat[athma$HumLivTDMR_H2H1_pvalue > 0.05 & p.adjust(athma$HumLivTDMR_H1H0_pvalue, method="BH") < 0.01], athma$HumLivTDMR_H1H0_pvalue[athma$HumLivTDMR_H2H1_pvalue > 0.05 & p.adjust(athma$HumLivTDMR_H1H0_pvalue, method="BH") < 0.01], pch=16, col=rgb(0,0,0,0.1), xlab="BS-seq tstat", ylab="Microarray p-value", log="y")
dev.off()
## still an enrichment for very low p-values
## The last 3 plots are used as sup figures!

## In conclusion, all large tstats are associated with low p-values, but not all low p-values associated with large t-stats
## 


## Other way to ask the question:
## What is the proportion of CpG sites with low p-values that fall into one of our DMRs?
summary(athma$HumLivTDMR_H2H1_pvalue > 0.05 & p.adjust(athma$HumLivTDMR_H1H0_pvalue, method="BH") < 0.001) ## 840 CpGs ## Note: In Athma's paper, the correction used is Q-value
summary(athma$HumLivTDMR_H2H1_pvalue > 0.05 & p.adjust(athma$HumLivTDMR_H1H0_pvalue, method="BH") < 0.01) ## 1225
dmrs <- read.table("./DMRs/tissues/Human_liverSpecific_DMRs.txt", h=T, sep="\t")
## dmrs <- read.table("./DMRs/tissues/Human_heart_liver_DMRs.txt", h=T, sep="\t") ## not better

athma_significant <- athma[athma$HumLivTDMR_H2H1_pvalue > 0.05 & p.adjust(athma$HumLivTDMR_H1H0_pvalue, method="BH") < 0.001, ]

inside_dmr <- apply(athma_significant, 1, function(x) {
  chr <- unlist(strsplit(x[38], split=":"))[1]
  start <- as.numeric(unlist(strsplit(x[38], split=":"))[2])
  return(nrow(dmrs[dmrs$chr == chr & dmrs$start <= start & dmrs$end >= start,]))
})
summary(as.factor(inside_dmr))
##  0   1 
##751  89 
## Pretty low (10.6%), meaning that many significant CpGs from Athma are not within a DMR
## As a control, the significant CpGs overlap only 1 lung-specific and 0 heart-specific DMRs
## The rate is lower for heart-specific and similar for kidney-specific

## FDR 1% instead of 0.1%
athma_significant <- athma[athma$HumLivTDMR_H2H1_pvalue > 0.05 & p.adjust(athma$HumLivTDMR_H1H0_pvalue, method="BH") < 0.01, ]
## ... same code as above
summary(as.factor(inside_dmr))
##   0    1
##1131   94
## 7.7%: still pretty low!

## Compare to random DMR file(s)
## make vector of overlaps with random DMR files (length 100)
## repeat the test
athma_significant <- athma[athma$HumLivTDMR_H2H1_pvalue > 0.05 & p.adjust(athma$HumLivTDMR_H1H0_pvalue, method="BH") < 0.01, ]
all_random <- vector(length=100)
for (i in 1:100){
  print(i)
  dmrs <- read.table(paste0("../annotation/DMRs/tissues/randomized/Human_liverSpecific_randomDMRs_", i, ".bed"), h=F, sep="\t")
  names(dmrs) <- c("chr", "start", "end")
  ## Be careful with coordinates! End is 1-based here
  dmrs$end <- dmrs$end-1
  
  inside_dmr <- apply(athma_significant, 1, function(x) {
    chr <- unlist(strsplit(x[38], split=":"))[1]
    start <- as.numeric(unlist(strsplit(x[38], split=":"))[2])
    return(nrow(dmrs[dmrs$chr == chr & dmrs$start <= start & dmrs$end >= start,]))
  })
  all_random[i] <- sum(inside_dmr != 0)
}
## control for CpG density
all_random_density <- vector(length=100)
for (i in 1:100){
  print(i)
  dmrs <- read.table(paste0("../annotation/DMRs/tissues/randomized_control_CpG_density/Human_liverSpecific_randomDMRs_", i, ".bed"), h=F, sep="\t")
  names(dmrs) <- c("chr", "start", "end")
  ## Be careful with coordinates! End is 1-based here
  dmrs$end <- dmrs$end-1
  
  inside_dmr <- apply(athma_significant, 1, function(x) {
    chr <- unlist(strsplit(x[38], split=":"))[1]
    start <- as.numeric(unlist(strsplit(x[38], split=":"))[2])
    return(nrow(dmrs[dmrs$chr == chr & dmrs$start <= start & dmrs$end >= start,]))
  })
  all_random_density[i] <- sum(inside_dmr != 0)
}

## Boxplot
pdf(file = "comparisonAthmaMicroarray_DMRs_liverSpecific_overlap.pdf", width = 3, height = 5)
## draw empty plot to be able to add backrgoudn grey lines below boxes
par(mar=c(5, 5, 2, 2))
plot(1:2, type="n", axes=F, xlab="", ylab="", ylim=c(0,100), yaxs="i", xlim=c(0,3))
## add vertical gray dashed lines to help reading
abline(v=c(1:2), col="gray90", lty=3)
## add boxplot
mp <- boxplot(all_random, all_random_density, pch=20, cex=0.8, notch=F, ylab="Number of significant CpGs overlapping DMRs", xlab="", xaxt="n", ylim=c(0,100), yaxs="i", whisklty=1, medlty=0, col="dimgray", border="dimgray", add=T, las=2)
## text 6% below x-axis 
axis(1, labels=NA, at=c(1:2))
# text(c(1:2), par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.06, srt = 45, adj = 1, labels = c("", "+ control for CpG density"), xpd = TRUE, cex=1)
text(2, par("usr")[3]-(par("usr")[4]-par("usr")[3])*0.1, adj = 0.5, labels = "Control for\nCpG density", xpd = TRUE, cex=0.8)
## Add real overlap values
points(c(1,2), c(94,94), col=pal[1], pch=18, cex=1.3)
dev.off()
## Added as sup figure in paper


## include 1kb before and after DMR
athma_significant <- athma[athma$HumLivTDMR_H2H1_pvalue > 0.05 & p.adjust(athma$HumLivTDMR_H1H0_pvalue, method="BH") < 0.01, ]
inside_dmr <- apply(athma_significant, 1, function(x) {
  chr <- unlist(strsplit(x[38], split=":"))[1]
  start <- as.numeric(unlist(strsplit(x[38], split=":"))[2])
  return(nrow(dmrs[dmrs$chr == chr & (dmrs$start-1000) <= start & (dmrs$end+1000) >= start,]))
})
summary(as.factor(inside_dmr))
## 0   1   2   3 
## 994 201  27   3 
## Of course they can overlap multiple DMRs
## 231/1224 = 18.9%

## Take both tests into account H1H0 or H2H1 and compare to all pairwise DMRs involving liver
athma_significant <- athma[p.adjust(athma$HumLivTDMR_H2H1_pvalue, method="BH") < 0.001 | p.adjust(athma$HumLivTDMR_H1H0_pvalue, method="BH") < 0.001, ] ## 2340
dmrs <- rbind(read.table("./DMRs/tissues/Human_heart_liver_DMRs.txt", h=T, sep="\t"), read.table("./DMRs/tissues/Human_liver_lung_DMRs.txt", h=T, sep="\t"), read.table("./DMRs/tissues/Human_liver_kidney_DMRs.txt", h=T, sep="\t")) ## 54204, some dmrs are redundant

inside_dmr <- apply(athma_significant, 1, function(x) {
  chr <- unlist(strsplit(x[38], split=":"))[1]
  start <- as.numeric(unlist(strsplit(x[38], split=":"))[2])
  return(nrow(dmrs[dmrs$chr == chr & dmrs$start <= start & dmrs$end >= start,]))
})
summary(inside_dmr > 0)
## 381/2340 = 0.1628205
## This is better


## For all our DMRs, look for Athma's sites p-values
dmrs <- read.table("./DMRs/tissues/Human_liverSpecific_DMRs.txt", h=T, sep="\t")
inside_dmr <- apply(athma, 1, function(x) {
  chr <- unlist(strsplit(x[38], split=":"))[1]
  start <- as.numeric(unlist(strsplit(x[38], split=":"))[2])
  return(nrow(dmrs[dmrs$chr == chr & dmrs$start <= start & dmrs$end >= start,]))
})
summary(as.factor(inside_dmr))
## 205 = 2.1%. This is lower than for Athma's significant sites only: good!

summary(p.adjust(athma$HumLivTDMR_H1H0_pvalue, method="BH")[inside_dmr != 0] < 0.001)
## 177/205 = 86.3%: almost all sites inside our DMRs are confirmed by Athma!
summary(p.adjust(athma$HumLivTDMR_H1H0_pvalue, method="BH")[inside_dmr != 0] < 0.01)
## 190/205 = 92.6%: This number used in paper
summary(athma$HumLivTDMR_H2H1_pvalue[inside_dmr != 0] > 0.05 & p.adjust(athma$HumLivTDMR_H1H0_pvalue, method="BH")[inside_dmr != 0] < 0.001)
## 89/205 = 43.4% also rejected H2
summary(athma$HumLivTDMR_H2H1_pvalue[inside_dmr != 0] > 0.05 & p.adjust(athma$HumLivTDMR_H1H0_pvalue, method="BH")[inside_dmr != 0] < 0.01)
## 94/205 = 45.9% also rejected H2. This number used in paper
summary(p.adjust(athma$HumLivTDMR_H2H1_pvalue, method="BH")[inside_dmr != 0] < 0.001 | p.adjust(athma$HumLivTDMR_H1H0_pvalue, method="BH")[inside_dmr != 0] < 0.001)
## 89.8% Both significant

## TO DO: Repeat for other tissues 
## athma$HumKidTDMR_H1H0_pvalue
## athma$HumHrtTDMR_H1H0_pvalue
## athma$ChiLivTDMR_H1H0_pvalue
## athma$ChiKidTDMR_H1H0_pvalue
## athma$ChiHrtTDMR_H1H0_pvalue

## TO DO: Species DMRs: where is the information? Doesn't seem to be calculated by Athma...

## TO DO: correlations on exact same samples used in both studies (C1 to C4 hearts were used in Athma's paper)
athma2 <- read.table("../Athma/CompleteDataSet_Auto.txt", h=T, sep="\t", check.names = F)
colnames(athma2)
## Here there is the methylation level for each sample
## 6 samples + 2 technical replicates per sample.
infoAthma <- as.factor(c(rep("HumLiv", times=12), rep("ChiLiv", times=12), rep("HumKid", times=12), rep("ChiKid", times=11), rep("HumHrt", times=12), rep("ChiHrt", times=12)))

tapply(as.numeric(athma2[1,14:84]), list(infoAthma), mean)
## This is similar to: 
athma[1,7:12]

## TO DO: for a human tissue, do the correlation on the whole set of CpG on the chip
##        Does smoothing strenghten the correlation? On the human+chimp CpGs this is not really the case

############### lower mean methylation in chimps compared to human? #######################
athma <- read.table("../Athma/table_S1.txt", h=T, sep="\t", check.names = F)
## Data for 9911 autosomal probes

pdf(file = "distributionMethylation_Athma.pdf", height = 6, width=6)
for (species in c("Hum", "Chi")){
  i <- 0
  j <- 1
  for (tissue in c("Hrt", "Kid", "Liv")){
    if (i == 0){
      plot(density(athma[, paste0(species, tissue, '.mean')]), col=pal[j], main=species, xlab="Fraction of methylated DNA", lwd=2, lty=1, xlim=c(-0.1, 1.1))
      i <- 1
    }
    else {
      lines(density(athma[, paste0(species, tissue, '.mean')]), col=pal[j], lwd=2, lty=1)
    }
    j <- j+1
  }
  legend("top", legend=c("Hrt", "Kid", "Liv"), lty=c(1,1,1), col=c(pal[1:3])) 
}
## Both species together
i <- 0
k <- 1
for (species in c("Hum", "Chi")){
  j <- 1
  for (tissue in c("Hrt", "Kid", "Liv")){
    if (i == 0){
      plot(density(athma[, paste0(species, tissue, '.mean')]), col=pal[j], main="", xlab="Fraction of methylated DNA", lwd=2, lty=k, xlim=c(-0.1, 1.1))
      i <- 1
    }
    else {
      lines(density(athma[, paste0(species, tissue, '.mean')]), col=pal[j], lwd=2, lty=k)
    }
    j <- j+1
  }
  k <- k+1
  legend("top", legend=c("Heart", "Kidney", "Liver", "Human", "Chimp"), lty=c(1,1,1,1,2), col=c(pal[1:3], "black", "black"), lwd=2) 
}
dev.off()
## We see a slight difference here too!

## Redo with complete dataset: athma2
dataAthma2 <- athma2[,14:84]

infoAthma <- data.frame(condition=as.factor(c(rep("HumLiv", times=12), rep("ChiLiv", times=12), rep("HumKid", times=12), rep("ChiKid", times=11), rep("HumHrt", times=12), rep("ChiHrt", times=12))))
infoAthma$tissue <- as.factor(c(rep("liver", times=24), rep("kidney", times=23), rep("heart", times=24)))
infoAthma$species <- as.factor(c(rep("Human", times=12), rep("Chimp", times=12), rep("Human", times=12), rep("Chimp", times=11), rep("Human", times=12), rep("Chimp", times=12)))
infoAthma$individual <- c(1,1,2,2,3,3,4,4,5,5,6,6,1,1,2,2,3,3,4,4,5,5,6,6,1,1,2,2,3,3,4,4,5,5,6,6,1,1,2,3,3,4,4,5,5,6,6,1,1,2,2,3,3,4,4,5,5,6,6,1,1,2,2,3,3,4,4,5,5,6,6)
row.names(infoAthma) <- colnames(dataAthma2)

pdf(file = "distributionMethylation_Athma_full_by_tissues.pdf", height = 6, width=18)
par(mfrow=c(1,3), ps = 14, cex = 1, cex.main = 1) ## Used to keep labels at Font 14pt
for (tissue in sort(unique(infoAthma$tissue))){
  i <- 0
  for (species in sort(unique(infoAthma$species))){
    for (sample in colnames(dataAthma2[, infoAthma$species == species & infoAthma$tissue == tissue])){
      cat(sample, "\n")
      if (i == 0){
        plot(density(dataAthma2[,sample]), col=pal[as.integer(as.factor(infoAthma[sample,]$species))+9], main=tissue, xlab=expression("Methylation percentage (" ~ beta ~ ")"), lwd=1.5, lty=infoAthma[sample,]$individual, xlim=c(-0.1,1.1), ylim=c(0,5.5))
        i <- 1
      }
      else {
        lines(density(dataAthma2[,sample]), col=pal[as.integer(as.factor(infoAthma[sample,]$species))+9], lwd=1.5, lty=infoAthma[sample,]$individual)
      }
    }
  }
  if (tissue==sort(unique(infoAthma$tissue))[1]){
    legend("topright", legend=c(as.character(sort(unique(infoAthma$species))), paste0("Individual ", sort(unique(infoAthma$individual))), "(2 technical replicates each)"), lty=c(1,1, 1:6, NA), lwd=1.5, col=c(pal[10:11], rep("black", 6), NA), bty = "n") 
  }
}
dev.off()
## Except kidney, we see a very similar pattern as in BS-seq data!
## -> Methylation structure is different in chimp and human. Probably related to CpG density, etc

#########################################
## Correlation with Irizarry's results ##
#########################################
## See ~/Methylation/Irizarry/

#############################
##     That all folks!     ##
#############################
