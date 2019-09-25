## Feb 28, 2014
## QC
## Main task: perform some PCA/clustering on reassembled smoothed samples
## This is done on the cluster (32g memory seems OK)
setwd("~/Methylation/bsseq/")
source("functions.R")

library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))

## source("http://www.bioconductor.org/biocLite.R")
## biocLite()
## biocLite("BiocUpgrade")
## biocLite("bsseq")
## biocLite("bsseq", lib="~/R/x86_64-redhat-linux-gnu-library/2.15/") ## on the cluster
library(bsseq) # version 0.10
sessionInfo()

###############
## Smoothing ##
###############
## correlation of smoothed values for a sample smoothed alone, or smoothed at all sites obesrevd in all samples
## Do it on sample with lowest coverage
load("smooth_data/smoothed_samples/H1K_Smoothed.RDa")
data.alone.fit <- data.fit
load("smooth_data/combined_samples/separated_samples/H1K_Smoothed.RDa")
data.fit.subset <- subsetByOverlaps(data.fit, granges(data.alone.fit))
tab.alone <- getMeth(data.alone.fit, type = "smooth")
tab.combined <- getMeth(data.fit.subset, type = "smooth")
summary(tab.alone[,1] == tab.combined[,1]) ## Not a single difference!

#########
## PCA ##
#########
load("smooth_data/combined_samples/combinedSmoothedCommonSites.RDa")
info <- pData(allData.fit.subset)

## Correct swapping of Rhesus individuals
## R3Li <-> R2Li, R3Lu <-> R2Lu
temp <- info["R2Li", ]
info["R2Li", ] <- info["R3Li", ]
info["R3Li", ] <- temp
row.names(info)[43] <- "R2Li2"
row.names(info)[39] <- "R3Li"
row.names(info)[43] <- "R2Li"

temp <- info["R2Lu", ]
info["R2Lu", ] <- info["R3Lu", ]
info["R3Lu", ] <- temp
row.names(info)[44] <- "R2Lu2"
row.names(info)[40] <- "R3Lu"
row.names(info)[44] <- "R2Lu"

sampleNames(allData.fit.subset) <- sampleNames(allData.fit.subset)[c(1:38, 43:44, 41:42, 39:40, 45:48)]
pData(allData.fit.subset) <- info

## perform PCA on coverage data ########################################################
tab <- getCoverage(allData.fit.subset)
row.names(tab) <- paste(seqnames(granges(allData.fit.subset)), start(granges(allData.fit.subset)), sep=":")
colnames(tab) <- sampleNames(allData.fit.subset)
dim(tab) ## 5,228,469 sites
print(head(tab))
## perform PCA 
pca1 <- prcomp(t(tab), scale = FALSE) ## turn off the argument "scale".
print(summary(pca1)) ## 33% variance explained by 3 first components
## TO DO scaled PCA too?

## loadings: coordinates of PCs in original coordinate system
loadings <- pca1$rotation
## scores: coordinates of scaled data (observations) projected onto the PCs.
scores <- pca1$x
print(head(scores))
rownames(scores) <- colnames(tab)

## biplot: PCs as axes, scores as dots, original variables as arrows. Angles between variables reflect if they are correlated. But there are too mnay variables, so here we just plot the scores
n <- 20
pdf(file = paste0("PCA_PCs1-", n+1, "_coverage.pdf"), width = 6, height = 6)
for (i in 1:n){
  ## colors for tissues
  ## symbols for species
  col.v <- pal[as.integer(info$Tissue)]
  pchs <- c(15,16,17)[as.integer(info$Species)]
  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=F,legend=F)
}
dev.off()
pdf(file = paste0("PCA_PCs1-", n+1, "_coverage_symbols.pdf"), width = 6, height = 6)
for (i in 1:n){
  ## colors for tissues
  ## symbols for species
  col.v <- pal[as.integer(info$Tissue)]
  pchs <- c(15,16,17)[as.integer(info$Species)]
  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=T,legend=T)
  legend("topright", c("Human","Chimpanzee", "Macaque"), pch=c(16,15,17), cex=1.3, col="black")
  legend("bottomright", c("Heart","Kidney", "Liver", "Lung"), pch=16, cex=1.3, col=pal[1:4])
}
dev.off()
## No obvious pattern but there is some tissue/species patterns:
## PC1: maybe liver+heart samples a bit more on right, lung+kidney slightly more on the left
## PC2: rhesus a bit more at bottom
## PC2+3: isolates macaque
## PC4: isolates human and chimp


## perform PCA on smooted methylation estimates ########################################
## ## get the smoothed methylation estimates
## tab <- getMeth(allData.fit.subset, type = "smooth")
## row.names(tab) <- paste(seqnames(granges(allData.fit.subset)), start(granges(allData.fit.subset)), sep=":")
## colnames(tab) <- sampleNames(allData.fit.subset)
## sds <- apply(tab, 1, sd) ## long!
## tab <- tab[sds > 0.01, ] ## filter rows with low variance. 
## tab <- na.omit(tab)
## dim(tab) ## 4,879,356 sites left
## print(head(tab))
## save(tab, file="PCA_combinedSmoothedCommonSites_table.RDa")
load("PCA_combinedSmoothedCommonSites_table.RDa")

## Data are not normally distributed: should we use another measure?
## For example M-values (Beta/1-Beta)
tabMvalues <- log10(tab / (1 - tab))
plot(density(tabMvalues[,1]))
## But the M-values distribution is still bimodal!

## Centering: substract the mean methylation on each row.
tabDiff <- t(apply(tab, 1, function(x){ return(x-mean(x)) }))
plot(density(tab[,1]))
## This is unimodal!
## But actually it is equivalent to using center=T in prcomp: see below

## perform PCA 
pca1 <- prcomp(t(tab), scale = FALSE) ## turn off the argument "scale".
## Here, unlike gene expression, the magnitude of %methylation change is varying only between 0 and 1. Scaling would equate a change of 10% methylation for 1 locus with a change of 50% methylation in another locus. Unscaling guarantees that the second contributes more.
print(summary(pca1)) ## 49.7% variance explained by 3 first components

## loadings: coordinates of PCs in original coordinate system
loadings <- pca1$rotation
## scores: coordinates of scaled data (observations) projected onto the PCs.
scores <- pca1$x
print(head(scores))
rownames(scores) <- colnames(tab)

## biplot: PCs as axes, scores as dots, original variables as arrows. Angles between variables reflect if they are correlated. But there are too mnay variables, so here we just plot the scores
n <- 20
pdf(file = paste0("PCA_PCs1-", n+1, "_smoothed.pdf"), width = 6, height = 6)
## pdf(file = paste0("PCA_PCs1-", n+1, "_smoothed_symbols.pdf"), width = 6, height = 6)
for (i in 1:n){
  ## colors for tissues
  ## symbols for species
  col.v <- pal[as.integer(info$Tissue)]
  pchs <- c(15,16,17)[as.integer(info$Species)]
  ##plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=T,legend=T)
  ##legend("topright", c("Human","Chimpanzee", "Macaque"), pch=c(16,15,17), cex=1.3, col="black")
  ##legend("bottomright", c("Heart","Kidney", "Liver", "Lung"), pch=16, cex=1.3, col=pal[1:4])

  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=F,legend=F)
}
dev.off()
## - The clustering is first by species, and then by tissue, opposite to RNA-seq.
## - liver seems to be the most different tissue: similar to Athma's results!
## - doubts on R3H since the mapping efficiency was a lot higher than the other macaque samples, but it groups very nicely with the other hearts
## - no PC correlating with read coverage: cool! (see below)
## - other technical factors?
## - Some of the deeper PCs seem to isolate individuals. For example, PC7 (Chimp3), PC13 (Rhesus1 and Rhesus4), PC18 (Chimp1 and Chimp2), PC20 (Chimp 4), PC21? (Human3 and Human4). Human individuals seem to be harder to separate than chimp and rhesus individuals. Probably there is less diversity of methylation signal in humans (similarly to lower genetic diversity)?
## - I am wondering if there is a mixing of individuals between Rhesus2 and Rhesus3: see PC7 and it's even clearer on PC16. This will be clarified by the SNP calling. 

## correlation of PCs to technical factors: see below

## Remove X chromosome
tabNoX <- tab[!grepl("^chrX", row.names(tab), perl=T),]
dim(tabNoX) ## 4,627,805 sites
pca1 <- prcomp(t(tabNoX), scale = FALSE) 
print(summary(pca1)) ## 51% variance explained by 3 first components
loadings <- pca1$rotation
scores <- pca1$x
rownames(scores) <- colnames(tab)
# plot
n <- 20
pdf(file = paste0("PCA_PCs1-", n+1, "_smoothed_noX.pdf"), width = 6, height = 6)
## pdf(file = paste0("PCA_PCs1-", n+1, "_smoothed_noX_symbols.pdf"), width = 6, height = 6)
for (i in 1:n){
  ## colors for tissues
  ## symbols for species
  col.v <- pal[as.integer(info$Tissue)]
  pchs <- c(15,16,17)[as.integer(info$Species)]
  ##plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=T,legend=T)
  ##legend("topright", c("Human","Chimpanzee", "Macaque"), pch=c(16,15,17), cex=1.3, col="black")
  ##legend("bottomright", c("Heart","Kidney", "Liver", "Lung"), pch=16, cex=1.3, col=pal[1:4])

  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=F,legend=F)
}
dev.off()

## Figure 2A
pdf(file = "PCA_PCs1-2_smoothed_noX_symbols.pdf", width = 6, height = 6)
col.v <- pal[as.integer(info$Tissue)]
col.v <- addalpha(col.v, c(0.65,1,0.4)[as.integer(info$Species)])
pchs <- c(15,16,17)[as.integer(info$Species)]
plot_scores(pca1, scores, 1, 2, cols=col.v, pch=pchs, points=T,legend=F)
## legend("topright", c("Human","Chimpanzee", "Macaque"), pch=c(16,15,17), cex=1.3, col=addalpha(rep("black", times=3), c(1,0.65,0.4)))
## legend("bottomright", c("Heart","Kidney", "Liver", "Lung"), pch=16, cex=1.3, col=pal[1:4])
dev.off()


## just include human and chimp samples:
tabHC <- tab[, grepl("^C", colnames(tab), perl=T) | grepl("^H", colnames(tab), perl=T)]
dim(tabHC)
print(head(tabHC))
pca1 <- prcomp(t(tabHC), scale = FALSE) ## turn off the argument "scale".
print(summary(pca1)) ## 50.5% variance explained by 3 first components
loadings <- pca1$rotation
scores <- pca1$x
rownames(scores) <- colnames(tabHC)
## plot
n <- 20
pdf(file = paste0("PCA_PCs1-", n+1, "_smoothed_human_chimp.pdf"), width = 6, height = 6)
for (i in 1:n){
  ## colors for tissues
  ## symbols for species
  col.v <- pal[as.integer(info$Tissue)]
  pchs <- c(15,16,17)[as.integer(info$Species)]
  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=F,legend=F)
}
dev.off()
## PCs 1,2,3: species and tissues
## PC4: lung
## PC5: C3
## PC11: C1 and C2
## PC12: C4?
## PC13: H4, H3

## TO DO: do it on the whole set of shared CpG sites between human and chimp (a lot more than 4M!)

## just include human samples: do we cluster individuals more easily?
tabH <- tab[, grepl("^H", colnames(tab), perl=T)]
dim(tabH)
print(head(tabH)) 
pca1 <- prcomp(t(tabH), scale = FALSE) ## turn off the argument "scale".
print(summary(pca1)) ## 58.1% variance explained by 3 first components
loadings <- pca1$rotation
scores <- pca1$x
rownames(scores) <- colnames(tabH)
## plot
n <- 15
pdf(file = paste0("PCA_PCs1-", n+1, "_smoothed_human.pdf"), width = 6, height = 6)
for (i in 1:n){
  ## colors for tissues
  ## symbols for species
  col.v <- pal[as.integer(info$Tissue)]
  pchs <- c(15,16,17)[as.integer(info$Species)]
  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=F,legend=F)
}
dev.off()
## tissues first
## then: unknown factors
## PC5,6,7: H3
## PC7,8: H2
## PC9: H4

## TO DO: do it on the whole set of CpG sites covered in human


## just include X chromosome: do we cluster male and female separately?
tabX <- tab[grepl("^chrX", row.names(tab), perl=T),]
dim(tabX) ## 251,551 sites
print(head(tabX))

pca1 <- prcomp(t(tabX), scale = FALSE) ## unscaled
print(summary(pca1)) ## 51% variance explained by 3 first components
loadings <- pca1$rotation
scores <- pca1$x
rownames(scores) <- colnames(tabX)
## plot
n <- 20
pdf(file = paste0("PCA_PCs1-", n+1, "_smoothed_chrX.pdf"), width = 6, height = 6)
for (i in 1:n){
  ## colors for tissues
  ## symbols for species
  col.v <- pal[as.integer(info$Tissue)]
  pchs <- c(15,16,17)[as.integer(info$Species)]
  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=F,legend=F)
}
dev.off()

## Figure S...
## 1 color for males and 1 for females
pdf(file = "PCA_PCs1-2_smoothed_chrX_symbols.pdf", width = 6, height = 6)
## symbols + transparency for species
## nothing for tissues
col.v <- rep(pal[10], times=48)
col.v[info$IndividualID == "R2" | info$IndividualID == "C3"] = pal[11]
col.v <- addalpha(col.v, c(0.65,1,0.4)[as.integer(info$Species)])
pchs <- c(15,16,17)[as.integer(info$Species)]
plot_scores(pca1, scores, 1, 2, cols=col.v, pch=pchs, points=T,legend=T)
legend("topright", c("Male","Female"), pch=16, cex=1.3, col=pal[10:11])
legend("bottomright", c("Human","Chimpanzee", "Macaque"), pch=c(16,15,17), cex=1.3, col=addalpha(rep("black", times=3), c(1,0.65,0.4)))
dev.off()


## perform PCA on raw methylation estimates ########################################
## tab.raw <- getMeth(allData.fit.subset, type = "raw")
## colnames(tab.raw) <- sampleNames(allData.fit.subset)
## row.names(tab.raw) <- paste(seqnames(granges(allData.fit.subset)), start(granges(allData.fit.subset)), sep=":")
## tab.raw <- tab.raw[rowSds(tab.raw) > 0.01, ] ## filter rows with low variance. 1,333,344 sites left
## tab.raw <- na.omit(tab.raw)
## dim(tab.raw)
## print(head(tab.raw))
## save(tab.raw, file="PCA_combinedRawCommonSites_table.RDa")

load("PCA_combinedRawCommonSites_table.RDa")
colnames(tab.raw)= sampleNames(allData.fit.subset)

## tab.raw <- tab.raw[!grepl("^chrX", row.names(tab.noX), perl=T),]
## only 3000 sites on chrX, let's keep the whole object
pca1 <- prcomp(t(tab.raw), scale = FALSE) 
print(summary(pca1)) ## 22.4% variance explained by 3 first components

## loadings: coordinates of PCs in original coordinate system
loadings <- pca1$rotation
## scores: coordinates of scaled data (observations) projected onto the PCs.
scores <- pca1$x
print(head(scores))
rownames(scores) <- colnames(tab.raw)

## plot
n <- 20
pdf(file = paste0("PCA_PCs1-", n+1, "_raw.pdf"), width = 6, height = 6)
for (i in 1:n){
  col.v <- pal[as.integer(info$Tissue)]
  pchs <- c(15,16,17)[as.integer(info$Species)]
  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=F,legend=F)
}
dev.off()
pdf(file = paste0("PCA_PCs1-", n+1, "_raw_symbols.pdf"), width = 6, height = 6)
for (i in 1:n){
  col.v <- pal[as.integer(info$Tissue)]
  pchs <- c(15,16,17)[as.integer(info$Species)]
  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=T,legend=T)
  legend("topright", c("Human","Chimpanzee", "Macaque"), pch=c(16,15,17), cex=1.3, col="black")
  legend("bottomright", c("Heart","Kidney", "Liver", "Lung"), pch=16, cex=1.3, col=pal[1:4])
}
dev.off()
## very nice separation too!
## PC1: species (rhesus)
## PC2: tissue (all, but esp. liver)
## PC3: tissue (all, but esp. kidney)
## PC4: 3 species
## PC5: lung
## PC6: H1K: something wrong? Lower read depth?
## PC10: rhesus individuals: R4, R2/R3 mixed?
## PC11: rhesus R1
## PC14: chimp C2
## PC16: chimp C1


## TO DO? subsample a smaller number of random shared sites to perform PCA
## TO DO? Do we introduce any bias toward human (since we map to human)? Do we make the species distinction (PC1) more clear because of this?
## TO DO: keep only sites with high coverage in several samples in the three species?


## correlation of PCs to technical factors ############################
Cov <- getCoverage(allData.fit.subset)
coverage <- colMeans(Cov)
apply(scores, 2, function(x){ cor(coverage, x, method="spearman") } )
apply(scores, 2, function(x){ cor.test(coverage, x, method="spearman", exact=F)$p.value } )
## No significant association with any PC: good!
## plot
coverage.std <- (coverage-min(coverage))/(max(coverage-min(coverage))) ## standardize coverage
## col.v <- rev(brewer.pal(11, "RdYlBu"))[round(coverage.std * 10) + 1]
col.v <- rev(colorRampPalette(brewer.pal(10,"RdBu"))(100))[round(coverage.std * 99) + 1]
n <- 20
## pdf(file = paste0("PCA_PCs1-", n+1, "_raw_coverage.pdf"), width = 6.5, height = 6)
pdf(file = paste0("PCA_PCs1-", n+1, "_smoothed_coverage.pdf"), width = 6.5, height = 6)
pchs <- c(15,16,17)[as.integer(info$Species)]
for (i in 1:n){
  print(i)
  layout(matrix(c(1,2), nrow=1), widths=c(5,1))
  par(mar=c(5, 4, 4, 0.5)) ##bottom, left, top, right
  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=F,legend=F)
  par(mar=c(5, 0.5, 4, 4), las=1) ##bottom, left, top, right
  image(t(as.matrix(1:100)), col=rev(colorRampPalette(brewer.pal(10,"RdBu"))(100)), axes=F)
  axis(4, at=seq(0, 1, length.out=10), labels=signif(seq(min(coverage), max(coverage), length.out=10), 2), col="white", col.ticks="black")
  mtext("Coverage", side=1, line=1, at=1)
}
dev.off()


## mean % methylation of samples
M <- getMeth(allData.fit.subset)
methylation <- colMeans(M)
apply(scores, 2, function(x){ cor(methylation, x, method="spearman") } )
apply(scores, 2, function(x){ cor.test(methylation, x, method="spearman", exact=F)$p.value } )
## Strongly correlated with PC1 (rho=-0.92)
## plot
methylation.std <- (methylation-min(methylation))/(max(methylation-min(methylation))) ## standardize methylation
## col.v <- rev(brewer.pal(11, "RdYlBu"))[round(methylation.std * 10) + 1]
col.v <- rev(colorRampPalette(brewer.pal(10,"RdBu"))(100))[round(methylation.std * 99) + 1]
n <- 20
pdf(file = paste0("PCA_PCs1-", n+1, "_smoothed_percent_methylation.pdf"), width = 6.5, height = 6)
pchs <- c(15,16,17)[as.integer(info$Species)]
for (i in 1:n){
  print(i)
  layout(matrix(c(1,2), nrow=1), widths=c(5,1))
  par(mar=c(5, 4, 4, 0.5)) ##bottom, left, top, right
  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=F,legend=F)
  par(mar=c(5, 0.5, 4, 4), las=1) ##bottom, left, top, right
  image(t(as.matrix(1:100)), col=rev(colorRampPalette(brewer.pal(10,"RdBu"))(100)), axes=F)
  axis(4, at=seq(0, 1, length.out=10), labels=signif(seq(min(methylation), max(methylation), length.out=10), 2), col="white", col.ticks="black")
  mtext("Methylation", side=1, line=1, at=1)
}
dev.off()
## Interesting: mean methylation estimates are lower in macaque than in chimp and human is the highest! This is explaining PC1 almost entirely. Is this only because of species differences? Not only: even inside each species there is a relation
pdf(file = paste0("PCA_PCs1_smoothed_percent_methylation.pdf"), width = 6, height = 6)
plot(methylation, scores[,1], pch=16, xlab="% methylation", ylab="PC1 scores", col=rep(pal[1:3], each=16))
dev.off()

## Prop of UMR, LMR, HMR (see analysis.R)
props <-  read.table("distributionMethylationClasses.txt", sep="\t", row.names=1, h=T)
apply(scores, 2, function(x){ cor(props$Low.methylation, x, method="spearman") } )
apply(scores, 2, function(x){ cor.test(props$Low.methylation, x, method="spearman", exact=F)$p.value } )
## UMR: PC2 (rho=0.59), PC5 (rho=0.51)
apply(scores, 2, function(x){ cor(props$Intermediate.methylation, x, method="spearman") } )
apply(scores, 2, function(x){ cor.test(props$Intermediate.methylation, x, method="spearman", exact=F)$p.value } )
##LMR: PC1 (rho=0.82)
apply(scores, 2, function(x){ cor(props$High.methylation, x, method="spearman") } )
apply(scores, 2, function(x){ cor.test(props$High.methylation, x, method="spearman", exact=F)$p.value } )
##HMR: PC1 (rho=-0.82)

## plot for LMR
methylation.std <- (props$Intermediate.methylation-min(props$Intermediate.methylation))/(max(props$Intermediate.methylation-min(props$Intermediate.methylation))) ## standardize props$Intermediate.methylation
## col.v <- rev(brewer.pal(11, "RdYlBu"))[round(methylation.std * 10) + 1]
col.v <- rev(colorRampPalette(brewer.pal(10,"RdBu"))(100))[round(methylation.std * 99) + 1]
n <- 20
pdf(file = paste0("PCA_PCs1-", n+1, "_smoothed_prop_LMR.pdf"), width = 6.5, height = 6)
pchs <- c(15,16,17)[as.integer(info$Species)]
for (i in 1:n){
  print(i)
  layout(matrix(c(1,2), nrow=1), widths=c(5,1))
  par(mar=c(5, 4, 4, 0.5)) ##bottom, left, top, right
  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=F,legend=F)
  par(mar=c(5, 0.5, 4, 4), las=1) ##bottom, left, top, right
  image(t(as.matrix(1:100)), col=rev(colorRampPalette(brewer.pal(10,"RdBu"))(100)), axes=F)
  axis(4, at=seq(0, 1, length.out=10), labels=signif(seq(min(props$Intermediate.methylation), max(props$Intermediate.methylation), length.out=10), 2), col="white", col.ticks="black")
  mtext("Prop. LMR", side=1, line=1, at=1)
}
dev.off()
## This is not very informative. The % LMR is tightly linked to the overall methylation levels, so it is higher in human vs. macaque, and so human has more LMRs...

## investigate if the %methylation difference between species is not artifactuals #########################
## It seems strange that there is overall higher methylation in human compared to chimp and macaque... This could be artifactual.

## let's filter out differentially methylated positions (in each tissue), and check if species is still a first PC
allDMRs <- vector()
allDMRs.filtered <- vector()
for (pair in list(c("Human", "Chimp"), c("Human", "Rhesus"), c("Chimp", "Rhesus"))){
  for (tissue in c("heart", "kidney", "liver", "lung")){
    cat("Species:", pair, " / Tissue:", tissue, "\n")
    load(paste0("DMRs/species/", pair[1], pair[2], "_", tissue, "_tstat.RDa"))

    ## find DMRs
    dmrs0 <- dmrFinder(allData.tstat, cutoff = c(-4.6, 4.6), maxGap = 300, stat = "tstat.corrected", verbose=T)
    dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1) ## filtered
    
    allDMRs <- sort(unique(c(allDMRs, paste(dmrs0$chr, dmrs0$start, sep=":"))))
    allDMRs.filtered <- sort(unique(c(allDMRs.filtered, paste(dmrs$chr, dmrs$start, sep=":"))))
  }
}

tab.noDMR <- tab[!row.names(tab) %in% allDMRs,]
dim(tab.noDMR) ## 4,696,414
tab.noDMR.filtered <- tab[!row.names(tab) %in% allDMRs.filtered,] ## less stringent filtering: useful?
dim(tab.noDMR.filtered) ## 4,778,157

## redo PCA
pca1 <- prcomp(t(tab.noDMR), scale = FALSE) 
print(summary(pca1)) ## 49.12% variance explained by 3 first components
loadings <- pca1$rotation
scores <- pca1$x
rownames(scores) <- colnames(tab.noDMR)

n <- 20
pdf(file = paste0("PCA_PCs1-", n+1, "_smoothed_noSpeciesDMR.pdf"), width = 6, height = 6)
for (i in 1:n){
  print(i)
  ## colors for tissues
  col.v <- pal[as.integer(info$Tissue)]
  ## symbols for species (not used)
  pchs <- c(15,16,17)[as.integer(info$Species)]
  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=F, legend=F)
}
dev.off()
## This is virtually identical to the PCA with differentially methylated positions included!
## The species segragation on PCA1 could be due to bias... or to many non-significant "real" differences between species
## TO DO? maybe signal is driven by sites next to DMRs: fiter these out!


## PCA on DMRs sites only: what is the species pattern?
tab.onlyDMR <- tab[row.names(tab) %in% allDMRs,]
dim(tab.onlyDMR) ## 182,942
tab.onlyDMR.filtered <- tab[!row.names(tab) %in% allDMRs.filtered,] ## less stringent filtering: useful?
dim(tab.onlyDMR.filtered) ## 4,778,157

## redo PCA
pca1 <- prcomp(t(tab.onlyDMR), scale = FALSE) 
print(summary(pca1)) ## 57.8% variance explained by 3 first components
loadings <- pca1$rotation
scores <- pca1$x
rownames(scores) <- colnames(tab.onlyDMR)

n <- 20
pdf(file = paste0("PCA_PCs1-", n+1, "_smoothed_onlySpeciesDMR.pdf"), width = 6, height = 6)
for (i in 1:n){
  print(i)
  ## colors for tissues
  col.v <- pal[as.integer(info$Tissue)]
  ## symbols for species (not used)
  pchs <- c(15,16,17)[as.integer(info$Species)]
  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=F, legend=F)
}
dev.off()
## interesting: the changes are not so profound compared to PCA_PCs1-21_smoothed.pdf
## PC1: rhesus vs. human+chimp
## PC2: tissues: esp. liver!
## PC3: species

## TO DO? PCA on same number of sites taken randomly as a control
## This is likely to be very similar to PCA_PCs1-21_smoothed.pdf (see below)

## Using only human-chimp DMRs:
allDMRs <- vector()
for (tissue in c("heart", "kidney", "liver", "lung")){
  cat("Tissue:", tissue, "\n")
  load(paste0("DMRss/pecies/HumanChimp_", tissue, "_tstat.RDa"))
  ## find DMRs
  dmrs0 <- dmrFinder(allData.tstat, cutoff = c(-4.6, 4.6), maxGap = 300, stat = "tstat.corrected", verbose=T)
  dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1) ## filtered
    
  allDMRs <- sort(unique(c(allDMRs, paste(dmrs0$chr, dmrs0$start, sep=":"))))
}

tab.onlyDMR <- tab[row.names(tab) %in% allDMRs,]
dim(tab.onlyDMR) ## 19920 only

## redo PCA
pca1 <- prcomp(t(tab.onlyDMR), scale = FALSE) 
print(summary(pca1)) ## 58% variance explained by 3 first components
loadings <- pca1$rotation
scores <- pca1$x
rownames(scores) <- colnames(tab.onlyDMR)

n <- 20
pdf(file = paste0("PCA_PCs1-", n+1, "_smoothed_onlyHumanChimpDMR.pdf"), width = 6, height = 6)
for (i in 1:n){
  print(i)
  ## colors for tissues
  col.v <- pal[as.integer(info$Tissue)]
  ## symbols for species (not used)
  pchs <- c(15,16,17)[as.integer(info$Species)]
  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=F, legend=F)
}
dev.off()
## PC1 and PC2: species, very strongly
## PC3: tissue


## Only Human-chimp DMRs in 1 tissue
allDMRs <- vector()
tissue <- "heart"
load(paste0("DMRs/species/HumanChimp_", tissue, "_tstat.RDa"))
## find DMRs
dmrs0 <- dmrFinder(allData.tstat, cutoff = c(-4.6, 4.6), maxGap = 300, stat = "tstat.corrected", verbose=T)
dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1) ## filtered
allDMRs <- sort(unique(c(allDMRs, paste(dmrs0$chr, dmrs0$start, sep=":"))))

tab.onlyDMR <- tab[row.names(tab) %in% allDMRs,]
dim(tab.onlyDMR) ## 7880 sites only

pca1 <- prcomp(t(tab.onlyDMR), scale = FALSE) 
print(summary(pca1)) ## 59.5% variance explained by 3 first components
loadings <- pca1$rotation
scores <- pca1$x
rownames(scores) <- colnames(tab.onlyDMR)

n <- 20
pdf(file = paste0("PCA_PCs1-", n+1, "_smoothed_onlyHumanChimpHeartDMR.pdf"), width = 6, height = 6)
for (i in 1:n){
  print(i)
  ## colors for tissues
  col.v <- pal[as.integer(info$Tissue)]
  ## symbols for species (not used)
  pchs <- c(15,16,17)[as.integer(info$Species)]
  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=F, legend=F)
}
dev.off()
## PC1 and PC2: species, very strongly. heart samples at the extreme between human and chimp
## PC3: tissue


## Species-specific DMRs
## Using only human-specific DMRs:
allDMRs <- vector()
for (tissue in c("heart", "kidney", "liver", "lung")){
  cat("Tissue:", tissue, "\n")
  load(paste0("DMRs/species/HumanSpecific_", tissue, "_tstat.RDa"))
  ## find DMRs
  dmrs0 <- dmrFinder(allData.tstat, cutoff = c(-4.6, 4.6), maxGap = 300, stat = "tstat.corrected", verbose=T)
  dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1) ## filtered
    
  allDMRs <- sort(unique(c(allDMRs, paste(dmrs0$chr, dmrs0$start, sep=":"))))
}

tab.onlyDMR <- tab[row.names(tab) %in% allDMRs,]
dim(tab.onlyDMR) ## 22862 sites

## redo PCA
pca1 <- prcomp(t(tab.onlyDMR), scale = FALSE) 
print(summary(pca1)) ## 65% variance explained by 3 first components
loadings <- pca1$rotation
scores <- pca1$x
rownames(scores) <- colnames(tab.onlyDMR)

n <- 20
pdf(file = paste0("PCA_PCs1-", n+1, "_smoothed_onlyHumanSpecificDMR.pdf"), width = 6, height = 6)
for (i in 1:n){
  print(i)
  ## colors for tissues
  col.v <- pal[as.integer(info$Tissue)]
  ## symbols for species (not used)
  pchs <- c(15,16,17)[as.integer(info$Species)]
  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=F, legend=F)
}
dev.off()
## very similar to Human-chimp DMRs only


## tissue DMRs (in human only)
allDMRs <- vector()
for (pair in list(c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "liver"))){
  cat("Testing tissue pair: ", pair, "\n")
  load(paste0("DMRs/tissues/Human_", pair[1], "_", pair[2], "_tstat.RDa"))

  ## find DMRs
  dmrs0 <- dmrFinder(allData.tstat, cutoff = c(-4.6, 4.6), maxGap = 300, stat = "tstat.corrected", verbose=T)
  allDMRs <- sort(unique(c(allDMRs, paste(dmrs0$chr, dmrs0$start, sep=":"))))
}
tab.onlyDMR <- tab[row.names(tab) %in% allDMRs,]
dim(tab.onlyDMR) ## 30899 sites

## redo PCA
pca1 <- prcomp(t(tab.onlyDMR), scale = FALSE) 
print(summary(pca1)) ## 54.7% variance explained by 3 first components
loadings <- pca1$rotation
scores <- pca1$x
rownames(scores) <- colnames(tab.onlyDMR)

n <- 20
pdf(file = paste0("PCA_PCs1-", n+1, "_smoothed_onlyHumanTissueDMR.pdf"), width = 6, height = 6)
for (i in 1:n){
  print(i)
  ## colors for tissues
  col.v <- pal[as.integer(info$Tissue)]
  ## symbols for species (not used)
  pchs <- c(15,16,17)[as.integer(info$Species)]
  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=F, legend=F)
}
dev.off()
## PC1 now has tissues!
## PC2: rhesus (+liver a bit)


## tissue-specific DMRs (any species): e.g., heart
allDMRs <- vector()
for (species in c("Human", "Chimp", "Rhesus")){
  cat("Testing species: ", species, "\n")
  load(paste0("DMRs/tissues/", species, "_heartSpecific_tstat.RDa"))

  ## find DMRs
  dmrs0 <- dmrFinder(allData.tstat, cutoff = c(-4.6, 4.6), maxGap = 300, stat = "tstat.corrected", verbose=T)
  allDMRs <- sort(unique(c(allDMRs, paste(dmrs0$chr, dmrs0$start, sep=":"))))
}
tab.onlyDMR <- tab[row.names(tab) %in% allDMRs,]
dim(tab.onlyDMR) ## 6996 sites

## redo PCA
pca1 <- prcomp(t(tab.onlyDMR), scale = FALSE) 
print(summary(pca1)) ## 58% variance explained by 3 first components
loadings <- pca1$rotation
scores <- pca1$x
rownames(scores) <- colnames(tab.onlyDMR)

n <- 20
pdf(file = paste0("PCA_PCs1-", n+1, "_smoothed_onlyHeartSpecificDMR.pdf"), width = 6, height = 6)
for (i in 1:n){
  print(i)
  ## colors for tissues
  col.v <- pal[as.integer(info$Tissue)]
  ## symbols for species (not used)
  pchs <- c(15,16,17)[as.integer(info$Species)]
  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=F, legend=F)
}
dev.off()
## PC1: tissues (and particularly heart!)
## PC2: rhesus

## TO DO? conserved tissue-specific DMRs

## TO DO: individual DMRs


################################################################
## let's filter out all positions were unmethylated reads are observed on both strands (not SNPs for sure), and check if species is still a first PC. The PCA will be biased because we tend to remove highly methylated positions, but we'll see if there is indeed a species component on first PC
## Matrix of evidence codes produced by script combine_all_evidence.gz

evidence <- read.table("combine_data_perl/evidence.gz", h=T, sep="\t")
row.names(evidence) <- paste(evidence[,1], evidence[,2], sep=":")
evidence <- evidence[,3:50]

## PCA on E sites
## see http://stackoverflow.com/questions/25667941/r-efficiently-grep-characters-in-rows-of-large-data-frame
rows <- which(
  rowSums(
    `dim<-`(grepl("E", as.matrix(evidence), fixed=TRUE), dim(evidence))
  ) > 0
)

subset_PCA(tab, row.names(tab) %in% row.names(evidence)[rows], 20, "PCA_PCs1-21_smoothed_erroneous_sites.pdf")
## 52826 sites tested
## 53.163% of variance explained by the 3 first components
## pattern simialr to exons etc
## Rhesus segregeating on PC1 and liver on PC2
## -> these sites do not bias the analysis

## filter NA
evidence.filtered <- na.omit(evidence) ## 1,354,302 sites left

## filter out E sites 
summary(apply(evidence.filtered, 1, function(x){ sum(grepl("E", as.matrix(x), perl=T)) > 0 })) ##10,507
evidence.filtered <- evidence.filtered[ apply(evidence.filtered, 1, function(x){ sum(grepl("E", as.matrix(x), perl=T)) == 0 }),]

## filter out any site with S evidence
summary(apply(evidence.filtered, 1, function(x){ sum(grepl("S", as.matrix(x), perl=T)) > 0 })) ##118
## not enough sites left!

## Actually we only need one U evidence per individual!
evidence.strand <- matrix(ncol=12, nrow=nrow(evidence.filtered))
for (i in 1:12){ ##12 individuals
  print(i)
  evidence.strand[,i] <- apply(evidence.filtered[,(i*4-3):(i*4)] == "U", 1, sum) > 0
}
## we want all individuals to have some unmethylated reads evicence on both strands
summary(apply(evidence.strand, 1, sum) == 12) 
##   Mode   FALSE    TRUE    NA's 
##logical 1119620  224175       0 
evidence.filtered <- evidence.filtered[apply(evidence.strand, 1, sum) == 12 ,]

## TO DO?: do S sites tend to be S on all samples of the same individuals: likely SNPs! 

## redo PCA using only filtered sites
## - on smoothed estimates
summary(row.names(tab) %in% row.names(evidence.filtered))
##   Mode   FALSE    TRUE    NA's 
##logical 4687774  191582       0 
## Only a small fraction of the sites are kept
tab.strand <- tab[row.names(tab) %in% row.names(evidence.filtered),]

## redo PCA
pca1 <- prcomp(t(tab.strand), scale = FALSE) 
print(summary(pca1)) ## 52.18% variance explained by 3 first components
loadings <- pca1$rotation
scores <- pca1$x
rownames(scores) <- colnames(tab.strand)

n <- 20
pdf(file = paste0("PCA_PCs1-", n+1, "_smoothed_both_strands_evidence.pdf"), width = 6, height = 6)
for (i in 1:n){
  print(i)
  ## colors for tissues
  col.v <- pal[as.integer(info$Tissue)]
  ## symbols for species (not used)
  pchs <- c(15,16,17)[as.integer(info$Species)]
  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=F, legend=F)
}
dev.off()
## This is very close to unfiltered tab. Rhesus still separates on PC1.
## This tends to show that methylation differences are not due to CpG -> TpG SNPs biasing the methylation estimates


## Redo PCA with same number of random sites to see how the smaller number of sites influences the results
pca1 <- prcomp(t(tab[sample(1:length(tab[,1]), length(tab.strand[,1])),]), scale = FALSE) 
print(summary(pca1)) 
loadings <- pca1$rotation
scores <- pca1$x
rownames(scores) <- colnames(tab.strand)
n <- 20
pdf(file = "temp.pdf", width = 6, height = 6)
for (i in 1:n){
  col.v <- pal[as.integer(info$Tissue)]
  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=F, legend=F)
}
dev.off()
## There is virtually no difference between using 5M sites and 191,582 sites
## The differences between PCA_PCs1-21_smoothed.pdf and PCA_PCs1-21_smoothed_both_strands_evidence.pdf are due to the choice of a biased set of sites


## Maybe smoothing is influenced by SNPs biased estimates?
## ->check with PCA on raw methylation estimates
tab.strand.raw <- tab.raw[row.names(tab.raw) %in% row.names(evidence.filtered),] ## 224175

## redo PCA
pca1 <- prcomp(t(tab.strand.raw), scale = FALSE) 
print(summary(pca1)) ## 24.4% variance explained by 3 first components
loadings <- pca1$rotation
scores <- pca1$x
rownames(scores) <- colnames(tab.strand.raw)

n <- 20
pdf(file = paste0("PCA_PCs1-", n+1, "_raw_both_strands_evidence.pdf"), width = 6, height = 6)
for (i in 1:n){
  print(i)
  ## colors for tissues
  col.v <- pal[as.integer(info$Tissue)]
  ## symbols for species (not used)
  pchs <- c(15,16,17)[as.integer(info$Species)]
  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=F, legend=F)
}
dev.off()
## PC1+2: rhesus and tissues(liver)
## Less importance of species (maybe?) when using this subset of positions... but it is still very present! This was also present when using teh smoothed data (some "liver vs. other tissues" signal on first PC)
## quite similar to PCA on all sites (only "rotation" on first PCs)



## Is there a correlation with age of the samples? ################################################
## This could be partly confounded with species
## See file Samples_used.xls
age <- c(rep(19.92, 4), rep(40.08, 4), rep(22.92, 4), rep(18.83, 4),  rep(70, 4),  rep(69, 4),  rep(39, 4),  rep(67, 4),  rep(7.25, 4),  rep(14, 4),  rep(6.42, 4),  rep(7.33, 4)) 
apply(scores, 2, function(x){ cor(age, x, method="spearman", use="pairwise.complete.obs") } )
apply(scores, 2, function(x){ cor.test(age, x, method="spearman", exact=F, use="pairwise.complete.obs")$p.value } )
## Strongly correlated with PC1 (rho=-0.81)
## plot
age.std <- (na.omit(age)-min(na.omit(age)))/(max(na.omit(age)-min(na.omit(age)))) ## standardize age
## col.v <- rev(brewer.pal(11, "RdYlBu"))[round(age.std * 10) + 1]
col.v <- rev(colorRampPalette(brewer.pal(10,"RdBu"))(100))[round(age.std * 99) + 1]
n <- 20
pdf(file = paste0("PCA_PCs1-", n+1, "_smoothed_age.pdf"), width = 6.5, height = 6)
pchs <- c(15,16,17)[as.integer(info$Species)]
for (i in 1:n){
  print(i)
  layout(matrix(c(1,2), nrow=1), widths=c(5,1))
  par(mar=c(5, 4, 4, 0.5)) ##bottom, left, top, right
  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=F,legend=F)
  par(mar=c(5, 0.5, 4, 4), las=1) ##bottom, left, top, right
  image(t(as.matrix(1:100)), col=rev(colorRampPalette(brewer.pal(10,"RdBu"))(100)), axes=F)
  axis(4, at=seq(0, 1, length.out=10), labels=signif(seq(min(na.omit(age)), max(na.omit(age)), length.out=10), 2), col="white", col.ticks="black")
  mtext("Age", side=1, line=1, at=1)
}
dev.off()
## Problem: it is not possible to compare age of human and age of other species, esp. macaque.
## See http://genomics.senescence.info/species/entry.php?species=Pan_troglodytes
##     http://genomics.senescence.info/species/entry.php?species=Macaca_mulatta
## plot versus PC1: 
pdf(file = paste0("PCA_PCs1_smoothed_age.pdf"), width = 6, height = 6)
plot(age, scores[,1], pch=16, xlab="Age of samples", ylab="PC1 scores", col=rep(pal[1:3], each=16))
dev.off()
## No correlation within species: this doesn't seem to be driving the correlation with % methylation (both between and within species)
pdf(file = "percent_methylation_age.pdf", width = 6, height = 6)
plot(age, methylation, pch=16, xlab="Age of samples", ylab="Mean % methylation", col=rep(pal[1:3], each=16))
dev.off()
## Older samples have higher methylation (/!\ this is probably confounded with species! It doesn't seem to be observed within species)

## We take the relative age compared to the max lifespan: 40y for macaque, 59.4y for chimp, 90y for human (see AnAge)
## "Compared to other species, of course, the maximum longevity of humans is based on a considerably larger sample. Therefore, it has been argued that, for comparative purposes, it is more adequate to use as human maximum longevity 90 or 100 years."
rel_age <- c(c(rep(19.92, 4), rep(40.08, 4), rep(22.92, 4), rep(18.83, 4))/59.4,  c(rep(70, 4),  rep(69, 4),  rep(39, 4),  rep(67, 4))/90,  c(rep(7.25, 4),  rep(14, 4),  rep(6.42, 4),  rep(7.33, 4))/40) 
apply(scores, 2, function(x){ cor(rel_age, x, method="spearman", use="pairwise.complete.obs") } )
apply(scores, 2, function(x){ cor.test(rel_age, x, method="spearman", exact=F, use="pairwise.complete.obs")$p.value } )
## Less strongly correlated with PC1 (rho=-0.74)
## plot
rel_age.std <- (na.omit(rel_age)-min(na.omit(rel_age)))/(max(na.omit(rel_age)-min(na.omit(rel_age)))) ## standardize rel_age
## col.v <- rev(brewer.pal(11, "RdYlBu"))[round(rel_age.std * 10) + 1]
col.v <- rev(colorRampPalette(brewer.pal(10,"RdBu"))(100))[round(rel_age.std * 99) + 1]
n <- 20
pdf(file = paste0("PCA_PCs1-", n+1, "_smoothed_relative_age.pdf"), width = 6.5, height = 6)
pchs <- c(15,16,17)[as.integer(info$Species)]
for (i in 1:n){
  print(i)
  layout(matrix(c(1,2), nrow=1), widths=c(5,1))
  par(mar=c(5, 4, 4, 0.5)) ##bottom, left, top, right
  plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=F,legend=F)
  par(mar=c(5, 0.5, 4, 4), las=1) ##bottom, left, top, right
  image(t(as.matrix(1:100)), col=rev(colorRampPalette(brewer.pal(10,"RdBu"))(100)), axes=F)
  axis(4, at=seq(0, 1, length.out=10), labels=signif(seq(min(na.omit(rel_age)), max(na.omit(rel_age)), length.out=10), 2), col="white", col.ticks="black")
  mtext("Age relative to max lifespan", side=1, line=1, at=1)
}
dev.off()
pdf(file = paste0("PCA_PCs1_smoothed_relative_age.pdf"), width = 6, height = 6)
plot(rel_age, scores[,1], pch=16, xlab="Age relative to max lifespan", ylab="PC1 scores", col=rep(pal[1:3], each=16))
dev.off()
pdf(file = "percent_methylation_relative_age.pdf", width = 6, height = 6)
plot(rel_age, methylation, pch=16, xlab="Age relative to max lifespan", ylab="Mean % methylation", col=rep(pal[1:3], each=16))
dev.off()
## No correlation within species: this doesn't seem to be driving the correlation with % methylation (both between and within species)

## TO DO? other technical factors: % duplicates / Mapping efficiency / conversion rate?
## They do not seem to be correlated, and they correspond to technical replicates: see below


## PCA of technical replicates #############################################
## Use raw methylation estimates, only at 5 million common sites
load("smooth_data/combined_samples/Methylation_technical.RDa")
## Problem: many NAs! only 5 rows without any NA :(

## Idea: remove samples with higher number of NAs
nas <- apply(Meth, 2, function(x){ sum(is.na(x)) })
sort(nas) ## from 489,023 to 4,991,069
summary(nas > length(Meth[,1])/2) ## 14 samples have more than half NAs, but removing only these leaves only 124 rows
Meth <- Meth[, nas < 1000000] ## 50 samples, 183383 rows
Meth <- Meth[, nas < 1500000] ## 75 samples, 18068 row
Meth <- Meth[, grepl("H\\dLi", colnames(Meth), perl=T)]## We want ot keep all human liver samples if there was some mix up of technical replicates. 1284318 left

Meth <- na.omit(Meth)
dim(Meth)
## TO DO? missing values imputation for ignored samples? pcaMethods package
## http://www.bioconductor.org/packages/release/bioc/vignettes/pcaMethods/inst/doc/missingValues.pdf

sds <- apply(Meth, 1, sd)
Meth <- Meth[sds > 0.01, ] ## filter rows with low variance. 
dim(Meth)
pca1 <- prcomp(t(Meth), scale = FALSE) ## turn off the argument "scale".
print(summary(pca1)) ## 50 samples: 20.4% variance explained by 3 first components
                     ## 75 samples: 17.2% variance explained by 3 first components
                     ## human livers only: 53% variance explained by 3 first components
loadings <- pca1$rotation
scores <- pca1$x
rownames(scores) <- colnames(Meth)
## plot
samples <- read.table("../raw_data/SamplesDirectories.txt", sep="\t", h=T)
samples <- unique(samples[samples$Type == "BS-seq", ])
row.names(samples) <- paste(samples$Flow.cell, samples$SampleID, sep="_")
samples <- samples[,c(4,5,6,7)]
samples <- samples[order(samples$SampleID),]
summary(names(nas) == row.names(samples))
## samples <- samples[nas < 1000000,]
## samples <- samples[nas < 1500000,]
samples <- samples[grepl("H\\dLi", samples$Condition, perl=T), ]

## n <- 20
## pdf(file = paste0("PCA_PCs1-", n+1, "_raw_50_technical_replicates.pdf"), width = 6, height = 6)
## pdf(file = paste0("PCA_PCs1-", n+1, "_raw_75_technical_replicates.pdf"), width = 6, height = 6)
for (n in 1:n){
  ## colors for tissues
  ## symbols for species
  col.v <- pal[as.integer(samples$Tissue)]
  pchs <- c(15,16,17)[as.integer(samples$Species)]
  plot_scores(pca1, scores, n, n+1, cols=col.v, pch=pchs, points=F,legend=F)
}
dev.off()

n <- 7
pdf(file = paste0("PCA_PCs1-", n+1, "_raw_human_liver_technical_replicates.pdf"), width = 6, height = 6)
for (n in 1:n){
  ## colors and symbols for different individuals
  col.v <- pal[as.integer(factor(samples$Condition))]
  pchs <- c(15,16,17,18)[as.integer(samples$Condition)]
  plot_scores(pca1, scores, n, n+1, cols=col.v, pch=pchs, points=F,legend=T)
  legend("bottomright", legend=unique(samples$Condition), pch=c(15,16,17,18), cex=1.3, col=pal[1:4])
}
dev.off()

## Expected large influence of coverage, but actually not!!! Very nice clustering of tehcnical replcicates together
## PC1: rhesus vs. human/chimp
## PC2: Liver
## PC3,4: tissues
## All technical replicates of samples for which there is some doubt (H1Li): group together

## TO DO? other technical factors: % duplicates / Mapping efficiency / conversion rate?
## They do not seem to be correlated, and have no obvious reason to be: not investigated for now.



################################################################
## PCA on some genomic features only:
## Using 15_features_sharedPositions.gz, we can produce PCAs for specific features (e.g., CGI shores)

## filter out chrX, which can introduce male/female separation
tabNoX <- tab[!grepl("^chrX", row.names(tab), perl=T),]
dim(tabNoX) ## 4,627,805 sites

## read the data 
features <- read.table("../annotation/15_features_sharedPositions.gz", h=T, sep="\t")
row.names(features) <- paste(features[,1], features[,2], sep=":")
features <- features[,3:17]
features <- features[row.names(features) %in% row.names(tab), ]
summary(row.names(tab) == row.names(features))
featuresNoX <- features[!grepl("^chrX", row.names(tab), perl=T), ]

## first perform some checks ###############
## we expect no overlap in these comparisons
summary(features$exon == 1 & features$intron == 1)
summary(features$exon == 1 & features$intergenic == 1)
summary(features$intron == 1 & features$intergenic == 1)
summary(features$intron == 1 & features$CDS == 1)
summary(features$"5_UTR" == 1 & features$CDS == 1) 
## Some overlap here! This is because a UTR from one transcript can be within a CDS from another transcript. Same thing for 3'UTR

summary(features$exon == 1 | features$intron == 1 | features$intergenic == 1)
##   Mode   FALSE    TRUE    NA's 
## logical  385213 4494143       0 
## This corresponds to exons and introns of non coding genes

summary(features$exon == 1 & features$conserved == 1)
##    Mode   FALSE    TRUE    NA's 
## logical 4493881  385475       0 
## Most conserved positions are in exons

summary(features$promoter == 1 & features$CpG_island == 1)
##    Mode   FALSE    TRUE    NA's 
## logical 4477377  401979       0 
summary(features$promoter_coding == 1 & features$CpG_island == 1)
##    Mode   FALSE    TRUE    NA's 
## logical 4513845  365511       0 
## 1/4 to 1/3 of promoters are in CpG islands


## PCA on promoters ######################################################
## 1378849 positions
subset_PCA(tabNoX, featuresNoX$promoter == 1, 20, "PCA_PCs1-21_smoothed_promoter.pdf")
## 53% variance explained by 3 first components
## Rhesus separation on PC1 seems stronger than on PCA with all sites
## Tissue separation cleaner on PC2
## Idem on PC3: human-chimp-rhesus separation is cleaner

## PCA on coding promoters 
subset_PCA(tabNoX, featuresNoX$promoter_coding == 1, 20, "PCA_PCs1-21_smoothed_promoter_coding.pdf")
## 53% of variance explained by the 3 first components
## identical to PCA on all promoters

## PCA on coding promoters 
subset_PCA(tabNoX, featuresNoX$promoter_proximal == 1, 20, "PCA_PCs1-21_smoothed_promoter_proximal.pdf")
## 55% of variance explained by the 3 first components
## identical to PCA on all promoters

## PCA on promoters (positions in conserved elements) 
subset_PCA(tabNoX, featuresNoX$promoter == 1 & featuresNoX$conserved == 1 , 20, "PCA_PCs1-21_smoothed_promoter_conserved.pdf")
## 53% of variance explained by the 3 first components
## Getting some tissue component on PC1

## PCA on promoters (positions within CpG island)
subset_PCA(tabNoX, featuresNoX$promoter == 1 & featuresNoX$CpG_island == 1 , 20, "PCA_PCs1-21_smoothed_promoter_CpG_island.pdf")
## 61% of variance explained by the 3 first components!
## PC1/2: species (first rhesus, then human-chimp)
## Very surpising pattern: C3 and R2/3 female individuals are separated on PC2!
## PC3 only" liver
## This is fun because it corresponds to almost all CpG island positions (401979/535981)
subset_PCA(tabNoX, featuresNoX$promoter_proximal == 1 & featuresNoX$CpG_island == 1 , 20, "PCA_PCs1-21_smoothed_promoter_proximal_CpG_island.pdf")
## 62% of variance explained by the 3 first components!

## PCA on promoters (positions NOT within CpG island)
subset_PCA(tabNoX, featuresNoX$promoter == 1 & featuresNoX$CpG_island == 0 , 20, "PCA_PCs1-21_smoothed_promoter_not_CpG_island.pdf")
## 51% of variance explained by the 3 first components
## very similar to all promoters (maybe even cleaner actually)

## PCA on promoters FAR from CpG island)
subset_PCA(tabNoX, featuresNoX$promoter == 1 & featuresNoX$CpG_island == 0 & featuresNoX$CpG_island_shore == 0 & featuresNoX$CpG_island_shelf== 0 , 20, "PCA_PCs1-21_smoothed_promoter_far_from_CpG_island.pdf")
## 51% of variance explained by the 3 first components
subset_PCA(tabNoX, featuresNoX$promoter_proximal == 1 & featuresNoX$CpG_island == 0 & featuresNoX$CpG_island_shore == 0 & featuresNoX$CpG_island_shelf== 0 , 20, "PCA_PCs1-21_smoothed_promoter_proximal_far_from_CpG_island.pdf")
## 54% of variance explained by the 3 first components

## PCA on promoters (positions in CpG island & conserved elements)
subset_PCA(tabNoX, featuresNoX$promoter == 1 & featuresNoX$conserved == 1 & featuresNoX$CpG_island == 1 , 20, "PCA_PCs1-21_smoothed_promoter_CpG_island_conserved.pdf")
## 58% of variance explained by the 3 first components
## very similar to pattern at all CpG islands
subset_PCA(tabNoX, featuresNoX$promoter_proximal == 1 & featuresNoX$conserved == 1 & featuresNoX$CpG_island == 1 , 20, "PCA_PCs1-21_smoothed_promoter_proximal_CpG_island_conserved.pdf")
## 57% of variance explained by the 3 first components

## PCA on promoters (positions in CpG island & not conserved elements)
subset_PCA(tabNoX, featuresNoX$promoter == 1 & featuresNoX$conserved == 0 & featuresNoX$CpG_island == 1 , 20, "PCA_PCs1-21_smoothed_promoter_CpG_island_not_conserved.pdf")
## 63% of variance explained by the 3 first components!!
## similar to promoter & all CpG islands
subset_PCA(tabNoX, featuresNoX$promoter_proximal == 1 & featuresNoX$conserved == 0 & featuresNoX$CpG_island == 1 , 20, "PCA_PCs1-21_smoothed_promoter_proximal_CpG_island_not_conserved.pdf")
## 64% of variance explained by the 3 first components!!

## See below: the partition intergenic/genic seems relevant

## PCA on exons ##########################################################
subset_PCA(tabNoX, featuresNoX$exon == 1, 20, "PCA_PCs1-21_smoothed_exon.pdf")
## 53% variance explained by 3 first components
## Quite similar to PCA on promoter
## Cleaner than PCA on all sites (but same pattern)

## PCA on first exons ######################################################
subset_PCA(tabNoX, featuresNoX$first_exon == 1, 20, "PCA_PCs1-21_smoothed_first_exon.pdf")
## 53% variance explained by 3 first components
## Quite similar to PCA on all exons:
## Maybe species separation a bit more strong and tissue separation a bit less string (but very small trend)

## PCA on repeats ######################################################
subset_PCA(tabNoX, featuresNoX$"repeat" == 1, 20, "PCA_PCs1-21_smoothed_repeats.pdf")
## 50% variance explained by 3 first components
## There is a clear tissue component to PC1 now (although rhesus still clearly separated)
## But actually not so many changes afterwards

## repeats & intergenic
subset_PCA(tabNoX, featuresNoX$"repeat" == 1 & featuresNoX$intergenic == 1, 20, "PCA_PCs1-21_smoothed_repeats_intergenic.pdf")
## 52% variance
## quite similar to full set of repeats (explians a bit more variance)

## repeats & genic (UTR+CDS)
subset_PCA(tabNoX, featuresNoX$"repeat" == 1 & featuresNoX$intergenic == 0, 20, "PCA_PCs1-21_smoothed_repeats_genic.pdf")
## 48% of variance explained by the 3 first components
## quite similar to PCA on all sites

## repeats & conserved
subset_PCA(tabNoX, featuresNoX$"repeat" == 1 & featuresNoX$conserved == 1, 20, "PCA_PCs1-21_smoothed_repeats_conserved.pdf")
## 49% of variance explained by the 3 first components
## female individuals C3 and R2 separated on PC3

## repeats & NOT conserved
subset_PCA(tabNoX, featuresNoX$"repeat" == 1 & featuresNoX$conserved == 0, 20, "PCA_PCs1-21_smoothed_repeats_not_conserved.pdf")
## 50% of variance explained by the 3 first components
## similar to whole set of repeats


## PCA on CpG islands ######################################################
subset_PCA(tabNoX, featuresNoX$CpG_island == 1, 20, "PCA_PCs1-21_smoothed_CpG_island.pdf")
## 61% variance explained by 3 first components !!
## A lot of variance explained: CpG islands are important ;)
## PC1 very clearly separates Rhesus
## PC2 separates tissues + also chimp, rhesus and human
## PC3 separates female individuals (PC7 on full PCA)
## PC4/5: separates kidney (PC3/4 in full PCA) 

## CpG island & promoter: see above. Surprising clustering of individuals on PC2
## CpG island & promoter & conserved: see above
## CpG island & promoter & not conserved: see above

## CpG island & conserved
subset_PCA(tabNoX, featuresNoX$CpG_island == 1 & featuresNoX$conserved == 1, 20, "PCA_PCs1-21_smoothed_CpG_island_conserved.pdf")
## 57% of variance explained by the 3 first components
## quite similar to whole set of CpG islands

## CpG island & NOT conserved
subset_PCA(tabNoX, featuresNoX$CpG_island == 1 & featuresNoX$conserved == 0, 20, "PCA_PCs1-21_smoothed_CpG_island_not_conserved.pdf")
## 62%% of variance explained by the 3 first components!
## individuals in PC2/3 (liver also)
## individual pattern less clear than on promoters+CpG_island(+not_conserved)

## CpG island & NOT promoter
subset_PCA(tabNoX, featuresNoX$CpG_island == 1 & featuresNoX$promoter == 0, 20, "PCA_PCs1-21_smoothed_CpG_island_not_promoter.pdf")
## 60% of variance explained by the 3 first components
## PC1: rhesus
## PC2: tissues
## PC3: 3 species
## PC4: tissues

## PCA on CDS ########################################################
subset_PCA(tabNoX, featuresNoX$CDS == 1, 20, "PCA_PCs1-21_smoothed_CDS.pdf")
## 53% variance explained by 3 first components
## very similar to PCA on exon

## PCA on intron ######################################################
subset_PCA(tabNoX, featuresNoX$intron == 1, 20, "PCA_PCs1-21_smoothed_intron.pdf")
## 50% variance explained by 3 first components
## better separating human and chimp on PC1 than PCA on exon
## better separation of heart on PC2

## PCA on intergenic regions ##########################################
subset_PCA(tabNoX, featuresNoX$intergenic == 1, 20, "PCA_PCs1-21_smoothed_intergenic.pdf")
## 52% of variance explained by 3 first components
## species + tissue (liver) on PC1
## PC2: liver (+ chimp)
## similar to repeats

## intergenic & conserved
subset_PCA(tabNoX, featuresNoX$intergenic == 1 & featuresNoX$conserved == 1, 20, "PCA_PCs1-21_smoothed_intergenic_conserved.pdf")
## 51% of variance explained by the 3 first components
## quite different than intergenic
## PC1: rhesus + chimp + human more separated

## intergenic & NOT conserved
subset_PCA(tabNoX, featuresNoX$intergenic == 1 & featuresNoX$conserved == 0, 20, "PCA_PCs1-21_smoothed_intergenic_not_conserved.pdf")
## 52% of variance explained by the 3 first components
## virtually identical to intergenic only (makes sense given the number of loci)

## intergenic and NOT promoter
subset_PCA(tabNoX, featuresNoX$intergenic == 1 & featuresNoX$promoter == 0, 20, "PCA_PCs1-21_smoothed_intergenic_not_promoter.pdf")
## 52% of variance explained by the 3 first components
## similar to intergenic, but even more amplified tissue trend on PC1

## intergenic and promoter
subset_PCA(tabNoX, featuresNoX$intergenic == 1 & featuresNoX$promoter == 1, 20, "PCA_PCs1-21_smoothed_intergenic_promoter.pdf")
## 55% of variance explained by the 3 first components
## Nice! PC1 and 2 are separting the 3 species
## slight tissue trend on PC2
## PC3: liver
## PC4" kidney

## PCA on conserved regions ###########################################
subset_PCA(tabNoX, featuresNoX$conserved == 1, 20, "PCA_PCs1-21_smoothed_conserved.pdf")
## % of variance explained by 3 first components
## seems quite similar to introns
## PC2: separating better heart from liver


######################
## heatmap plotting ##
######################
## filter out chrX, which can introduce male/female separation
tabNoX <- tab[!grepl("^chrX", row.names(tab), perl=T),]
dim(tabNoX) ## 4,627,805 sites

## Use smoothed methylation data in variable "tab" (see above PCA)
library(gplots)
## colors <- rev(colorRampPalette(brewer.pal(10,"RdBu"))(100)) 
colors <- colorRampPalette(c(brewer.pal(9, "Blues")[1],brewer.pal(9, "Blues")[9]))(100)

## Clustering using all sites: too long!
## subset_heatmap(tabNoX, NA, colors, info, "heatmap_smoothed_spearman.pdf")

## Use pearson correlation
cors <- cor(tabNoX, method="pearson", use="pairwise.complete.obs") ## long, but feasible. 
labels <- paste(info$Species, info$Tissue, sep=" ")
pdf(file = "heatmap_smoothed_pearson.pdf", width = 12, height = 8)
heatmap.2( cors, scale="none", col = colors, margins = c(12, 12), trace='none', denscol="white", labCol=labels, ColSideColors=pal[as.integer(as.factor(info$Species))], RowSideColors=pal[10:13][as.integer(as.factor(info$Tissue))], cexCol = 1.5)
dev.off()
## However pearson correlation is not very adapted to these data (strongly bimodal)

## randomly subset 10000 sites
subset_heatmap_random_sites(tabNoX, !is.na(row.names(tabNoX)), colors, info, "heatmap_smoothed_10000_random.pdf", 10000)
## - Good and robust!
## - Only strange pattern with H1Li, clustering with Human heart and lung? On PCA it is also distant from the others. But no mixing of technical replicates, the tissue quality must be lower.
## - rhesus samples cluster out from (human and chimp)
## - within human and chimp, it is more complex: liver human+chimp group out, but heart, lung and kidney group by species
reuse_subset_heatmap_random_sites(tabNoX, "heatmap_smoothed_10000_random.txt", colors, info, "heatmap_smoothed_10000_random2.pdf")
## Figure 2B


## randomly subset 100000 sites
subset_heatmap_random_sites(tabNoX, !is.na(row.names(tabNoX)), colors, info, "heatmap_smoothed_100000_random.pdf", 100000)
## very similar to 10,000 sites

## top 1% of variable sites
sds <- apply(tabNoX, 1, sd) ## long!
subset <- sds > quantile(sds, probs=0.99) ## 46279
subset_heatmap(tabNoX, subset, colors, info, "heatmap_smoothed_top_sd.pdf")
## - very clean pattern. Here human and chimp kidney do not cluster together, but with other tissues in each species 
## Interesting: there are negative correlation values + very small variation in each condition: we probably enrich for DMRs

## more highly covered sites
Cov <- getCoverage(allData.fit.subset)
row.names(Cov) <- paste(seqnames(granges(allData.fit.subset)), start(granges(allData.fit.subset)), sep=":")
colnames(Cov) <- row.names(info)
Cov <- Cov[row.names(tab), ]
CovNoX <- Cov[row.names(tabNoX), ]

## coverage >= 4X in 4 samples in each species
## keepLoci <- which(rowSums(Cov[, grepl("^C", colnames(Cov), perl=T)] >= 4) >= 4 & rowSums(Cov[, grepl("^R", colnames(Cov), perl=T)] >= 4) >= 4 & rowSums(Cov[, grepl("^H", colnames(Cov), perl=T)] >= 4) >= 4)
## Almost all loci fulfill this condition (3.8M)

## coverage > 6X in 24 samples
keepLoci <- which(rowSums(CovNoX > 6) >= 24)
length(keepLoci) ## 75614
subset_heatmap(tabNoX, keepLoci, colors, info, "heatmap_smoothed_high_coverage.pdf")

## select top coverage sites
covs <- apply(CovNoX, 1, mean) ## long!
subset <- covs > quantile(covs, probs=0.99) ## 45935
subset_heatmap(tabNoX, subset, colors, info, "heatmap_smoothed_top_coverage.pdf")

## select only sites in particular genomic features ####################
features <- read.table("../annotation/15_features_sharedPositions.gz", h=T, sep="\t", check.names=F)
row.names(features) <- paste(features[,1], features[,2], sep=":")
features <- features[,3:17]
features <- features[row.names(features) %in% row.names(tab), ]
summary(row.names(tab) == row.names(features))
featuresNoX <- features[!grepl("^chrX", row.names(tab), perl=T), ]

## promoter
sum(featuresNoX$promoter == 1) ## 1312821
subset_heatmap_random_sites(tabNoX, featuresNoX$promoter == 1, colors, info, "heatmap_smoothed_promoter_10000_random.pdf", 10000)
## - kidney clustering in
## - no more problem with H1Li

## promoter proximal
sum(featuresNoX$promoter_proximal == 1) ## 400316
subset_heatmap_random_sites(tabNoX, featuresNoX$promoter_proximal == 1, colors, info, "heatmap_smoothed_promoter_proximal_10000_random.pdf", 10000)
## similar to promoter
reuse_subset_heatmap_random_sites(tabNoX, "heatmap_smoothed_promoter_proximal_10000_random.txt", colors, info, "heatmap_smoothed_promoter_proximal_10000_random.pdf")
## Used as Figure S...


## CGI
sum(featuresNoX$CpG_island == 1) ## 508307
subset_heatmap_random_sites(tabNoX, featuresNoX$CpG_island == 1, colors, info, "heatmap_smoothed_CpG_island_10000_random.pdf", 10000)
## see also "heatmap_smoothed_CpG_island_10000_random2.pdf"
## clustering by species first, and then by tissue (a bit noisy)
## livers show good correlations, even across species (even if not always clustering together)

## CGI shore
sum(featuresNoX$CpG_island_shore == 1) ## 598024
subset_heatmap_random_sites(tabNoX, featuresNoX$CpG_island_shore == 1, colors, info, "heatmap_smoothed_CpG_island_shore_10000_random.pdf", 10000)
## - kidney clustering in
## - no more problem with H1Li

## CGI shelf
sum(featuresNoX$CpG_island_shelf == 1) ## 233936
subset_heatmap_random_sites(tabNoX, featuresNoX$CpG_island_shelf == 1, colors, info, "heatmap_smoothed_CpG_island_shelf_10000_random.pdf", 10000)
## not so stable:
## - same as global pattern
## - or clustering by tissue for human and chimp samples

## repeat
sum(featuresNoX$"repeat" == 1 & featuresNoX$intergenic == 1) ## 662762
subset_heatmap_random_sites(tabNoX, featuresNoX$"repeat" == 1 & featuresNoX$intergenic == 1, colors, info, "heatmap_smoothed_repeat_10000_random.pdf", 10000)
## same as global pattern

## CDS
sum(featuresNoX$CDS == 1) ## 471812
subset_heatmap_random_sites(tabNoX, featuresNoX$CDS == 1, colors, info, "heatmap_smoothed_CDS_10000_random.pdf", 10000)
## human and chimp tissue cluster together

## CGI & promoter
sum(featuresNoX$CpG_island == 1 & featuresNoX$promoter == 1) ## 378407
subset_heatmap_random_sites(tabNoX, featuresNoX$CpG_island == 1 & featuresNoX$promoter == 1, colors, info, "heatmap_smoothed_promoter_CpG_island_10000_random.pdf", 10000)
sum(featuresNoX$CpG_island == 1 & featuresNoX$promoter_proximal == 1) ## 186247
subset_heatmap_random_sites(tabNoX, featuresNoX$CpG_island == 1 & featuresNoX$promoter_proximal == 1, colors, info, "heatmap_smoothed_promoter_proximal_CpG_island_10000_random.pdf", 10000)
## clustering by species, messy within each species
reuse_subset_heatmap_random_sites(tabNoX, "heatmap_smoothed_promoter_proximal_CpG_island_10000_random.txt", colors, info, "heatmap_smoothed_promoter_proximal_CpG_island_10000_random.pdf")
## Figure S...

## CGI shore & promoter
sum(featuresNoX$CpG_island_shore == 1 & featuresNoX$promoter == 1) ## 407421
subset_heatmap_random_sites(tabNoX, featuresNoX$CpG_island_shore == 1 & featuresNoX$promoter == 1, colors, info, "heatmap_smoothed_promoter_CpG_island_shore_10000_random.pdf", 10000)
sum(featuresNoX$CpG_island_shore == 1 & featuresNoX$promoter_proximal == 1) ## 95259
subset_heatmap_random_sites(tabNoX, featuresNoX$CpG_island_shore == 1 & featuresNoX$promoter_proximal == 1, colors, info, "heatmap_smoothed_promoter_proximal_CpG_island_shore_10000_random.pdf", 10000)
## rhesus out, then liver, then chimp and human
## -> ~ similar to promoter

## promoter & conserved
sum(featuresNoX$conserved == 1 & featuresNoX$promoter == 1) ## 269947
subset_heatmap_random_sites(tabNoX, featuresNoX$conserved == 1 & featuresNoX$promoter == 1, colors, info, "heatmap_smoothed_promoter_conserved_10000_random.pdf", 10000)
sum(featuresNoX$conserved == 1 & featuresNoX$promoter_proximal == 1) ## 113035
subset_heatmap_random_sites(tabNoX, featuresNoX$conserved == 1 & featuresNoX$promoter_proximal == 1, colors, info, "heatmap_smoothed_promoter_proximal_conserved_10000_random.pdf", 10000)
## ~similar to promoters

## promoter & CGI & conserved
sum(featuresNoX$conserved == 1 & featuresNoX$promoter == 1 & featuresNoX$CpG_island == 1 ) ## 113182
subset_heatmap_random_sites(tabNoX, featuresNoX$conserved == 1 & featuresNoX$promoter == 1 & featuresNoX$CpG_island == 1, colors, info, "heatmap_smoothed_promoter_CpG_island_conserved_10000_random.pdf", 10000)
sum(featuresNoX$conserved == 1 & featuresNoX$promoter_proximal == 1 & featuresNoX$CpG_island == 1) ## 54314
subset_heatmap_random_sites(tabNoX, featuresNoX$conserved == 1 & featuresNoX$promoter_proximal == 1 & featuresNoX$CpG_island == 1, colors, info, "heatmap_smoothed_promoter_proximal_CpG_island_conserved_10000_random.pdf", 10000)
## ~ similar to promoters

## promoter & !CGI & !CGI shore & !CGI shelf 
sum(featuresNoX$CpG_island == 0 & featuresNoX$CpG_island_shore == 0 & featuresNoX$CpG_island_shore == 0 & featuresNoX$promoter == 1) ## 583303
subset_heatmap_random_sites(tabNoX, featuresNoX$CpG_island == 0 & featuresNoX$CpG_island_shore == 0 & featuresNoX$CpG_island_shore == 0 & featuresNoX$promoter == 1, colors, info, "heatmap_smoothed_promoter_far_from_CpG_island_10000_random.pdf", 10000)
sum(featuresNoX$CpG_island == 0 & featuresNoX$CpG_island_shore == 0 & featuresNoX$CpG_island_shore == 0 & featuresNoX$promoter_proximal == 1) ## 140224
subset_heatmap_random_sites(tabNoX, featuresNoX$CpG_island == 0 & featuresNoX$CpG_island_shore == 0 & featuresNoX$CpG_island_shore == 0 & featuresNoX$promoter_proximal == 1, colors, info, "heatmap_smoothed_promoter_proximal_far_from_CpG_island_10000_random.pdf", 10000)
reuse_subset_heatmap_random_sites(tabNoX, "heatmap_smoothed_promoter_proximal_far_from_CpG_island_10000_random.txt", colors, info, "heatmap_smoothed_promoter_proximal_far_from_CpG_island_10000_random.pdf")
## Figure S...

## ...

