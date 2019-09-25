## Jul 18, 2014
## There is an error when calculating the Chimp heart-liver contrast. Try to fix this manually in this script

library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))

## source("http://www.bioconductor.org/biocLite.R")
## biocLite()
## biocLite("BiocUpgrade")
## biocLite("bsseq")
## biocLite("bsseq", lib="~/R/x86_64-redhat-linux-gnu-library/2.15/") ## on the cluster
library(bsseq) # version 0.10
sessionInfo()

species <- "Chimp"
pair <- c("heart", "liver")

load("smooth_data/combined_samples/pData.RDa")
keep <- pData$Species == species
pData <- pData[keep,]
  
load(paste0("DMRs/tissues/", species, "_pairwiseTissues.RDa"))
  
keep <- pData$Tissue == pair[1] | pData$Tissue == pair[2] 
  
## Only keep sites with at least 2 samples with coverage of 2 in each tissue
keepLoci <- which(rowSums(getCoverage(allData.fit[, pData$Tissue == pair[1]]) >= 2) >= 2 & rowSums(getCoverage(allData.fit[, pData$Tissue == pair[2]]) >= 2) >= 2)

## remove loci with NA smoothed values
keepLoci <- keepLoci[!is.na(apply(getMeth(allData.fit[keepLoci, keep], type = "smooth"), 1, sum))]

## remove loci with mean coverage > 10X (likely to be mapped to repeated regions)
keepLoci <- keepLoci[rowMeans(getCoverage(allData.fit[keepLoci, keep])) <= 10]
cat("Tested sites: ", format(length(keepLoci), big.mark=","), "\n")

allData.fit.subset <- allData.fit[keepLoci, keep]

summary(getMeth(allData.fit.subset, type = "smooth"))

## compute t-statistics
## allData.tstat <- BSmooth.tstat(allData.fit.subset,
##                                group1 = rownames(pData[pData$Tissue == pair[1],]),
##                                group2 = rownames(pData[pData$Tissue  == pair[2],]),
##                                estimate.var = "same",
##                                local.correct = TRUE, ## large-scale (low-frequency) mean correction. This is especially important when large-scale methylation differences between 2 conditions (e.g., cancer and normals).
##                                verbose = TRUE)

## This bugging!
## Taking BSmooth.tstat code and going step by step to catch the problem
group1 <- rownames(pData[pData$Tissue == pair[1],])
group2 <- rownames(pData[pData$Tissue  == pair[2],])
estimate.var <- "same"
local.correct <- TRUE
qSd <- 0.75
k <- 101
mc.cores <- 1
compute.correction <- function(idx, qSd = 0.75) {
  xx <- start(BSseq)[idx]
  yy <- tstat[idx]
  suppressWarnings({
    drange <- diff(range(xx, na.rm = TRUE))
  })
  if (drange <= 25000)
    return(yy)
  tstat.function <- approxfun(xx, yy)
  xx.reg <- seq(from = min(xx), to = max(xx), by = 2000)
  yy.reg <- tstat.function(xx.reg)
  fit <- locfit(yy.reg ~ lp(xx.reg, h = 25000, deg = 2,
                            nn = 0), family = "huber", maxk = 50000)
  correction <- predict(fit, newdata = data.frame(xx.reg = xx))
  yy - correction
}
smoothSd <- function(Sds, k) {
  k0 <- floor(k/2)
  if (all(is.na(Sds)))
    return(Sds)
  thresSD <- pmax(Sds, quantile(Sds, qSd, na.rm = TRUE),
                  na.rm = TRUE)
  addSD <- rep(median(Sds, na.rm = TRUE), k0)
  sSds <- as.vector(runmean(Rle(c(addSD, thresSD, addSD)),
                            k = k))
  sSds
}

maxGap <- allData.fit.subset@parameters$maxGap
clusterIdx <- bsseq:::makeClusters(allData.fit.subset, maxGap = maxGap)
## computing stats within groups
allPs <- getMeth(allData.fit.subset, type = "smooth", what = "perBase", confint = FALSE)
group1.means <- rowMeans(allPs[, group1, drop = FALSE], na.rm = TRUE)
group2.means <- rowMeans(allPs[, group2, drop = FALSE], na.rm = TRUE)
## computing stats across groups
rawSds <- sqrt(((length(group1) - 1) * rowVars(allPs[, group1, drop = FALSE]) + (length(group2) - 1) * rowVars(allPs[, group2, drop = FALSE]))/(length(group1) + length(group2) - 2))
smoothSds <- do.call(c, mclapply(clusterIdx, function(idx) {
  smoothSd(rawSds[idx], k = k)
}, mc.cores = mc.cores))
scale <- sqrt(1/length(group1) + 1/length(group2))
tstat.sd <- smoothSds * scale
tstat <- (group1.means - group2.means)/tstat.sd
is.na(tstat)[tstat.sd == 0] <- TRUE
## local correction
tstat.corrected <- do.call(c, mclapply(clusterIdx, compute.correction, qSd = qSd, mc.cores = mc.cores))
## Pb here! Take the function compute correction step by step (with cluster 1)
xx <- start(allData.fit.subset)[clusterIdx[[1]]]
yy <- tstat[clusterIdx[[1]]]
drange <- diff(range(xx, na.rm = TRUE))
if (drange <= 25000)
  return(yy)
tstat.function <- approxfun(xx, yy)
xx.reg <- seq(from = min(xx), to = max(xx), by = 2000)
yy.reg <- tstat.function(xx.reg)
library(locfit) ## added: doesn't work without this
fit <- locfit(yy.reg ~ lp(xx.reg, h = 25000, deg = 2, nn = 0), family = "huber", maxk = 50000)
correction <- predict(fit, newdata = data.frame(xx.reg = xx))
yy - correction

compute.correction <- function(idx, BSseq, qSd = 0.75) { ## added BSseq object in arguments
  library(locfit) ## added
  xx <- start(BSseq)[idx]
  yy <- tstat[idx]
  suppressWarnings({
    drange <- diff(range(xx, na.rm = TRUE))
  })
  if (drange <= 25000)
    return(yy)
  tstat.function <- approxfun(xx, yy)
  xx.reg <- seq(from = min(xx), to = max(xx), by = 2000)
  yy.reg <- tstat.function(xx.reg)
  fit <- locfit(yy.reg ~ lp(xx.reg, h = 25000, deg = 2, nn = 0), family = "huber", maxk = 50000)
  correction <- predict(fit, newdata = data.frame(xx.reg = xx))
  yy - correction
}
tstat.corrected <- do.call(c, lapply(clusterIdx, compute.correction, BSseq = allData.fit.subset, qSd = qSd)) ## we also replaced mclapply by lapply here (but this was not the reason for bug)
## -> now we recreated the bug!
## Error in preplot.locfit.raw(object, newdata, where, what, band) : 
##   NA/NaN/Inf in foreign function call (arg 2)
## In addition: Warning message:
## In lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth,  :
##   max_nr reduction problem

## try to identify which cluster is bugging
for (i in 1:length(clusterIdx)){
  print(i)
  compute.correction(clusterIdx[[i]], BSseq = allData.fit.subset, qSd = qSd)
}
## cluster 8!
summary(clusterIdx[[8]])
compute.correction(clusterIdx[[8]], BSseq = allData.fit.subset, qSd = qSd)
for (i in 9:length(clusterIdx)){
  print(i)
  compute.correction(clusterIdx[[i]], BSseq = allData.fit.subset, qSd = qSd)
}
## other clusters are OK

## investigate deeper
xx <- start(allData.fit.subset)[clusterIdx[[8]]] ## full object has 640336 positions
xx <- start(allData.fit.subset)[clusterIdx[[8]][1:640263]]
yy <- tstat[clusterIdx[[8]][1:640263]]
drange <- diff(range(xx, na.rm = TRUE))
tstat.function <- approxfun(xx, yy)
xx.reg <- seq(from = min(xx), to = max(xx), by = 2000)
yy.reg <- tstat.function(xx.reg)
fit <- locfit(yy.reg ~ lp(xx.reg, h = 25000, deg = 2, nn = 0), family = "huber", maxk = 50000)
correction <- predict(na.omit(fit), newdata = data.frame(xx.reg = xx))
tstat.corrected.8 <- yy - correction
## problem with position 640264. But this works with cluster 8 truncated to 640263 first positions!

## clusterIdx are positions of CpG sites:
sum(unlist(lapply( clusterIdx, length)))
length(allData.fit.subset)

clusterIdx[[8]][640264:length(clusterIdx[[8]])]
allData.fit.subset2 <- allData.fit.subset[c(1:clusterIdx[[8]][640263], (clusterIdx[[8]][length(clusterIdx[[8]])]+1):length(allData.fit.subset)), ]

clusterIdx2 <- clusterIdx
clusterIdx2[[8]] <- clusterIdx2[[8]][1:640263]
tstat.corrected <- do.call(c, lapply(clusterIdx2, compute.correction, BSseq = allData.fit.subset2, qSd = qSd))
## Doesn't work... why???
## next step if it worked:
## stats <- cbind(rawSds, tstat.sd, group2.means, group1.means, tstat, tstat.corrected)
## colnames(stats) <- c("rawSds", "tstat.sd", "group2.means", "group1.means", "tstat", "tstat.corrected")
## parameters <- c(allData.fit.subset@parameters, list(tstatText = sprintf("BSmooth.tstat (local.correct = %s, maxGap = %d)", local.correct, maxGap), group1 = group1, group2 = group2, k = k, qSd = qSd, local.correct = local.correct, maxGap = maxGap))
## allData.tstat <- BSseqTstat(gr = granges(allData.fit.subset), stats = stats, parameters = parameters)
    
## try from beginning: should work now
allData.tstat <- BSmooth.tstat(allData.fit.subset2,
                               group1 = rownames(pData[pData$Tissue == pair[1],]),
                               group2 = rownames(pData[pData$Tissue  == pair[2],]),
                               estimate.var = "same",
                               local.correct = TRUE, ## large-scale (low-frequency) mean correction. This is especially important when large-scale methylation differences between 2 conditions (e.g., cancer and normals).
                               verbose = TRUE)

save(allData.tstat, file = paste0("DMRs/tissues/", species, "_", pair[1], "_", pair[2],"_tstat.RDa"))
pdf(file = paste0("DMRs/tissues/", species, "_", pair[1], "_", pair[2],"_histTstats.pdf"), width = 6, height = 6)
plot(allData.tstat)
dev.off()

## TO DO: look at object at problematic positions
allData.fit.subset[clusterIdx[[8]][640264:length(clusterIdx[[8]])], ]


## find DMRs
## with local correction: small-scale DMRs
dmrs0 <- dmrFinder(allData.tstat, cutoff = c(-4.6, 4.6), maxGap = 300, stat = "tstat.corrected", verbose=T)

## without local correction: large-scale DMRs (such as in cancer).
## /!\ the statistic may not be symmetrical.
## For cutoff, see Hansen et al. 2014
dmrs0.uncorrected <- dmrFinder(allData.tstat, cutoff = c(-2, 2), maxGap = 300, stat = "tstat", verbose=T)

## TO DO: adjust cutoff? comparable between species?
dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
cat("DMRs found: ", nrow(dmrs), "\n")
print(head(dmrs, n = 5))
write.table(dmrs, file= paste0("DMRs/tissues/", species, "_", pair[1], "_", pair[2],"_DMRs.txt"), sep="\t", quote=F, row.names=F, col.names=T)

dmrs.uncorrected <- subset(dmrs0.uncorrected, width >= 5000)
cat("DMRs (uncorrected t-stat) found: ", nrow(dmrs.uncorrected), "\n")
print(head(dmrs.uncorrected, n = 5))
write.table(dmrs.uncorrected, file= paste0("DMRs/tissues/", species, "_", pair[1], "_", pair[2],"_DMRs.uncorrected.txt"), sep="\t", quote=F, row.names=F, col.names=T)

## Plotting: 2 tissues of interest included only
n <- 100
if (nrow(dmrs) < n){
  n <- nrow(dmrs)
}
pData$col <- rep(c(pal[1], pal[2], pal[3], pal[4]), times = 4)
pData(allData.fit) <- pData
keepLoci <- which(rowSums(getCoverage(allData.fit[, keep])) >= 1)
pdf(file = paste0("DMRs/tissues/", species, "_", pair[1], "_", pair[2], "_top", n, "DMRs.pdf"), width = 10, height = 5)
plotManyRegions(allData.fit[keepLoci, keep], dmrs[1:n,], extend = 5000, addRegions = dmrs)
dev.off()

n <- 100
if (nrow(dmrs.uncorrected) < n){
  n <- nrow(dmrs.uncorrected)
}
if (n > 0){
  pdf(file = paste0("DMRs/tissues/", species, "_", pair[1], "_", pair[2], "_top", n, "DMRs.uncorrected.pdf"), width = 10, height = 5)
  plotManyRegions(allData.fit[keepLoci, keep], dmrs.uncorrected[1:n,], extend = 5000, addRegions = dmrs.uncorrected)
  dev.off()
}

print(warnings())
