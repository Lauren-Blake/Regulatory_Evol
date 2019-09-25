library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))

counts <- read.table("SamplesDirectories+rescued+counts.txt", h=T, sep="\t",comment.char="")
counts[,10] <- counts[,10]/4 ## when using wc there are 4 lines per read in the fastq files
names(counts)[10] <- "Read.count"

#############
## RNA-seq ##
#############
## plot the extent of # reads (all samples in all lanes)

sum_all_samples <- counts[counts$Type=="RNA-seq", ]
sum_all_samples <- sum_all_samples[order(sum_all_samples$Read.count), ]

par(mar=c(8, 4, 4, 2) + 0.1) # increase margins
mp <- barplot(log10(sum_all_samples$Read.count), ylab="log10(number of sequenced reads)", xlab="", xaxt="n")
# text(mp, -0.2, srt = 45, adj = 1, labels = paste(sum_all_samples$SampleId, " (Lane ", sum_all_samples$Lane, ", Flow-cell ", sum_all_samples$Flow.cell,")", sep=""), xpd = TRUE, cex=0.7)
mp <- barplot(log10(sum_all_samples$Read.count), ylab="log10(number of sequenced reads)", col=ifelse(grepl("Rescued", sum_all_samples$Directory), pal[3], pal[4]), xlab="", xaxt="n", border=NA)

## check that we have 48 different samples
aggregate(sum_all_samples$Condition, list(sum_all_samples$Condition), length)

## gather by sample
sum_samples <- aggregate(sum_all_samples$Read.count, list(sum_all_samples$Condition), sum) ## 48 samples
num_samples <- aggregate(sum_all_samples$Read.count, list(sum_all_samples$Condition), length)
num_samples <- num_samples[order(sum_samples[,2]), ]
sum_samples <- sum_samples[order(sum_samples[,2]), ]

sum_samples_no_rescue <- aggregate(sum_all_samples$Read.count[!grepl("Rescued", sum_all_samples$Directory)], list(sum_all_samples$Condition[!grepl("Rescued", sum_all_samples$Directory)]), sum) ## 48 samples
sum_samples_no_rescue <- sum_samples_no_rescue[order(match(sum_samples_no_rescue[,1],sum_samples[,1])),]
summary(sum_samples_no_rescue[,1] == sum_samples[,1])

par(mar=c(8, 4, 4, 2) + 0.1) # increase margins
## mp <- barplot(log10(sum_samples[,2]), ylab="Number of sequenced reads", xlab="", xaxt="n")
mp <- barplot(sum_samples[,2], ylab="Number of sequenced reads", xlab="", xaxt="n", col="white")
mp <- barplot(sum_samples_no_rescue[,2], add=T, col="darkgrey", axes=F)
text(mp, -2000000, srt = 45, adj = 1, labels = sum_samples[,1], xpd = TRUE, cex=0.8)
## adds the number of technical replicates
## text(mp, 3000000, srt = 0, labels = num_samples[,2], xpd = TRUE, cex=0.7)
legend("topleft", c("Demultiplexed reads", "Rescued reads"), fill = c("darkgrey", "white"), bty="n")
## What's special with C1Li? It was the sample that worked great in on the master mixes in teh flow cells that didn't work 
## For flow cell 1, in some lanes, almost all adapters have a N in them: the reads were almost all rescued
dev2bitmap("barplot_RNA-seq_rescue_rates.pdf", type="pdfwrite", method="pdf", w=10, h=6)

## arbitrary threshold of 25M reads
abline(h=2.5e7, lty=2, col="blue")
summary(sum_samples[,2] < 2.5e7)


## What is the ratio of rescued over total for all samples?
temp <- sum_all_samples[order(sum_all_samples$Flow.cell, sum_all_samples$Lane, sum_all_samples$Condition, sum_all_samples$Directory ),]
summary(temp[!grepl("Rescued", temp$Directory),5] == temp[grepl("Rescued", temp$Directory),5])
summary(temp[!grepl("Rescued", temp$Directory),4] == temp[grepl("Rescued", temp$Directory),4])
hist(temp[grepl("Rescued", temp$Directory),]$Read.count / (temp[!grepl("Rescued", temp$Directory),]$Read.count + temp[grepl("Rescued", temp$Directory),]$Read.count), breaks = 50, main="", xlab="Percentage of rescued reads")
## Some technical replicates have almost all their reads rescued!

summary(sum_samples[,2])
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 36050000 48160000 55510000 55360000 60640000 83860000

## gather by lane
sum_lanes <- aggregate(sum_all_samples$Read.count, list(sum_all_samples$Flow.cell, sum_all_samples$Lane), sum) ## 18 lanes

sum_lanes_no_rescue <- aggregate(sum_all_samples$Read.count[!grepl("Rescued", sum_all_samples$Directory)], list(sum_all_samples$Flow.cell[!grepl("Rescued", sum_all_samples$Directory)], sum_all_samples$Lane[!grepl("Rescued", sum_all_samples$Directory)]), sum)
sum_lanes_no_rescue <- sum_lanes_no_rescue[order(sum_lanes[,3]),]
sum_lanes <- sum_lanes[order(sum_lanes[,3]), ]
summary(sum_lanes_no_rescue[,1:2] == sum_lanes[,1:2])

par(mar=c(8, 4, 4, 2) + 0.1) # increase margins
mp <- barplot(sum_lanes[,3], ylab="Number of sequenced reads", xlab="", xaxt="n", col="white")
mp <- barplot(sum_lanes_no_rescue[,3], add=T, col=ifelse(sum_lanes[,2] %in% c(1,3,5,7), pal[1], pal[2]))
text(mp, -1000000, srt = 45, adj = 1, labels = paste("Lane ", sum_lanes[,2], ", Flow-cell ", sum_lanes[,1], sep=""), xpd = TRUE, cex=0.8)
legend("topleft", c("Master mix 1", "Master mix 2"), fill = c(pal[1], pal[2]), bty="n")
(sum_lanes[,3]-sum_lanes_no_rescue[,3])/sum_lanes[,3]
## saved from 7% to 66% of reads

## ## gather by species (3 species)
## sum_species <- aggregate(sum_all_samples$Read.count, list(sum_all_samples$Species), sum) 
## sum_species <- sum_species[order(sum_species[,2]), ]
## par(mar=c(8, 4, 4, 2) + 0.1) # increase margins
## mp <- barplot(sum_species[,2], ylab="Number of sequenced reads", xlab="", xaxt="n", col="darkgrey")
## text(mp, -10000000, srt = 45, adj = 1, labels = sum_species[,1], xpd = TRUE, cex=0.7)

## ## gather by tissue
## sum_tissue <- aggregate(sum_all_samples$Read.count, list(sum_all_samples$Tissue), sum) 
## sum_tissue <- sum_tissue[order(sum_tissue[,2]), ]
## par(mar=c(8, 4, 4, 2) + 0.1) # increase margins
## mp <- barplot(sum_tissue[,2], ylab="Number of sequenced reads", xlab="", xaxt="n", col="darkgrey")
## text(mp, -10000000, srt = 45, adj = 1, labels = sum_tissue[,1], xpd = TRUE, cex=0.7)

## ## gather by index
## mean_index <- aggregate(sum_all_samples$Read.count, list(sum_all_samples$Index), mean) ## 24 indexes
## par(mar=c(8, 4, 4, 2) + 0.1) # increase margins
## mp <- barplot(sort(mean_index[,2]), ylab="Number of reads", xlab="", xaxt="n")
## text(mp, -300000, srt = 45, adj = 1, labels = mean_index[order(mean_index[,2]),1], xpd = TRUE, cex=1)
## ## adds the number of technical replicates for each index
## num_index <- aggregate(sum_all_samples$Read.count, list(sum_all_samples$Index), length)
## ## text(mp, 100000, srt = 0, labels = num_index[order(mean_index[,2]),2], xpd = TRUE, cex=1) ## 10 for all

#####################
## RNA-seq mapping ##
#####################
tab <- read.table("../tophat/reports.txt", h=T, sep="\t", row.names=1) 

## samples ordered by number of reads mapped
tab <- tab[order(tab$Total.number.of.reads.mapped),]
par(mar=c(4, 4, 4, 2) + 0.1) # increase margins
mp <- barplot(tab$Total.number.of.reads.processed, ylab="Number of reads", xlab="", xaxt="n", col="white")
mp <- barplot(tab$Total.number.of.reads.mapped, add=T, col="darkgrey")
text(mp, -2000000, srt = 45, adj = 1, labels = row.names(tab), xpd = TRUE, cex=0.8)
legend("topleft", c("Unmapped", "Mapped"), fill = c("white", "darkgrey"), bty="n")
dev2bitmap("barplot_RNA-seq_mapping_rates.pdf", type="pdfwrite", method="pdf", w=10, h=6)

## same thing but ordered by total number of reads sequenced
tab <- tab[order(tab$Total.number.of.reads.processed),]
par(mar=c(4, 4, 4, 2) + 0.1) # increase margins
mp <- barplot(tab$Total.number.of.reads.processed, ylab="Number of reads", xlab="", xaxt="n", col="white")
mp <- barplot(tab$Total.number.of.reads.mapped, add=T, col="darkgrey")
text(mp, -2000000, srt = 45, adj = 1, labels = row.names(tab), xpd = TRUE, cex=0.8)
legend("topleft", c("Unmapped", "Mapped"), fill = c("white", "darkgrey"), bty="n")

## Percentage of reads mapped
tab <- tab[order(tab$Percentage.of.reads.mapped),]
par(mar=c(4, 4, 4, 2) + 0.1) # increase margins
mp <- barplot(tab$Percentage.of.reads.mapped, ylab="Percentage of mapped reads", xlab="", xaxt="n", col="darkgrey", ylim=c(0,100))
text(mp, -5, srt = 45, adj = 1, labels = row.names(tab), xpd = TRUE, cex=0.8)
## lower mappability on heart samples: why? see below

summary(tab$Total.number.of.reads.mapped)
##    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 28540000 40430000 47050000 46510000 51090000 72810000

summary(tab$Percentage.of.reads.mapped)
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  71.04   79.59   85.91   84.83   88.73   94.06

## Number of junctions
tab <- tab[order(tab$Number.of.junctions),]
par(mar=c(4, 4, 4, 2) + 0.1) # increase margins
mp <- barplot(tab$Number.of.junctions, ylab="Number of junctions discovered", xlab="", xaxt="n", col="darkgrey")
text(mp, -5000, srt = 45, adj = 1, labels = row.names(tab), xpd = TRUE, cex=0.8)

## ratio number junctions / number reads mapped
plot(tab$Number.of.junctions, tab$Total.number.of.reads.mapped, pch=16)

########################
## percentage mapping ##
########################
## More advanced stats on mapability: we noticed large differences in mapping % between samples: why is it the case?
stats <- read.table("~/work/methylation/stats_samples_RNA-seq.txt", h=T, sep="\t",comment.char="")
## stats_samples_RNA-seq.txt is exported from the supTable.xls excel sheet (first copy values only to a blank sheet to get rid of formatting)
## this table is quite similar to tab (above), but here we added more external infos (RIN, ...)
summary(tab[order(row.names(tab)),]$Number.of.junctions == stats$Number.of.junctions)
row.names(stats) <- stats[,1]
stats <- stats[stats[,1]!="H1H",] ## H1H is probably a liver sample

## differences between species
boxplot(stats$Percentage.of.reads.mapped ~ stats$Species, col=pal[1:3], main="", ylab="Mapping efficiency (%)", notch=F, pch=16, boxwex=0.7)

## differences between tissues
boxplot(stats$Percentage.of.reads.mapped ~ stats$Tissue, col=pal[1:4], main="", ylab="Mapping efficiency (%)", notch=F, pch=16, boxwex=0.7)
dev2bitmap("boxplot_RNA-seq_mapping_efficiency.pdf", type="pdfwrite", method="pdf", w=6, h=6)

summary(lm(stats$Percentage.of.reads.mapped ~ stats$Tissue)) ## p-value: 7.645e-11
t.test(stats$Percentage.of.reads.mapped[stats$Tissue == "heart"], stats$Percentage.of.reads.mapped[stats$Tissue == "kidney"]) ##  p-value = 4.345e-05
t.test(stats$Percentage.of.reads.mapped[stats$Tissue == "liver"], stats$Percentage.of.reads.mapped[stats$Tissue == "kidney"]) ##  p-value = 0.1039
t.test(stats$Percentage.of.reads.mapped[stats$Tissue == "liver"], stats$Percentage.of.reads.mapped[stats$Tissue == "lung"]) ## p-value = 0.05721
t.test(stats$Percentage.of.reads.mapped[stats$Tissue == "kidney"], stats$Percentage.of.reads.mapped[stats$Tissue == "lung"]) ## p-value = 0.001054

## split by species
boxplot(stats$Percentage.of.reads.mapped ~ interaction(stats$Tissue, stats$Species), col=pal[1:4], main="", ylab="Mapping efficiency (%)", notch=F, pch=16, boxwex=0.7, names=rep(c("heart", "kidney", "liver", "lung"), 3), at =c(1,2,3,4, 6,7,8,9, 11,12,13,14))
mtext(unique(stats$Species)[1], side = 1, line = 3, at=2.5)
mtext(unique(stats$Species)[2], side = 1, line = 3, at=7.5)
mtext(unique(stats$Species)[3], side = 1, line = 3, at=12.5)

## influence of RIN score?
plot(stats$Percentage.of.reads.mapped, stats$RIN, pch=16, xlab="% mapped reads", ylab="RIN score")
cor.test(stats$Percentage.of.reads.mapped, stats$RIN, method="spearman") ## rho =0.068, p=0.66
boxplot(stats$RIN ~ stats$Tissue, col=pal[1:4], main="", ylab="RIN score", notch=F, pch=16, boxwex=0.7)
## probably not: no big influence of tissue

## differences of RIN scores between species?
boxplot(stats[-17,]$RIN ~ stats[-17,]$Species, pch=16, xlab="", ylab="RIN score", col=pal[10:12], notch=T, boxwex=0.7)
dev2bitmap("boxplot_RNA-seq_RIN_species.pdf", type="pdfwrite", method="pdf", w=6, h=6)
## reflects lower quality of tissues samples in human, leading to less DE genes in human


## influence of concentration of RNA after extraction
plot(stats$Percentage.of.reads.mapped, stats$RNA.extraction.concentration..ng.uL., pch=16, xlab="% mapped reads", ylab="RNA extration concentration", log="y")
cor.test(stats$Percentage.of.reads.mapped, stats$RNA.extraction.concentration..ng.uL., method="spearman") ## rho=0.21, p=0.16
boxplot(stats$RNA.extraction.concentration..ng.uL. ~ stats$Tissue, col=pal[1:4], main="", ylab="RNA concentration", notch=F, pch=16, boxwex=0.7, log="y")
## lung is as low as heart!

## influence of concentration of RNA-seq library
plot(stats$Percentage.of.reads.mapped, stats$Concentration..ng.uL., pch=16, xlab="% mapped reads", ylab="RNA-seq library concentration")
cor.test(stats$Percentage.of.reads.mapped, stats$Concentration..ng.uL., method="spearman") ## rho=-0.33, p=0.025
boxplot(stats$Concentration..ng.uL. ~ stats$Tissue, col=pal[1:4], main="", ylab="library concentration", notch=F, pch=16, boxwex=0.7, log="y")
summary(lm(stats$Concentration..ng.uL. ~ stats$Tissue)) ## NS: p-value: 0.28
## Not same ranking for kidney and liver

## library complexity: because of PCR over-amplification, or because of differences in complexity of expressed transcriptome
## Preseq: gives number of unique reads for 28M subsampled reads
complexity <- read.table("~/clusterhome/Methylation/preseq/complexity_all.txt", h=F, sep="\t",comment.char="")
row.names(complexity) <- gsub(".txt:28000000", "", unlist(lapply(strsplit(as.character(complexity[,1]), "_"), function(x){ return(x[[2]])})))
complexity <- complexity[row.names(complexity)!="H1H",] ## H1H is probably a liver sample
summary(row.names(complexity)==stats$Sample.ID)
plot(stats$Percentage.of.reads.mapped, complexity[,2])
cor.test(stats$Percentage.of.reads.mapped, complexity[,2], method="spearman") ## rho=0.28, p=0.055
boxplot(complexity[,2] ~ stats$Tissue, col=pal[1:4], main="", ylab="complexity (preseq)", notch=F, pch=16, boxwex=0.7, log="y")
## much lower for liver: correlates well with cumulative distribution curves (see ../limma/analysis.R)

## other way of estimating transcriptome complexity:
## TMM normalization factors from edgeR (based on exons read counts)
## low value means that much of the sequencing effort is going to very highly expressed genes, leading to  under-sampling of the majority of other genes: low complexity
## tmm <- read.table("~/clusterhome/Methylation/limma/allSamples_libSize_normFactors_exons.txt", h=T, sep="\t",comment.char="", row.names=1)
## summary(row.names(tmm)==stats$sample.ID)
## cor.test(stats$Percentage.of.reads.mapped, tmm$norm.factors, method="spearman") ## rho=0.16, p=0.27
## boxplot(tmm$norm.factors ~ stats$Tissue, col=pal[1:4], main="", ylab="TMM", notch=F, pch=16, boxwex=0.7, log="y")
## This measure is very similar to preseq estimates
## plot(tmm$norm.factors, complexity[,2])
## cor.test(tmm$norm.factors, complexity[,2], method="spearman") ## rho=0.86, p< 2.2e-16

## other way of estimating transcriptome complexity:
## number of exons expressed?
## tab <- read.table("../orthoExon/counts.txt", sep="\t", h=T)
## rpkm <- read.table("../orthoExon/RPKM_normalized.txt", sep="\t", h=T)
## ## Compare histograms of log(rpkm): different for different tissues?
## par(mfrow=c(6,8))
## sapply(colnames(rpkm[,-c(1,2)]), function(x){
##   par(mar=c(2,1,1,1))
##   mp <- hist(log10(rpkm[,x]), breaks=100, main=x, xlab="", ylab="", xlim=c(-2,4))
##   abline(v=0, lty=2, col="red")
## })
## ## difference between tissues are not obvious. /!\ we don't see rpkm=0 (~one fourth of the points) in those plots
## lowlyExpressed <- sapply(colnames(rpkm[,-c(1,2)]), function(x){sum(rpkm[,x] < 1)})
## plot(stats$Percentage.of.reads.mapped, lowlyExpressed, type="n")
## text(stats$Percentage.of.reads.mapped, lowlyExpressed, names(lowlyExpressed))
## cor.test(stats$Percentage.of.reads.mapped, lowlyExpressed, method="spearman") ## rho=-0.32, p=0.030
## ## samples with many lowly expressed genes have a lower % mapping
## ## works with different cutoffs (0, 1, 10). RPKM = 1 seems to be a good cutoff to differentiate expressed/not expressed genes (see literature)
## boxplot(lowlyExpressed ~ stats$Tissue, col=pal[1:4], main="", ylab="Number of genes not expressed", notch=F, pch=16, boxwex=0.7, log="y")
## ## much lower for liver: correlates well with cumulative distribution curves (see ../limma/analysis.R)
## ## this measure is highly correlated to TMM and preseq measures
## plot(tmm$norm.factors, lowlyExpressed)
## cor.test(tmm$norm.factors, lowlyExpressed, method="spearman") ## rho=-0.82, p<2.2e-16
## cor.test(complexity[,2], lowlyExpressed, method="spearman") ## rho=-0.77, p<2.2e-16

## ## Lower complexity: only a few very highly expressed genes. Could these generate a lot of new isoforms (difficult to map)?
## ## we use as a proxy, the number of junctions, divided by the effective library size
## plot(stats$Percentage.of.reads.mapped, stats$Number.of.junctions/(tmm[,2]*tmm[,3]))
## cor.test(stats$Percentage.of.reads.mapped, stats$Number.of.junctions/(tmm[,2]*tmm[,3]), method="spearman") ## -0.21, p=0.14
## ## Libraries with more junctions discovered (normalized by sequencing effort) tend to have lower mappability
## plot(stats$Number.of.junctions/(tmm[,2]*tmm[,3]), lowlyExpressed, type="n")
## text(stats$Number.of.junctions/(tmm[,2]*tmm[,3]), lowlyExpressed, names(lowlyExpressed))
## cor.test(lowlyExpressed, stats$Number.of.junctions/(tmm[,2]*tmm[,3]), method="spearman") ## rho=0.62, p=3.6e-06
## ## Problem: this comparison is biased by the TMM correction factor (correlated with lowlyExpressed)?
## ## When we don't use thisfactor, the relation is opposite !!!
## plot(stats$Number.of.junctions/stats$Total.number.of.reads.mapped, lowlyExpressed, type="n")
## text(stats$Number.of.junctions/stats$Total.number.of.reads.mapped, lowlyExpressed, names(lowlyExpressed))
## cor.test(lowlyExpressed, stats$Number.of.junctions/stats$Total.number.of.reads.mapped, method="spearman") ## rho=-0.49, p=0.0005043
## boxplot(stats$Number.of.junctions/stats$Total.number.of.reads.mapped ~ stats$Tissue, col=pal[1:4], main="", ylab="Number of junctions", notch=F, pch=16, boxwex=0.7, log="y")

## ## TO DO: how to get a better proxy of complexity of splicing: see soumillon et al.
## Idea: Use % of reads overlapping a junction?
## plot(stats$Percentage.of.reads.mapped, stats$Percentage.of.mapped.reads.overlapping.a.junction)
## boxplot(stats$Percentage.of.mapped.reads.overlapping.a.junction ~ stats$Tissue, col=pal[1:4], main="", ylab="Percentage of reads overlapping a junction", notch=F, pch=16, boxwex=0.7)
## Liver higher than the others: artifact?

## Test of all potential remaining variables:
## Number of reads sequenced
boxplot(stats$Total.number.of.reads.sequenced ~ stats$Tissue, col=pal[1:4], main="", ylab="", notch=F, pch=16, boxwex=0.7) ## no significant difference

## Trimming
boxplot(stats$Percentage.of.bps.trimmed..adapters. ~ stats$Tissue, col=pal[1:4], main="", ylab="", notch=F, pch=16, boxwex=0.7, log="y")
cor.test(stats$Percentage.of.bps.trimmed..adapters., stats$Percentage.of.reads.mapped, method="spearman") ## rho=0.47, p=0.001 ## more trimming results in better mapping (probably just spurious correlation)
boxplot(stats$Number.of.reads.shorter.than.20bp.removed ~ stats$Tissue, col=pal[1:4], main="", ylab="", notch=F, pch=16, boxwex=0.7, log="y")

## Number of reads mapped
boxplot(stats$Total.number.of.reads.mapped ~ stats$Tissue, col=pal[1:4], main="", ylab="", notch=F, pch=16, boxwex=0.7) ## same ranking of tissues (logic!), but differences not signficant
cor.test(stats$Total.number.of.reads.mapped, stats$Percentage.of.reads.mapped, method="spearman") ## rho=0.29, p=0.05

## fragment size
boxplot(stats$Fragments.size..bp. ~ stats$Tissue, col=pal[1:4], main="", ylab="", notch=F, pch=16, boxwex=0.7)
## differences not signficant


## Mean quality scores of unmapped reads: lower for heart samples?
quality <- read.table("../FASTQC/quality_RNA-seq_unmapped.txt", h=F, sep="\t")
plot(c(min(quality[,2]):max(quality[,2])), c(min(quality[,2]):max(quality[,2])), ylim=c(min(quality[,3]), max(quality[,3])), type="n", ylab="Number of reads", xlab="Mean quality score")
for (sample in unique(quality[,1])){
  cat(sample, "\t", "\n")
  color <- "black"
  if(grepl("H$", sample, perl=T)){ color <- pal[1] }
  points(quality[quality[,1] == sample, 2], quality[quality[,1] == sample, 3], type="l", col=color)
}
## normalize by number of reads
plot(c(min(quality[,2]):max(quality[,2])), c(min(quality[,2]):max(quality[,2])), ylim=c(0, 0.5), type="n", ylab="Number of reads", xlab="Mean quality score")
for (sample in unique(quality[,1])){
  cat(sample, "\t", (quality[quality[,1] == sample, 3]/sum(quality[quality[,1] == sample, 3]))[quality[quality[,1] == sample, 2]==38], "\n")
  color <- "black"
  line <- 1
  if(grepl("H$", sample, perl=T)){ color <- pal[1] }
  if(sample %in% c("C1Li", "C4Lu", "R4Li", "C3H", "R3K")){ line <- 2 }
  points(quality[quality[,1] == sample, 2], quality[quality[,1] == sample, 3]/sum(quality[quality[,1] == sample, 3]), type="l", col=color, lty=line)
}
## 5 samples are peculiar: C1Li, C4Lu, R4Li, C3H and R3K. Only 1 heart
stats[c("C1Li", "C4Lu", "R4Li", "C3H", "R3K"), ]
## These are all the ones that were sequenced only in the lab, not the core facility!
## Do we observe the same pattern for mapped reads?
quality <- read.table("../FASTQC/quality_RNA-seq_mapped.txt", h=F, sep="\t")
plot(c(min(quality[,2]):max(quality[,2])), c(min(quality[,2]):max(quality[,2])), ylim=c(0, 0.7), type="n", ylab="Number of reads", xlab="Mean quality score")
for (sample in unique(quality[,1])){
  cat(sample, "\t", (quality[quality[,1] == sample, 3]/sum(quality[quality[,1] == sample, 3]))[quality[quality[,1] == sample, 2]==38], "\n")
  color <- "black"
  line <- 1
  if(grepl("H$", sample, perl=T)){ color <- pal[1] }
  if(sample %in% c("C1Li", "C4Lu", "R4Li", "C3H", "R3K")){ line <- 2 }
  points(quality[quality[,1] == sample, 2], quality[quality[,1] == sample, 3]/sum(quality[quality[,1] == sample, 3]), type="l", col=color, lty=line)
}
## The same 5 samples have peculiar profiles (less dramatic)


## Are the lower mappings due to the expression of more paralogous genes that are not mapped because we forced unique mapping in tophat?
## Remappng with multiple hits allowed
tab2 <- read.table("../tophat/tophat_multiple_hits/reports.txt", h=T, sep="\t", row.names=1) 
tab2 <- tab2[row.names(tab2)!="H1H",]
tab2 <- tab2[row.names(stats),]
summary(row.names(tab2) == row.names(stats))
stats[,25] <- tab2$Percentage.of.reads.mapped
names(stats)[25] <- "Percentage.of.reads.mapped.multiple"

## differences between species
boxplot(stats$Percentage.of.reads.mapped.multiple ~ stats$Species, col=pal[1:3], main="", ylab="Mapping efficiency (%)", notch=F, pch=16, boxwex=0.7)
## differences between tissues
boxplot(stats$Percentage.of.reads.mapped.multiple ~ stats$Tissue, col=pal[1:4], main="", ylab="Mapping efficiency (%)", notch=F, pch=16, boxwex=0.7) ## similar!
dev2bitmap("boxplot_RNA-seq_mapping_efficiency_multiple.pdf", type="pdfwrite", method="pdf", w=6, h=6)
summary(lm(stats$Percentage.of.reads.mapped.multiple ~ stats$Tissue)) ## p-value: 1.556e-10
## split by species
boxplot(stats$Percentage.of.reads.mapped.multiple ~ interaction(stats$Tissue, stats$Species), col=pal[1:4], main="", ylab="Mapping efficiency (%)", notch=F, pch=16, boxwex=0.7, names=rep(c("heart", "kidney", "liver", "lung"), 3), at =c(1,2,3,4, 6,7,8,9, 11,12,13,14))
mtext(unique(stats$Species)[1], side = 1, line = 3, at=2.5)
mtext(unique(stats$Species)[2], side = 1, line = 3, at=7.5)
mtext(unique(stats$Species)[3], side = 1, line = 3, at=12.5)

## correlation multiple / uniquely mappings
cor.test(stats$Percentage.of.reads.mapped.multiple, stats$Percentage.of.reads.mapped, method="spearman")
##  rho=0.9986124 p-value < 2.2e-16
plot(stats$Percentage.of.reads.mapped.multiple, stats$Percentage.of.reads.mapped, pch=16, col=ifelse(grepl("H$", stats[,1], perl=T), pal[1], "black"))


## Mapping on Ensembl genomes (more complete, include MT and Y chromosomes, and patches)
tab2 <- read.table("../tophat/tophat_Ensembl_genomes/reports.txt", h=T, sep="\t", row.names=1) 
tab2 <- tab2[row.names(tab2)!="H1H",]
tab2 <- tab2[row.names(stats),]
summary(row.names(tab2) == row.names(stats))
stats[,26] <- tab2$Percentage.of.reads.mapped
names(stats)[26] <- "Percentage.of.reads.mapped.Ensembl"

## correlation Ensembl / UCSC
cor.test(stats$Percentage.of.reads.mapped.Ensembl, stats$Percentage.of.reads.mapped, method="spearman")
##  rho=0.9986124 p-value < 2.2e-16
plot(stats$Percentage.of.reads.mapped.Ensembl, stats$Percentage.of.reads.mapped, pch=16, col=ifelse(grepl("H$", stats[,1], perl=T), pal[1], "black"), xlim=c(70,100), ylim=c(70,100))
abline(a=0, b=1)

## Percentage of reads mapped
tab2 <- tab2[order(tab2$Percentage.of.reads.mapped),]
par(mar=c(4, 4, 4, 2) + 0.1) # increase margins
mp <- barplot(tab2$Percentage.of.reads.mapped, ylab="Percentage of mapped reads", xlab="", xaxt="n", col="darkgrey", ylim=c(0,100))
text(mp, -5, srt = 45, adj = 1, labels = row.names(tab2), xpd = TRUE, cex=0.8)
summary(tab2$Percentage.of.reads.mapped)
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  88.33   92.53   93.26   93.43   94.19   96.75 
## very high for all samples!

## differences between species
boxplot(stats$Percentage.of.reads.mapped.Ensembl ~ stats$Species, col=pal[1:3], main="", ylab="Mapping efficiency (%)", notch=F, pch=16, boxwex=0.7)
## Differences are still here:
## - quality of genome assemblies
## - polymorphism
## - ...

## differences between tissues
boxplot(stats$Percentage.of.reads.mapped.Ensembl ~ stats$Tissue, col=pal[1:4], main="", ylab="Mapping efficiency (%)", notch=F, pch=16, boxwex=0.7)
dev2bitmap("boxplot_RNA-seq_mapping_efficiency_Ensembl.pdf", type="pdfwrite", method="pdf", w=6, h=6)
summary(lm(stats$Percentage.of.reads.mapped.Ensembl ~ stats$Tissue)) ## p-value: 0.9813
## split by species
boxplot(stats$Percentage.of.reads.mapped.Ensembl ~ interaction(stats$Tissue, stats$Species), col=pal[1:4], main="", ylab="Mapping efficiency (%)", notch=F, pch=16, boxwex=0.7, names=rep(c("heart", "kidney", "liver", "lung"), 3), at =c(1,2,3,4, 6,7,8,9, 11,12,13,14))
mtext(unique(stats$Species)[1], side = 1, line = 3, at=2.5)
mtext(unique(stats$Species)[2], side = 1, line = 3, at=7.5)
mtext(unique(stats$Species)[3], side = 1, line = 3, at=12.5)

## Are differences explained by technical factors?
cor.test(stats$Percentage.of.reads.mapped.Ensembl, stats$RIN, method="spearman") ## rho=-0.11 NS
cor.test(stats$Percentage.of.reads.mapped.Ensembl, stats$Percentage.of.bps.trimmed..adapters., method="spearman") ## rho=0.15 NS
cor.test(stats$Percentage.of.reads.mapped.Ensembl, stats$RNA.extraction.concentration..ng.uL., method="spearman") ## rho=0.41 p=0.0046
plot(stats$Percentage.of.reads.mapped.Ensembl, stats$RNA.extraction.concentration..ng.uL., pch=16, log="y")
cor.test(stats$Percentage.of.reads.mapped.Ensembl, stats$Concentration..ng.uL., method="spearman") ## rho=-0.38 p=0.008232
plot(stats$Percentage.of.reads.mapped.Ensembl, stats$Concentration..ng.uL., pch=16, log="y")
## Maybe some influence of the technical complexity of the libraries
## Anyway, after removing the species effects, the differences between samples should be very small

############
## BS-seq ##
############

## plot the extent of # reads (all samples in all lanes)

sum_all_samples <- counts[counts$Type=="BS-seq", ]
sum_all_samples <- sum_all_samples[order(sum_all_samples$ReadCount), ]

par(mar=c(12, 4, 4, 2) + 0.1) # increase margins
mp <- barplot(sum_all_samples$ReadCount, ylab="Number of sequenced reads", xlab="", xaxt="n")
text(mp, -5000000, srt = 45, adj = 1, labels = paste(sum_all_samples$SampleID, " (Lane ", sum_all_samples$Lane, ", Flow-cell ", sum_all_samples$Flow.cell,")", sep=""), xpd = TRUE, cex=0.7)

## aggregate by condition: how many sequenced samples?
aggregate(sum_all_samples$Condition, list(sum_all_samples$Condition), length)
## aggregate by condition: how many technical replicates? (some libraries were sequenced several times)
aggregate(sum_all_samples$SampleID, list(sum_all_samples$Condition), function(x){ return(length(unique(x))) })

## gather by sample
sum_samples <- aggregate(sum_all_samples$ReadCount, list(sum_all_samples$Condition), sum) ## 48 samples
## number of times we tried to sequence each sample (all technical replicates mixed)
num_seq <- aggregate(sum_all_samples$ReadCount, list(sum_all_samples$Condition), length)
num_seq <- num_seq[order(sum_samples[,2]), ]
## number of different technical replicates sequenced
num_samples <- aggregate(sum_all_samples$SampleID, list(sum_all_samples$Condition), function(x){ length(unique(x))})
num_samples <- num_samples[order(sum_samples[,2]), ]
sum_samples <- sum_samples[order(sum_samples[,2]), ]

par(mar=c(8, 4, 4, 2) + 0.1) # increase margins
mp <- barplot(sum_samples[,2], ylab="Number of sequenced reads", xlab="", xaxt="n", col="lightgrey")
text(mp, -10000000, srt = 45, adj = 1, labels = sum_samples[,1], xpd = TRUE, cex=0.8)
## adds the number of technical replicates
text(mp, 8000000, srt = 0, labels = num_samples[,2], xpd = TRUE, cex=0.7)
## ## arbitrary threshold of 4X = 12Gb = 240M reads. Not very meaningful on sequenced reads... See below with mapped/deduplicated reads
## abline(h=240e6, lty=2, col="blue")
## legend("top", "4X coverage", col="blue", lty=2, bty="n")

summary(sum_samples[,2] > 24e7) ## 44 samples
sum_samples[sum_samples[,2] > 24e7,1]
## save a list of all samples for which no more sequencing needs to be done, so that we can smooth them and perform QC (PCA)
## write.table(sum_samples[sum_samples[,2] > 24e7,1], file="../bsseq/samples_to_smooth.txt", col.names=F, sep="\t", row.names=F, quote=F)


## number of mapped reads
mapped <- read.table("../bismark/checked_samples.txt", h=T, sep="\t",comment.char="")
## remove commas
for (i in c(4,5,7,9)){
  mapped[,i] <- as.numeric(gsub(",","", mapped[,i]))
}
## sort by number of deduplicated reads
mapped <- mapped[order(mapped$Number.of.reads.after.deduplication),]
## sort by number of mapped reads
## mapped <- mapped[order(mapped$Number.of.mapped.reads),]
## sort by number of sequenced reads
## mapped <- mapped[order(mapped$Number.of.sequenced.reads),]

par(mar=c(8, 4, 4, 2) + 0.1) # increase margins
## mp <- barplot(mapped$Number.of.sequenced.reads, ylab="Number of sequenced reads", xlab="", xaxt="n")
mp <- barplot(as.matrix(rbind(mapped$Number.of.sequenced.reads, mapped$Number.of.mapped.reads, mapped$Number.of.reads.after.deduplication)), ylab="Number of reads", xlab="", xaxt="n", beside=T, col=pal[c(2,10,22)])
text(mp[2,], -5000000, srt = 45, adj = 1, labels = paste(mapped[,2], " (", mapped[,1], ")", sep=""), xpd = TRUE, cex=0.5)
legend("topleft", c("Sequenced", "Mapped", "After removing duplicates"), fill=pal[c(2,10,22)], bty="n")

## same plot by grouped by sample
sum_sequenced <- aggregate(mapped$Number.of.sequenced.reads, list(mapped$Condition), sum)
sum_mapped <- aggregate(mapped$Number.of.mapped.reads, list(mapped$Condition), sum)
sum_dedup <- aggregate(mapped$Number.of.reads.after.deduplication, list(mapped$Condition), sum)
## number of different technical replicates sequenced
num_samples <- aggregate(mapped$Number.of.sequenced.reads, list(mapped$Condition), function(x){ length(unique(x))})

## sort samples:
sum_sequenced <- sum_sequenced[order(sum_dedup[,2]),]
sum_mapped <- sum_mapped[order(sum_dedup[,2]),] 
num_samples <- num_samples[order(sum_dedup[,2]),] 
sum_dedup <- sum_dedup[order(sum_dedup[,2]),]

par(mar=c(8, 4, 4, 2) + 0.1) # increase margins
mp <- barplot(as.matrix(rbind(sum_sequenced[,2], sum_mapped[,2], sum_dedup[,2])), ylab="Number of reads", xlab="", xaxt="n", beside=T, col=pal[c(2,10,22)])
text(mp[2,], -8000000, srt = 45, adj = 1, labels = sum_sequenced[,1], xpd = TRUE, cex=1)
## adds the number of technical replicates
text(mp[2,], sum_sequenced[,2]+7000000, srt = 0, labels = num_samples[,2], xpd = TRUE, cex=1)
legend("topleft", c("Sequenced", "Mapped", "After removing duplicates"), lty=c(0, 0, 0), col=c("black", "black", "black"), bty="n", pch = c(22, 22, 22), pt.bg = c(pal[c(2,10,22)]), pt.cex = 2)
## TO DO: move legend up with Illustator?
dev2bitmap("barplot_BS-seq_mapping_rates.pdf", type="pdfwrite", method="pdf", w=10, h=6)

## arbitrary threshold of 4X = 12Gb = 240M reads
abline(h=240e6, lty=2, col=pal[11])
legend("topleft", c("Sequenced", "Mapped", "After removing duplicates", "~4X coverage"), lty=c(0, 0, 0, 2), col=c("black", "black", "black", pal[11]), bty="n", pch = c(22, 22, 22, NA), pt.bg = c(pal[c(2,10,22)], NA), pt.cex = 2)


## barplot of number of sites with >=2X coverage 
numberSites <- read.table("~/clusterhome/Methylation/bsseq/number_sites_used.txt", h=F, row.names=1)
numberSites <- numberSites[order(numberSites[,1]),, drop = FALSE]
par(mar=c(8, 4, 4, 2) + 0.1) # increase margins2
mp <- barplot(numberSites[,1], ylab="Number of sites covered >=2X", xlab="", xaxt="n", beside=T, , col="lightgrey")
text(mp, -500000, srt = 45, adj = 1, labels = row.names(numberSites), xpd = TRUE, cex=1)

## TO DO: barplot of mean coverage (see analysis.R)
load("../bsseq/smooth_data/combined_samples/combinedSmoothedCommonSites.RDa")
round(colMeans(getCoverage(allData.fit.subset)), 2)
summary(colMeans(getCoverage(allData.fit.subset)))
#...
