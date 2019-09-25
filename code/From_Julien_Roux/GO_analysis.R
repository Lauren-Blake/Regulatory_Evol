## Mar 24, 2016
## R commands can be run locally on macBook
source("functions.R")

library(RColorBrewer)
## display.brewer.all()
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))
## pal <- brewer.pal(9, "Set1")
## pal <- colorRampPalette(pal)(23)
## pie(rep(1,23), col=pal)

## sample info
samples <- read.table("../raw_data/SamplesDirectories.txt", sep="\t", h=T)
samples <- unique(samples[samples$Type == "RNA-seq" & samples$Flow.cell == "RNA-seq_1" , c(5,6,7,4)]) ## we din't need the whole file, just one flow cell with all the samples
samples <- unique(samples)
samples[,5] <- as.factor(unlist(lapply(strsplit(as.character(samples$Condition), split=""), function(x) { return(x[2]) })))
names(samples)[5] <- "IndividualID"
samples[,6] <- as.factor(unlist(lapply(strsplit(as.character(samples$SampleID), split="\\s+"), function(x) { return(x[2]) })))
names(samples)[6] <- "Individual"
samples <- samples[order(samples[,1]),]
row.names(samples) <- samples[,1]
samples <- samples[-17,]

#########################
## GO enrichment tests ##
#########################

## first we need the GO mapping: download it with BiomaRt
## library("biomaRt")
## listMarts(host='apr2013.archive.ensembl.org')[,1]
## ensembl71 <- useMart(host='apr2013.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
## listFilters(ensembl71)[grep("GO", listFilters(ensembl71)[,2]),]
## listAttributes(ensembl71)[grep("GO", listAttributes(ensembl71)[,2]),]
## GO <- getBM(attributes=c('go_id', 'name_1006', 'namespace_1003'), mart = ensembl71)
## save(GO, file="./GO_enrichment/GO.RDa")
load("GO_enrichment/GO.RDa")
dim(GO)
## 13,587 terms + root

## Ensembl2GO <- getBM(attributes=c('ensembl_gene_id', 'go_id'), filters = 'with_go_id', values = TRUE, mart = ensembl71)
## save(Ensembl2GO, file="GO_enrichment/Ensembl2GO.RDa")
load("GO_enrichment/Ensembl2GO.RDa")
dim(Ensembl2GO)
## 243,392 gene to term relations (20,695 genes involved)

gene2GO <- tapply(as.character(Ensembl2GO[,2]), as.character(Ensembl2GO[,1]), unique)
library(topGO)
GO2gene <- inverseList(gene2GO)

## topGO #################################################################
## make function to test all 3 ontologies (fisher + elim + weight), calculate FDR and return table of significant terms
## library(qvalue)
testAll <- function(geneList, gene2GO, outfile, nodesize){
  ## build all topGO objects
  GOdataBP <- new("topGOdata",
                  ontology = "BP",
                  allGenes = geneList, 
                  nodeSize = nodesize,
                  annot = annFUN.gene2GO,
                  gene2GO = gene2GO
                  )
  GOdataMF <- new("topGOdata",
                  ontology = "MF",
                  allGenes = geneList, 
                  nodeSize = nodesize,
                  annot = annFUN.gene2GO,
                  gene2GO = gene2GO
                  )
  GOdataCC <- new("topGOdata",
                  ontology = "CC",
                  allGenes = geneList, 
                  nodeSize = nodesize,
                  annot = annFUN.gene2GO,
                  gene2GO = gene2GO
                  )

  ##   ## test enrichement (elim)
  ##   resultsUpElimBP <- runTest(GOdataBP, algorithm = "elim", statistic = "fisher")
  ##   resultsUpElimMF <- runTest(GOdataMF, algorithm = "elim", statistic = "fisher")
  ##   resultsUpElimCC <- runTest(GOdataCC, algorithm = "elim", statistic = "fisher")

  ## test enrichement (weight)
  resultsUpWeightBP <- runTest(GOdataBP, algorithm = "weight", statistic = "fisher")
  resultsUpWeightMF <- runTest(GOdataMF, algorithm = "weight", statistic = "fisher")
  resultsUpWeightCC <- runTest(GOdataCC, algorithm = "weight", statistic = "fisher")

  ## Summary table:
  ## - only output weight scores
  ## - retrieve scores for all categories (length(score(resultsUp)))
  ## - display 100 characters for the name of the GO category
  ## - append the 3 ontologies
  allResUp <- GenTable(GOdataBP, "P-value" = resultsUpWeightBP, topNodes = length(score(resultsUpWeightBP)), numChar=100)
  allResUp <- rbind(allResUp, GenTable(GOdataMF, "P-value" = resultsUpWeightMF, topNodes = length(score(resultsUpWeightMF)), numChar=100))
  allResUp <- rbind(allResUp, GenTable(GOdataCC, "P-value" = resultsUpWeightCC, topNodes = length(score(resultsUpWeightCC)), numChar=100))
  allResUp$"P-value"[allResUp$"P-value"=="< 1e-30"] <- "1e-30"
  allResUp$"P-value" <- as.numeric(allResUp$"P-value")

  ## add column with name of ontology
  allResUp <- cbind(allResUp, c(rep("BP", length(score(resultsUpWeightBP))),rep("MF", length(score(resultsUpWeightMF))),rep("CC", length(score(resultsUpWeightCC)))))
  names(allResUp)[7] <- "Ontology"
  allResUp <- allResUp[,c(1,7,2,3,4,5,6)]
  ## sort by weight p-value
  allResUp <- allResUp[order(allResUp$"P-value"),]
  ## plot histogram of p-values
  hist(allResUp$"P-value", breaks=100, main="Histogram of p-values")

  ## FDR correction: on whole set of weight p-values
  FDR <- p.adjust(p=allResUp$"P-value", method = "fdr")
  allResUp <- cbind(allResUp, FDR)
    
  ## export output file
  if(!is.na(outfile)){
    write.table(allResUp, file=outfile, col.names=T, sep="\t", row.names=F, quote=F)
  }
  return(allResUp)
}

################################# Genes next to shared DMRs##############################
shared <- read.table("DMRs/species/HumanChimp_liver_DMRs_tissueShared_genes.txt", sep="\t", h=F)

## The universe: all genes next to DMRs
universe <- unique(shared[,1])
## The set of interest: all genes next to shared DMRs (TRUE)
geneList <- shared[shared[,2] == TRUE, 1]
## TO DO: we could also check that those genes are never seen next to tissue-specific DMRs
geneList <- as.factor(as.numeric(universe %in% geneList))
names(geneList) <- universe
summary(geneList)
##   0    1 
## 5332 1770 
results <- testAll(geneList, gene2GO, "GO_enrichment/HumanChimp_liver_DMRs_tissueShared_genes.txt", 10)

################################# Genes next to tDMRs / conserved tDMRs ##############################
## We will take the genes overlapping the promoter of a 1-to-1 gene
overlap <- read.table("../annotation/DMRs/tissues/Human_heart_liver_DMRs_conserved_features.gz", sep="\t", h=T)

## Foreground: all genes next to DMRs
## Extract all gene IDs in overlap$promoter
geneList <- unique(unlist(strsplit(as.character(overlap$promoter[!is.na(overlap$promoter)]), split=",")))
## Background: all 1-to-1 orth genes
universe <- unique(read.table("../orthoExon/metaOrthoExon.v2.hg19-panTro3-rheMac2.10_05_2011.txt", sep="\t", h=T)[,1])

geneList <- as.factor(as.numeric(as.character(universe) %in% geneList))
names(geneList) <- universe
summary(geneList)
results <- testAll(geneList, gene2GO, "GO_enrichment/Human_heart_liver_DMRs_genes.txt", 10)

## Now looking at enrichment of only genes next to conserved tDMRs 
geneList <- unique(unlist(strsplit(as.character(overlap$promoter_conserved_panTro3.rheMac2[!is.na(overlap$promoter_conserved_panTro3.rheMac2)]), split=",")))
geneList <- as.factor(as.numeric(as.character(universe) %in% geneList))
names(geneList) <- universe
summary(geneList)
results <- testAll(geneList, gene2GO, "GO_enrichment/Human_heart_liver_conserved_DMRs_genes.txt", 10)

## Same thing but taking all genes next to tDMRs as universe
geneList <- unique(unlist(strsplit(as.character(overlap$promoter_conserved_panTro3.rheMac2[!is.na(overlap$promoter_conserved_panTro3.rheMac2)]), split=",")))
universe <- unique(unlist(strsplit(as.character(overlap$promoter[!is.na(overlap$promoter)]), split=",")))
geneList <- as.factor(as.numeric(as.character(universe) %in% geneList))
names(geneList) <- universe
summary(geneList)
results <- testAll(geneList, gene2GO, "GO_enrichment/Human_heart_liver_conserved_DMRs_genes_universe_genes_next_to_DMRs.txt", 10)

## Repeat for all tissue-specific tDMRs 
universe <- unique(read.table("../orthoExon/metaOrthoExon.v2.hg19-panTro3-rheMac2.10_05_2011.txt", sep="\t", h=T)[,1])
tissues <- c("heartSpecific", "liverSpecific", "kidneySpecific", "lungSpecific")
for (tissue in tissues){
  cat("Testing tissue: ", tissue, "\n")
  ## Reading file of DMR overlapping promoters
  overlap <- read.table(paste0("../annotation/DMRs/tissues/Human_", tissue, "_DMRs_conserved_features.gz"), sep="\t", h=T)

  ## Foreground: all genes next to DMRs
  ## Extract all gene IDs in overlap$promoter
  genes <- unique(unlist(strsplit(as.character(overlap$promoter[!is.na(overlap$promoter)]), split=",")))
  ## Background: all 1-to-1 orth genes
  geneList <- as.factor(as.numeric(as.character(universe) %in% genes))
  names(geneList) <- universe
  summary(geneList)
  results <- testAll(geneList, gene2GO, paste0("GO_enrichment/genes_next_to_Human_", tissue, "_DMRs.txt"), 10)

  ## Now looking at enrichment of only genes next to conserved tDMRs 
  genesCons <- unique(unlist(strsplit(as.character(overlap$promoter_conserved_panTro3.rheMac2[!is.na(overlap$promoter_conserved_panTro3.rheMac2)]), split=",")))
  geneList <- as.factor(as.numeric(as.character(universe) %in% genesCons))
  names(geneList) <- universe
  summary(geneList)
  results <- testAll(geneList, gene2GO, paste0("GO_enrichment/genes_next_to_conserved_Human_", tissue, "_DMRs.txt"), 10)

  ## Same thing but taking all genes next to tDMRs as universe
  geneList <- as.factor(as.numeric(as.character(genes) %in% genesCons))
  names(geneList) <- genes
  summary(geneList)
  results <- testAll(geneList, gene2GO, paste0("GO_enrichment/genes_next_to_conserved_Human_", tissue, "_DMRs_universe_genes_next_to_Human_", tissue, "_DMRs.txt"), 10)
}

## TO DO do not look at genes wiht tDMRs overlapping promoter, but look at genes closest to any tDMrs (record distance too)
