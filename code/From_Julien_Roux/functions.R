## May 20, 2015
## Julien Roux
## File of (custom-made) functions used in ~/Methylation/bsseq scripts

#########
## PCA ##
#########

## plotting function
plot_scores <- function(pca, scores, n, m, cols, points=F, pchs, legend=F){
  xmin <- min(scores[,n]) - (max(scores[,n]) - min(scores[,n]))*0.05
  if (legend == T){ ## let some room (35%) for a legend
    xmax <- max(scores[,n]) + (max(scores[,n]) - min(scores[,n]))*0.35
  }
  else {
    xmax <- max(scores[,n]) + (max(scores[,n]) - min(scores[,n]))*0.05
  }
  ymin <- min(scores[,m]) - (max(scores[,m]) - min(scores[,m]))*0.05
  ymax <- max(scores[,m]) + (max(scores[,m]) - min(scores[,m]))*0.05
  
  plot(scores[,n], scores[,m], xlab=paste("PC", n, ": ", round(summary(pca)$importance[2,n],3)*100, "% variance explained", sep=""), ylab=paste("PC", m, ": ", round(summary(pca)$importance[2,m],3)*100, "% variance explained", sep=""), xlim=c(xmin, xmax), ylim=c(ymin, ymax), type="n")

  if (points == F){
    text(scores[,n],scores[,m], labels=rownames(scores), col=cols, cex=1)
  }
  else {
    points(scores[,n],scores[,m], col=cols, pch=pchs, cex=1.3)
  }
}

## perform a PCA on a subset of a table
subset_PCA <- function(tab, condition, n, pdfname){
  print("Lauching PCA...")
  print(paste0(length(tab[condition, 1]), " sites tested"))
  pca1 <- prcomp(t(tab[condition, ]), scale = FALSE) 
  print(paste0(summary(pca1)$importance[3,3]*100, "% of variance explained by the 3 first components"))
  loadings <- pca1$rotation
  scores <- pca1$x
  rownames(scores) <- colnames(tab)

  print("Plotting PCA...")
  pdf(file = pdfname, width = 6, height = 6)
  for (i in 1:n){
    ## colors for tissues
    col.v <- pal[as.integer(info$Tissue)]
    ## symbols for species (not used)
    pchs <- c(15,16,17)[as.integer(info$Species)]
    plot_scores(pca1, scores, i, i+1, cols=col.v, pch=pchs, points=F, legend=F)
  }
  dev.off()
  print("Done")
}

##############
## Heatmaps ##
##############
library(gplots)
library(RColorBrewer)
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))

## perform a heatmap on correlation matrix of a all or a subset of a methylation sites
subset_heatmap <- function(tab, condition, colors, info, pdfname){
  print(paste0("Calculating spearman correlation matrix using ", length(tab[condition, 1]), " sites..."))
  cors <- cor(tab[condition,], method="spearman", use="pairwise.complete.obs") 

  print("Plotting heatmap...")
  pdf(file = pdfname, width = 12, height = 8)
  
  labels <- paste(info$Species, info$Tissue, sep=" ")
  heatmap.2( cors, scale="none", col = colors, margins = c(12, 12), trace='none', denscol="white", labCol=labels, ColSideColors=pal[as.integer(as.factor(info$Species))+9], RowSideColors=pal[as.integer(as.factor(info$Tissue))], cexCol = 1.5)
  dev.off()
  print("Done")
}

subset_heatmap_random_sites <- function(tab, condition, colors, info, pdfname, n=10000){
  print(paste0("Calculating spearman correlation matrix using ", n, " random sites out of ", length(tab[condition, 1]), " sites..."))
  subset <- sample(row.names(tab[condition,]), n, replace=FALSE)
  write.table(subset, file=gsub("pdf", "txt", pdfname), sep="\t", quote=F, row.names=F, col.names=F)

  cors <- cor(tab[subset,], method="spearman", use="pairwise.complete.obs") 

  print("Plotting heatmap...")
  pdf(file = pdfname, width = 12, height = 8)
  
  labels <- paste(info$Species, info$Tissue, sep=" ")
  heatmap.2( cors, scale="none", col = colors, margins = c(12, 12), trace='none', denscol="white", labCol=labels, ColSideColors=pal[as.integer(as.factor(info$Species))+9], RowSideColors=pal[as.integer(as.factor(info$Tissue))], cexCol = 1.5)
  dev.off()
  print("Done")
}

reuse_subset_heatmap_random_sites <- function(tab, filename, colors, info, pdfname){
  print(paste0("Reading the random sites saved in ", filename))
  subset <- read.table(filename, sep="\t", h=F)[,1]
  print("Calculating correlation matrix")
  cors <- cor(tab[as.character(subset),], method="spearman", use="pairwise.complete.obs") 

  print("Plotting heatmap...")
  pdf(file = pdfname, width = 12, height = 8)
  
  labels <- paste(info$Species, info$Tissue, sep=" ")
  heatmap.2( cors, scale="none", col = colors, margins = c(12, 12), trace='none', denscol="white", labCol=labels, ColSideColors=pal[as.integer(as.factor(info$Species))+9], RowSideColors=pal[as.integer(as.factor(info$Tissue))], cexCol = 1.5)
  dev.off()
  print("Done")
}


plot_dmrs <- function(contrast){
  ## Load a DMR file and record coordinates, for example: 
  dmrs <- read.table(paste0("./DMRs/", contrast, "_DMRs.txt"), sep="\t", h=T)
  ## create GRanges object
  library("GenomicRanges")
  dmrsGr <- GRanges(dmrs$chr, IRanges(dmrs$start, dmrs$end))

  if (file.exists(paste0("DMRs/", contrast, "_DMRs_mean_Methylation.RDa"))){
    print("RDa file found") 
    load(paste0("DMRs/", contrast, "_DMRs_mean_Methylation.RDa"))
  } else {
    print("RDa file does not exist") 
    
    ## Get mean methylation at DMR positions for all samples
    meanMeth <- matrix(nrow=length(dmrs[,1]), ncol=48)
    colnames(meanMeth) <- unique(samples$Condition)
    for (sample in unique(samples$Condition)){
      cat("  ", sample, "\n")
      ## load sample and calculate the mean methylation at each DMR location
      load(paste0("./smooth_data/smoothed_samples/", sample, "_Smoothed.RDa"))
      meanMeth[,sample] <- getMeth(data.fit, regions = dmrsGr, type = "smooth", what = "perRegion")
    }
    ## Save object for potential future use
    save(meanMeth, file=paste0("DMRs/", contrast, "_DMRs_mean_Methylation.RDa"))
  }    
  
  library(gplots)
  colors <- colorRampPalette(c(brewer.pal(9, "Blues")[1],brewer.pal(9, "Blues")[9]))(100)
  ## Create breaks manually to span full range of methylation values
  palette.breaks <- seq(0, 1, 0.01)

  ## plot heatmap
  print("Plotting heatmap...")
  pdf(file = paste0("DMRs/", contrast, "_DMRs_heatmap.pdf"), width = 12, height = 8)
  h <- heatmap.2( na.omit(meanMeth), scale="none", col = colors, breaks = palette.breaks, margins = c(12, 12), trace='none', denscol="white", labCol=paste(samples$Species, samples$Tissue, sep=" "), labRow=NA, ColSideColors=pal[as.integer(as.factor(samples$Species))+9], cexCol = 1.5, dendrogram='col', key=TRUE)
  dev.off()  
}

plot_dmrs_correlation <- function(contrast){
  ## Load a DMR file and record coordinates, for example: 
  dmrs <- read.table(paste0("./DMRs/", contrast, "_DMRs.txt"), sep="\t", h=T)
  ## create GRanges object
  library("GenomicRanges")
  dmrsGr <- GRanges(dmrs$chr, IRanges(dmrs$start, dmrs$end))
    
  if (file.exists(paste0("DMRs/", contrast, "_DMRs_mean_Methylation.RDa"))){
    print("RDa file found") 
    load(paste0("DMRs/", contrast, "_DMRs_mean_Methylation.RDa"))
  } else {
    print("RDa file does not exist") 
    
    ## Get mean methylation at DMR positions for all samples
    meanMeth <- matrix(nrow=length(dmrs[,1]), ncol=48)
    colnames(meanMeth) <- unique(samples$Condition)
    for (sample in unique(samples$Condition)){
      cat("  ", sample, "\n")
      ## load sample and calculate the mean methylation at each DMR location
      load(paste0("./smooth_data/smoothed_samples/", sample, "_Smoothed.RDa"))
      meanMeth[,sample] <- getMeth(data.fit, regions = dmrsGr, type = "smooth", what = "perRegion")
    }
    ## Save object for potential future use
    save(meanMeth, file=paste0("DMRs/", contrast, "_DMRs_mean_Methylation.RDa"))
  }
  
  ## calculate correlation matrix
  cors <- cor(meanMeth, method="spearman", use="pairwise.complete.obs") 

  ## plot heatmap
  print("Plotting heatmap...")
  pdf(file = paste0("DMRs/", contrast, "_DMRs_heatmap_correlation.pdf"), width = 12, height = 8)
  h <- heatmap.2( cors, scale="none", col = colors, margins = c(12, 12), trace='none', denscol="white", labCol=paste(samples$Species, samples$Tissue, sep=" "), labRow=samples$Condition, ColSideColors=pal[as.integer(as.factor(samples$Species))+9], RowSideColors=pal[as.integer(as.factor(samples$Tissue))], cexCol = 1.5)
  dev.off()
  print("Done")
}

plot_dmrs_correlation_conserved <- function(contrast, conservationDegree){
  ## Load a DMR file and record coordinates, for example: 
  dmrs <- read.table(paste0("./DMRs/", contrast, "_DMRs.txt"), sep="\t", h=T)
  ## create GRanges object
  library("GenomicRanges")
  dmrsGr <- GRanges(dmrs$chr, IRanges(dmrs$start, dmrs$end))
    
  if (file.exists(paste0("DMRs/", contrast, "_DMRs_mean_Methylation.RDa"))){
    print("RDa file found") 
    load(paste0("DMRs/", contrast, "_DMRs_mean_Methylation.RDa"))
  } else {
    print("RDa file does not exist") 
    
    ## Get mean methylation at DMR positions for all samples
    meanMeth <- matrix(nrow=length(dmrs[,1]), ncol=48)
    colnames(meanMeth) <- unique(samples$Condition)
    for (sample in unique(samples$Condition)){
      cat("  ", sample, "\n")
      ## load sample and calculate the mean methylation at each DMR location
      load(paste0("./smooth_data/smoothed_samples/", sample, "_Smoothed.RDa"))
      meanMeth[,sample] <- getMeth(data.fit, regions = dmrsGr, type = "smooth", what = "perRegion")
    }
    ## Save object for potential future use
    save(meanMeth, file=paste0("DMRs/", contrast, "_DMRs_mean_Methylation.RDa"))
  }

  ## Load DMRs conservation
  cons <- read.table(paste0("DMRs/", contrast, "_DMRs_conservation.txt"), h=T)
  ## Subset data 
  meanMeth <- meanMeth[cons$conservation == conservationDegree,]
  
  ## calculate correlation matrix
  cors <- cor(meanMeth, method="spearman", use="pairwise.complete.obs") 

  ## plot heatmap
  print("Plotting heatmap...")
  pdf(file = paste0("DMRs/", contrast, "_DMRs_heatmap_correlation_conserved_", conservationDegree, ".pdf"), width = 12, height = 8)
  h <- heatmap.2( cors, scale="none", col = colors, margins = c(12, 12), trace='none', denscol="white", labCol=paste(samples$Species, samples$Tissue, sep=" "), labRow=samples$Condition, ColSideColors=pal[as.integer(as.factor(samples$Species))+9], RowSideColors=pal[as.integer(as.factor(samples$Tissue))], cexCol = 1.5)
  dev.off()
  print("Done")
}
## We usually see a lot more tissue clustering. Makes sense since we look at regions wiht tissue differences conserved acrosss species. Still this is a positive control showing that there is not too much batch effects across species
## TO DO: mean methylation at minimal overlapping region between tDMRs of different species
##        This is a lot more complex to implement, and I don't think the result will be more informative than taking the mean methylation over whole DMR length (if any difference, the species effects would be even smaller)

##############
## Graphics ##
##############
## Taken from https://github.com/mylesmharrison/colorRampPaletteAlpha/blob/master/colorRampPaletteAlpha.R
## See http://www.everydayanalytics.ca/2014/03/colorRampPalette-alpha-in-R.html
addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

## add density plots in margins of figure
## Adapted from function found on http://sas-and-r.blogspot.ch/2012/09/example-103-enhanced-scatterplot-with.html
scatterhist <- function(x, y, xlab = "", ylab = "", xsize=1, cleanup=TRUE,...){
  # save the old graphics settings-- they may be needed
  def.par <- par(no.readonly = TRUE)
  
  zones <- matrix(c(0,4,0,
                    1,5,3,
                    0,2,0), ncol = 3, byrow = TRUE)
  layout(zones, widths=c(0.5,4,1), heights = c(5,16,2))
  ## layout.show(n=5)

  
  # tuning to plot histograms nicely
  xhist <- hist(x, plot = FALSE, breaks=20)
  yhist <- hist(y, plot = FALSE, breaks=20)
  top <- max(c(xhist$counts, yhist$counts))
  
  # for all titles: 
  # drop the axis titles and omit boxes, set up margins
  par(xaxt="n", yaxt="n", bty="n",  mar = c(.3,2,.3,0) +.05)
  # fig 1
  # drop the axis titles, set up margins
  plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1))
  text(0,0,paste(ylab), cex=1.5, srt=90)

  # fig 2
  plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1))
  text(0,0,paste(xlab), cex=1.5)
  
  # fig 3, the first histogram, needs different margins
  # no margin on the left
  par(mar = c(2,0,0,1)) 
  barplot(yhist$counts, axes = FALSE, xlim = c(0, top), space = 0, horiz = TRUE, col="white")

  # fig 4, other histogram needs no margin on the bottom
  par(mar = c(0,2,1,0)) ##bottom, left, top, and right
  barplot(xhist$counts, axes = FALSE, ylim = c(0, top), space = 0, col="white")
  
  # fig 5, finally, the scatterplot-- needs regular axes, different margins
  par(mar = c(2,2,.5,.5), xaxt="s", yaxt="s", bty="o")
  # this color allows transparency & overplotting-- useful if a lot of points
  ## plot(x, y , pch=19, col="#00000022", cex=xsize, ...)
  # Take what I use usually instead:
  plot(x, y , pch=16, col=rgb(0,0,0,0.1), cex=xsize, ...)

  # reset the graphics, if desired 
  if(cleanup) {
    par(def.par)
  }

  ## TO DO:
  ## - how to stick histograms to box of scatterplot?
  ## - same function with density plot instead of histogram
  ## - there is a subtle shift if we compare the histogram to the scatterplot: the histograms to not cover the full span of the scatterplot data: how to fix this: xlim chnaged to 0.95*xlim? 
}

##########
## DMRs ##
##########

## species DMRs
plot_speciesDMRs_characteristics2 <- function(baseName, DMRcolumn, legend=T){
  library(RColorBrewer)
  ## display.brewer.all()
  pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))
  
  ## species DMRs, grouped by tissue
  speciesPairs <- list(c("Human", "Chimp"), c("Human", "Rhesus"), c("Chimp", "Rhesus"), c("Human", "Specific"), c("Chimp", "Specific"))
  allTissues <- c("heart", "kidney", "liver", "lung")
  
  pdf(file = paste0(baseName, ".pdf"), width = 6, height = 5)
  i <- 1
  allLegendText <- vector()
  for (tissue in allTissues){
    j <- 1
    for (pair in speciesPairs){
      cat("Testing species pair: ", pair, "in", tissue, "\n")

      ## Reading the DMR file
      dmrs <- read.table(paste0("./DMRs/species/", pair[1], pair[2], "_", tissue, "_DMRs.txt"), sep="\t", h=T)   
      cat("  ", length(dmrs[,1]), "DMRs read\n")
      
      if (i == 1 & j == 1){
        ## In case, take 50% larger y-axis limits
        ylim <- c(0, max(density(dmrs[,DMRcolumn], adjust=1)$y)*1.5)
        plot(density(dmrs[,DMRcolumn], adjust=1), col=pal[i], main="", xlab=paste0("DMR ", DMRcolumn), lwd=1, lty=j, ylim=ylim)
      }
      else {
        lines(density(dmrs[,DMRcolumn], adjust=1), col=pal[i], lwd=1, lty=j)
      }
      j <- j + 1
    }
    i <- i+1
    ## legend
    legendText <- paste0(do.call(rbind, speciesPairs)[,1], " vs. ", do.call(rbind, speciesPairs)[,2], " in ", tissue)
    legendText <- gsub(" vs. Specific", "-specific", legendText, perl=TRUE)
    allLegendText <- c(allLegendText, legendText)
  }
  if (legend == TRUE){
    legend("topright", legend=allLegendText, lty=rep(c(1:5), 4), col=rep(pal[1:4], each=5)) 
  }
  dev.off()
}

plot_speciesDMRs_characteristics <- function(baseName, DMRcolumn){
  ## species DMRs, grouped by tissue
  speciesPairs <- list(c("Human", "Chimp"), c("Human", "Rhesus"), c("Chimp", "Rhesus"), c("Human", "Specific"), c("Chimp", "Specific"))
  allTissues <- c("heart", "kidney", "liver", "lung")
  
  pdf(file = paste0(baseName, "_by_tissue.pdf"), width = 6, height = 5)
  for (tissue in allTissues){
    i <- 1
    for (pair in speciesPairs){
      cat("Testing species pair: ", pair, "in", tissue, "\n")

      ## We can start by reading the output file: no need to recalculate all DMRs
      dmrs <- read.table(paste0("~/clusterhome/Methylation/bsseq/DMRs/species/", pair[1], pair[2], "_", tissue, "_DMRs.txt"), sep="\t", h=T)   
      cat("  ", length(dmrs[,1]), "DMRs read\n")
      
      if (i == 1){
        ylim <- c(0, max(density(dmrs[,DMRcolumn], adjust=1)$y)*1.5)
        plot(density(dmrs[,DMRcolumn], adjust=1), col=pal[i], main="", xlab=paste0("DMR ", DMRcolumn), lwd=1, lty=i, ylim=ylim)
        ## plot legend
        legendText <- paste0(do.call(rbind, speciesPairs)[,1], " vs. ", do.call(rbind, speciesPairs)[,2], " in ", tissue)
        legendText <- gsub(" vs. Specific", "-specific", legendText, perl=TRUE)

        legend("topright", legend=legendText, lty=c(1:5), col=pal[1:5]) 
      }
      else {
        lines(density(dmrs[,DMRcolumn], adjust=1), col=pal[i], lwd=1, lty=i)
      }
      i <- i+1
    }
  }
  dev.off()

  ## species DMRs, grouped by contrast
  pdf(file = paste0(baseName, "_by_species_contrast.pdf"), width = 6, height = 5)
  for (pair in speciesPairs){
    i <- 1
    for (tissue in allTissues){
      cat("Testing species pair: ", pair, "in", tissue, "\n")

      ## We can start by reading the output file: no need to recalculate all DMRs
      dmrs <- read.table(paste0("~/clusterhome/Methylation/bsseq/DMRs/species/", pair[1], pair[2], "_", tissue, "_DMRs.txt"), sep="\t", h=T)   
      cat("  ", length(dmrs[,1]), "DMRs read\n")
      
      if (i == 1){
        ylim <- c(0, max(density(dmrs[,DMRcolumn], adjust=1)$y)*1.1)
        plot(density(dmrs[,DMRcolumn], adjust=1), col=pal[i], main="", xlab=paste0("DMR ", DMRcolumn), lwd=1, lty=i, ylim=ylim)
        ## plot legend
        legendText <- paste0(pair[1], " vs. ", pair[2], " in ", allTissues)
        legendText <- gsub(" vs. Specific", "-specific", legendText, perl=TRUE)

        legend("topright", legend=legendText, lty=c(1:5), col=pal[1:5]) 
      }
      else {
        lines(density(dmrs[,DMRcolumn], adjust=1), col=pal[i], lwd=1, lty=i)
      }
      i <- i+1
    }
  }
  dev.off()
}

## tissue DMRs
plot_tissueDMRs_characteristics2 <- function(baseName, DMRcolumn, legend=T, scaling){
  library(RColorBrewer)
  ## display.brewer.all()
  pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))
  ##  lineTypes <- c(1:6, "1F", "F1", "4C88C488", "12345678")
  lineTypes <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "1F", "F1", "4C88C488", "12345678")
  ## 1F": dash length 1, gap length F (15)
  ## "F1": dash length F (15), gap length 1
  ## "4C88C488": dash (4), gap (C=12), dash (8), gap (8), dash (C=12), ...
  ## "12345678": dash (1), gap (2), dash (3), gap (4), ...
  ## http://stackoverflow.com/questions/25788945/how-to-define-more-line-types-for-graphs-in-r
  ## Problem not all seem to be workign with Illustrator (8 characters ones) 
  lineTypes <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "1A", "F1", "5181B1", "113151")

  ## tissue DMRs, grouped by species
  allSpecies <- c("Chimp", "Human", "Rhesus")
  tissuePairs <- list(c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "liver"), c("heart", "Specific"), c("liver", "Specific"), c("kidney", "Specific"), c("lung", "Specific"))
  
  pdf(file = paste0(baseName, ".pdf"), width = 6, height = 5)
  i <- 1
  allLegendText <- vector()
  for (species in allSpecies){
    j <- 1
    for (pair in tissuePairs){
      cat("Testing tissue pair: ", pair, "in", species, "\n")

      ## Reading DMR file
      if (pair[2] == "Specific"){
        dmrs <- read.table(paste0("./DMRs/tissues/", species, "_", pair[1], pair[2], "_DMRs.txt"), sep="\t", h=T)
      }
      else {
        dmrs <- read.table(paste0("./DMRs/tissues/", species, "_", pair[1], "_", pair[2], "_DMRs.txt"), sep="\t", h=T)   
      }
      cat("  ", length(dmrs[,1]), "DMRs read\n")
      
      if (i == 1 & j == 1){
        ## In case, take larger y-axis limits
        ylim <- c(0, max(density(dmrs[,DMRcolumn], adjust=1)$y)*scaling)
        plot(density(dmrs[,DMRcolumn], adjust=1), col=pal[i+9], main="", xlab=paste0("DMR ", DMRcolumn), lwd=1, lty=lineTypes[j], ylim=ylim)
      }
      else {
        lines(density(dmrs[,DMRcolumn], adjust=1), col=pal[i+9], lwd=1, lty=lineTypes[j])
      }
      j <- j + 1
    }
    i <- i + 1
    ## legend
    legendText <- paste0(do.call(rbind, tissuePairs)[,1], " vs. ", do.call(rbind, tissuePairs)[,2], " in ", species)
    legendText <- gsub(" vs. Specific", "-specific", legendText, perl=TRUE)
    allLegendText <- c(allLegendText, legendText)
  }
  if (legend == TRUE){
    ## print(allLegendText)
    legend("topright", legend=allLegendText, lty=rep(lineTypes, 3), col=rep(pal[10:12], each=10)) 
  }
  dev.off()
}

plot_tissueDMRs_characteristics <- function(baseName, DMRcolumn){
  ## tissue DMRs, grouped by species
  allSpecies <- c("Human", "Chimp", "Rhesus")
  tissuePairs <- list(c("heart", "lung"), c("heart", "kidney"), c("liver", "lung"), c("liver", "kidney"), c("lung", "kidney"), c("heart", "liver"), c("heart", "Specific"), c("liver", "Specific"), c("kidney", "Specific"), c("lung", "Specific"))

  pdf(file = paste0(baseName, "_by_species.pdf"), width = 6, height = 5)
  for (species in allSpecies){
    i <- 1
    for (pair in tissuePairs){
      cat("Testing tissue pair: ", pair, "in", species, "\n")

      ## We can start by reading the output file: no need to recalculate all DMRs
      if (pair[2] == "Specific"){
        dmrs <- read.table(paste0("~/clusterhome/Methylation/bsseq/DMRs/tissues/", species, "_", pair[1], pair[2], "_DMRs.txt"), sep="\t", h=T)
      }
      else {
        dmrs <- read.table(paste0("~/clusterhome/Methylation/bsseq/DMRs/tissues/", species, "_", pair[1], "_", pair[2], "_DMRs.txt"), sep="\t", h=T)   
      }
      cat("  ", length(dmrs[,1]), "DMRs read\n")
      
      if (i == 1){
        ylim <- c(0, max(density(dmrs[,DMRcolumn], adjust=1)$y)*1.5)
        plot(density(dmrs[,DMRcolumn], adjust=1), col=pal[i], main="", xlab=paste0("DMR ", DMRcolumn), lwd=1, lty=i, ylim=ylim)
        ## plot legend
        legendText <- paste0(do.call(rbind, tissuePairs)[,1], " vs. ", do.call(rbind, tissuePairs)[,2], " in ", species)
        legendText <- gsub(" vs. Specific", "-specific", legendText, perl=TRUE)

        legend("topright", legend=legendText, lty=c(1:10), col=pal[1:10]) 
      }
      else {
        lines(density(dmrs[,DMRcolumn], adjust=1), col=pal[i], lwd=1, lty=i)
      }
      i <- i+1
    }
  }
  dev.off()

  ## tissues DMRs, grouped by contrast
  pdf(file = paste0(baseName, "_by_tissue_contrast.pdf"), width = 6, height = 5)
  for (pair in tissuePairs){
    i <- 1
    for (species in allSpecies){
      cat("Testing tissue pair: ", pair, "in", species, "\n")

      ## We can start by reading the output file: no need to recalculate all DMRs
      if (pair[2] == "Specific"){
        dmrs <- read.table(paste0("~/clusterhome/Methylation/bsseq/DMRs/tissues/", species, "_", pair[1], pair[2], "_DMRs.txt"), sep="\t", h=T)
      }
      else {
        dmrs <- read.table(paste0("~/clusterhome/Methylation/bsseq/DMRs/tissues/", species, "_", pair[1], "_", pair[2], "_DMRs.txt"), sep="\t", h=T)   
      }
      cat("  ", length(dmrs[,1]), "DMRs read\n")
      
      if (i == 1){
        ylim <- c(0, max(density(dmrs[,DMRcolumn], adjust=1)$y)*1.1)
        plot(density(dmrs[,DMRcolumn], adjust=1), col=pal[i], main="", xlab=paste0("DMR ", DMRcolumn), lwd=1, lty=i, ylim=ylim)
        ## plot legend
        legendText <- paste0(pair[1], " vs. ", pair[2], " in ", allSpecies)
        legendText <- gsub(" vs. Specific", "-specific", legendText, perl=TRUE)

        legend("topright", legend=legendText, lty=c(1:10), col=pal[1:10]) 
      }
      else {
        lines(density(dmrs[,DMRcolumn], adjust=1), col=pal[i], lwd=1, lty=i)
      }
      i <- i+1
    }
  }
  dev.off()
}

plot_individualDMRs_characteristics <- function(baseName, DMRcolumn){
  ## species DMRs, grouped by individual
  allSpecies <- c("Human", "Chimp", "Rhesus")
  individualPairs <- list(c("1", "2"), c("1", "3"), c("1", "4"), c("2", "3"), c("2", "4"), c("3", "4"))
                          
  pdf(file = paste0(baseName, "_by_species.pdf"), width = 6, height = 5)
  for (species in allSpecies){
    i <- 1
    for (pair in individualPairs){
      cat("Testing individual pair: ", pair, "in", species, "\n")
    ## baseName <- paste0(species, "_", pair[1], "_", pair[2])

      
      ## We can start by reading the output file: no need to recalculate all DMRs
      if (pair[2] == "Specific"){
        dmrs <- read.table(paste0("~/clusterhome/Methylation/bsseq/DMRs/individuals/", species, "_", pair[1], pair[2], "_DMRs.txt"), sep="\t", h=T)
      }
      else {
        dmrs <- read.table(paste0("~/clusterhome/Methylation/bsseq/DMRs/individuals/", species, "_", pair[1], "_", pair[2], "_DMRs.txt"), sep="\t", h=T)   
      }
      cat("  ", length(dmrs[,1]), "DMRs read\n")
      
      if (i == 1){
        ylim <- c(0, max(density(dmrs[,DMRcolumn], adjust=1)$y)*1.5)
        plot(density(dmrs[,DMRcolumn], adjust=1), col=pal[i], main="", xlab=paste0("DMR ", DMRcolumn), lwd=1, lty=i, ylim=ylim)
        ## plot legend
        legendText <- paste0(do.call(rbind, individualPairs)[,1], " vs. ", do.call(rbind, individualPairs)[,2], " in ", species)
        legendText <- gsub(" vs. Specific", "-specific", legendText, perl=TRUE)

        legend("topright", legend=legendText, lty=c(1:6), col=pal[1:6]) 
      }
      else {
        lines(density(dmrs[,DMRcolumn], adjust=1), col=pal[i], lwd=1, lty=i)
      }
      i <- i+1
    }
  }
  dev.off()

  ## individuals DMRs, grouped by contrast: no meaningful since individuals differ between species
}
