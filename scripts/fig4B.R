library(DESeq2)
library(tidyr)
library(ggplot2)
library(patchwork)
library(VennDiagram) 
library(data.table)
library(sva)
library(patchwork)
rm(list=ls())

### Genes ####

countdata <- read.csv(file = "dat/pro_raw_counts.txt",
                      sep="\t",header=TRUE,skip=1)

rownames(countdata) <- countdata$Geneid

countdata <- countdata[,-c(1:6)]
metadata <- read.csv(file = "dat/metadata.csv",sep = ",", header=TRUE)
metadata$samplecols <- colnames(countdata)

metadata$mutant <- factor(metadata$mutant,levels = c("wt","sof1","sof2","sof3","dbm"))

# Ensure columns are all present in metadata (and any subsets, if added in)
countdata <- countdata[,metadata$samplecols]

# Pre-filter low counts
countdata <- countdata[rowSums(countdata >= 10) >= 2, ]

dds <- DESeqDataSetFromMatrix(countData = countdata, colData = metadata,
                              design =~mutant*treatment)

dds <- DESeq(dds)

generate_plot = function(interaction="coef",coef){

  # WT results
  if (interaction=="coef") {
    res <- lfcShrink(dds, coef = coef, type = "ashr")
  }
  else{
    res <- lfcShrink(dds, contrast = list(c("treatment_e2_vs_dmso",coef)), type = "ashr")
    
  }
  resdata <- as.data.frame(res)
  resdata$refseq <- rownames(resdata)
  resdata <- separate_wider_delim(resdata, cols = refseq, delim = ".", names = c("refid", "iso"))
  refseq_to_cid <- read.table(file = "dat/refseq_to_common_id.txt",
                              sep="\t",header=FALSE)
  resdata_id <- merge(resdata, refseq_to_cid,by.x="refid",by.y="V1")
  resdata_id <- resdata_id[!(duplicated(resdata_id$V2)),]
  upgenes <- resdata_id[resdata_id$padj < 0.1 & resdata_id$log2FoldChange > 0,]
  downgenes <- resdata_id[resdata_id$padj < 0.1 & resdata_id$log2FoldChange < 0,]
  res_plotting <- resdata_id
  res_plotting$reg <- ifelse(res_plotting$log2FoldChange <= 0,"down","up")
  size_bound <- 2
  padj_bound <- 1e-10
  
  # Bound data for plotting purposes.
  res_plotting$significant <- ifelse(res_plotting$padj > 0.1, FALSE,TRUE)
  res_plotting$outlier <- ifelse(abs(res_plotting$log2FoldChange) > size_bound | 
                                   res_plotting$padj < padj_bound, "outlier", "normal")
  
  res_plotting$log2FoldChange[res_plotting$log2FoldChange > size_bound] <- size_bound
  res_plotting$log2FoldChange[res_plotting$log2FoldChange < -size_bound] <- -size_bound
  
  res_plotting$padj[res_plotting$padj < padj_bound] <- padj_bound
  
  res_plotting$SignificanceRegulation <- with(res_plotting, 
                                              ifelse(significant & reg == "up", "Upregulated",
                                                     ifelse(significant & reg == "down", "Downregulated", "Not Significant")))
  
  volcanoplot <- ggplot(data = res_plotting, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(size = 1.5, alpha = 1.5, aes(color = SignificanceRegulation, shape = outlier)) +
    scale_color_manual(
      values = c("Not Significant" = "darkgrey",
                 "Upregulated" = "red",
                 "Downregulated" = "blue"), 
      drop = FALSE # Prevents ggplot from automatically dropping unused levels
    ) +
    xlab(NULL) +
    ylab(" ") +
    xlim(-size_bound,size_bound) +
    theme_classic(base_size = 14) +
    theme(legend.position = "none")
  
  barchart <- ggplot(data = na.omit(res_plotting[res_plotting$padj < 0.1,]), aes(x=reg)) +
    geom_bar(color="red",fill="red") +
    xlab("") +
    ylab("Significant Genes") +
    theme_classic(base_size=8) 
  
chart <- volcanoplot + patchwork::inset_element(barchart, .8, .6, 1, 1)
return(chart)
}

wt <- generate_plot(interaction="coef",coef="treatment_e2_vs_dmso")
sof3 <- generate_plot(interaction="contrast",coef="mutantsof3.treatmente2")
dbm <- generate_plot(interaction="contrast",coef="mutantdbm.treatmente2")

wt / sof3 / dbm
