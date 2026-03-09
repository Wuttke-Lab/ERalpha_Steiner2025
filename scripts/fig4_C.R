# Fig 4C

library(DESeq2)
library(tidyr)
library(ggplot2)
library(patchwork)
library(VennDiagram) 
library(data.table)
library(sva)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)
rm(list=ls())
### Genes ####



countdata <- read.csv(file = "dat/pro_raw_counts.txt",
                      sep="\t",header=TRUE,skip=1)

rownames(countdata) <- countdata$Geneid
countdata <- countdata[,-c(1:6)]
metadata <- read.csv(file = "dat/metadata.csv",sep = ",", header=TRUE)
metadata$samplecols <- colnames(countdata)
metadata$mutant <- factor(metadata$mutant,levels = c("wt","sof1","sof2","sof3","dbm"))
countdata <- countdata[,metadata$samplecols]
countdata <- countdata[rowSums(countdata >= 10) >= 2, ]

dds <- DESeqDataSetFromMatrix(countData = countdata, colData = metadata,
                              design =~mutant*treatment)

dds <- DESeq(dds)
rlog <- rlog(dds,blind = TRUE)

results_to_fetch <- c("treatment_e2_vs_dmso",
                      "mutantsof1.treatmente2", 
                      "mutantsof2.treatmente2", 
                      "mutantsof3.treatmente2", 
                      "mutantdbm.treatmente2")

all_sig_genes <- data.frame()
for (result in results_to_fetch) {
  res <- results(dds, name = result)
  resdata <- as.data.frame(res)
  resdata$result <- result
  to_keep <- na.omit(resdata[resdata$padj <= 0.1,])
  all_sig_genes <- rbind(all_sig_genes, to_keep)
}
all_sig_genes$rownames <- rownames(all_sig_genes)

results_singlediff <- c("treatment_e2_vs_dmso",
                        "mutantsof1.treatmente2",
                        "mutantsof2.treatmente2",
                        "mutantsof3.treatmente2", 
                        "mutantdbm.treatmente2")

heatmap_results <- list()
for (result in results_singlediff) {
  if (result=="treatment_e2_vs_dmso") {
    res <- results(dds, name = result)
    resdata <- as.data.frame(res)
    to_keep <- resdata[rownames(resdata) %in% rownames(all_sig_genes),]
    to_keep <- to_keep %>%
      select(log2FoldChange)
    colnames(to_keep) <- result
    heatmap_results[[result]] <- to_keep
  } else{
    res <- results(dds, contrast = list(c("treatment_e2_vs_dmso",result)))
    resdata <- as.data.frame(res)
    to_keep <- resdata[rownames(resdata) %in% rownames(all_sig_genes),]
    to_keep <- to_keep %>%
      select(log2FoldChange)
    colnames(to_keep) <- result
    heatmap_results[[result]] <- to_keep
  }
}

heatmap_df <- as.matrix(as.data.frame(heatmap_results))
my_palette <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)

# Define color breaks centered around 0
max_val <- max(abs(heatmap_df), na.rm = TRUE)
breaks <- seq(-3, 3, length.out = 101)

# Perform k-means clustering manually
# Perform k-means clustering
set.seed(123)
# Compute and plot wss for k = 2 to k = 15.
k.max <- 15
data <- heatmap_df
wss <- sapply(1:k.max, 
              function(k){kmeans(data, k, nstart=50,iter.max = 15 )$tot.withinss})
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")


k <- 4
kmeans_result <- kmeans(heatmap_df, centers = k)

# Create annotation dataframe
annotation_row <- data.frame(Cluster = factor(kmeans_result$cluster))
rownames(annotation_row) <- rownames(heatmap_df)

# Reorder matrix and annotation by cluster
ordering <- order(annotation_row$Cluster)
heatmap_df_ordered <- heatmap_df[ordering, ]
annotation_row_ordered <- annotation_row[ordering, , drop = FALSE]

annotation_colors <- list(
  Cluster = setNames(RColorBrewer::brewer.pal(k, "Set3"), levels(annotation_row$Cluster))
)

# Can also define the desired cluster order. This is mostly aesthetic
desired_order <- c(4, 1, 3, 2)

# Reorder the rows based on cluster order
heatmap_df_reordered <- do.call(rbind, lapply(desired_order, function(cluster) {
  heatmap_df_ordered[annotation_row_ordered$Cluster == cluster, , drop = FALSE]
}))

annotation_row_reordered <- do.call(rbind, lapply(desired_order, function(cluster) {
  annotation_row_ordered[annotation_row_ordered$Cluster == cluster, , drop = FALSE]
}))

# Now plot with the reordered data
Fig4C <- pheatmap(
  heatmap_df_reordered,
  color = my_palette,
  breaks = breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  annotation_row = annotation_row_reordered,
  annotation_colors = annotation_colors
)

View(annotation_row)
heatmap_df_reordered
