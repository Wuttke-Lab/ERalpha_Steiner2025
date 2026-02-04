
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

### Genes ####

countdata <- read.csv(file = "dat/counts.txt",
                      sep="\t",header=TRUE,skip=1)

rownames(countdata) <- countdata$Geneid
countdata <- countdata[,-c(1:6)]
metadata <- read.csv(file = "dat/metdata.csv",sep = ",", header=TRUE)
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
breaks <- seq(-2, 2, length.out = 101)

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

# Optional: define custom cluster colors
annotation_colors <- list(
  Cluster = setNames(RColorBrewer::brewer.pal(k, "Set3"), levels(annotation_row$Cluster))
)

# Plot heatmap with grouped clusters
pheatmap(
  heatmap_df_ordered,
  color = my_palette,
  breaks = breaks,
  cluster_rows = FALSE,                # we already sorted the rows
  cluster_cols = FALSE,
  show_rownames = FALSE,
  annotation_row = annotation_row_ordered,
  annotation_colors = annotation_colors
)

# Extract gene names for cluster 2
annotation_row <- annotation_row %>%
  rownames_to_column() %>%
  mutate(gene_id = sub("\\.\\d+$", "", rowname))

refseq.df <- read.csv("dat/refseq_to_common_id.txt",
                      header=FALSE,sep="\t") 

clustered_genes <- merge(annotation_row, refseq.df, by.x="gene_id",by.y="V1")
# Read MSigDB GMT file
gmt_lines <- readLines("/dat/h.all.v2024.1.Hs.symbols.gmt.txt")

# Parse into a named list: gene set name -> gene vector
gmt_list <- strsplit(gmt_lines, "\t") %>%
  set_names(map_chr(., 1)) %>%
  map(~ .x[-c(1:2)])  # remove name and description

e2_genes <- gmt_list[["HALLMARK_ESTROGEN_RESPONSE_LATE"]]  # or whichever is relevant
# Universe = all genes assigned to any cluster

# All genes assigned to any cluster
all_genes <- unique(clustered_genes$V2)

# Function to test enrichment
check_enrichment <- function(cluster_id, gene_df, e2_genes, universe) {
  cluster_genes <- gene_df %>% filter(Cluster == cluster_id) %>% pull(V2)
  
  # Contingency table
  in_set <- cluster_genes %in% e2_genes
  not_in_set <- !cluster_genes %in% e2_genes
  
  in_not_cluster <- setdiff(universe, cluster_genes)
  in_set_other <- in_not_cluster %in% e2_genes
  not_in_set_other <- !in_not_cluster %in% e2_genes
  
  mat <- matrix(
    c(sum(in_set), sum(not_in_set), sum(in_set_other), sum(not_in_set_other)),
    nrow = 2,
    byrow = TRUE
  )
  
  fisher <- fisher.test(mat)
  
  tibble(
    Cluster = as.character(cluster_id),
    OddsRatio = fisher$estimate,
    p_value = fisher$p.value,
    Overlap = sum(in_set),
    ClusterSize = length(cluster_genes)
  )
}

# Apply to all clusters (assumes clusters are labeled as 1, 2, 3, 4)
results_df <- bind_rows(
  lapply(1:4, function(k) check_enrichment(k, clustered_genes, e2_genes, all_genes))
)
ggplot(results_df, aes(x = Cluster, y = OddsRatio)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = paste0("p=", signif(p_value, 2))), 
            vjust = -0.5, size = 4) +
  labs(title = "Enrichment in Gene Set",
       y = "Fisher's Odds Ratio",
       x = "Cluster") +
  theme_classic(base_size=20)

