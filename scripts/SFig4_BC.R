library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)
library(DESeq2)
library(tidyr)
library(VennDiagram) 
library(data.table)
library(sva)
library(pheatmap)
library(RColorBrewer)

rm(list=ls())

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
  to_keep <- na.omit(resdata[resdata$padj < 0.1,])
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
SFig4A <- pheatmap(
  heatmap_df_reordered,
  color = my_palette,
  breaks = breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  annotation_row = annotation_row_reordered,
  annotation_colors = annotation_colors
)
SFig4A

# Extract gene names for the second displayed cluster
annotation_row <- annotation_row %>%
  rownames_to_column() %>%
  mutate(gene_id = sub("\\.\\d+$", "", rowname))

refseq.df <- read.csv("dat/refseq_to_common_id.txt",
                      header=FALSE,sep="\t") 

clustered_genes <- merge(annotation_row, refseq.df, by.x="gene_id",by.y="V1")



countdata <- read.csv(
  file = "dat/pro_raw_counts.txt",
  sep = "\t", header = TRUE, skip = 1
)

active_esr1_sites <- read.table(
  file = "dat/active_esr1_sites.bed",
  sep = "\t", header = FALSE
)

# ESR1 sites GRanges
gr_esr1 <- GRanges(
  seqnames = active_esr1_sites$V1,
  ranges = IRanges(start = active_esr1_sites$V2, end = active_esr1_sites$V3),
  strand = active_esr1_sites$V6,
  score = active_esr1_sites$V5,
  id = active_esr1_sites$V4
)

get_promoters_for_cluster <- function(cluster_num, clustered_genes, countdata) {
  
  cluster_genes <- clustered_genes[clustered_genes$Cluster == cluster_num, ]
  
  cluster_annot <- countdata[countdata$Geneid %in% cluster_genes$rowname, ] %>%
    select(Chr, Start, End, Geneid, Length, Strand) %>%
    mutate(
      Start = suppressWarnings(as.numeric(Start)),
      End   = suppressWarnings(as.numeric(End))
    ) %>%
    filter(
      !is.na(Chr),
      !is.na(Start),
      !is.na(End),
      Start > 0,
      End >= Start
    )
  
  if (nrow(cluster_annot) == 0) {
    stop(paste("No valid genomic coordinates for cluster", cluster_num))
  }
  
  gr_genes <- GRanges(
    seqnames = cluster_annot$Chr,
    ranges   = IRanges(start = cluster_annot$Start, end = cluster_annot$End),
    strand   = cluster_annot$Strand,
    gene_id  = cluster_annot$Geneid
  )
  
  promoters(gr_genes, upstream = 2000, downstream = 200)
}

find_esr1_overlaps <- function(gr_esr1, promoter_gr, cluster_num) {
  hits <- findOverlaps(gr_esr1, promoter_gr)
  df <- data.frame(
    ESR1_chr = seqnames(gr_esr1)[queryHits(hits)],
    ESR1_start = start(gr_esr1)[queryHits(hits)],
    ESR1_end = end(gr_esr1)[queryHits(hits)],
    ESR1_id = gr_esr1$id[queryHits(hits)],
    ESR1_score = gr_esr1$score[queryHits(hits)],
    Gene = promoter_gr$gene_id[subjectHits(hits)],
    Strand = strand(promoter_gr)[subjectHits(hits)],
    cluster = as.character(cluster_num)
  )
  return(df)
}

compute_gc_content <- function(promoter_gr, cluster_num) {
  seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, promoter_gr)
  gc_vals <- letterFrequency(seqs, letters = "GC", as.prob = TRUE) * 100
  data.frame(GC_content = gc_vals, cluster = as.character(cluster_num))
}

clusters <- 1:4
overlaps_all <- list()
gc_all <- list()

for (c in clusters) {
  cat("Processing cluster", c, "...\n")
  
  promoter_gr <- get_promoters_for_cluster(c, clustered_genes, countdata)
  
  overlaps_df <- find_esr1_overlaps(gr_esr1, promoter_gr, c)
  overlaps_all[[as.character(c)]] <- overlaps_df
  
  # Compute GC content
  gc_df <- compute_gc_content(promoter_gr, c)
  gc_all[[as.character(c)]] <- gc_df
}

# Combine all clusters
overlaps_tot <- do.call(rbind, overlaps_all)
gc_df_tot <- do.call(rbind, gc_all)

# Cluster order/labeling is arbitrary, but it's manually set
# here to keep the ordering consistent with other figures


desired_order <- c(4, 1, 3, 2)
overlaps_tot$cluster <- factor(overlaps_tot$cluster, levels=desired_order)
gc_df_tot$cluster <- factor(gc_df_tot$cluster, levels=desired_order)

SFig4B <- ggplot(overlaps_tot, aes(x = cluster, y = ESR1_score)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  theme_classic() +
  labs(x = "Cluster", y = "ESR1 motif score") +
  stat_compare_means(method = "anova", label = "p.format") +
  scale_x_discrete(labels = c("1", "2","3", "4")) +
  scale_fill_brewer(palette = "Set2")

SFig4C <- ggplot(gc_df_tot, aes(x = cluster, y = G.C)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.3) +
  theme_classic() +
  labs(x = "Cluster", y = "Promoter GC content") +
  scale_x_discrete(labels = c("1", "2","3", "4")) +
  stat_compare_means(method = "anova", label.x = 1, label.y = max(gc_df_tot$G.C))

SFig4B | SFig4C

