# Promoter Pause Ratio comparison

library(ggplot2)
library(ggpubr)
library(DESeq2)
library(data.table)
rm(list=ls())
color_scale=c("red","skyblue","purple","darkblue")


gene_body <- read.csv("dat/genebody_counts_250down_tes.antisense.txt",
                      sep="\t",header = TRUE,skip=1)

proms <- read.csv("dat/pause_counts_50up_250down.antisense.txt",
                  sep="\t",header=TRUE,skip=1)

unique(proms$Geneid == gene_body$Geneid)

metadata <- read.csv(file = "dat/metadata.csv",
                     sep = ",", header=TRUE)

metadata$clean_name <- paste0(metadata$mutant,"_",metadata$treatment)

rownames(gene_body) <- gene_body$Geneid
rownames(proms) <- proms$Geneid

length_prom <- as.numeric(proms$End) - as.numeric(proms$Start)

length_gb <- as.numeric(gene_body$End) - as.numeric(gene_body$Start)

proms <- proms %>% select(-c(Geneid,Chr,Start,End,Strand,Length))
gene_body <- gene_body %>% select(-c(Geneid,Chr,Start,End,Strand,Length))

for (col in colnames(proms)){
  proms[[col]] <- as.numeric(proms[[col]])/length_prom
}
for (col in colnames(gene_body)){
  gene_body[[col]] <- as.numeric(gene_body[[col]])/length_gb
}

pause_ratio <- log2(((proms)/(gene_body)) + 1)
colnames(pause_ratio) <- metadata$clean_name

pause_ratio$gene_name <- rownames(pause_ratio)

pause_ratio
pause_long <- pause_ratio %>% pivot_longer(cols = -gene_name,names_to = "clean_name",values_to = "pause_ratio")

pause_joined <- pause_long %>%
  left_join(metadata, by = "clean_name")

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
res <- lfcShrink(dds, coef = "treatment_e2_vs_dmso", type = "ashr")
resdata <- as.data.frame(res)


upgenes <- rownames(na.omit(resdata[resdata$padj < 0.1 & resdata$log2FoldChange > 0,]))
all_genes <- pause_joined$gene_name

# clustered genes from heatmap_genes_cluster

pause_test <- pause_joined %>%
  filter(gene_name %in% upgenes)

pause_summary <- pause_test %>%
  group_by(gene_name, mutant, treatment) %>%
  summarise(mean_pause = mean(pause_ratio))

pause_compare_to_wt <- pause_test %>%
  group_by(gene_name, treatment) %>%
  # join WT values to all mutants within each gene/treatment
  left_join(
    pause_test %>% 
      filter(mutant == "wt") %>%
      select(gene_name, treatment, wt_pause_ratio = pause_ratio),
    by = c("gene_name", "treatment")
  ) %>%
  # compute comparison
  mutate(
    fold_change = pause_ratio / wt_pause_ratio,
    diff_from_wt = pause_ratio - wt_pause_ratio
  )

# Average across replicates if not done

pause_mean <- pause_compare_to_wt %>%
  group_by(gene_name, treatment, mutant) %>%
  summarize(mean_pause = mean(pause_ratio, na.rm = TRUE), .groups = "drop") %>%
  filter(treatment == "e2", gene_name %in% upgenes)

# Join WT pause values for each gene
pause_joined <- pause_mean %>%
  filter(treatment == "e2") %>%
  left_join(
    pause_mean %>%
      filter(treatment == "e2", mutant == "wt") %>%
      select(gene_name, wt_pause = mean_pause),
    by = "gene_name"
  ) %>%
  filter(mutant != "wt")  # only keep mutants for plotting

# Plot WT vs each mutant in separate facets

ggplot(pause_joined, aes(x = mutant, y = mean_pause)) +
  # Lines connecting WT to each mutant per gene
  geom_segment(
    aes(
      x = "wt", xend = mutant,
      y = wt_pause, yend = mean_pause,
      group = gene_name
    ),
    color = "gray70", alpha = 0.25
  ) +
  # WT points
  geom_point(
    data = pause_joined %>% select(gene_name, wt_pause) %>% distinct(),
    aes(x = "wt", y = wt_pause),
    color = "black", size = 1
  ) +
  geom_boxplot(
    data = pause_joined %>% select(gene_name, wt_pause) %>% distinct(),
    aes(x = "wt", y = wt_pause),
    color = "darkgrey",
    alpha = 0.2,
    position = position_nudge(x = 0.2),
    width=0.2,
    outlier.shape = NA
  ) +
  # Mutant points
  geom_point(aes(color = mutant), size = 1) +
  # Boxplots for mutant distributions
  geom_boxplot(
    aes(x = mutant, y = mean_pause,color=mutant),
    alpha = 0.2,
    width=0.2,
    position = position_nudge(x = -0.2),
    outlier.shape = NA
  ) +
  scale_color_manual(values=color_scale) +
  facet_wrap(~ mutant, scales = "free_x") +
  labs(
    x = "WT vs mutant",
    y = "Mean pause ratio"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )

# PERMUTATION TEST

n_genes <- length(upgenes)  
n_perm <- 5000             

# Function to compute mean mutant vs WT difference for a set of genes
mean_diff <- function(genes, mutant_of_interest) {
  mean(pause_joined$mean_pause[pause_joined$gene_name %in% genes & pause_joined$mutant == mutant_of_interest]) -
    mean(pause_joined$wt_pause[pause_joined$gene_name %in% genes])
}

# Compute observed difference for each mutant
mutants <- unique(pause_joined$mutant)
obs_diff <- sapply(mutants, function(m) mean_diff(upgenes, m))

# Generate null distribution via permutation
perm_diff <- replicate(n_perm, {
  sampled_genes <- sample(all_genes, n_genes)  # sample same number of genes
  sapply(mutants, function(m) mean_diff(sampled_genes, m))
})

# Compute permutation p-values
perm_p <- sapply(seq_along(mutants), function(i) {
  mean(abs(perm_diff[i, ]) >= abs(obs_diff[i]))
})

# Combine results into a table
perm_results <- data.frame(
  mutant = mutants,
  obs_diff = obs_diff,
  perm_p = perm_p
)
perm_results

