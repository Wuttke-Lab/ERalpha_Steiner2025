library(enrichR)
library(dplyr)
library(fgsea)
library(ggplot2)
library(dplyr)
library(igraph)
library(ggraph)
library(tidygraph)

### GSEA ####

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

resdata <- as.data.frame(res)
resdata$refseq <- rownames(resdata)


resdata <- separate_wider_delim(resdata, cols = refseq, delim = ".", names = c("refid", "iso"))
resdata

refseq_to_cid <- read.table(file = "dat/refseq_to_common_id.txt",
                            sep="\t",header=FALSE)

resdata_id <- merge(resdata, refseq_to_cid,by.x="refid",by.y="V1")
resdata_id <- resdata_id[!(duplicated(resdata_id$V2)),]


resdata_id <- na.omit(resdata_id)


gsea_vec <- -log(resdata_id$pvalue + 1e-20) * sign(resdata_id$log2FoldChange)
names(gsea_vec) <- resdata_id$V2


library(fgsea)


GSEA = function(gene_list, GO_file, pval) {
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  myGO = fgsea::gmtPathways(GO_file)
  
  fgRes <- fgsea::fgsea(pathways = myGO,
                        stats = gene_list,
                        minSize=10, ## minimum gene set size
                        maxSize=400 ## maximum gene set size
                        ) %>% 
    as.data.frame() %>% 
    dplyr::filter(padj < !!pval) %>% 
    arrange(desc(NES))
  message(paste("Number of signficant gene sets =", nrow(fgRes)))
  
  message("Collapsing Pathways -----")
  concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
                                      pathways = myGO,
                                      stats = gene_list)
  fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  message(paste("Number of gene sets after collapsing =", nrow(fgRes)))
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 20),
                  tail(fgRes, n = 20 ))
  
  total_up = sum(fgRes$Enrichment == "Up-regulated")
  total_down = sum(fgRes$Enrichment == "Down-regulated")
  header = paste0("Top 10 (Total pathways: Up=", total_up,", Down=",    total_down, ")")
  
  colos = setNames(c("firebrick2", "dodgerblue2"),
                   c("Up-regulated", "Down-regulated"))
  
  g1 = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_point( aes(fill = Enrichment, size = size), shape=21) +
    scale_fill_manual(values = colos ) +
    scale_size_continuous(range = c(2,10)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=header) +
    theme_classic(base_size = 14)
  
  output = list("Results" = fgRes, "Plot" = g1)
  return(output)
}

GO_file = "/dat/h.all.v2024.1.Hs.symbols.gmt.txt"

res = GSEA(gsea_vec, GO_file, pval = 0.05)
res$Plot

pathway_to_plot = "HALLMARK_ESTROGEN_RESPONSE_EARLY" 

plot_single_enrichment = function(gene_list, GO_file, pathway_name) {
  library(fgsea)
  
  # Load and check gene list
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  
  # Load pathways
  myGO = fgsea::gmtPathways(GO_file)
  
  # Check if pathway exists
  if (!(pathway_name %in% names(myGO))) {
    stop(paste("Pathway", pathway_name, "not found in GO file"))
  }
  
  # Plot running enrichment score
  p = plotEnrichment(myGO[[pathway_name]], gene_list) +
    labs(title = pathway_name) +
    ylim(-1,1) +
    theme_classic(base_size = 16)
  
  return(p)
}

plot_single_enrichment(sort(gsea_vec,decreasing=TRUE), GO_file, pathway_to_plot)
