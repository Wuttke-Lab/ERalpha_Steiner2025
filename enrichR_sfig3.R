library(enrichR)
library(dplyr)
library(fgsea)
library(ggplot2)
library(dplyr)
library(igraph)
library(ggraph)
library(tidygraph)

### ENRICHR ####
setEnrichrSite("Enrichr")
websiteLive <- TRUE
dbs <- listEnrichrDbs()

if (is.null(dbs)) websiteLive <- FALSE

if (websiteLive) head(dbs)

if (websiteLive) {
  enriched <- enrichr(na.omit(cluster2_genes), dbs$libraryName)
}

go_molecular <- enriched$GO_Molecular_Function_2025
go_molecular <- go_molecular[go_molecular$Adjusted.P.value<0.05,]
go_molecular <- go_molecular[order((-go_molecular$P.value)),]
go_molecular$P.value <- as.numeric(as.character(go_molecular$P.value))

go_molecular$Term <- as.factor(go_molecular$Term)
go_molecular$Term <- factor(go_molecular$Term, levels=as.character(go_molecular$Term))

ggplot(data=go_molecular,aes(x= -log10(P.value),y=Term)) +
  geom_bar(stat="identity",color="black",fill="skyblue") + 
  theme_classic()


### GSEA ####

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
