# Residence Time Plotting

library(ggsignif)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(data.table)
library(rstatix)
library(patchwork)
rm(list=ls())

color_scale=c("darkgrey","skyblue","purple","darkblue","red")

residence_fits <- read.table("dat/model_fits_summary_full.csv",sep=",",header=TRUE) 

residence_fits <- residence_fits[residence_fits$Molecule.Name=="Estradiol",]
residence_fits <- residence_fits[residence_fits$ess_bulk > 1000,]


residence_fits$Cell.Line.ID <- factor(residence_fits$Cell.Line.ID,levels = (c("WT (U2OS)",
                                                                              "SOF1 (U2OS)",
                                                                              "SOF2 (U2OS)",
                                                                              "SOF3 (U2OS)",
                                                                              "DBM (U2OS)"
                                                                                          )))

pairwise_comparisons <- list(
  c("WT (U2OS)", "SOF1 (U2OS)"),
  c("WT (U2OS)", "SOF2 (U2OS)"),
  c("WT (U2OS)", "SOF3 (U2OS)"),
  c("WT (U2OS)", "DBM (U2OS)")
)

# Select parameter of interest
residence_fits_fbound <- residence_fits %>%
  filter(parameter == "weights[1]")

stat_df <- residence_fits_fbound %>%
  wilcox_test(
    mean ~ Cell.Line.ID,
    ref.group = "WT (U2OS)"
  ) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

stat_df

FigE <- ggplot(residence_fits_fbound, aes(x = Cell.Line.ID, y = mean, fill=Cell.Line.ID)) +
  scale_fill_manual(values=color_scale) +
  geom_boxplot(outlier.shape = NA) +
  ylab("f_slow") +
  xlab(NULL) +
  theme_classic()

# Select parameter of interest
residence_fits_inv <- residence_fits %>%
  filter(parameter == "koff_fast") %>%
  mutate(
    residence_time = 1 / mean
  )

wt_ref <- residence_fits_inv %>%
  filter(Cell.Line.ID == "WT (U2OS)") %>%
  summarise(wt_mean = mean(residence_time)) %>%
  pull(wt_mean)

residence_fits_inv <- residence_fits_inv %>%
  mutate(
    fold_change = residence_time / wt_ref
  )

stat_df <- residence_fits_inv %>%
  wilcox_test(
    fold_change ~ Cell.Line.ID,
    ref.group = "WT (U2OS)"
  ) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

stat_df

FigF <- ggplot(residence_fits_inv, aes(x = Cell.Line.ID, y = fold_change, fill=Cell.Line.ID)) +
  scale_fill_manual(values=color_scale) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  ylab("Transient Residence time fold-change (vs WT)") +
  xlab(NULL) +
  theme_classic()

# Select parameter of interest
residence_fits_inv <- residence_fits %>%
  filter(parameter == "koff_slow") %>%
  mutate(
    residence_time = 1 / mean
  )

wt_ref <- residence_fits_inv %>%
  filter(Cell.Line.ID == "WT (U2OS)") %>%
  summarise(wt_mean = mean(residence_time)) %>%
  pull(wt_mean)

residence_fits_inv <- residence_fits_inv %>%
  mutate(
    fold_change = residence_time / wt_ref
  )

stat_df <- residence_fits_inv %>%
  wilcox_test(
    fold_change ~ Cell.Line.ID,
    ref.group = "WT (U2OS)"
  ) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

FigG <- ggplot(residence_fits_inv, aes(x = Cell.Line.ID, y = fold_change, fill = Cell.Line.ID)) +
  scale_fill_manual(values=color_scale) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  ylab("Staable Residence time fold-change (vs WT)") +
  xlab(NULL) +
  theme_classic()


# Sample fbound and productive fraction together, multiply it for total productive fraction
residence_fits <- read.table("dat/model_fits_summary_full.csv",sep=",",header=TRUE) 

# Poorly fitted models can be selected out by filtering on effective sample size
residence_fits <- residence_fits[residence_fits$ess_bulk > 1000,]

fbound <- read.table("dat/all_filtered_fov_dfs.csv", sep=",",header=TRUE)

# Function to estimate P(Productive) with uncertainty for a single cell line
estimate_productive_probability <- function(cell_line, residence_data, bound_data, molecule="Estradiol", 
                                            concentration=100, n_samples = 1000) {
  
  
  prod_vals <- residence_data$mean[residence_data$Cell.Line.ID == cell_line & 
                                     residence_data$Molecule.Name==molecule & 
                                     residence_data$Molecule.Concentration==concentration & 
                                     residence_data$parameter=="weights[1]"]
  bound_vals <- bound_data$fraction_bound[bound_data$user_metadata.Cell.Line.ID == cell_line & 
                                            bound_data$Molecule.Name==molecule & 
                                            bound_data$Molecule.Concentration==concentration]
  
  if (length(prod_vals) < 1 || length(bound_vals) < 1) return(NULL)
  
  resampled_prod <- sample(prod_vals, size = n_samples, replace = TRUE)
  resampled_bound <- sample(bound_vals, size = n_samples, replace = TRUE)
  p_prod <- resampled_prod * resampled_bound
  
  data.frame(Cell.Line.ID = cell_line,
             P_productive = p_prod)
}

# Identify common cell lines across both datasets
shared_cell_lines <- intersect(unique(residence_fits$Cell.Line.ID),
                               unique(fbound$user_metadata.Cell.Line.ID))


# Estimate for all shared cell lines
results <- do.call(rbind, lapply(shared_cell_lines, function(cl) {
  estimate_productive_probability(cl, residence_fits, fbound)
}))

results$Cell.Line.ID <- factor(results$Cell.Line.ID, levels = c("WT (U2OS)",
                                                                "SOF1 (U2OS)",
                                                                "SOF2 (U2OS)",
                                                                "SOF3 (U2OS)",
                                                                "DBM (U2OS)"))


wt_id <- "WT (U2OS)"
boot_df <- do.call(rbind, lapply(shared_cell_lines, function(cl) {
  estimate_productive_probability(cl, residence_fits, fbound)
}))

# Perform pairwise comparisons to WT
wt_boot <- boot_df %>% filter(Cell.Line.ID == wt_id) %>% pull(P_productive)

comparison_results <- lapply(shared_cell_lines[shared_cell_lines != wt_id], function(mutant_id) {
  mutant_boot <- boot_df %>% filter(Cell.Line.ID == mutant_id) %>% pull(P_productive)
  diffs <- mutant_boot - wt_boot
  prob_lt <- mean(mutant_boot > wt_boot)
  data.frame(Mutant = mutant_id, P_value = prob_lt)
})

sig_df <- bind_rows(comparison_results)

# Add WT with empty significance to match
sig_df <- bind_rows(sig_df, data.frame(Mutant = wt_id, P_value = NA))

# Merge significance labels into boot_df for plotting
boot_df <- boot_df %>%
  left_join(sig_df, by = c("Cell.Line.ID" = "Mutant")) %>%
  mutate(Cell.Line.ID = factor(Cell.Line.ID, levels = c("WT (U2OS)",
                                                        "SOF1 (U2OS)",
                                                        "SOF2 (U2OS)",
                                                        "SOF3 (U2OS)",
                                                        "DBM (U2OS)")))

FigH <- ggplot(boot_df, aes(x = Cell.Line.ID, y = P_productive, fill = Cell.Line.ID)) +
  scale_fill_manual(values=color_scale) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic() +
  ylab("Estimated P(Productive)") +
  xlab("Cell Line") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

(FigE | FigF) /( FigG | FigH)
