# Residence Time Plotting

library(ggsignif)
library(ggplot2)
library(dplyr)
library(ggpubr)

residence_fits_eikon <- read.table("dat/model_fits_summary_full.csv",sep=",",header=TRUE) 


residence_fits_eikon <- residence_fits_eikon[residence_fits_eikon$Molecule.Name=="Estradiol",]
residence_fits_eikon <- residence_fits_eikon[residence_fits_eikon$Cell.Line.ID %like% "U2OS",]

# Fits that failed will have a low effective sample size. Filter these.
residence_fits_eikon <- residence_fits_eikon[residence_fits_eikon$ess_bulk > 1000,]


residence_fits_eikon$Cell.Line.ID <- factor(residence_fits_eikon$Cell.Line.ID,levels = (c("WT (U2OS)",
                                                                                          "SOF1 (U2OS)",
                                                                                          "SOF2 (U2OS)",
                                                                                          "SOF3 (U2OS)",
                                                                                          "DBM (U2OS)")))
pairwise_comparisons <- list(
  c("WT (U2OS)", "SOF1 (U2OS)"),
  c("WT (U2OS)", "SOF2 (U2OS)"),
  c("WT (U2OS)", "SOF3 (U2OS)"),
  c("WT (U2OS)", "DBM (U2OS)")
)

residence_fits_eikon

ggplot(data=residence_fits_eikon[residence_fits_eikon$parameter=="weights[0]",],aes(x=Cell.Line.ID,y=mean)) +
  geom_boxplot(fill="skyblue") +
  #geom_jitter() +
  ylab("fast_fraction") +
  geom_signif(comparisons = pairwise_comparisons, 
              map_signif_level=TRUE,y_position = c(0.89,0.9,0.91,0.92)) +
  theme_classic()


fbound_eikon <- read.table("dat/all_filtered_fov_dfs.csv", sep=",",header=TRUE)


# --- Parameters ---
n_samples <- 10000
wt_id <- "WT (U2OS)"

# --- Function to bootstrap P(Productive) ---
estimate_productive_probability <- function(cell_line, residence_data, bound_data, molecule="Estradiol", 
                                            concentration=100, n_samples = 10000) {
  
  
  prod_vals <- residence_data$mean[residence_data$Cell.Line.ID == cell_line & residence_data$Molecule.Name==molecule &  residence_data$Molecule.Concentration==concentration & residence_data$parameter=="weights[1]"]
  bound_vals <- bound_data$fraction_bound[bound_data$user_metadata.Cell.Line.ID == cell_line & bound_data$Molecule.Name==molecule & bound_data$Molecule.Concentration==concentration]
  
  if (length(prod_vals) < 1 || length(bound_vals) < 1) return(NULL)
  
  resampled_prod <- sample(prod_vals, size = n_samples, replace = TRUE)
  resampled_bound <- sample(bound_vals, size = n_samples, replace = TRUE)
  p_prod <- resampled_prod * resampled_bound
  
  #data.frame(Cell.Line.ID = cell_line,
  #           P_productive_mean = mean(p_prod),
  #           P_productive_lower = quantile(p_prod, 0.025),
  #           P_productive_upper = quantile(p_prod, 0.975))
  
  data.frame(Cell.Line.ID = cell_line,
             P_productive = p_prod)
  
}

# --- Bootstrap all shared cell lines ---
shared_cell_lines <- intersect(unique(residence_fits_eikon$Cell.Line.ID),
                               unique(fbound_eikon$user_metadata.Cell.Line.ID))

# Bootstrap results as a long data frame

boot_df <- do.call(rbind, lapply(shared_cell_lines, function(cl) {
  estimate_productive_probability(cl, residence_fits_eikon, fbound_eikon)
}))
boot_df
# --- Perform pairwise comparisons to WT ---
wt_boot <- boot_df %>% filter(Cell.Line.ID == wt_id) %>% pull(P_productive)
diffs <- mutant_boot - wt_boot

comparison_results <- lapply(shared_cell_lines[shared_cell_lines != wt_id], function(mutant_id) {
  mutant_boot <- boot_df %>% filter(Cell.Line.ID == mutant_id) %>% pull(P_productive)
  diffs <- mutant_boot - wt_boot
  p_val <- mean(abs(diffs) >= abs(mean(diffs)))
  prob_lt <- mean(mutant_boot > wt_boot)
  
  # Determine significance asterisk
  sig_label <- case_when(
    p_val < 0.001 ~ "***",
    p_val < 0.01 ~ "**",
    p_val <= 0.05 ~ "*",
    p_val > 0.05 ~ "N.S.",
    TRUE ~ ""
  )
  
  data.frame(Mutant = mutant_id, P_value = p_val, Significance = sig_label)
})

sig_df <- bind_rows(comparison_results)

sig_df

# Add WT with empty significance to match
sig_df <- bind_rows(sig_df, data.frame(Mutant = wt_id, P_value = NA, Significance = ""))

# Merge significance labels into boot_df for plotting
boot_df <- boot_df %>%
  left_join(sig_df, by = c("Cell.Line.ID" = "Mutant")) %>%
  mutate(Cell.Line.ID = factor(Cell.Line.ID, levels = c("WT (U2OS)",
                                                                 "SOF1 (U2OS)",
                                                                 "SOF2 (U2OS)",
                                                                 "SOF3 (U2OS)",
                                                                 "DBM (U2OS)")))

# --- Plot ---
ggplot(boot_df, aes(x = Cell.Line.ID, y = P_productive, fill = Cell.Line.ID)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7,fill="steelblue") +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2, fill = "white") +
  theme_classic() +
  ylab("Estimated P(Productive)") +
  xlab("Cell Line") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")




wt_boot <- boot_df %>% filter(Cell.Line.ID == wt_id) %>% pull(P_productive)
mutant_boot <- boot_df %>% filter(Cell.Line.ID == "SOF2 (U2OS)") %>% pull(P_productive)

obs_diff <- mean(wt_boot) - mean(mutant_boot)
combined <- c(wt_boot, mutant_boot)
n_wt <- length(wt_boot)
set.seed(42)  # for reproducibility
n_perm <- 10000
perm_diffs <- numeric(n_perm)

for (i in 1:n_perm) {
  permuted <- sample(combined)
  perm_wt <- permuted[1:n_wt]
  perm_mut <- permuted[(n_wt + 1):length(permuted)]
  perm_diffs[i] <- mean(perm_wt) - mean(perm_mut)
}

p_value <- mean(abs(perm_diffs) >= abs(obs_diff))

hist(perm_diffs, breaks = 50, col = "lightgray",
     main = sprintf("Permutation Test\np = %.4f", p_value),
     xlab = "Difference in Means")
abline(v = obs_diff, col = "red", lwd = 2, lty = 2)
abline(v = -obs_diff, col = "red", lwd = 2, lty = 2)


