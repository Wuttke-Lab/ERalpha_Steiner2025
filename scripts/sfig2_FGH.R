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

color_scale=c("darkgrey","skyblue","purple","darkblue","red", "darkgrey")


residence_fits <- read.table("dat/model_fits_summary_full.csv",sep=",",header=TRUE) 

residence_fits <- residence_fits[residence_fits$Molecule.Name=="DMSO",]
# Poorly fitted models can be selected out based on low effective sample size
residence_fits <- residence_fits[residence_fits$ess_bulk > 1000,]


residence_fits$Cell.Line.ID <- factor(residence_fits$Cell.Line.ID,levels = (c("WT (U2OS)",
                                                                              "SOF1 (U2OS)",
                                                                              "SOF2 (U2OS)",
                                                                              "SOF3 (U2OS)",
                                                                              "DBM (U2OS)",
                                                                              "H2B")))

residence_fits
pairwise_comparisons <- list(
  c("WT (U2OS)", "SOF1 (U2OS)"),
  c("WT (U2OS)", "SOF2 (U2OS)"),
  c("WT (U2OS)", "SOF3 (U2OS)"),
  c("WT (U2OS)", "H2B"),
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

SFigF <- ggplot(residence_fits_fbound, aes(x = Cell.Line.ID, y = mean, fill=Cell.Line.ID)) +
  scale_fill_manual(values = color_scale) +
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

SFigG <- ggplot(residence_fits_inv, aes(x = Cell.Line.ID, y = fold_change, fill=Cell.Line.ID)) +
  scale_fill_manual(values = color_scale) +
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

stat_df

SFigH <- ggplot(residence_fits_inv, aes(x = Cell.Line.ID, y = fold_change, fill=Cell.Line.ID)) +
  scale_fill_manual(values = color_scale) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  ylab("Stable Residence time fold-change (vs WT)") +
  xlab(NULL) +
  theme_classic()

SFigF | SFigG | SFigH

