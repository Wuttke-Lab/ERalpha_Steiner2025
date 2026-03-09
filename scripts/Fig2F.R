# Linear Models of kinetic parameters and DNA/RNA affinity
library(tidyverse)
library(dplyr)
library(ggplot2)
library(tidyr)
library(broom)
library(patchwork)
rm(list=ls())
R.version

affinity_df <- read.csv("dat/kds_binding_affinity.csv",
                        sep=",",
                        header=TRUE)
affinity_df <- affinity_df %>%
  select(-c(kd,log_kd)) %>%
  distinct(mutant, mean, .keep_all = TRUE)
  
fbound_df <- read.table("dat/all_filtered_fov_dfs.csv",
                        sep=",", 
                        header=TRUE)

fits <- read.csv("dat/fbound_fits_byconc.csv",header=TRUE,sep=",")

fbound_df <- fbound_df %>%
  filter(Molecule.Concentration == 100) %>%
  filter(Molecule.Name == "Estradiol") %>%
  select(c(user_metadata.Cell.Line.ID,Molecule.Name,fraction_bound,Molecule.Concentration))

# Merge with fits to zero out the baseline
fbound_df <- merge(fbound_df, fits,by="user_metadata.Cell.Line.ID")
fbound_df$d_fbound <- fbound_df$fraction_bound - fbound_df$Ymin

fbound_df <- merge(fbound_df,affinity_df,by.x="user_metadata.Cell.Line.ID",by.y="mutant")
mutants <- unique(affinity_df$mutant)

#View(fbound_df)

reg_stats <- bind_rows(
  # DNA regression
  lm(d_fbound ~ mean, data = fbound_df[fbound_df$nucleic_acid=="DNA",]) %>%
    tidy() %>%
    filter(term == "mean") %>%
    mutate(
      predictor = "DNA",
      r2 = summary(lm(d_fbound ~ mean, data = fbound_df[fbound_df$nucleic_acid=="DNA",]))$r.squared
    ),
  
  # RNA regression
  lm(d_fbound ~ mean, data = fbound_df[fbound_df$nucleic_acid=="RNA",]) %>%
    tidy() %>%
    filter(term == "mean") %>%
    mutate(
      predictor = "RNA",
      r2 = summary(lm(d_fbound ~ mean, data = fbound_df[fbound_df$nucleic_acid=="RNA",]))$r.squared
    )
) %>%
  transmute(
    predictor,
    label = sprintf(
      "R² = %.2f\np = %.2g",
      r2, p.value
    )
  )
fbound_df

summary_df <- fbound_df %>%
  group_by(nucleic_acid, user_metadata.Cell.Line.ID) %>%
  summarise(mutant_mean = mean(d_fbound),n = n(),
            mutant_sd = sd(d_fbound),
            mutant_se = mutant_sd/sqrt(n),
            kd = mean(mean))

Fig2F_left <- ggplot(data=fbound_df[fbound_df$nucleic_acid=="DNA",],aes(x=mean,y=d_fbound)) +
  geom_errorbar(data=summary_df[summary_df$nucleic_acid=="DNA",],aes(ymin=mutant_mean-mutant_sd, y=mutant_mean, x=kd,
                                                                     ymax=mutant_mean+mutant_sd), width=0.1,
                position=position_dodge(0.05)) +
  geom_point(data=summary_df[summary_df$nucleic_acid=="DNA",],
             aes(x=kd,y=mutant_mean,color=user_metadata.Cell.Line.ID),size=1) +
  geom_smooth(method="lm") +
  geom_text(
    data = reg_stats[reg_stats$predictor=="DNA",],
    aes(
      x = -Inf, y = 0.4,
      label = label
    ),
    hjust = -0.05, vjust = -0.1,
    size = 4,
    inherit.aes = FALSE
  ) +
  theme_classic(base_size=18) +
  theme(legend.position = "none")


Fig2F_right <- ggplot(data=fbound_df[fbound_df$nucleic_acid=="RNA",],aes(x=mean,y=d_fbound)) +
  geom_errorbar(data=summary_df[summary_df$nucleic_acid=="RNA",],aes(ymin=mutant_mean-mutant_sd, y=mutant_mean, x=kd,
                                                                     ymax=mutant_mean+mutant_sd), width=0.1,
                position=position_dodge(0.05)) +
  geom_point(data=summary_df[summary_df$nucleic_acid=="RNA",],
             aes(x=kd,y=mutant_mean,color=user_metadata.Cell.Line.ID),size=1) +
  geom_smooth(method="lm") +
  geom_text(
    data = reg_stats[reg_stats$predictor=="RNA",],
    aes(
      x = -Inf, y = 0.4,
      label = label
    ),
    hjust = -0.05, vjust = -0.1,
    size = 4,
    inherit.aes = FALSE
  ) +
  theme_classic(base_size=18)

Fig2F_left | Fig2F_right

