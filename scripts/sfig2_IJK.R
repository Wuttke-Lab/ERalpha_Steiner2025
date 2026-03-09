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
  
residence_imaging_df <- read.table("dat/koff_fast_as_fc_to_wt.csv",
                                   sep=",",
                                   header=TRUE) 
residence_imaging_df <- merge(residence_imaging_df, affinity_df,by.x="Cell.Line.ID",by.y="mutant")


reg_stats <- bind_rows(
  # DNA regression
  lm(fold_change ~ mean.y, data = residence_imaging_df[residence_imaging_df$nucleic_acid=="DNA",]) %>%
    tidy() %>%
    filter(term == "mean.y") %>%
    mutate(
      predictor = "DNA",
      r2 = summary(lm(fold_change ~ mean.y, data = residence_imaging_df[residence_imaging_df$nucleic_acid=="DNA",]))$r.squared
    ),
  
  # RNA regression
  lm(fold_change ~ mean.y, data = residence_imaging_df[residence_imaging_df$nucleic_acid=="RNA",]) %>%
    tidy() %>%
    filter(term == "mean.y") %>%
    mutate(
      predictor = "RNA",
      r2 = summary(lm(fold_change ~ mean.y, data = residence_imaging_df[residence_imaging_df$nucleic_acid=="RNA",]))$r.squared
    )
) %>%
  transmute(
    predictor,
    label = sprintf(
      "R² = %.2f\np = %.2g",
      r2, p.value
    )
  )

summary_df <- residence_imaging_df %>%
  group_by(nucleic_acid, Cell.Line.ID) %>%
  summarise(mutant_mean = mean(fold_change),
            kd = mean(mean.y), n = n(),
            mutant_sd = sd(fold_change),
            mutant_se = mutant_sd/sqrt(n))


FigI_left <- ggplot(data=residence_imaging_df[residence_imaging_df$nucleic_acid=="DNA",],aes(x=mean.y,y=fold_change)) +
  geom_point(data=summary_df[summary_df$nucleic_acid=="DNA",],aes(x=kd,y=mutant_mean,color=Cell.Line.ID),size=3) +
  geom_errorbar(data=summary_df[summary_df$nucleic_acid=="DNA",],aes(ymin=mutant_mean-mutant_se, y=mutant_mean, x=kd,
                                                                     ymax=mutant_mean+mutant_se), width=0.1,
                position=position_dodge(0.05)) +
  geom_smooth(method="lm") +
  geom_text(
    data = reg_stats[reg_stats$predictor=="DNA",],
    aes(
      x = -Inf, y = 0.89,
      label = label
    ),
    hjust = -0.05, vjust = -0.1,
    size = 4,
    inherit.aes = FALSE
  ) +
  theme_classic(base_size=18) +
  theme(legend.position = "none")

FigI_right <- ggplot(data=residence_imaging_df[residence_imaging_df$nucleic_acid=="RNA",],aes(x=mean.y,y=fold_change)) +
  geom_point(data=summary_df[summary_df$nucleic_acid=="RNA",],aes(x=kd,y=mutant_mean,color=Cell.Line.ID),size=3) +
  geom_errorbar(data=summary_df[summary_df$nucleic_acid=="RNA",],aes(ymin=mutant_mean-mutant_se, y=mutant_mean, x=kd,
                                                                     ymax=mutant_mean+mutant_se), width=0.1,
                position=position_dodge(0.05)) +
  geom_smooth(method="lm") +
  geom_text(
    data = reg_stats[reg_stats$predictor=="RNA",],
    aes(
      x = -Inf, y = 0.89,
      label = label
    ),
    hjust = -0.05, vjust = -0.1,
    size = 4,
    inherit.aes = FALSE
  ) +
  theme_classic(base_size=18)

FigI_left | FigI_right

### koff_slow
rm(list=ls())
affinity_df <- read.csv("dat/kds_binding_affinity.csv",
                        sep=",",
                        header=TRUE)
affinity_df <- affinity_df %>%
  select(-c(kd,log_kd)) %>%
  distinct(mutant, mean, .keep_all = TRUE)

residence_imaging_df <- read.table("dat/koff_slow_as_fc_to_wt.csv",
                                   sep=",",
                                   header=TRUE) 

residence_imaging_df <- merge(residence_imaging_df, affinity_df,by.x="Cell.Line.ID",by.y="mutant")


reg_stats <- bind_rows(
  # DNA regression
  lm(fold_change ~ mean.y, data = residence_imaging_df[residence_imaging_df$nucleic_acid=="DNA",]) %>%
    tidy() %>%
    filter(term == "mean.y") %>%
    mutate(
      predictor = "DNA",
      r2 = summary(lm(fold_change ~ mean.y, data = residence_imaging_df[residence_imaging_df$nucleic_acid=="DNA",]))$r.squared
    ),
  
  # RNA regression
  lm(fold_change ~ mean.y, data = residence_imaging_df[residence_imaging_df$nucleic_acid=="RNA",]) %>%
    tidy() %>%
    filter(term == "mean.y") %>%
    mutate(
      predictor = "RNA",
      r2 = summary(lm(fold_change ~ mean.y, data = residence_imaging_df[residence_imaging_df$nucleic_acid=="RNA",]))$r.squared
    )
) %>%
  transmute(
    predictor,
    label = sprintf(
      "R² = %.2f\np = %.2g",
      r2, p.value
    )
  )

summary_df <- residence_imaging_df %>%
  group_by(nucleic_acid, Cell.Line.ID) %>%
  summarise(mutant_mean = mean(fold_change),
            kd = mean(mean.y), n = n(),
            mutant_sd = sd(fold_change),
            mutant_se = mutant_sd/sqrt(n))

FigJ_left <- ggplot(data=residence_imaging_df[residence_imaging_df$nucleic_acid=="DNA",],aes(x=mean.y,y=fold_change)) +
  geom_point(data=summary_df[summary_df$nucleic_acid=="DNA",],aes(x=kd,y=mutant_mean,color=Cell.Line.ID),size=3) +
  geom_errorbar(data=summary_df[summary_df$nucleic_acid=="DNA",],aes(ymin=mutant_mean-mutant_se, y=mutant_mean, x=kd,
                                                                     ymax=mutant_mean+mutant_se), width=0.1,
                position=position_dodge(0.05)) +
  geom_smooth(method="lm") +
  geom_text(
    data = reg_stats[reg_stats$predictor=="DNA",],
    aes(
      x = -Inf, y = 0.89,
      label = label
    ),
    hjust = -0.05, vjust = -0.1,
    size = 4,
    inherit.aes = FALSE
  ) +
  theme_classic(base_size=18) +
  theme(legend.position = "none")

FigJ_right <- ggplot(data=residence_imaging_df[residence_imaging_df$nucleic_acid=="RNA",],aes(x=mean.y,y=fold_change)) +
  geom_point(data=summary_df[summary_df$nucleic_acid=="RNA",],aes(x=kd,y=mutant_mean,color=Cell.Line.ID),size=3) +
  geom_errorbar(data=summary_df[summary_df$nucleic_acid=="RNA",],aes(ymin=mutant_mean-mutant_se, y=mutant_mean, x=kd,
                                                                     ymax=mutant_mean+mutant_se), width=0.1,
                position=position_dodge(0.05)) +
  geom_smooth(method="lm") +
  geom_text(
    data = reg_stats[reg_stats$predictor=="RNA",],
    aes(
      x = -Inf, y = 0.89,
      label = label
    ),
    hjust = -0.05, vjust = -0.1,
    size = 4,
    inherit.aes = FALSE
  ) +
  theme_classic(base_size=18)
FigJ_left | FigJ_right

### fslow
rm(list=ls())
affinity_df <- read.csv("dat/kds_binding_affinity.csv",
                        sep=",",
                        header=TRUE)

affinity_df <- affinity_df %>%
  select(-c(kd,log_kd)) %>%
  distinct(mutant, mean, .keep_all = TRUE)

residence_imaging_df <- read.table("dat/model_fits_summary_full.csv",
                                   sep=",",
                                   header=TRUE) 

residence_imaging_df <- merge(residence_imaging_df, affinity_df,by.x="Cell.Line.ID",by.y="mutant")

residence_imaging_df <- residence_imaging_df %>%
  filter(Molecule.Name=="Estradiol") %>%
  filter(parameter=="weights[1]")

reg_stats <- bind_rows(
  # DNA regression
  lm(mean.x ~ mean.y, data = residence_imaging_df[residence_imaging_df$nucleic_acid=="DNA",]) %>%
    tidy() %>%
    filter(term == "mean.y") %>%
    mutate(
      predictor = "DNA",
      r2 = summary(lm(mean.x ~ mean.y, data = residence_imaging_df[residence_imaging_df$nucleic_acid=="DNA",]))$r.squared
    ),
  
  # RNA regression
  lm(mean.x ~ mean.y, data = residence_imaging_df[residence_imaging_df$nucleic_acid=="RNA",]) %>%
    tidy() %>%
    filter(term == "mean.y") %>%
    mutate(
      predictor = "RNA",
      r2 = summary(lm(mean.x ~ mean.y, data = residence_imaging_df[residence_imaging_df$nucleic_acid=="RNA",]))$r.squared
    )
) %>%
  transmute(
    predictor,
    label = sprintf(
      "R² = %.2f\np = %.2g",
      r2, p.value
    )
  )

summary_df <- residence_imaging_df %>%
  group_by(nucleic_acid, Cell.Line.ID) %>%
  summarise(mutant_mean = mean(mean.x),
            kd = mean(mean.y), n = n(),
            mutant_sd = sd(mean.x),
            mutant_se = mutant_sd/sqrt(n))

FigK_left <- ggplot(data=residence_imaging_df[residence_imaging_df$nucleic_acid=="DNA",],aes(x=mean.y,y=mean.x)) +
  geom_point(data=summary_df[summary_df$nucleic_acid=="DNA",],aes(x=kd,y=mutant_mean,color=Cell.Line.ID),size=3) +
  geom_errorbar(data=summary_df[summary_df$nucleic_acid=="DNA",],aes(ymin=mutant_mean-mutant_se, y=mutant_mean, x=kd,
                                                                     ymax=mutant_mean+mutant_se), width=0.1,
                position=position_dodge(0.05)) +
  geom_smooth(method="lm") +
  geom_text(
    data = reg_stats[reg_stats$predictor=="DNA",],
    aes(
      x = -Inf, y = 0.14,
      label = label
    ),
    hjust = -0.05, vjust = -0.1,
    size = 4,
    inherit.aes = FALSE
  ) +
  theme_classic(base_size=18) +
  theme(legend.position = "none")

FigK_right <- ggplot(data=residence_imaging_df[residence_imaging_df$nucleic_acid=="RNA",],aes(x=mean.y,y=mean.x)) +
  geom_point(data=summary_df[summary_df$nucleic_acid=="RNA",],aes(x=kd,y=mutant_mean,color=Cell.Line.ID),size=3) +
  geom_errorbar(data=summary_df[summary_df$nucleic_acid=="RNA",],aes(ymin=mutant_mean-mutant_se, y=mutant_mean, x=kd,
                                                                     ymax=mutant_mean+mutant_se), width=0.1,
                position=position_dodge(0.05)) +
  geom_smooth(method="lm") +
  geom_text(
    data = reg_stats[reg_stats$predictor=="RNA",],
    aes(
      x = -Inf, y = 0.14,
      label = label
    ),
    hjust = -0.05, vjust = -0.1,
    size = 4,
    inherit.aes = FALSE
  ) +
  theme_classic(base_size=18)
FigK_left | FigK_right
