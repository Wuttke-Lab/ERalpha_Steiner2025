# Fbound by concentration line chart
library(ggplot2)
library(tidyverse)
library(dplyr)
library(purrr)
library(rootSolve)

kinetic_df <- read.table("/Users/samuelhunter/ERalpha_Steiner2025/dat/all_filtered_fov_dfs.csv",
                         sep=",", header=TRUE)

kinetic_df <- kinetic_df[kinetic_df$`Molecule.Name`=="Estradiol",]
colnames(kinetic_df)
epsilon <- 1e-10
kinetic_df$logConc <- log10(kinetic_df$Molecule.Concentration + epsilon)


summary_kinetic <- kinetic_df %>%
  group_by(`user_metadata.Cell.Line.ID`,logConc) %>%
  summarise_at(vars(fraction_bound),list(frac = mean))

# Define 4-parameter logistic function
fourPL <- function(x, Ymin, Ymax, logEC50, n) {
  Ymin + (Ymax - Ymin) / (1 + 10^((logEC50 - log10(x)) * n))
}

# Fit per group
View(kinetic_df)
fits <- kinetic_df %>%
  group_by(user_metadata.Cell.Line.ID) %>%
  group_modify(~ {
    tryCatch({
      nls_fit <- nls(
        fraction_bound ~ fourPL(Molecule.Concentration, Ymin, Ymax, logEC50, n),
        data = .x,
        start = list(
          Ymin = min(.x$fraction_bound, na.rm = TRUE),
          Ymax = max(.x$fraction_bound, na.rm = TRUE),
          logEC50 = 0,
          n = 1
        )
      )
      tibble(
        logEC50 = coef(nls_fit)["logEC50"],
        n = coef(nls_fit)["n"],
        Ymin = coef(nls_fit)["Ymin"],
        Ymax = coef(nls_fit)["Ymax"]
      )
    }, error = function(e) tibble())
  })

fits

# Create smooth curves per group
df_pred <- kinetic_df %>%
  left_join(fits, by = "user_metadata.Cell.Line.ID") %>%
  group_by(user_metadata.Cell.Line.ID) %>%
  do({
    grid <- tibble(Molecule.Concentration = 10^seq(-3, 2, length.out = 200))
    params <- distinct(., Ymin, Ymax, logEC50, n)
    grid$fraction_pred <- with(params,
                               fourPL(grid$Molecule.Concentration, Ymin, Ymax, logEC50, n)
    )
    grid
  })

ggplot(kinetic_df, aes(Molecule.Concentration, fraction_bound, color = user_metadata.Cell.Line.ID)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_line(data = df_pred, aes(y = fraction_pred), linewidth=1) +
  scale_x_log10(breaks=c(0.001, 0.01, 0.1, 1, 10, 100)) +
  theme_classic()

# Function to compute EC90 given model parameters
compute_EC90 <- function(Ymin, Ymax, logEC50, n) {
  # target fraction = 90% of response range
  target <- Ymin + 0.9 * (Ymax - Ymin)
  
  # equation to solve for x
  f <- function(x) {
    Ymin + (Ymax - Ymin) / (1 + 10^((logEC50 - log10(x)) * n)) - target
  }
  
  # search in a broad range (1e-4 to 1e4)
  out <- tryCatch(
    uniroot(f, interval = c(1e-4, 1e4))$root,
    error = function(e) NA
  )
  
  return(out)
}

# Compute EC50 and EC90 per mutant
ec_table <- fits %>%
  mutate(
    EC50 = 10^(logEC50),
    EC90 = pmap_dbl(
      list(Ymin, Ymax, logEC50, n),
      compute_EC90
    )
  )

print(ec_table)
