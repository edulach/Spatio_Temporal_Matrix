#### Set the path of files ####
base_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
results_dir <- file.path(base_dir, "results_test")
setwd(base_dir)


#1. Required packages 
library(dplyr)
library(ggplot2)
library(readxl)

#2. Frobenius Beta
data_beta <- read.csv(file.path(results_dir, "Beta_combined_100sims.csv"))

frobenius_beta <- data_beta %>%
  group_by(Simulation) %>%
  summarise(
    frobenius_error = sqrt(sum((Estimate - True)^2)),
    .groups = "drop"
  ) %>%
  mutate(Model = "L")

frobenius_beta

#3. Frobenius Sigma

# Determine matrix size from true Sigma
p <- nrow(Sigma)
paste0("Sigma.", rep(1:p, each = p), ".", rep(1:p, times = p), ".")

data_Sigma <- read.csv(file.path(results_dir, "Sigma_combined_100sims.csv"))

sigma_cols <- as.vector(outer(1:p, 1:p, function(i, j) paste0("Sigma.", i, ".", j, ".")))

data_Sigma %>%
  select(-Simulation) %>%
  summarise(across(everything(), mean))
# Compute Frobenius error for each simulation
frobenius_Sigma <- sapply(1:nrow(data_Sigma), function(i) {
  
  estimated_matrix <- matrix(
    as.numeric(data_Sigma[i, sigma_cols]),
    nrow = p, ncol = p, byrow = TRUE
  )
  
  norm(estimated_matrix - Sigma, "F")
})

frobenius_Sigma



##############################################################

