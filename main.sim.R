#### Set the path of files ####
base_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(base_dir)
source(file.path(base_dir, "codes", "Functions SpatioFile.R"))
#1. Required packages 
library(matrixNormal)
library(MixMatrix)
library(mvtnorm)
library(MomTrunc)
library(mnormt)
library(Matrix)

#2. Create output directory
output_dir <- "results_test" 
if (!dir.exists(output_dir)) dir.create(output_dir)

#3. Fix parameters and matrices

p <- 3      # number of variables
L <- 5      # number of spatial locations
T <- 12    # temporal observations per location
r <- L * T
s <- 3      # number of covariates per observation

Beta<- matrix(data = c(1,1.4,2,1,1.2,1,2,1,1.2), p, s, byrow = T)
phi <- 1.2
rho <- 0.7
sigma2 <- 1.1
corr <- "exponential"
Sigma <- matrix(c(
  1, 0.8, 0.64,
  0.8, 1, 0.8,
  0.64, 0.8, 1
), nrow = 3, byrow = TRUE)

#4. Storage results
all_beta <- list()
all_sigma <- list()
all_scalars <- list()
all_logliks <- list()
all_BIC <- numeric(100)

set.seed(123)  # For reproducibility

#5. Simulation
nsim<-10

for (sim in 1:nsim) {
  cat("Running simulation", sim, "...\n")
  
  # Simulate covariates
  X <- matrix(rnorm(s * r, 5, 1), s, r)
  M <- Beta %*% X
  
  # Random coordinates per sim
  r1 <- sample(seq(1, 10, length = 300), L)
  r2 <- sample(seq(1, 10, length = 300), L)
  coords <- cbind(r1, r2)
  
  # Simulate data
  Matrixdata <- GeraMatrixFull(M = M, coords = coords, phi = phi, rho = rho,
                               sigma2 = sigma2, Sigma = Sigma, corr = corr)
  dados <- Matrixdata$Y
  
  # Fit model
  Res_var <- ML.MatrixRegreSPatioT(dados, X, coords,
                                   precision = 1e-9, MaxIter = 50, corr = corr)
  
  # Save in memory
  all_beta[[sim]] <- as.vector(Res_var$Beta)
  all_sigma[[sim]] <- as.vector(Res_var$Sigma)
  all_scalars[[sim]] <- c(phi_est = Res_var$phi, rho_est = Res_var$rho, sigma2_est = Res_var$sigma2)
  all_logliks[[sim]] <- Res_var$loglik
  all_BIC[sim] <- Res_var$BIC
}

# Convert results to data frames
beta_df <- do.call(rbind, lapply(seq_along(all_beta), function(i) {
  data.frame(Simulation = i,
             Parameter = paste0("Beta[", seq_along(all_beta[[i]]), "]"),
             Estimate = all_beta[[i]],
             True = as.vector(Beta))
}))

sigma_df <- do.call(rbind, lapply(seq_along(all_sigma), function(i) {
  data.frame(Simulation = i, matrix(all_sigma[[i]], nrow = 1))
}))
colnames(sigma_df)[-1] <- paste0("Sigma[", rep(1:p, each = p), ",", rep(1:p, times = p), "]")

scalars_df <- do.call(rbind, lapply(seq_along(all_scalars), function(i) {
  data.frame(Simulation = i,
             Parameter = c("phi", "rho", "sigma2"),
             Estimate = all_scalars[[i]],
             True = c(phi, rho, sigma2))
}))

loglik_df <- do.call(rbind, lapply(seq_along(all_logliks), function(i) {
  data.frame(Simulation = i,
             Iteration = seq_along(all_logliks[[i]]),
             LogLik = all_logliks[[i]])
}))

BIC_df <- data.frame(Simulation = 1:100, BIC = all_BIC)

# Export results
write.csv(beta_df, file = file.path(output_dir, "Beta_combined_100sims.csv"), row.names = FALSE)
write.csv(sigma_df, file = file.path(output_dir, "Sigma_combined_100sims.csv"), row.names = FALSE)
write.csv(scalars_df, file = file.path(output_dir, "Scalars_combined_100sims.csv"), row.names = FALSE)
write.csv(loglik_df, file = file.path(output_dir, "Loglik_combined_100sims.csv"), row.names = FALSE)
write.csv(BIC_df, file = file.path(output_dir, "BIC_combined_100sims.csv"), row.names = FALSE)

cat("✅ Finished all 100 simulations. Results saved in:\n", normalizePath(output_dir), "\n")
