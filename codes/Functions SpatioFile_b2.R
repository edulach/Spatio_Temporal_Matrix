# library(matrixNormal) must be installed
# library(MixMatrix) must be installed
# library(MomTrunc) must be installed
# library(mnormt) must be installed
# library(Matrix) must be installed


################################################################################
### Covariances
################################################################################

# Temporal AR(1)
ar1_cor <- function(n, rho){
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) -
                    (1:n - 1))
  rho^exponent
}

ar1_inv <- function(n, rho) {
  mat <- matrix(0, n, n)
  diag(mat) <- 1 + rho^2
  mat[1, 1] <- 1
  mat[n, n] <- 1
  for (i in 1:(n-1)) {
    mat[i, i+1] <- -rho
    mat[i+1, i] <- -rho
  }
  mat <- mat / (1 - rho^2)
  return(mat)
}


# Spatio

spatial_corr_matrix <- function(coords, corr_fun = "exponential", phi, nu = NULL) {
  # Ensure coordinates are in matrix form.
  coords <- as.matrix(coords)
  # Compute the Euclidean distance matrix.
  D <- as.matrix(dist(coords))
  
  if (tolower(corr_fun) == "exponential") {
    # Exponential: rho(d) = exp(-d/phi)
    Corr <- exp(-D / phi)
    
  } else if (tolower(corr_fun) == "gaussian") {
    # Gaussian: rho(d) = exp(-(d/phi)^2)
    Corr <- exp(-(D / phi)^2)
  } else if (tolower(corr_fun) == "cubic") {
    
    r <- D / phi
    Corr <- (1 - 7 * r^2 + 8.75 * r^3 - 3.5 * r^5 + 0.75 * r^7)
    Corr[r >= 1] <- 0
  } else if (tolower(corr_fun) == "spherical") {
    r <- D / phi
    Corr <- (1 - 1.5 * r + 0.5 * r^3)
    Corr[r >= 1] <- 0
    
  } else if (tolower(corr_fun) == "matern") {
    # Mat?rn: requires a smoothness parameter nu.
    if (is.null(nu)) {
      stop("For matern correlation, please provide the nu parameter (smoothness).")
    }
    # Mat?rn correlation:
    #   rho(d) = (2^(1-nu)/gamma(nu)) * (d/phi)^nu * besselK(d/phi, nu),
    # with the convention that rho(0) = 1.
    r <- D / phi
    Corr <- ifelse(r == 0, 1, (2^(1 - nu) / gamma(nu)) * (r^nu) * besselK(r, nu))
    
  }  else {
    stop("Unknown corr_fun. Please use 'exponential', 'gaussian', or 'matern'.")
  }
  
  return(Corr)
}


#### Likelihood complete data ####


loglikelSpatiorho<-function(rho, PsiS, Sigma, dados, mu, TT){
  PsiT<-ar1_cor(TT,rho)
  Psi<-kronecker(PsiS, PsiT)
  Psi <- (Psi + t(Psi)) / 2
  
  # Check positive definiteness
  eigenvals <- eigen(Psi, symmetric = TRUE, only.values = TRUE)$values
  if(min(eigenvals) <= 0) {
    return(-1e10)  # Return large negative instead of -Inf
  }
  
  like<-matrixNormal::dmatnorm(dados, mu, Sigma, Psi, tol = .Machine$double.eps^0.5, log = TRUE)
  # Handle potential -Inf
  if(!is.finite(like)) return(-1e10)
  return(like)
}



loglikelSpatiophi<-function(phi, PsiT, Sigma, sigma2, dados, coords, mu, corr="exponential",nu=NULL){
  PsiS<-spatial_corr_matrix(coords, corr_fun = corr, phi,nu)*sigma2
  Psi<-kronecker(PsiS, PsiT)
  Psi <- (Psi + t(Psi)) / 2
  
  # Add small nugget for numerical stability
  diag(Psi) <- diag(Psi) + 1e-6
  
  # Verificar definição positiva antes de calcular likelihood
  ev <- eigen(Psi, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= 0) return(-1e10)

  like <- matrixNormal::dmatnorm(dados, mu, Sigma, Psi,
                                  tol = .Machine$double.eps^0.5, log = TRUE)
  if (!is.finite(like)) return(-1e10)
  return(like)
}


################################################################################
### Samples Generator
################################################################################


GeraMatrixFull<-function(M=M, coords, phi, rho, sigma2, Sigma, corr="exponential"){
  # coords: dimension 2xL. Same coordinates for all "p" observed variables
  # phi: Spatial Parameter
  # sigma2: Spatial Parameter
  # rho: Temporal Parameter
  
  L<-nrow(coords)
  p<-nrow(M)
  r<-ncol(M) # r=TL
  T<-r/L
  
  PsiT <- ar1_cor (T, rho = rho)
  PsiS<- spatial_corr_matrix(coords, corr_fun = corr, phi,nu)*sigma2
  Psi<- kronecker(PsiS, PsiT)
  
  X.or <- LaplacesDemon::rmatrixnorm(M=M, U=Sigma, V=Psi)
  return(list(Y=X.or))
}



################################################################################
#### Maximum Likelihood for complete data : Spatio-temporal ####################
#### Beta full matrix and M=Beta*X
################################################################################

ML.MatrixRegreSPatioT<-function(dados, X, coords, precision=0.0000001, MaxIter=50, corr="exponential",nu=NULL,  phi_range = c(0.0001, 2000), phi_init = NULL){
  
  ## Coords: Matrix of coordinates of dim L*2
  ## Dados: A matrix of dimension p*r
  ## X: Matrix of dimentsion sxr: r=TT*L, x : number of covar observed for each data point
  ## Beta of dimension p*s
  s<-nrow(X)
  L<-nrow(coords)
  p<-nrow(dados)
  r<-ncol(dados)
  TT<-r/L
  
  loglik<-numeric()
  criterio<-1
  count<-0
  
  # initial values
  
  Beta<-dados%*%t(X)%*%solve(X%*%t(X))
  
 
   # inicial values
  rho<- 0.5
  phi <- if (is.null(phi_init)) median(as.matrix(dist(coords))[as.matrix(dist(coords)) > 0]) else phi_init
  sigma2<- 0.5
  #construct matrix 
  PsiT<- ar1_cor (TT, rho = rho)
  PsiScorr<-spatial_corr_matrix(coords, corr_fun = corr, phi,nu)
  PsiS<- PsiScorr*sigma2
  Sigma<- diag(p)
  Psi<- kronecker(PsiS, PsiT)
  Psicorr<-kronecker(PsiScorr, PsiT)
  #Psicorr <- (Psicorr + t(Psicorr)) / 2  # symmetry
  
  while(criterio > precision){
    
    count <- count + 1
    ##update beta
    
    PsiS_inv  <- solve(PsiS)
    PsiT_inv     <- ar1_inv(TT, rho)
    PsiScorr_inv <- solve(PsiScorr) 
    
    Psi_inv   <- kronecker(PsiS_inv, PsiT_inv)
    Psicorr_inv  <- kronecker(PsiScorr_inv, PsiT_inv)
    Beta <- dados %*% Psi_inv %*% t(X) %*% solve(X %*% Psi_inv %*% t(X))
    
    mu<-Beta%*%X
    
    #update sigma2 
    PsiScorr_inv <- solve(PsiScorr)
    Psicorr_inv  <- kronecker(PsiScorr_inv, PsiT_inv)   # reuse PsiT_inv from above
    
    Aux <- solve(Sigma) %*% (dados - mu) %*% Psicorr_inv %*% t(dados - mu)
    
    sigma2_new <- 1/(p*r) * sum(diag(Aux))
    sigma2 <- max(min(sigma2_new, 1e4), 1e-6)  # truncar entre 1e-6 e 1e4
    print(sigma2)

    # update Sigma
    Sigma <- 1/r * (dados - mu) %*% Psi_inv %*% t(dados - mu)
    Sigma <- (Sigma + t(Sigma)) / 2

    # Forçar definição positiva
    ev    <- eigen(Sigma, symmetric = TRUE)
    lam   <- pmax(ev$values, 1e-6)
    Sigma <- ev$vectors %*% diag(lam) %*% t(ev$vectors)

    # Identificabilidade — reescalar pela diagonal (1,1)
    Sigma <- Sigma / Sigma[1,1]


   # Forçar definição positiva via correção espectral
    ev  <- eigen(Sigma, symmetric = TRUE)
    lam <- pmax(ev$values, 1e-6)
    Sigma <- ev$vectors %*% diag(lam) %*% t(ev$vectors)
    Sigma <- Sigma + diag(p) * max(0, -min(eigen(Sigma, only.values=TRUE)$values) + 1e-6)
    #update rho
    PsiS<- PsiScorr*sigma2
    #PsiS <- (PsiS + t(PsiS)) / 2  # symmetry
    Erho<-optimize(loglikelSpatiorho, c(0.0001,0.999), tol = 0.0000000001, maximum = TRUE, PsiS, Sigma, dados, mu, TT)
    rho<-Erho$maximum
    #update phi
    PsiT<- ar1_cor (TT, rho = rho)
    #PsiT <- (PsiT + t(PsiT)) / 2  # symmetry
    Ephi<-optimize(loglikelSpatiophi, phi_range, tol = 0.0000000001, maximum = TRUE, PsiT, Sigma, sigma2, dados, coords, mu,  corr=corr, nu=nu)
    phi<-Ephi$maximum
    
    #loop until log-likelihood falls below precision
    loglik[count]<-Ephi$objective
    PsiScorr<-spatial_corr_matrix(coords, corr_fun = corr, phi,nu)
    PsiS<- PsiScorr*sigma2
    Psi<-kronecker(PsiS, PsiT)
    #Psi <- (Psi + t(Psi)) / 2  # symmetry
    Psicorr<-kronecker(PsiScorr, PsiT)
    #Psicorr <- (Psicorr + t(Psicorr)) / 2  # sy  mmetry
    if (count > 2) {
  if (is.na(loglik[count]) || is.na(loglik[count-1]) || 
      !is.finite(loglik[count]) || !is.finite(loglik[count-1])) {
    criterio <- 0.000000000000001
  } else {
    criterio <- abs(1 - loglik[count-1] / loglik[count])
  }
}

if (count == MaxIter) {
  criterio <- 0.000000000000001
}
    
  }
  # Beta, rho, phi, sigma2, Sigma
  #compute criteria
  npar <- (p*s)+1+1+1+p*(p+1)/2
  BIC   <- -2*loglik[count] + npar*log(p*TT*L)    # to be minimized
  
  obj.out <- list(Beta=Beta, rho=rho, phi=phi, sigma2=sigma2, Sigma=Sigma, loglik=loglik, BIC=BIC, iter = count)
  
  class(obj.out) <- "ML.MatrixRegreSPatioT"
  return(obj.out)
  
}


#posterior_df <- as.data.frame(stan_exp_exp_model)
#stan_exp_exp_model["beta0"]
#mean(posterior_df$nugget)
