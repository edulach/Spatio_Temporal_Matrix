
################################################################################
### Covariances
################################################################################

# Temporal AR(1)
ar1_cor <- function(n, rho){
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) -
                    (1:n - 1))
  rho^exponent
}


# Spatio

spatial_corr_matrix <- function(coords, corr_fun = "exponential", phi, nu = NULL) {
  # Ensure coordinates are in matrix form.
  coords <- as.matrix(coords)
  # Compute the Euclidean distance matrix.
  D <- as.matrix(dist(coords))
  #xtable(as.matrix(dist(coords)))
  #xtable(as.matrix(dist(coords)), align = "r", digits = 3)
  #D <- D/ max(D, na.rm = TRUE)
  
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
    if (is.null(nu)) {
      stop("For matern correlation, please provide the nu parameter (smoothness).")
    }
    r <- D / phi
    Corr <- (1 / (2^(nu - 1) * gamma(nu))) * (r^nu) * besselK(r, nu)
    Corr[is.na(Corr) | is.infinite(Corr)] <- 1
  } else {
    stop("Unknown corr_fun. Please use 'exponential', 'gaussian', or 'matern'.")
  }
  
  # Ensure symmetry and diagonal = 1
  Corr <- (Corr + t(Corr)) / 2 ##symmetry
  diag(Corr) <- 1
  return(Corr)
}


#### Likelihood complete data ####

loglikel<-function(dados, mu, Sigma, PsiT, PsiS){
  ### Sigma: variables dependence
  ## PsiT: temporal dependence
  ## PsiS: spatial dependence
  Psi<- kronecker(PsiS, PsiT)
  like<- matrixNormal::dmatnorm(dados, mu, Sigma, Psi, tol = .Machine$double.eps^0.5, log = TRUE)
  return(like)
}


loglikelSpatiorho<-function(rho, PsiS, Sigma, dados, mu, TT){
  PsiT<-ar1_cor(TT,rho)
  Psi<-kronecker(PsiS, PsiT)
  Psi <- (Psi + t(Psi)) / 2
  like<-matrixNormal::dmatnorm(dados, mu, Sigma, Psi, tol = .Machine$double.eps^0.5, log = TRUE)
  return(like)
}


loglikelSpatiophi<-function(phi, PsiT, Sigma, sigma2, dados, coords, mu, corr="exponential",nu = NULL){
  #phi <- exp(theta)
  PsiS<-spatial_corr_matrix(coords, corr_fun = corr, phi, nu = nu)*sigma2
  # Add small nugget for stability
  #epsilon <- 1e-6
  #PsiS <- (PsiS + t(PsiS)) / 2
  #diag(PsiS) <- 1 + epsilon
  # Same for temporal matrix
  #PsiT <- (PsiT + t(PsiT)) / 2
  #diag(PsiT) <- 1 + epsilon
  Psi<-kronecker(PsiS, PsiT)
  Psi <- (Psi + t(Psi)) / 2
  like<-matrixNormal::dmatnorm(dados, mu, Sigma, Psi, tol = .Machine$double.eps^0.5, log = TRUE)
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
  
  X.or <- matrixNormal::rmatrixnorm(M=M, U=Sigma, V=Psi)
  return(list(Y=X.or))
}



################################################################################
#### Maximum Likelihood for complete data : Spatio-temporal ####################
#### Beta full matrix and M=Beta*X
################################################################################

ML.MatrixRegreSPatioT<-function(dados, X, coords, precision=0.0000001, MaxIter=100, corr="exponential",nu = NULL){
  
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
  #  Beta[2, 3] <- 0
  #  Beta[2, 2] <- 0 
  rho<- 0.5
  phi<- 0.1 #*min(dist(coords))  # Better init
  sigma2<- 1
  #corr="gaussian"
  
  #construct matrix 
  PsiT<- ar1_cor (TT, rho = rho)
  PsiScorr<-spatial_corr_matrix(coords, corr_fun = corr, phi,nu)
  PsiS<- PsiScorr*sigma2
  Sigma<- diag(p)
  Sigma[1, 1] <- 1  # Fix reference variance
  Psi<- kronecker(PsiS, PsiT)
  Psicorr<-kronecker(PsiScorr, PsiT)
  #Psicorr <- (Psicorr + t(Psicorr)) / 2  # symmetry
  
  
  
  while(criterio > precision){
    
    count <- count + 1
    ##update beta
    Beta<-dados%*%solve(Psi)%*%t(X)%*%solve(X%*%solve(Psi)%*%t(X))
    #Beta[2, 3] <- 0
    # Beta[2, 2] <- 0 
    mu<-Beta%*%X
    
    #update sigma2 
    Aux<-solve(Sigma)%*%(dados-mu)%*%solve(Psicorr)%*%t(dados-mu)
    sigma2<-1/(p*r)*sum(diag(Aux))
    print(sigma2)
    #update Sigma
    Sigma<-1/r*(dados-mu)%*%solve(Psi)%*%t(dados-mu)
    Sigma[1,1]<-1
    #Sigma <- (Sigma + t(Sigma)) / 2 # symmetry
    #update rho
    PsiS<- PsiScorr*sigma2
    #PsiS <- (PsiS + t(PsiS)) / 2  # symmetry
    Erho<-optimize(loglikelSpatiorho, c(0.0001,0.999), tol = 1e-12, maximum = TRUE, PsiS, Sigma, dados, mu, TT)
    rho<-Erho$maximum
    #update phi
    PsiT<- ar1_cor (TT, rho = rho)
    #PsiT <- (PsiT + t(PsiT)) / 2  # symmetry
    Ephi<-optimize(loglikelSpatiophi, c(0,2000), tol = 1e-12, maximum = TRUE, PsiT, Sigma, sigma2, dados, coords, mu,  corr=corr,nu = nu)
    phi<-Ephi$maximum
    
    
    
    
    #loop until log-likelihood falls below precision
    loglik[count]<-Ephi$objective
    PsiScorr<-spatial_corr_matrix(coords, corr_fun = corr, phi,nu)
    PsiS<- PsiScorr*sigma2
    Psi<-kronecker(PsiS, PsiT)
    #Psi <- (Psi + t(Psi)) / 2  # symmetry
    Psicorr<-kronecker(PsiScorr, PsiT)
    #Psicorr <- (Psicorr + t(Psicorr)) / 2  # sy  mmetry
    if (count>2){
      # at<- (loglik[count]-loglik[count-1])/(loglik[count-1]-loglik[count-2])
      #criterio<-(loglik[count]-loglik[count-1])/(1-at)
      criterio<-abs(1-loglik[count-1]/loglik[count])
    }
    
    if (count==MaxIter){
      criterio <- 1e-12
    }
    
    
  }
  # Beta, rho, phi, sigma2, Sigma
  #compute criteria
  npar <- (p*s)+1+1+1+p*(p+1)/2
  #npar <- (p*s)+1+1+1+p*(p+1)/2
  BIC   <- -2*loglik[count] + npar*log(p*TT*L)    # to be minimized
  
  obj.out <- list(Beta=Beta, rho=rho, phi=phi, sigma2=sigma2, Sigma=Sigma,nu=nu, loglik=loglik, BIC=BIC, iter = count)
  
  class(obj.out) <- "ML.MatrixRegreSPatioT"
  return(obj.out)
  
}



