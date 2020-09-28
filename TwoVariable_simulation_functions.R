# Section 4. Comparison on simulated data
# In mixedCCA R package, there is no option whether we want to use ML or MLBD.
# Here, the following functions for ML.
# most of the parts are copied from mixedCCA R package but removed the "cutoff" part.


source("Functions/bridgeInv.R") # same functions as in mixedCCA R package.
load("Functions/sysdata.rda") # the same saved precomputed values as in mixedCCA R package.

# K: Kendall's tau matrix.
# zratio: a column vector of zero proportion values.
fromKtoR_ml_nocutoff <- function(K, zratio = NULL, type = "trunc", tol = 1e-3) {
  K <- as.matrix(K)
  d1 <- nrow(K)
  # adjust for 1 dimensional case
  if (d1 == 1){return(as.matrix(1))}
  
  if (type == "continuous") {
    hatR <- sin(pi/2 * K)
  } else { # if the type is either "trunc" or "binary"
    upperR <- c(upper.tri(K))
    hatRupper <- rep(NA, sum(upperR))
    Kupper <- c(K[upperR])
    
    # based on the data type, select bridgeInv and cutoff functions.
    bridgeInv <- bridgeInv_select(type1 = type, type2 = type)
    
    #### Form zratio vectors
    zratio1vec = rep(zratio, d1)[upperR]
    zratio2vec = rep(zratio, each = d1)[upperR]
    
    # much faster multi-linear interpolation part using saved ipol function.
    hatRupper <- bridgeInv(Kupper, zratio1 = zratio1vec, zratio2 = zratio2vec)
    
    # Get upperR into hatR
    hatR <- matrix(0, d1, d1)
    hatR[upperR] <- hatRupper
    hatR <- hatR + t(hatR)
    diag(hatR) <- 1
  }
  
  return(hatR)
}

# K12: Kendall's tau matrix.
# zratio1: a vector of zero proportion values for row variables. The length should match with nrow of K12.
# zratio2: a vector of zero proportion values for column variables. The length should match with ncol of K12.
fromKtoR_ml_mixed_nocutoff <- function(K12, zratio1 = NULL, zratio2 = NULL, type1 = "trunc", type2 = "continuous", tol = 1e-3) {
  
  K12 <- as.matrix(K12)
  d1 <- nrow(K12)
  d2 <- ncol(K12)
  
  if (type1 == "continuous" & type2 == "continuous") {
    hatR <- sin(pi/2 * K12)
  } else {
    # if the case is either CT, TC, TT, BC, BB or TB.
    
    # based on the data type, select bridgeInv and cutoff functions.
    bridgeInv <- bridgeInv_select(type1 = type1, type2 = type2)
    
    zratio1vec = rep(zratio1, d2)
    zratio2vec = rep(zratio2, each = d1)
    
    # much faster multi-linear interpolation part using saved ipol function.
    hatR <- matrix(bridgeInv(c(K12), zratio1 = zratio1vec, zratio2 = zratio2vec), d1, d2)
  }
  return(hatR)
}


estimateR_mlonly <- function(X, type = "trunc", use.nearPD = TRUE, nu = 0.01, tol = 1e-3, verbose = FALSE){
  X <- as.matrix(X)
  p <- ncol(X)
  
  # shrinkage method
  if(nu < 0 | nu > 1){
    stop("nu must be be between 0 and 1.")
  }
  
  if (!(type %in% c("continuous", "binary","trunc"))){
    stop("Unrecognized type of data. Should be one of continuous, binary or trunc.")
  }
  
  if (type == "continuous"){
    if (any(is.na(X))){
      # If there are any missing measurements, use slower function
      K <- cor(X, method = "kendall", use = "pairwise.complete.obs")
    }else{
      K <- pcaPP::cor.fk(X)
    }
    R <- sin(pi/2 * K)
  } else {
    zratio <- colMeans(X == 0)
    if (type == "trunc"){
      # checking data type
      if(sum(X < 0) > 0) {
        stop("The data of truncated type contains negative values.")
      }
      # checking proportion of zero values
      if(sum(zratio) == 0){
        message("The data does not contain zeros. Consider changing the type to \"continuous\".")
      }
      if (sum(zratio == 1) > 0){
        stop("There are variables in the data that have only zeros. Filter those     variables before continuing. \n")
      }
    }else{
      # checking data type
      if(sum(!(X %in% c(0, 1))) > 0) {
        stop("The data is not \"binary\".")
      }
      if (sum(zratio == 1) > 0 | sum(zratio == 0) > 0){
        stop("There are binary variables in the data that have only zeros or only ones. Filter those variables before continuing. \n")
      }
    }
    K <- Kendall_matrix(X)
    
    R <- fromKtoR_ml_nocutoff(K, zratio = zratio, type = type, tol = tol)
  }
  
  # nearPD to make it semi pos-definite
  if (use.nearPD == TRUE){
    if (min(eigen(R)$values) < 0) {
      if(verbose){
        message(" minimum eigenvalue of correlation estimator is ", min(eigen(R)$values), "\n nearPD is used")
      }
      R <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
    }
  }
  
  # Shrinkage adjustment by nu
  R.final <- (1 - nu) * R + nu * diag(p)
  
  ### To keep the correct column names for each matrices
  if(length(colnames(X)) == p){
    colnames(R.final) <- rownames(R.final) <- c(colnames(X))
  }
  
  return(list(type = type, R = R.final))
}

estimateR_mixed_mlonly <- function(X1, X2, type1 = "trunc", type2 = "continuous", use.nearPD = TRUE, nu = 0.01, tol = 1e-3, verbose = FALSE){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  
  if (nrow(X1) != nrow(X2)){ # Check of they have the same sample size.
    stop ("X1 and X2 must have the same sample size.")
  }
  
  # shrinkage method
  if(nu < 0 | nu > 1){
    stop("nu must be be between 0 and 1.")
  }
  
  p1 <- ncol(X1); p2 <- ncol(X2)
  
  if (sum(c(type1, type2) %in% c("continuous", "binary", "trunc")) != 2){
    stop("Unrecognised type of variables. Should be one of continuous, binary or trunc.")
  }
  
  zratio1 <- colMeans(X1 == 0)
  if (type1 == "trunc"){
    if(sum(X1 < 0) > 0) {
      stop("The data X1 contains negative values.")
    }
    if(sum(zratio1) == 0){
      message("The data X1 does not contain zeros. Consider changing the type to \"continuous\".")
    }
    if (sum(zratio1 == 1) > 0){
      warning("There are truncated variables in the data that have only zeros.\n")
    }
  }
  if (type1 == "binary"){
    if(sum(!(X1 %in% c(0, 1))) > 0) {
      stop("The data X1 is not \"binary\".")
    }
    if (sum(zratio1 == 1) > 0 | sum(zratio1 == 0) > 0){
      warning("There are binary variables in the data that have only zeros or only ones.\n")
    }
  }
  
  zratio2 <- colMeans(X2 == 0)
  if (type2 == "trunc"){
    if(sum(X2 < 0) > 0) {
      stop("The data X2 contains negative values.")
    }
    if(sum(zratio2) == 0){
      message("The data X2 does not contain zeros. Consider changing the type to \"continuous\".")
    }
    if (sum(zratio2 == 1) > 0){
      warning("There are truncated variables in the data that have only zeros.\n")
    }
  }
  if (type2 == "binary"){
    if(sum(!(X2 %in% c(0, 1))) > 0) {
      stop("The data X2 is not \"binary\".")
    }
    if (sum(zratio2 == 1) > 0 | sum(zratio2 == 0) > 0){
      warning("There are binary variables in the data that have only zeros or only ones.\n")
    }
  }
  
  if (p1 == 1 & p2 == 1){
    # This is just pairwise correlation
    k12 = KendallTau(X1, X2)
    r12 = fromKtoR_ml_mixed_nocutoff(k12, zratio1 = zratio1, zratio2 = zratio2, type1 = type1, type2 = type2)
    return(list(type = c(type1, type2), R1 = 1, R2 = 1, R12 = r12, R = matrix(c(1, r12, r12, 1), 2, 2)))
  }
  
  if (type1 == type2) {
    ################### both datasets are of the same type. CC, TT or BB case.
    
    Xcomb <- cbind(X1, X2)
    R.final <- estimateR_mlonly(Xcomb, type = type1, method = method, use.nearPD = use.nearPD, nu = nu, tol = tol)$R
    R1 <- R.final[1:p1, 1:p1]
    R2 <- R.final[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)]
    R12 <- R.final[1:p1, (p1 + 1):(p1 + p2)]
    
  } else {
    ################### datasets are of different type
    # Start with 1st dataset
    if (p1 == 1){
      R1 <- 1
    }else if (type1 == "continuous"){
      if (any(is.na(X1))){
        K1 <- cor(X1, method = "kendall", use = "pairwise.complete.obs")
      }else{
        K1 <- pcaPP::cor.fk(X1)
      }
      R1 <- sin(pi/2 * K1)
    }else{
      K1 <- Kendall_matrix(X1)
      R1 <- fromKtoR_ml_nocutoff(K1, zratio = zratio1, type = type1, tol = tol)
    }
    # Continue with 2nd dataset
    if (p2 == 1){
      R2 <- 1
    }else if (type2 == "continuous"){
      if (any(is.na(X2))){
        K2 <- cor(X2, method = "kendall", use = "pairwise.complete.obs")
      }else{
        K2 <- pcaPP::cor.fk(X2)
      }
      R2 <- sin(pi/2 * K2)
    }else{
      K2 <- Kendall_matrix(X2)
      R2 <- fromKtoR_ml_nocutoff(K2, zratio = zratio2, type = type2, tol = tol)
    }
    # Do cross-product
    K12 <- Kendall_matrix(X1, X2)
    
    R12 <- fromKtoR_ml_mixed_nocutoff(K12, zratio1 = zratio1, zratio2 = zratio2, type1 = type1, type2 = type2, tol = tol)
    
    Rall <- rbind(cbind(R1, R12), cbind(t(R12), R2))
    
    if (use.nearPD == TRUE){
      if(min(eigen(Rall)$values) < 0) {
        if(verbose) {
          message(" minimum eigenvalue of correlation estimator is ", min(eigen(Rall)$values), "\n nearPD is used")
        }
        Rall <- as.matrix(Matrix::nearPD(Rall, corr = TRUE)$mat)
      }
    }
    
    # Shrinkage step based on nu
    R.final <- (1 - nu) * Rall + nu * diag(p1 + p2)
    
    ### To keep the column names in R according to column names that are originally supplied in each matrix
    if(length(colnames(X1)) == p1 & length(colnames(X2)) == p2){
      colnames(R.final) <- rownames(R.final) <- c(colnames(X1), colnames(X2))
    } else if(length(colnames(X1)) != p1 & length(colnames(X2)) == p2){
      colnames(R.final) <- rownames(R.final) <- rep(NA, p1 + p2)
      colnames(R.final)[(p1 + 1):(p1 + p2)] <- rownames(R.final)[(p1 + 1):(p1 + p2)] <- colnames(X2)
    } else if(length(colnames(X1)) == p1 & length(colnames(X2)) != p2){
      colnames(R.final) <- rownames(R.final) <- rep(NA, p1 + p2)
      colnames(R.final)[1:p1] <- rownames(R.final)[1:p1] <- colnames(X1)
    }
    
    # For convenience, split the R matrices
    R1 <- R.final[1:p1, 1:p1]
    R2 <- R.final[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)]
    R12 <- R.final[1:p1, (p1 + 1):(p1 + p2)]
  }
  
  return(list(type = c(type1, type2), R1 = R1, R2 = R2, R12 = R12, R = R.final))
}
