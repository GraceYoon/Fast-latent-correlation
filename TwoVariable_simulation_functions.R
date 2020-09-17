# Chapter 4. Comparison on simulated data
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
  # K is d1 by d1 square matrix. zratio should be a vector of length d1. If the type is "continuous", no input is necessary for zratio.
  
  if (type == "continuous") {
    hatR <- sin(pi/2 * K)
  } else { # if the type is either "trunc" or "binary"
    
    if (is.null(zratio)){ stop ("The input for zratio is required for \"trunc\" and \"binary\" types.") }
    if (length(zratio)!=d1){ stop ("The length of zratio must match with the number of columns in K.") }
    
    # based on the data type, select bridgeInv and cutoff functions.
    bridgeInv <- bridgeInv_select(type1 = type, type2 = type)
    
    # much faster multi-linear interpolation part using saved ipol function.
    hatR <- bridgeInv(K, zratio1 = zratio, zratio2 = zratio)
    
    # make sure the diagonal elements are all one and they are symmetric.
    diag(hatR) <- 1
    hatR <- (hatR + t(hatR))/2 # even without this step, hatR is very close to symmetric but not exactly. (symmetric within the error 1e-5)

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
  # K12 is d1 by d2 matrix.
  # zratio1 should be a vector of length d1. zratio2 should be a vector of length d2.
  # If the both types are "continuous", no input is necessary for zratio1 and zratio2.
  
  if (type1 == "continuous" & type2 == "continuous") {
    hatR <- sin(pi/2 * K12)
  } else {
    # if the case is either CT, TC, TT, BC, BB or TB.
    if (type1 != "continuous"){
      if (is.null(zratio1)){ stop ("The input for zratio1 is required for \"trunc\" and \"binary\" types.") }
      if (length(zratio1)!=d1){ stop ("The length of zratio1 must match with the number of rows in K12.") }
    }
    if (type2 != "continuous"){
      if (is.null(zratio2)){ stop ("The input for zratio2 is required for \"trunc\" and \"binary\" types.") }
      if (length(zratio2)!=d2){ stop ("The length of zratio2 must match with the number of columns in K12.") }
    }
    
    # based on the data type, select bridgeInv and cutoff functions.
    bridgeInv <- bridgeInv_select(type1 = type1, type2 = type2)
  
    # much faster multi-linear interpolation part using saved ipol function.
    hatR <- bridgeInv(K12, zratio1 = zratio1, zratio2 = zratio2)
    # Here, there is no step to make sure about the diagonal elements and symmetrize.
  }
  return(hatR)
}


estimateR_mlonly <- function(X, type = "trunc", method = "mlonly", use.nearPD = TRUE, nu = 0.01, tol = 1e-3, verbose = FALSE){
  X <- as.matrix(X)
  
  n <- nrow(X)
  p <- ncol(X)
  
  if (!(type %in% c("continuous", "binary","trunc"))){
    stop("Unrecognized type of data. Should be one of continuous, binary or trunc.")
  }
  ind.NaN <- NULL # initialization of the variable. The default is assuming there is no NA or NaN.
  
  if (type == "continuous"){
    K <- pcaPP::cor.fk(X)
    R <- sin(pi/2 * K)
  } else if (type == "trunc"){
    # checking data type
    if(sum(X < 0) > 0) {
      stop("The data contains negative values.")
    }
    # checking proportion of zero values
    zratio <- colMeans(X == 0)
    if(sum(zratio) == 0){
      message("The data does not contain zeros. Consider changing the type to \"continuous\".")
    }
    if (sum(zratio == 1) > 0){
      warning("There are variables in the data that have only zeros.\n")
    }
    
    K <- Kendall_matrix(X)
    if(sum(is.na(K)) > 0){
      warning("There are NaN values in Kendall's tau matrix.\n")
      ind.NaN <- which(colSums(is.na(K)) == (p-1))
      K <- K[-ind.NaN, -ind.NaN]
      zratio <- zratio[-ind.NaN]
    }
    
    if(method == "mlonly"){
      R <- fromKtoR_ml_nocutoff(K, zratio = zratio, type = type, tol = tol)
    } 
    
  } else if (type == "binary"){
    # checking data type
    if(sum(!(X %in% c(0, 1))) > 0) {
      stop("The data is not \"binary\".")
    }
    # checking proportion of zero values
    zratio <- colMeans(X == 0)
    if (sum(zratio == 1) > 0 | sum(zratio == 0) > 0){
      warning("There are binary variables in the data that have only zeros or only ones.\n")
    }
    
    K <- Kendall_matrix(X)
    if(sum(is.na(K)) > 0){
      warning("There are NaN values in Kendall's tau matrix.\n")
      ind.NaN <- which(colSums(is.na(K)) == (p-1))
      K <- K[-ind.NaN, -ind.NaN]
      zratio <- zratio[-ind.NaN]
    }
    
    if(method == "mlonly"){
      R <- fromKtoR_ml_nocutoff(K, zratio = zratio, type = type, tol = tol)
    } 
  }
  
  # nearPD
  if ( use.nearPD == TRUE & min(eigen(R)$values) < 0 ) {
    if( verbose ){
      message(" minimum eigenvalue of correlation estimator is ", min(eigen(R)$values), "\n nearPD is used")
    }
    R <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
  }
  # shrinkage method
  if(nu < 0 | nu > 1){
    stop("nu must be be between 0 and 1.")
  }
  
  if (length(ind.NaN) > 0){
    R <- (1 - nu)*R + nu*diag(p - length(ind.NaN))
    R.final <- matrix(NaN, nrow = p, ncol = p)
    R.final[-ind.NaN, -ind.NaN] <- R
  } else if (length(ind.NaN) == 0){
    R.final <- (1 - nu)*R + nu*diag(p)
  }
  
  ### To keep the correct column names for each matrices
  if(length(colnames(X)) == p){
    colnames(R.final) <- rownames(R.final) <- c(colnames(X))
  }
  
  return(list(type = type, R = R.final))
}

estimateR_mixed_mlonly <- function(X1, X2, type1 = "trunc", type2 = "continuous", method = "mlonly", use.nearPD = TRUE, nu = 0.01, tol = 1e-3, verbose = FALSE){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  
  if (nrow(X1) != nrow(X2)){ # Check of they have the same sample size.
    stop ("X1 and X2 must have the same sample size.")
  }
  
  ind.NaN <- NULL  # initialization of the variable. The default is assuming there is no NA or NaN.
  
  n <- length(X1)
  p1 <- ncol(X1); p2 <- ncol(X2)
  
  if (sum(c(type1, type2) %in% c("continuous", "binary", "trunc")) != 2){
    stop("Unrecognised type of variables. Should be one of continuous, binary or trunc.")
  }
  
  if (type1 == "trunc"){
    if(sum(X1 < 0) > 0) {
      stop("The data X1 contains negative values.")
    }
    zratio1 <- colMeans(X1 == 0)
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
    zratio1 <- colMeans(X1 == 0)
    if (sum(zratio1 == 1) > 0 | sum(zratio1 == 0) > 0){
      warning("There are binary variables in the data that have only zeros or only ones.\n")
    }
  }
  
  if (type2 == "trunc"){
    if(sum(X2 < 0) > 0) {
      stop("The data X2 contains negative values.")
    }
    zratio2 <- colMeans(X2 == 0)
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
    zratio2 <- colMeans(X2 == 0)
    if (sum(zratio2 == 1) > 0 | sum(zratio2 == 0) > 0){
      warning("There are binary variables in the data that have only zeros or only ones.\n")
    }
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
    if (type1 == "continuous"){
      ################### These are CT or CB case.
      K1 <- pcaPP::cor.fk(X1)
      R1 <- sin(pi/2 * K1)
      
      zratio2 <- colMeans(X2 == 0)
      # check Kendall's tau calculation.
      K2 <- Kendall_matrix(X2)
      K12 <- Kendall_matrix(X1, X2)
      if(sum(is.na(K2)) + sum(is.na(K12)) > 0){
        warning("There are NaN values in Kendall's tau matrix.\n")
        ind.NaN <- which(colSums(is.na(K2)) == (p2-1))
        K12 <- K12[, -ind.NaN]
        K2 <- K2[-ind.NaN, -ind.NaN]
        zratio2 <- zratio2[-ind.NaN]
        ind.NaN <- ind.NaN + p1
      }
      
      if(method == "mlonly"){
        R2 <- fromKtoR_ml_nocutoff(K2, zratio = zratio2, type = type2, tol = tol)
        R12 <- fromKtoR_ml_mixed_nocutoff(K12, zratio2 = zratio2, type1 = type1, type2 = type2, tol = tol)
      } 
      Rall <- rbind(cbind(R1, R12), cbind(t(R12), R2))
      
    } else if (type2 == "continuous"){
      ################### These are TC or BC case.
      K2 <- pcaPP::cor.fk(X2)
      R2 <- sin(pi/2 * K2)
      
      zratio1 <- colMeans(X1 == 0)
      # check Kendall's tau calculation.
      K1 <- Kendall_matrix(X1)
      K12 <- Kendall_matrix(X1, X2)
      if(sum(is.na(K1)) + sum(is.na(K12)) > 0){
        warning("There are NaN values in Kendall's tau matrix.\n")
        ind.NaN <- which(colSums(is.na(K1)) == (p1-1))
        K12 <- K12[-ind.NaN, ]
        K1 <- K1[-ind.NaN, -ind.NaN]
        zratio1 <- zratio1[-ind.NaN]
      }
      
      if(method == "mlonly"){
        R1 <- fromKtoR_ml_nocutoff(K1, zratio = zratio1, type = type1, tol = tol)
        R12 <- fromKtoR_ml_mixed_nocutoff(K12, zratio1 = zratio1, type1 = type1, type2 = type2, tol = tol)
      } 
      Rall <- rbind(cbind(R1, R12), cbind(t(R12), R2))
      
    } else {
      ################### These are TB or BT case.
      
      zratio1 <- colMeans(X1 == 0)
      zratio2 <- colMeans(X2 == 0)
      
      # check Kendall's tau calculation.
      K1 <- Kendall_matrix(X1)
      K2 <- Kendall_matrix(X2)
      K12 <- Kendall_matrix(X1, X2)
      if(sum(is.na(K1)) + sum(is.na(K2)) + sum(is.na(K12)) > 0){
        warning("There are NaN values in Kendall's tau matrix.\n")
        ind.NaN1 <- which(colSums(is.na(K1)) == (p1-1))
        ind.NaN2 <- which(colSums(is.na(K2)) == (p2-1))
        K12 <- K12[-ind.NaN1, -ind.NaN2]
        K1 <- K1[-ind.NaN1, -ind.NaN1]
        K2 <- K2[-ind.NaN2, -ind.NaN2]
        zratio1 <- zratio1[-ind.NaN1]
        zratio2 <- zratio2[-ind.NaN2]
        ind.NaN <- c(ind.NaN1, ind.NaN2+p1)
      }
      
      if(method == "mlonly"){
        R1 <- fromKtoR_ml_nocutoff(K1, zratio = zratio1, type = type1, tol = tol)
        R2 <- fromKtoR_ml_nocutoff(K2, zratio = zratio2, type = type2, tol = tol)
        R12 <- fromKtoR_ml_mixed_nocutoff(K12, zratio1 = zratio1, zratio2 = zratio2, type1 = type1, type2 = type2, tol = tol)
      } 
      
      Rall <- rbind(cbind(R1, R12), cbind(t(R12), R2))
    }
    
    if ( use.nearPD == TRUE & min(eigen(Rall)$values) < 0 ) {
      if( verbose ){
        message(" minimum eigenvalue of correlation estimator is ", min(eigen(Rall)$values), "\n nearPD is used")
      }
      Rall <- as.matrix(Matrix::nearPD(Rall, corr = TRUE)$mat)
    }
    # shrinkage method
    if(nu < 0 | nu > 1){
      stop("nu must be be between 0 and 1.")
    }
    
    if (length(ind.NaN) > 0){
      Rall <- (1 - nu)*Rall + nu*diag(p1 + p2 - length(ind.NaN))
      R.final <- matrix(NaN, nrow = (p1 + p2), ncol = (p1 + p2))
      R.final[-ind.NaN, -ind.NaN] <- Rall
    } else if (length(ind.NaN) == 0){
      R.final <- (1 - nu)*Rall + nu*diag(p1 + p2)
    }
    
    ### To keep the correct column names for each matrices
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