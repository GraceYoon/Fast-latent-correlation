# functions to estimate rank-based correlation value between two variables


# load precomputed interpolation functions and source codes from mixedCCA R package in ``fastapprox'' branch.
load("~/Dropbox/mixedCCA/R/sysdata.rda")
source("~/Dropbox/mixedCCA/R/KendallTau.R") # need to install "pcaPP" package.
source("~/Dropbox/mixedCCA/R/bridge.R") # need to install "fMultivar" package.
source("~/Dropbox/mixedCCA/R/fromKtoR.R")
source("~/Dropbox/mixedCCA/R/bridgeInv.R")
source("~/Dropbox/mixedCCA/R/fromKtoR_ipol.R")


#############################################################################

# wrapper function using ORG method
estimatedR_mixed_org <- function(X1, X2, type1, type2, tol = 1e-6){
  if(type1 == "trunc"){
    if(type2 == "continuous"){ ### TC case
      R12 <- estimatedR_mixed_org_tc(X1, X2)
    } else if(type2 == "trunc"){ ### TT case
      Xcomb <- as.matrix(cbind(X1, X2))
      R12 <- estimatedR_mixed_org_tt(Xcomb = cbind(X1, X2))[1, 2]
    } else if(type2 == "binary"){ ### TB case
      R12 <- estimatedR_mixed_org_tb(X1, X2)
    }
  } else if(type1 == "binary"){
    if(type2 == "continuous"){ ### BC case
      R12 <- estimatedR_mixed_org_bc(X1, X2)
    } else if(type2 == "binary"){ ### BB case
      Xcomb <- as.matrix(cbind(X1, X2))
      R12 <- estimatedR_mixed_org_bb(Xcomb = cbind(X1, X2))[1, 2]
    }
  }
  return(R12)
}


# wrapper function using ML method
estimatedR_mixed_ml <- function(X1, X2, type1, type2){
  if(type1 == "trunc"){
    if(type2 == "continuous"){
      R12 <- estimatedR_mixed_ml_tc(X1, X2)
    } else if(type2 == "trunc"){
      Xcomb <- as.matrix(cbind(X1, X2))
      R12 <- estimatedR_mixed_ml_tt(Xcomb = cbind(X1, X2))[1, 2]
    } else if(type2 == "binary"){
      R12 <- estimatedR_mixed_ml_tb(X1, X2)
    }
  } else if(type1 == "binary"){
    if(type2 == "continuous"){
      R12 <- estimatedR_mixed_ml_bc(X1, X2)
    } else if(type2 == "binary"){
      Xcomb <- as.matrix(cbind(X1, X2))
      R12 <- estimatedR_mixed_ml_bb(Xcomb = cbind(X1, X2))[1, 2]
    }
  }
  return(R12)
}


# wrapper function using MLBD method
estimatedR_mixed_mlbd <- function(X1, X2, type1, type2, tol = 1e-6){
  if(type1 == "trunc"){
    if(type2 == "continuous"){
      R12 <- estimatedR_mixed_mlbd_tc(X1, X2)
    } else if(type2 == "trunc"){
      Xcomb <- as.matrix(cbind(X1, X2))
      R12 <- estimatedR_mixed_mlbd_tt(Xcomb = cbind(X1, X2))[1, 2]
    } else if(type2 == "binary"){
      R12 <- estimatedR_mixed_mlbd_tb(X1, X2)
    }
  } else if(type1 == "binary"){
    if(type2 == "continuous"){
      R12 <- estimatedR_mixed_mlbd_bc(X1, X2)
    } else if(type2 == "binary"){
      Xcomb <- as.matrix(cbind(X1, X2))
      R12 <- estimatedR_mixed_mlbd_bb(Xcomb = cbind(X1, X2))[1, 2]
    }
  }
  return(R12)
}

#############################################################################

# three types of functions for TC case

# ORG method
estimatedR_mixed_org_tc <- function(X1, X2, type1 = "trunc", type2 = "continuous", tol = 1e-6){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  n <- length(X1)
  p1 <- ncol(X1); p2 <- ncol(X2)
  # Only pasting TC case.
  zratio1 <- colMeans(X1 == 0)
  R12 <- fromKtoR_mixed(Kendall_matrix(X1, X2), zratio1 = zratio1, type1 = type1, type2 = type2, tol = tol)
  return(R12)
}
# ML method
estimatedR_mixed_ml_tc <- function(X1, X2, type1 = "trunc", type2 = "continuous"){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  n <- length(X1)
  p1 <- ncol(X1); p2 <- ncol(X2)
  
  # Only pasting TC case.
  zratio1 <- colMeans(X1 == 0)
  R12 <- bridgeInv_tc(Kendall_matrix(X1, X2), zratio1 = zratio1)
  return(R12)
}
# MLBD method
estimatedR_mixed_mlbd_tc <- function(X1, X2, type1 = "trunc", type2 = "continuous", tol = 1e-6){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  n <- length(X1)
  p1 <- ncol(X1); p2 <- ncol(X2)
  
  # Only pasting TC case.
  zratio1 <- colMeans(X1 == 0)
  K12 <- Kendall_matrix(X1, X2)
  if (abs(K12) > 0.9*(1-zratio1^2)){
    warning("Kendall tau= ", round(K12, 2), " where the truncation rate = ", zratio1, "\n")
    R12 <- fromKtoR_mixed(K12, zratio1 = zratio1, type1 = type1, type2 = type2, tol = tol)
  } else {
    R12 <- bridgeInv_tc(K12, zratio1 = zratio1)
  }
  return(R12)
}



#############################################################################

# three types of functions for TT case

# ORG method
estimatedR_mixed_org_tt <- function(Xcomb, type = "trunc", tol = 1e-6){
  Xcomb <- as.matrix(Xcomb)
  n <- nrow(Xcomb)
  p <- ncol(Xcomb)
  
  R <- matrix(1, p, p) # do not calculate diagonal part.
  zratio <- colMeans(Xcomb == 0)
  # only calculating off-diagonal part.
  R[1, 2] <- R[2, 1] <- fromKtoR_mixed(Kendall_matrix(Xcomb)[1, 2], zratio1 = zratio[1], zratio2 = zratio[2], type1 = type, type2 = type, tol = tol)
  return(R)
}
# ML method
estimatedR_mixed_ml_tt <- function(Xcomb, type = "trunc"){
  Xcomb <- as.matrix(Xcomb)
  n <- nrow(Xcomb)
  p <- ncol(Xcomb)
  R <- matrix(1, p, p) # do not calculate diagonal part.
  zratio <- colMeans(Xcomb == 0)
  # only calculating off-diagonal part.
  R[1, 2] <- R[2, 1] <- bridgeInv_tt(Kendall_matrix(Xcomb)[1, 2], zratio1 = zratio[1], zratio2 = zratio[2])
  return(R)
}
# MLBD method
estimatedR_mixed_mlbd_tt <- function(Xcomb, type = "trunc", tol = 1e-6){
  Xcomb <- as.matrix(Xcomb)
  n <- nrow(Xcomb)
  p <- ncol(Xcomb)
  R <- matrix(1, p, p) # do not calculate diagonal part.
  zratio <- colMeans(Xcomb == 0)
  # only calculating off-diagonal part.
  K12 <- Kendall_matrix(Xcomb)[1, 2]
  if (abs(K12) > 0.9*(1-max(zratio)^2)){
    warning("exceeded.", K12, "with trate = ", zratio, "\n")
    R[1, 2] <- R[2, 1] <- fromKtoR_mixed(Kendall_matrix(Xcomb)[1, 2], zratio1 = zratio[1], zratio2 = zratio[2], type1 = type, type2 = type, tol = tol)
  } else {
    R[1, 2] <- R[2, 1] <- bridgeInv_tt(K12, zratio1 = zratio[1], zratio2 = zratio[2])
  }
  return(R)
}


#############################################################################

# three types of functions for TB case

# ORG method
estimatedR_mixed_org_tb <- function(X1, X2, type1 = "trunc", type2 = "binary", tol = 1e-6){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  n <- length(X1)
  p1 <- ncol(X1); p2 <- ncol(X2)
  zratio1 <- colMeans(X1 == 0)
  zratio2 <- colMeans(X2 == 0)
  R12 <- fromKtoR_mixed(Kendall_matrix(X1, X2), zratio1 = zratio1, zratio2 = zratio2, type1 = type1, type2 = type2, tol = tol)
  return(R12)
}
# ML method
estimatedR_mixed_ml_tb <- function(X1, X2, type1 = "trunc", type2 = "binary"){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  n <- length(X1)
  p1 <- ncol(X1); p2 <- ncol(X2)
  
  zratio1 <- colMeans(X1 == 0)
  zratio2 <- colMeans(X2 == 0)
  R12 <- bridgeInv_tb(Kendall_matrix(X1, X2), zratio1 = zratio1, zratio2 = zratio2)
  return(R12)
}
# MLBD method
estimatedR_mixed_mlbd_tb <- function(X1, X2, type1 = "trunc", type2 = "binary", tol = 1e-6){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  n <- length(X1)
  p1 <- ncol(X1); p2 <- ncol(X2)
  
  zratio1 <- colMeans(X1 == 0)
  zratio2 <- colMeans(X2 == 0)
  cutoff <- 0.9*min(1-zratio1^2, 2*zratio2*(1-zratio2))
  K12 <- Kendall_matrix(X1, X2)
  if (abs(K12) > cutoff){
    warning("Kendall tau= ", round(K12, 2), " where the truncation rate = ", zratio1, "\n")
    R12 <- fromKtoR_mixed(K12, zratio1 = zratio1, zratio2 = zratio2, type1 = type1, type2 = type2, tol = tol)
  } else {
    R12 <- bridgeInv_tb(K12, zratio1 = zratio1, zratio2 = zratio2)
  }
  return(R12)
}


#############################################################################

# three types of functions for BC case

# ORG method
estimatedR_mixed_org_bc <- function(X1, X2, type1 = "binary", type2 = "continuous", tol = 1e-6){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  n <- length(X1)
  p1 <- ncol(X1); p2 <- ncol(X2)
  # Only pasting BC case.
  zratio1 <- colMeans(X1 == 0)
  R12 <- fromKtoR_mixed(Kendall_matrix(X1, X2), zratio1 = zratio1, type1 = type1, type2 = type2, tol = tol)
  return(R12)
}
# ML method
estimatedR_mixed_ml_bc <- function(X1, X2, type1 = "binary", type2 = "continuous"){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  n <- length(X1)
  p1 <- ncol(X1); p2 <- ncol(X2)
  
  # Only pasting BC case.
  zratio1 <- colMeans(X1 == 0)
  R12 <- bridgeInv_bc(Kendall_matrix(X1, X2), zratio1 = zratio1)
  return(R12)
}
# MLBD method
estimatedR_mixed_mlbd_bc <- function(X1, X2, type1 = "binary", type2 = "continuous", tol = 1e-6){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  n <- length(X1)
  p1 <- ncol(X1); p2 <- ncol(X2)
  
  # Only pasting BC case.
  zratio1 <- colMeans(X1 == 0)
  cutoff <- 2*zratio1*(1 - zratio1)
  K12 <- Kendall_matrix(X1, X2)
  if (abs(K12) > cutoff){
    warning("Kendall tau= ", round(K12, 2), " where the truncation rate = ", zratio1, "\n")
    R12 <- fromKtoR_mixed(K12, zratio1 = zratio1, type1 = type1, type2 = type2, tol = tol)
  } else {
    R12 <- bridgeInv_bc(K12, zratio1 = zratio1)
  }
  return(R12)
}


#############################################################################

# three types of functions for BB case

# ORG method
estimatedR_mixed_org_bb <- function(Xcomb, type = "binary", tol = 1e-6){
  Xcomb <- as.matrix(Xcomb)
  n <- nrow(Xcomb)
  p <- ncol(Xcomb)
  
  R <- matrix(1, p, p) # do not calculate diagonal part.
  zratio <- colMeans(Xcomb == 0)
  # only calculating off-diagonal part.
  R[1, 2] <- R[2, 1] <- fromKtoR_mixed(Kendall_matrix(Xcomb)[1, 2], zratio1 = zratio[1], zratio2 = zratio[2], type1 = type, type2 = type, tol = tol)
  return(R)
}
# ML method
estimatedR_mixed_ml_bb <- function(Xcomb, type = "binary"){
  Xcomb <- as.matrix(Xcomb)
  n <- nrow(Xcomb)
  p <- ncol(Xcomb)
  R <- matrix(1, p, p) # do not calculate diagonal part.
  zratio <- colMeans(Xcomb == 0)
  # only calculating off-diagonal part.
  R[1, 2] <- R[2, 1] <- bridgeInv_bb(Kendall_matrix(Xcomb)[1, 2], zratio1 = zratio[1], zratio2 = zratio[2])
  return(R)
}

# MLBD method
estimatedR_mixed_mlbd_bb <- function(Xcomb, type = "binary", tol = 1e-6){
  Xcomb <- as.matrix(Xcomb)
  n <- nrow(Xcomb)
  p <- ncol(Xcomb)
  R <- matrix(1, p, p) # do not calculate diagonal part.
  zratio <- colMeans(Xcomb == 0)
  # only calculating off-diagonal part.
  K12 <- Kendall_matrix(Xcomb)[1, 2]
  if (abs(K12) > 0.9*2*min(zratio)*(1-max(zratio))){
    warning("exceeded.", K12, "with trate = ", zratio, "\n")
    R[1, 2] <- R[2, 1] <- fromKtoR_mixed(Kendall_matrix(Xcomb)[1, 2], zratio1 = zratio[1], zratio2 = zratio[2], type1 = type, type2 = type, tol = tol)
  } else {
    R[1, 2] <- R[2, 1] <- bridgeInv_bb(K12, zratio1 = zratio[1], zratio2 = zratio[2])
  }
  return(R)
}