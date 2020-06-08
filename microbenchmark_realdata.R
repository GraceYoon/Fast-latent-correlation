# check computation time using microbenchmark


rm(list=ls())

library(microbenchmark) # time the computation
library(SPRING) # to use SPRING function and the data "QMP" saved in SPRING package.
library(chebpol) # interpolation is based on this package.


# install the R package from the working branch ``fastapprox'' using the code below.
devtools::install_github("irinagain/mixedCCA", ref = "fastapprox", force = TRUE)
library(mixedCCA)


###### QMP dataset in SPRING package
data(QMP) # dimension 106 by 91

##### check the computation time of calculating "Kendall's tau" matrix only.
microbenchmark(Kendall_QMP <- Kendall_matrix(QMP), times = 100)
# Unit: milliseconds
#                               expr      min       lq     mean   median       uq      max neval
# Kendall_QMP <- Kendall_matrix(QMP) 957.8031 993.3659 1103.528 1012.281 1039.986 2164.093   100
##### Median: 1012 milliseconds = 1.012 seconds
save(Kendall_QMP, file = "Data/Kendallmatrix_QMP.Rda")

##### check the whole computation time of using MLBD method
microbenchmark(Kcor_QMP_fast <- estimateR(QMP, type = "trunc", method = "approx", tol = 1e-6, verbose = TRUE)$R, times = 100)
# Unit: seconds
# expr
# Kcor_QMP_fast <- estimateR(QMP, type = "trunc", method = "approx", tol = 1e-06, verbose = TRUE)$R
#      min       lq     mean   median       uq     max  neval
# 1.163961 1.221555 1.281554 1.255546 1.297741 2.54252    100
save(Kcor_QMP_fast, file = "Data/Kcor_QMP_fast.Rda")

##### check the whole computation time of using ORG method
microbenchmark(Kcor_QMP_org <- estimateR(QMP, type = "trunc", method = "original", tol = 1e-6, verbose = TRUE)$R, times = 10)
# Unit: seconds
#                                                                                             expr
# Kcor_QMP_org <- estimateR(QMP, type = "trunc", method = "original", tol = 1e-06, verbose = TRUE)$R
#     min       lq    mean   median       uq     max neval
# 52.97246 53.01858 53.1843 53.11137 53.27196 53.6568   10
save(Kcor_QMP_org, file = "Data/Kcor_QMP_org.Rda")

###### amgut dataset in SPRING package
load("Data/amgutpruned.rdata") # pruned amgut data is saved.
# 6482 by 481 matrix.

##### check the computation time of calculating "Kendall's tau" matrix only.
microbenchmark(Kendall_amgut <- Kendall_matrix(amgutpruned), times = 10)
# Unit: seconds
#       expr
# Kendall_amgut <- Kendall_matrix(amgutpruned)
#    min       lq     mean   median       uq      max neval
# 332.53 332.8834 334.8921 333.0824 333.3579 352.0477    10
save(Kendall_amgut, file = "Data/Kendallmatrix_amgut.Rda")

##### check the whole computation time of using MLBD method
microbenchmark(Kcor_amgut_fast <- estimateR(amgutpruned, type = "trunc", method = "approx", tol = 1e-6, verbose = TRUE)$R, times = 10)
# Unit: seconds
#     expr
# Kcor_amgut_fast <- estimateR(amgutpruned, type = "trunc", method = "approx", tol = 1e-06, verbose = TRUE)$R
#      min       lq     mean   median       uq      max neval
# 341.1044 341.4393 341.7057 341.5378 342.0269 342.6305    10
save(Kcor_amgut_fast, file = "Data/Kcor_amgut_fast.Rda")

##### check the whole computation time of using ORG method
microbenchmark(Kcor_amgut_org <- estimateR(amgutpruned, type = "trunc", method = "original", tol = 1e-6, verbose = TRUE)$R, times = 1)
# Unit: seconds
#   expr
# Kcor_amgut_org <- estimateR(amgutpruned, type = "trunc", method = "original", tol = 1e-06, verbose = TRUE)$R
#      min       lq     mean   median       uq      max neval
# 3102.334 3102.334 3102.334 3102.334 3102.334 3102.334     1
save(Kcor_amgut_org, file = "Data/Kcor_amgut_org.Rda")

###### TCGA dataset in mixedCCA package
load("Data/matchedTCGA.Rdata")
# X1: 500 by 981 continuous, X2: 500 by 431 zero-inflated

##### check the computation time of calculating "Kendall's tau" matrix only.
microbenchmark(Kendall_TCGA <- Kendall_matrix(cbind(matchedTCGA$X1, matchedTCGA$X2)), times = 10)
# Unit: seconds
#   expr
# Kendall_TCGA <- Kendall_matrix(cbind(matchedTCGA$X1, matchedTCGA$X2))
#      min       lq     mean   median       uq      max neval
# 332.9015 333.1105 333.5827 333.4394 333.8567 334.9375    10
save(Kendall_TCGA, file = "Data/Kendallmatrix_tcga")

##### check the whole computation time of using MLBD method
microbenchmark(Kcor_TCGA_fast <- estimateR_mixed(X1 = matchedTCGA$X1, X2 = matchedTCGA$X2, type1 = "continuous", type2 = "trunc", method = "approx", tol = 1e-6, verbose = TRUE)$R, times = 10)
# Unit: seconds
# expr
# Kcor_TCGA_fast <- estimateR_mixed(X1 = matchedTCGA$X1, X2 = matchedTCGA$X2, type1 = "continuous", type2 = "trunc", method = "approx", tol = 1e-06, verbose = TRUE)$R
#      min       lq     mean   median       uq       max neval
# 497.9776 498.3583 499.4291 499.6474 500.1961  500.9116    10
save(Kcor_TCGA_fast, file = "Data/Kcor_TCGA_fast.rdata")

##### check the whole computation time of using ORG method
microbenchmark(Kcor_TCGA_org <- estimateR_mixed(X1 = matchedTCGA$X1, X2 = matchedTCGA$X2, type1 = "continuous", type2 = "trunc", method = "original", tol = 1e-6, verbose = TRUE)$R, times = 1)
# Unit: seconds
#                                                                                  expr
# Kcor_TCGA_org <- estimateR_mixed(X1 = matchedTCGA$X1, X2 = matchedTCGA$X2, type1 = "continuous", type2 = "trunc", method = "original", tol = 1e-06, verbose = TRUE)$R
#      min       lq     mean   median       uq      max neval
# 2490.727 2490.727 2490.727 2490.727 2490.727 2490.727     1
save(Kcor_TCGA_org, file = "Data/Kcor_TCGA_org.rdata")


###### QMP dataset in SPRING package
source("SPRING_v2.R") # revised based on updated "estimateR" function with method = "approx".

n = dim(QMP)[1]
p = dim(QMP)[2]
rep.num = 50 # the repetition number of subsampling
nlam = 50 # the number of lambda sequence values
lambda.min.ratio = 1e-2 # the ratio of lambda.min over lambda.max
thresh = 0.1 # threshold for StARS criterion
subsample.ratio = 0.8 # subsample size ratio over the total samples
nc = 2 # numbere of cores for subsampling in parallel mode
seed = 10010 # seed for subsampling


### SPRING with ORG. The lambda sequence is generated based on estimated rank-based correlation.
ptm <- proc.time()
fit.spring <- SPRING::SPRING(QMP, quantitative = TRUE, lambdaseq = "data-specific", lambda.min.ratio = lambda.min.ratio, nlambda = nlam, subsample.ratio = subsample.ratio, seed = seed, ncores = nc, rep.num = rep.num)
proc.time() - ptm # 1353 seconds = 23 minutes, # 1910.154 seconds = 31.83 minutes.
save(fit.spring, file = "Data/QMP_SPRING_org_ddlam.rdata")


### SPRING with ORG. The lambda sequence is the default lambda sequence from 0.6 to 0.006 of length 50.
ptm <- proc.time()
fit.spring_deflam <- SPRING::SPRING(QMP, quantitative = TRUE, lambda.min.ratio = lambda.min.ratio, nlambda = nlam, thresh = thresh, subsample.ratio = subsample.ratio, seed = seed, ncores = nc, rep.num = rep.num)
proc.time() - ptm # 1697.584 seconds = 28.23 minutes
save(fit.spring_deflam, file = "Data/QMP_SPRING_org_defaultlam.rdata")


### SPRING with MLBD. The lambda sequence is generated based on estimated rank-based correlation.
ptm <- proc.time()
fit.spring_fast <- SPRING_v2(QMP, quantitative = TRUE, lambdaseq = "data-specific", lambda.min.ratio = lambda.min.ratio, nlambda = nlam, thresh = thresh, subsample.ratio = subsample.ratio, seed = seed, ncores = nc, rep.num = rep.num)
proc.time() - ptm
# user  system elapsed 
# 151.984  51.655 105.131 
save(fit.spring_fast, file = "Data/QMP_SPRING_fast_ddlam.rdata")


### SPRING with MLBD. The lambda sequence is the default lambda sequence from 0.6 to 0.006 of length 50.
ptm <- proc.time()
fit.spring_fast_deflam <- SPRING_v2(QMP, quantitative = TRUE, lambda.min.ratio = lambda.min.ratio, nlambda = nlam, thresh = thresh, subsample.ratio = subsample.ratio, seed = seed, ncores = nc, rep.num = rep.num)
proc.time() - ptm
# user  system elapsed 
# 150.475  47.295 101.142 
save(fit.spring_fast_deflam, file = "Data/QMP_SPRING_fast_defaultlam.rdata")

