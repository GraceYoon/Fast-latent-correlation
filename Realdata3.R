# check real data results with new CRAN package

# to check accuracy & computation time of multilinear approximation

rm(list = ls())

# devtools::install_github("GraceYoon/SPRING", force = TRUE)
# devtools::install_github("zdk123/SpiecEasi", force = TRUE)

library(microbenchmark)
library(SPRING) # for QMP data
library(mixedCCA) # from CRAN

#################################################
####### for QMP data in SPRING package

microbenchmark(Kendall_QMP <- Kendall_matrix(QMP), times = 100)

microbenchmark(Kcor_QMP_org <- estimateR(QMP, type = "trunc", method = "original", tol = 1e-6, verbose = TRUE)$R, times= 10)

microbenchmark(Kcor_QMP_fast <- estimateR(QMP, type = "trunc", method = "approx", tol = 1e-6, verbose = TRUE)$R, times= 100)

max(abs(c(Kcor_QMP_fast - Kcor_QMP_org)))
mean(abs(c(Kcor_QMP_fast - Kcor_QMP_org)))

#################################################
###### for AGP data
load("Data/amgutpruned.rdata")

microbenchmark(Kendall_amgut <- Kendall_matrix(amgutpruned), times = 10)

microbenchmark(Kcor_amgut_org <- estimateR(amgutpruned, type = "trunc", method = "original", tol = 1e-6, verbose = TRUE)$R, times = 1)

microbenchmark(Kcor_amgut_fast <- estimateR(amgutpruned, type = "trunc", method = "approx", tol = 1e-6, verbose = TRUE)$R, times = 10)

max(abs(c(Kcor_amgut_fast - Kcor_amgut_org)))
mean(abs(c(Kcor_amgut_fast - Kcor_amgut_org)))

#################################################
###### for TCGA-BRCA data
load("Data/matchedTCGA.Rdata")

microbenchmark(Kendall_TCGA <- Kendall_matrix(cbind(matchedTCGA$X1, matchedTCGA$X2)), times = 10)

microbenchmark(Kcor_TCGA_org <- estimateR_mixed(X1 = matchedTCGA$X1, X2 = matchedTCGA$X2, type1 = "continuous", type2 = "trunc", method = "original", tol = 1e-6, verbose = TRUE)$R, times = 1)

microbenchmark(Kcor_TCGA_fast <- estimateR_mixed(X1 = matchedTCGA$X1, X2 = matchedTCGA$X2, type1 = "continuous", type2 = "trunc", method = "approx", tol = 1e-6, verbose = TRUE)$R, times = 10)

max(abs(c(Kcor_TCGA_fast - Kcor_TCGA_org)))
mean(abs(c(Kcor_TCGA_fast - Kcor_TCGA_org)))


save(Kcor_QMP_org, Kcor_QMP_fast, Kcor_amgut_org, Kcor_amgut_fast, Kcor_TCGA_org, Kcor_TCGA_fast, Kendall_QMP, Kendall_amgut, Kendall_TCGA, file = "Data/Realdata3.Rda")
