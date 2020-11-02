# Compare the results (accuracy and computation time) on 3 real datasets
# WARNING: this will take a couple of hours to run as original method takes a lot of time on real data

rm(list = ls())

# Install the required packages first
# devtools::install_github("irinagain/mixedCCA", force = TRUE) # version 1.4.1
# devtools::install_github("GraceYoon/SPRING", force = TRUE)

library(microbenchmark)
library(SPRING) # for QMP data
library(mixedCCA)

# set up to save results
CompRes <- matrix(NA, nrow = 3, ncol = 3)
rownames(CompRes) <- c("ORG", "MLBD", "Kendall")
colnames(CompRes) <- c("QMP", "AGP", "TCGA-BRCA")

AccRes <- matrix(NA, nrow = 2, ncol = 3)
rownames(AccRes) <- c("Max AE", "Mean AE")
colnames(AccRes) <- c("QMP", "AGP", "TCGA-BRCA")


#################################################
####### for QMP data in SPRING package

colLoc <- 1 # QMP data result needs to be saved in 1st column


CompRes[1, colLoc] <- print(microbenchmark(Kcor_QMP_org <- estimateR(QMP, type = "trunc", method = "original", tol = 1e-6, verbose = TRUE)$R, times= 10), unit = "s")[, 5]

CompRes[2, colLoc] <- print(microbenchmark(Kcor_QMP_fast <- estimateR(QMP, type = "trunc", method = "approx", tol = 1e-6, verbose = TRUE)$R, times= 100), unit = "s")[, 5]

CompRes[3, colLoc] <- print(microbenchmark(Kendall_QMP <- Kendall_matrix(QMP), times = 100), unit = "s")[, 5]


AccRes[1, colLoc] <- max(abs(c(Kcor_QMP_fast - Kcor_QMP_org)))
AccRes[2, colLoc] <- mean(abs(c(Kcor_QMP_fast - Kcor_QMP_org)))

#################################################
###### for AGP data
load("Data/amgutpruned.rdata")

colLoc <- 2 # AGP data result needs to be saved in 1st column

CompRes[1, colLoc] <- print(microbenchmark(Kcor_amgut_org <- estimateR(amgutpruned, type = "trunc", method = "original", tol = 1e-6, verbose = TRUE)$R, times = 1), unit = "s")[, 5]

CompRes[2, colLoc] <- print(microbenchmark(Kcor_amgut_fast <- estimateR(amgutpruned, type = "trunc", method = "approx", tol = 1e-6, verbose = TRUE)$R, times = 10), unit = "s")[, 5]

CompRes[3, colLoc] <- print(microbenchmark(Kendall_amgut <- Kendall_matrix(amgutpruned), times = 10), unit = "s")[, 5]

AccRes[1, colLoc] <- max(abs(c(Kcor_amgut_fast - Kcor_amgut_org)))
AccRes[2, colLoc] <- mean(abs(c(Kcor_amgut_fast - Kcor_amgut_org)))


#################################################
###### for TCGA-BRCA data
load("Data/matchedTCGA.Rdata")

colLoc <- 3 # TCGA-BRCA data result needs to be saved in 3rd column

CompRes[1, colLoc] <- print(microbenchmark(Kcor_TCGA_org <- estimateR_mixed(X1 = matchedTCGA$X1, X2 = matchedTCGA$X2, type1 = "continuous", type2 = "trunc", method = "original", tol = 1e-6, verbose = TRUE)$R, times = 1), unit = "s")[, 5]

CompRes[2, colLoc] <- print(microbenchmark(Kcor_TCGA_fast <- estimateR_mixed(X1 = matchedTCGA$X1, X2 = matchedTCGA$X2, type1 = "continuous", type2 = "trunc", method = "approx", tol = 1e-6, verbose = TRUE)$R, times = 10), unit = "s")[, 5]


KendallTC <- function(X1, X2){
   K1 <- pcaPP::cor.fk(X1)
   K2 <- Kendall_matrix(X2)
   K12 <- Kendall_matrix(X1, X2)
   return(rbind(cbind(K1, K12), cbind(t(K12), K2)))
}

CompRes[3, colLoc] <- print(microbenchmark(Kendall_TCGA <- KendallTC(matchedTCGA$X1, matchedTCGA$X2), times = 10), unit = "s")[, 5]


AccRes[1, colLoc] <- max(abs(c(Kcor_TCGA_fast - Kcor_TCGA_org)))
AccRes[2, colLoc] <- mean(abs(c(Kcor_TCGA_fast - Kcor_TCGA_org)))


save(Kcor_QMP_org, Kcor_QMP_fast, Kcor_amgut_org, Kcor_amgut_fast, Kcor_TCGA_org, Kcor_TCGA_fast, Kendall_QMP, Kendall_amgut, Kendall_TCGA, CompRes, AccRes, file = "Data/Realdata3_tab.Rda")
