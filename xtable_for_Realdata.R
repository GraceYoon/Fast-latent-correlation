# to extract saved results on real data application to table format.

rm(list = ls())


# load saved computation results (Realdata: QMP, AGP: American Gut Project, TCGA-BRCA)
load("Data/Realdata3_tab.Rda")

CompRes1 <- CompRes # Computation time for ORG, MLBD and Kendall
AccRes1 <- AccRes # Accuracy in terms of Maximum Absolute Error (AE) and Mean AE


# load saved computation results (SPRING: Semi-Parametric Rank-based approach for INference in Graphical model.)
load("Data/Realdata-SPRING_tab.Rda")

# dd: data driven lambda sequence (from 0.01*sigma_max to sigma_max), sigma_max is largest off-diagonal element in estimated correlation values. 
# df: default lambda sequence (from 0.006 to 0.6)
# both sequcnes are in log scale.

CompRes2 <- CompRes
AccRes2 <- AccRes


library(xtable)
xtable(cbind(CompRes1, rbind(CompRes2, NA)))


digitmat <- matrix(c(rep(0, 2), rep(c(4, 6), 3), rep(c(4, 7), 2)), nrow = 2)
xtable(cbind(AccRes1, AccRes2), digits = digitmat)

       