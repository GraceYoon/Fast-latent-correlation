# to draw the dimension of variables versus the computation time.
# want to offset of Kendall's tau calculation.

rm(list = ls())

library(microbenchmark)
library(mixedCCA)
load("Data/amgutpruned.rdata") # 6482 by 481 matrix.

plist <- c(2, 5, 10, 20, 50, 100, 200, 300, 400, 481)

for (dim in c(100, nrow(amgutpruned))){
  
  time_K <- time_org <- time_v2 <- matrix(NA, nrow = length(plist), ncol = 6) # the ncol is from microbenchmark print format.
  # min, Q1, mean, median, Q3, max -> 6 columns are needed.
  Kendall_sub <- Kcor_sub_org <- Kcor_sub_v2 <- list()
  
  for (i in 1:length(plist)){
    p <- plist[i]
    subdat <- amgutpruned[1:dim, 1:p]
    
    time_K[i, ] <- unlist(print(microbenchmark(Kendall_sub[[i]] <- mixedCCA::Kendall_matrix(subdat), times = 10), unit = "s")[1, 2:7]) # to use fixed unit: "seconds"
    
    time_org[i, ] <- unlist(print(microbenchmark(Kcor_sub_org[[i]] <- mixedCCA::estimateR(subdat, type = "trunc", method = "original", tol = 1e-6, verbose = TRUE)$R, times = 2), unit = "s")[1, 2:7]) # to use fixed unit: "seconds"
    
    time_v2[i, ] <- unlist(print(microbenchmark(Kcor_sub_v2[[i]] <- mixedCCA::estimateR(subdat, type = "trunc", method = "approx", tol = 1e-6, verbose = TRUE)$R, times = 10), unit = "s")[1, 2:7]) # to use fixed unit: "seconds"
    
    cat("Done with n = ", dim, " p = ", p, "\n\n\n")
  }
  
  save(time_K, time_org, time_v2, file = paste0("Data/RunTimePlot_v2_range_", dim, ".Rda"))
}




