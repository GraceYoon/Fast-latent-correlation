# Section 4. to compare the computation time of latent correlation between two simulated variables


# to create a table of computation time.
rm(list = ls())

library(xtable)
options("scipen"=100, "digits"=4) # not to use expoential format



##### median of median values across 10 microbenchmark results.


tseq <- NULL
for (type in c("TC", "TT", "BC", "BB", "TB")){
  load(paste0("Data/TwoSim_", type, "_rep10.Rda"))
  
  tseq <- cbind(tseq,
                c(median(df_comptime$medianTime[df_comptime$method == "org"]),
                  median(df_comptime$medianTime[df_comptime$method == "ipol"]),
                  median(df_comptime$medianTime[df_comptime$method == "ipol_UB"]))
  )
  
  
}

colnames(tseq) <- c("TC", "TT", "BC", "BB", "TB")
rownames(tseq) <- c("ORG", "ML", "MLBD")
xtable(tseq)

# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrr}
# \hline
# & TC & TT & BC & BB & TB \\ 
# \hline
# ORG & 3767.28 & 24255.43 & 2516.40 & 2894.14 & 2177.13 \\ 
# ML & 350.72 & 454.88 & 352.73 & 446.21 & 452.80 \\ 
# MLBD & 362.72 & 479.92 & 368.92 & 483.67 & 493.77 \\ 
# \hline
# \end{tabular}
# \end{table}



# to create a figure of maximum absolute error.

rm(list = ls())

library(ggplot2)
library(tidyverse)
library(RColorBrewer)

##### Maximum Absolute Error
for (type in c("TC", "TT", "BC", "BB", "TB")){
  load(paste0("Data/TwoSim_", type, "_rep10.Rda"))
  g <- df_accuracy %>% ggplot(aes(x = as.factor(TruncRate), y = MaxAD, fill = method)) + geom_bar(stat = "identity", position = "dodge2", aes(fill = method, color = method, alpha = LatentR), size = 0.1) + scale_fill_discrete(name = "Method", label = c("ML", "MLBD")) + guides(color = F) + scale_alpha_continuous(breaks = seq(0.05, 0.91, length.out = 9)) + xlab("Zero proportion") + ylab("Maximum Absolute Error") 
  
  pdf(file = paste0("Figures/", type, "_MaxAD.pdf"), width=8, height=4)
  print(g)
  dev.off()
}


##### Mean Absolute Error
for (type in c("TC", "TT", "BC", "BB", "TB")){
  load(paste0("Data/TwoSim_", type, "_rep10.Rda"))
  g <- df_accuracy %>% ggplot(aes(x = as.factor(TruncRate), y = MeanAD, fill = method)) + geom_bar(stat = "identity", position = "dodge2", aes(fill = method, color = method, alpha = LatentR), size = 0.1) + scale_fill_discrete(name = "Method", label = c("ML", "MLBD")) + guides(color = F) + scale_alpha_continuous(breaks = seq(0.05, 0.91, length.out = 9)) + xlab("Zero proportion") + ylab("Mean Absolute Error") 
  
  pdf(file = paste0("Figures/", type, "_MeanAD.pdf"), width=8, height=4)
  print(g)
  dev.off()
}



