# Section 4. to compare the computation time of latent correlation between two simulated variables


# to create a table of computation time.
rm(list = ls())

library(xtable)
options("scipen"=100, "digits"=4) # not to use expoential format


##### 1 repetition result.
tseq <- NULL
for (type in c("TC", "TT", "BC", "BB", "TB")){
  load(paste0("Data/TwoSimVariable", type, "_symmZR.Rda"))
  
  tseq <- cbind(tseq,
                c(median(df_comptime$medianTime[df_comptime$method == "org"]),
                  median(df_comptime$medianTime[df_comptime$method == "ipol"]),
                  median(df_comptime$medianTime[df_comptime$method == "ipol_UB"]))
  )
  
  
}

colnames(tseq) <- c("TC", "TT", "BC", "BB", "TB")
rownames(tseq) <- c("ORG", "ML", "MLBD")
xtable(tseq)

# % latex table generated in R 4.0.2 by xtable 1.8-4 package
# % Thu Sep 17 00:22:18 2020
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrr}
# \hline
# & TC & TT & BC & BB & TB \\ 
# \hline
# ORG & 3731.44 & 20604.44 & 2793.95 & 3042.49 & 2408.03 \\ 
# ML & 793.00 & 759.23 & 790.37 & 766.62 & 983.85 \\ 
# MLBD & 877.77 & 804.12 & 891.44 & 828.17 & 1119.90 \\ 
# \hline
# \end{tabular}
# \end{table}


##### 10 repetition result.
tseq <- NULL
for (type in c("TC", "TT", "BC", "BB", "TB")){
  load(paste0("Data/TwoSimVariable", type, "_symmZR_rep10.Rda"))
  
  tseq <- cbind(tseq,
                c(median(df_comptime$medianTime[df_comptime$method == "org"]),
                  median(df_comptime$medianTime[df_comptime$method == "ipol"]),
                  median(df_comptime$medianTime[df_comptime$method == "ipol_UB"]))
  )
  
  
}

colnames(tseq) <- c("TC", "TT", "BC", "BB", "TB")
rownames(tseq) <- c("ORG", "ML", "MLBD")
xtable(tseq)

# % latex table generated in R 4.0.2 by xtable 1.8-4 package
# % Thu Sep 17 00:22:29 2020
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrr}
# \hline
# & TC & TT & BC & BB & TB \\ 
# \hline
# ORG & 3481.29 & 19261.48 & 2597.56 & 2771.11 & 2182.01 \\ 
# ML & 741.91 & 710.98 & 743.59 & 717.45 & 915.82 \\ 
# MLBD & 809.21 & 745.77 & 820.21 & 761.62 & 1034.39 \\ 
# \hline
# \end{tabular}
# \end{table}




# to create a figure of maximum absolute error.

rm(list = ls())

library(ggplot2)
library(tidyverse)
library(RColorBrewer)


for (type in c("TC", "TT", "BC", "BB", "TB")){
  load(paste0("Data/TwoSimVariable", type, "_symmZR.Rda"))
  g <- df_accuracy %>% ggplot(aes(x = TruncRate, y = MaxAD, group = interaction(as.factor(LatentR), method))) + geom_bar(stat = "identity", width = 0.09, position = "dodge", aes(fill = method, color = method, alpha = LatentR), size = 0.1) + scale_x_continuous(breaks = unique(df_accuracy$TruncRate))  + scale_fill_discrete(name = "Method", label = c("ML", "MLBD")) + guides(color = F) + scale_alpha_continuous(breaks = seq(0.05, 0.91, length.out = 9)) + xlab("Zero proportion") + ylab("Maximum Absolute Error") 
  
  pdf(file = paste0("Figures/", type, "_symmZR.pdf"), width=8, height=4)
  print(g)
  dev.off()
}





