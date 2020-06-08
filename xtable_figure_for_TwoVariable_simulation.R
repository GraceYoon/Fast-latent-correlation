# to compare the computation time of latent correlation between two simulated variables


# to create a table of computation time.
rm(list = ls())

library(xtable)
options("scipen"=100, "digits"=4) # not to use expoential format

tseq <- NULL
for (type in c("TC", "TT", "BC", "BB", "TB")){
  load(paste0("Data/TwoSimVariable", type, ".Rda"))

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
# ORG & 5790.27 & 21315.62 & 2830.30 & 4507.15 & 6068.79 \\ 
# ML & 457.09 & 573.89 & 439.93 & 523.94 & 582.67 \\ 
# MLBD & 476.82 & 613.58 & 457.38 & 557.49 & 621.40 \\ 
# \hline
# \end{tabular}
# \end{table}




# to create a figure of maximum absolute error.

rm(list = ls())

library(ggplot2)
library(tidyverse)
library(RColorBrewer)


for (type in c("TC", "TT", "BC", "BB", "TB")){
  load(paste0("Data/TwoSimVariable", type, ".Rda"))
  g <- df_accuracy %>% ggplot(aes(x = TruncRate, y = MaxAD, group = interaction(as.factor(LatentR), method))) + geom_bar(stat = "identity", width = 0.09, position = "dodge", aes(fill = method, color = method, alpha = LatentR), size = 0.1) + scale_x_continuous(breaks = unique(df_accuracy$TruncRate))  + scale_fill_discrete(name = "Method", label = c("ML", "MLBD")) + guides(color = F) + scale_alpha_continuous(breaks = seq(0.05, 0.91, length.out = 9)) + xlab("Zero proportion") + ylab("Maximum Absolute Error") 

  pdf(file = paste0("Figures/", type, "_MaxAE_color.pdf"), width=8, height=4)
  print(g)
  dev.off()
}









