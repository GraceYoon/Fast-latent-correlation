### log10 dim vs log10 computation time
### to compare ORG, MLBD and Kendall only
### Figure 4 in manuscript.

rm(list=ls())

library(ggplot2) # for figure
library(tidyverse) # for %>%

###################################################
# to extract ggplot default colors.
# n is the number of colors you need.
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
###################################################

# for different dimension size.
plist <- c(2, 5, 10, 20, 50, 100, 200, 300, 400, 481)

# for each dimension.
for (ndim in c(100, 6482)){
  load(paste0("Data/RunTimePlot_v2_range_", ndim, ".Rda"))
  
  # 1 min. 2 Q1. 3 mean. 4 median. 4 Q3. 5 max.
  df_tmp <- data.frame(timemed = c(time_K[, 4], time_org[, 4], time_v2[, 4]), timemin = c(time_K[, 1], time_org[, 1], time_v2[, 1]), timemax = c(time_K[, 6], time_org[, 6], time_v2[, 6]), timemean = c(time_K[, 3], time_org[, 3], time_v2[, 3]), method = rep(c("Kendall", "ORG", "ML"), rep(length(plist), 3)), dim = rep(plist, 3))
  
  #df_tmp$dim <- as.factor(df_tmp$dim)
  assign(paste0("dftime_", ndim), df_tmp)
}

rm(df_tmp)


# combine two results with different sample size
dfall <- rbind.data.frame(data.frame(dftime_100, sampsize = 100),
                          data.frame(dftime_6482, sampsize = 6482))

sampsize.labs <- c("n = 100", "n = 6482")
names(sampsize.labs) <- c(100, 6482)
dfall$sampsize <- as.factor(dfall$sampsize)
dfall$method <- factor(dfall$method, levels = c("ORG", "ML", "Kendall"))





# to dodge position slightly for each method and sample size
### -0.02 for ORG 100, ML 100, ML 6482
### 0.02 for ORG 6482, KL100, KL6482
dfall$pd <- rep(-0.01, nrow(dfall))

dfall$pd[which(dfall$method == "Kendall" & dfall$sampsize == 100)] <- 0.01
dfall$pd[which(dfall$method == "Kendall" & dfall$sampsize == 6482)] <- 0.01

dfall$pd[which(dfall$method == "ORG" & dfall$sampsize == 100)] <- 0.01
dfall$pd[which(dfall$method == "ORG" & dfall$sampsize == 6482)] <- -0.01

dfall$pd[dfall$method == "ORG" & dfall$sampsize == 100] # -0.01
dfall$pd[dfall$method == "ORG" & dfall$sampsize == 6482] # 0.01
dfall$pd[dfall$method == "ML" & dfall$sampsize == 100] # -0.01
dfall$pd[dfall$method == "ML" & dfall$sampsize == 6482] # -0.01
dfall$pd[dfall$method == "Kendall" & dfall$sampsize == 100] # 0.01
dfall$pd[dfall$method == "Kendall" & dfall$sampsize == 6482] # 0.01




pdf(file = "Figures/RunTimePlot_v2_loglog.pdf", height = 6, width = 8)

dfall %>% filter(dfall$dim > 10) %>% # smaller dimension less than 20 is omitted.
  ggplot(aes(x = log10(dim) + 1.5*pd, y = log10(timemed), fill = method, group = interaction(method, sampsize))) + ylab("log10 of Computation Time (in seconds)") + geom_errorbar(aes(ymin = log10(timemin), ymax = log10(timemax), colour = method), width = 0.05) + geom_line(aes(colour = method, linetype = sampsize)) + geom_point(aes(colour = method)) + scale_colour_manual("Method", values = gg_color_hue(3), labels = c("ORG", "MLBD", "Kendall")) + scale_fill_manual("Method", values = gg_color_hue(3), labels = c("ORG", "MLBD", "Kendall")) + scale_linetype_discrete("Sample Size") + scale_x_continuous(name = "log10 of Dimension", breaks = log10(plist[-c(1:3)]), labels = paste0("log10 (", plist[-c(1:3)], ")")) + theme(text = element_text(size = 18), axis.text.x = element_text(angle = 45, vjust = 1.1, hjust = 1))

dev.off()

