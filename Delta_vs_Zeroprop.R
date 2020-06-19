# draw a plot zero proportion versus Delta

rm(list = ls())

library(ggplot2)

zp <- seq(0.001, 0.999, by = 0.001)
Delta <- qnorm(zp)
df <- data.frame(ZeroProportion = c(0, zp, 1), Delta = c(-max(abs(Delta)), Delta, max(abs(Delta))))

pdf(file = "Figures/Delta_zeroprop.pdf", height = 6, width = 6)

ggplot(df, aes(x = ZeroProportion, y = Delta)) + geom_path() + xlab(expression(paste("Zero proportion ",pi[0]))) + ylab(expression(Delta)) + theme(plot.margin = margin(t = 0.5, r = 1, b = 0.5, l = 0.5, unit = "cm"), text = element_text(size = 18), axis.title.y = element_text(angle = 0, vjust = 0.5), axis.text.x = element_text(angle = 45)) + scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) + scale_y_continuous(breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  geom_segment(aes(x = 0.1, y = -max(abs(Delta)), xend = 0.1, yend = qnorm(0.1)), colour = "red", linetype = "dotted") + geom_segment(aes(x = 0, y = qnorm(0.1), xend = 0.1, yend = qnorm(0.1)), colour = "red", linetype = "dotted") +
  geom_segment(aes(x = 0.99, y = -max(abs(Delta)), xend = 0.99, yend = qnorm(0.99)), colour = "red", linetype = "dotted") + geom_segment(aes(x = 0, y = qnorm(0.99), xend = 0.99, yend = qnorm(0.99)), colour = "red", linetype = "dotted") +
  geom_segment(aes(x = 0.95, y = -max(abs(Delta)), xend = 0.95, yend = qnorm(0.95)), colour = "red", linetype = "dotted") + geom_segment(aes(x = 0, y = qnorm(0.95), xend = 0.95, yend = qnorm(0.95)), colour = "red", linetype = "dotted") + 
  # add texts to x-axis ## needed to adjust to avoid overlaps
  geom_text(x = 0.93, y = -3.4, hjust = 0, size = 5, label = "0.99", angle = 45) + 
  geom_text(x = 0.88, y = -3.4, hjust = 0, size = 5, label = "0.95", angle = 45) + 
  geom_text(x = 0.07, y = -3.37, hjust = 0, size = 5, label = "0.1", angle = 45) + 
  # add texts to y-axis
  geom_text(x = -0.1, y = qnorm(0.1), hjust = 0, size = 5, label = as.character(round(qnorm(0.1), 2)), angle = 0) + 
  geom_text(x = -0.08, y = qnorm(0.95), hjust = 0, size = 5, label = as.character(round(qnorm(0.95), 2)), angle = 0) + 
  geom_text(x = -0.08, y = qnorm(0.99), hjust = 0, size = 5, label = as.character(round(qnorm(0.99), 2)), angle = 0) + coord_cartesian(xlim = c(0, 1), ylim = c(-3.09, 3.09), clip = "off", expand = FALSE)

dev.off()


