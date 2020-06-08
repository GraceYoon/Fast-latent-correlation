# 3D graph of bridge inverse function for TC case
# Used in Figure 1 in the manuscript.

# x-axis: Kendall's tau (tau)
# y-axis: Zero proportion (z1)
# z-axis: Rank-based correlation = BridgeTCinverse(tau, z1)

rm(list=ls())


# load the data for TC case example
load("Data/tcInv_equigrid.Rda") # precomputed results. contained gridTCinv and z1.
# there are precomputed values in list for 199 tau values for each 100 zero proportion value.

# change the list format of the precomputed results into matrix format.
# now each row represents tau values and each column represents zero proportion level.
TCvalue <- matrix(unlist(gridTCinv), ncol = length(z1), byrow = FALSE)


# using plotly package, we can create 3d figure.
library(plotly)

plot_ly(z = t(TCvalue), colorbar = list(tickmode = "array", tickvals = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1), x = 0.82, y = 0.75)) %>% # tickvals and x and y are for legend ticks and location.
  add_surface(opacity = 0.9) %>% # for 3d figure, we need to add surface.
  layout(scene = list(xaxis = list(tickfont = list(size = 13), title = list(text = "Kendall's tau", font = list(size = 15)), tickmode = "array", nticks = 5, tickvals = c(0, 50, 100, 150, 200), ticktext = c("-1", "-0.5", "0", "0.5", "1")), # for x-axis. (tau)
                      yaxis = list(tickfont = list(size = 13), title = list(text = "Zero Proportion", font = list(size = 15)), tickmode = "array", nticks = 5, tickvals = c(1, 21, 41, 61, 81), ticktext = c("0", "0.2", "0.4", "0.6", "0.8")), # for y-axis. (zero proportion)
                      zaxis = list(tickfont = list(size = 13), title = list(text = "Latent correlation", font = list(size = 15))), legend = list(title = "correlation") # for z-axis (latent correlation, r)
))
