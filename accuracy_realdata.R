# to check accuracy of multilinear approximation

rm(list = ls())

#################################################
####### for QMP data
load("Data/Kcor_QMP_org.Rda")
load("Data/Kcor_QMP_fast.Rda")

max(abs(c(Kcor_QMP_fast - Kcor_org)))
mean(abs(c(Kcor_QMP_fast - Kcor_org)))
# > max(abs(c(Kcor_QMP_fast - Kcor_QMP_org)))
# [1] 0.0006314209
# > mean(abs(c(Kcor_QMP_fast - Kcor_QMP_org)))
# [1] 8.354191e-05

#################################################
###### for AGP data
load("Data/Kcor_amgut_org.Rda")
load("Data/Kcor_amgut_fast.Rda")

max(abs(c(Kcor_amgut_fast - Kcor_amgut_org)))
mean(abs(c(Kcor_amgut_fast - Kcor_amgut_org)))
# > max(abs(c(Kcor_amgut_fast - Kcor_amgut_org)))
# [1] 0.000573721
# > mean(abs(c(Kcor_amgut_fast - Kcor_amgut_org)))
# [1] 7.99355e-05

#################################################
###### for TCGA-BRCA data
load("Data/Kcor_TCGA_fast.Rda")
load("Data/Kcor_TCGA_org.Rda")

max(abs(c(Kcor_TCGA_fast - Kcor_TCGA_org)))
mean(abs(c(Kcor_TCGA_fast - Kcor_TCGA_org)))
# > max(abs(c(Kcor_TCGA_fast - Kcor_TCGA_org)))
# [1] 0.000487254
# > mean(abs(c(Kcor_TCGA_fast - Kcor_TCGA_org)))
# [1] 1.16256e-05


#################################################
###### for sparse graphical network estimation on QMP data
###### need to find estimated partial correlation coefficient based on estimated optimal sparsity level.

###### With lambda sequence generated from estimated latent correlation
load("Data/QMP_SPRING_org_ddlam.rdata")
load("Data/QMP_SPRING_fast_ddlam.rdata")

# For StARS criterion
thresh = 0.1
# According to the threshold, find opt.index and save adjacency matrix and coefficient among the path.
opt.K_fast <- max(which(fit.spring_fast$output$stars$summary<thresh))
opt.K_org <- max(which(fit.spring$output$stars$summary<thresh))

adj.K_fast <- fit.spring_fast$fit$est$path[[opt.K_fast]]; adj.K2 <- as.matrix(adj.K_fast)
adj.K_org <- fit.spring$fit$est$path[[opt.K_org]]; adj.K <- as.matrix(adj.K_org)

Coef.K_fast <- SpiecEasi::symBeta(fit.spring_fast$output$est$beta[[opt.K_fast]], mode='maxabs') # isSymmetric(Coef.K2) # TRUE
Coef.K_fast <- as.matrix(Coef.K_fast)
Coef.K_org <- SpiecEasi::symBeta(fit.spring$output$est$beta[[opt.K_org]], mode='maxabs') # isSymmetric(Coef.K) # TRUE
Coef.K_org <- as.matrix(Coef.K_org)

##### maximum absolute error 
max(abs(Coef.K_org-Coef.K_fast))
# [1] 0.00133607

##### mean absolute error 
mean(abs(Coef.K_org-Coef.K_fast))
# [1] 9.571935e-06


###### With a default lambda sequence from 0.6 to 0.006 of length 50.
load("Data/QMP_SPRING_org_defaultlam.rdata")
load("Data/QMP_SPRING_fast_defaultlam.rdata")

# For StARS criterion
thresh = 0.1
# According to the threshold, find opt.index and save adjacency matrix and coefficient among the path.
opt.K_fast <- max(which(fit.spring_fast_deflam$output$stars$summary<thresh))
opt.K_org <- max(which(fit.spring_deflam$output$stars$summary<thresh))

adj.K_fast <- fit.spring_fast_deflam$fit$est$path[[opt.K_fast]]; adj.K2 <- as.matrix(adj.K_fast)
adj.K_org <- fit.spring_deflam$fit$est$path[[opt.K_org]]; adj.K <- as.matrix(adj.K_org)

Coef.K_fast <- SpiecEasi::symBeta(fit.spring_fast_deflam$output$est$beta[[opt.K_fast]], mode='maxabs') # isSymmetric(Coef.K2) # TRUE
Coef.K_fast <- as.matrix(Coef.K_fast)
Coef.K_org <- SpiecEasi::symBeta(fit.spring_deflam$output$est$beta[[opt.K_org]], mode='maxabs') # isSymmetric(Coef.K) # TRUE
Coef.K_org <- as.matrix(Coef.K_org)

##### maximum absolute error 
max(abs(Coef.K_org-Coef.K_fast))
# [1] 0.001073247

##### mean absolute error 
mean(abs(Coef.K_org-Coef.K_fast))
# [1] 7.96922e-06



