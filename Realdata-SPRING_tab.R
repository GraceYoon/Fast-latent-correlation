# check real data results
# to check accuracy & computation time of multilinear approximation

rm(list = ls())

# devtools::install_github("irinagain/mixedCCA", ref = "TBcutoff", force = TRUE)
# devtools::install_github("GraceYoon/SPRING", force = TRUE)

library(microbenchmark)
library(SPRING) # for QMP data
library(mixedCCA)

# set up to save results
CompRes <- matrix(NA, nrow = 2, ncol = 2)
rownames(CompRes) <- c("ORG", "MLBD")
colnames(CompRes) <- c("SPRING-dd", "SPRING-df") # dd: data-driven, df: default

AccRes <- matrix(NA, nrow = 2, ncol = 2)
rownames(AccRes) <- c("Max AE", "Mean AE")
colnames(AccRes) <- c("SPRING-dd", "SPRING-df") # dd: data-driven, df: default

#################################################
###### for sparse graphical network estimation on QMP data
###### need to find estimated partial correlation coefficient based on estimated optimal sparsity level.

###### With lambda sequence generated from estimated latent correlation
n = dim(QMP)[1]
p = dim(QMP)[2]
rep.num = 50 # the repetition number of subsampling
nlam = 50 # the number of lambda sequence values
lambda.min.ratio = 1e-2 # the ratio of lambda.min over lambda.max
thresh = 0.1 # threshold for StARS criterion
subsample.ratio = 0.8 # subsample size ratio over the total samples
nc = 2 # numbere of cores for subsampling in parallel mode
seed = 10010 # seed for subsampling


colLoc <- 1 # Data-specific lambda sequence result will be saved in 1st column

ptm <- proc.time()
fit.spring <- SPRING(QMP, quantitative = TRUE, Rmethod = "original", lambdaseq = "data-specific", lambda.min.ratio = lambda.min.ratio, nlambda = nlam, thresh = thresh, subsample.ratio = subsample.ratio, seed = seed, ncores = nc, rep.num = rep.num)
CompRes[1, colLoc] <- (proc.time() - ptm)[3]

ptm <- proc.time()
fit.spring_fast <- SPRING(QMP, quantitative = TRUE, Rmethod = "approx", lambdaseq = "data-specific", lambda.min.ratio = lambda.min.ratio, nlambda = nlam, thresh = thresh, subsample.ratio = subsample.ratio, seed = seed, ncores = nc, rep.num = rep.num)
CompRes[2, colLoc] <- (proc.time() - ptm)[3]



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
AccRes[1, colLoc] <- max(abs(Coef.K_org-Coef.K_fast))

##### mean absolute error 
AccRes[2, colLoc] <- mean(abs(Coef.K_org-Coef.K_fast))





###### With a default lambda sequence from 0.6 to 0.006 of length 50.
colLoc <- 2 # default lambda sequence result will be saved in 2nd column

ptm <- proc.time()
fit.spring_deflam <- SPRING(QMP, quantitative = TRUE, Rmethod = "original", lambda.min.ratio = lambda.min.ratio, nlambda = nlam, thresh = thresh, subsample.ratio = subsample.ratio, seed = seed, ncores = nc, rep.num = rep.num)
CompRes[1, colLoc] <- (proc.time() - ptm)[3]

ptm <- proc.time()
fit.spring_fast_deflam <- SPRING(QMP, quantitative = TRUE, Rmethod = "approx", lambda.min.ratio = lambda.min.ratio, nlambda = nlam, thresh = thresh, subsample.ratio = subsample.ratio, seed = seed, ncores = nc, rep.num = rep.num)
CompRes[2, colLoc] <- (proc.time() - ptm)[3]

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
AccRes[1, colLoc] <- max(abs(Coef.K_org-Coef.K_fast))

##### mean absolute error 
AccRes[2, colLoc] <- mean(abs(Coef.K_org-Coef.K_fast))




save(fit.spring_deflam, fit.spring, fit.spring_fast_deflam, fit.spring_fast, CompRes, AccRes, file = "Data/Realdata-SPRING_tab.Rda")