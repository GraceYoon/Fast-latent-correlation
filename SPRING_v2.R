##### updated SPRING code to use updated "estimateR" function in mixedCCA.
##### plan to update SPRING package in github in the near future. ("https://github.com/GraceYoon/SPRING/")

hugeKmb_v2 <- function(data, lambda, type = "trunc", sym = "or", verbose = TRUE, verboseR = TRUE, Rmethod = "approx", tol = 1e-6) {
  S    <- estimateR(data, type = type, method = Rmethod, tol = tol, verbose = verboseR)$R
  est  <- huge::huge.mb(S, lambda, sym = sym, verbose = verbose)
  est
}


SPRING_v2 <- function(data, quantitative = FALSE, method = "mb", lambda.min.ratio = 1e-2, nlambda = 20, lambdaseq = exp(seq(log(0.6), log(0.6*lambda.min.ratio), length.out = nlambda)), seed = 10010, ncores = 1, thresh = 0.1, subsample.ratio = 0.8, rep.num = 20, Rtol = 1e-6, verbose = TRUE, verboseR = FALSE, Rmethod = "approx"){
  
  if (any(data < 0)) {
    stop("Negative values are detected, but either quantitative or compositional counts are expected.\n")
  }
  p <- ncol(data)
  if (quantitative){
    if (max(rowSums(data)) <= 1 | isTRUE(all.equal(max(rowSums(data)), 1))){
      warning("The input data is normalized, but quantitative count data is expected.\n")
    }
    qdat <- data
  } else {
    qdat <- SPRING::mclr(data)
  }
  rm(data)
  gc()
  
  if(is.character(lambdaseq)){
    if(lambdaseq == "data-specific"){
      Kcor <- estimateR(qdat, type = "trunc", method = Rmethod, tol = Rtol, verbose = verboseR)$R
      # generate lambda sequence
      lambda.max <- max(max(Kcor-diag(p)), -min(Kcor-diag(p)))
      lambda.min <- lambda.min.ratio * lambda.max
      lambdaseq <- exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
    } else {
      stop("The input for lambdaseq is not correct.\n")
    }
  }
  
  if(method == "mb"){
    fun <- hugeKmb_v2
  }
  
  out1.K_count <- pulsar::pulsar(qdat, fun = fun, fargs = list(lambda = lambdaseq, tol = Rtol, verbose = verbose, verboseR = verboseR), rep.num = rep.num, criterion = 'stars', seed = seed, ncores = ncores, thresh = thresh, subsample.ratio = subsample.ratio)
  
  fit1.K_count <- pulsar::refit(out1.K_count)
  
  return(list(output = out1.K_count, fit = fit1.K_count, lambdaseq = lambdaseq))
}
