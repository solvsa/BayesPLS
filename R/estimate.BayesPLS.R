estimate.BayesPLS <- function(obj, start=1, stop=NULL, thin=1, probs = c(0.025, 0.975)){
  if(is.null(stop)){
    stop <- length(obj$theta$solu)
  }
  useiter <- seq(start, stop, by=thin)
  if(is.null(useiter)){stop("There are no values to estimate from\n")}
  
  betas <- obj$betas
  theta <- obj$theta$solu
  if(obj$scale){ 
    sdY <- attr(obj$Y, "scaled:scale")
    sdX <- attr(obj$X, "scaled:scale")
    betas <- betas*sdY/sdX
    theta <- sdY^2*theta
  }
  meanY <- attr(obj$Y, "scaled:center")
  meanX <- attr(obj$X, "scaled:center")
  
  betahat <- apply(betas[useiter,],2,mean)
  betaquants <- apply(betas[useiter,],2, quantile, probs = probs)
  beta0s <- meanY - betas%*%meanX
  beta0 <- mean(beta0s[useiter])
  beta0quants <- quantile(beta0s[useiter], probs = probs)
  sigma.sq <- mean(theta)
  sigma.sq.quants <- quantile(theta, probs=probs)
  
  res <- list()
  res$coefficients <- betahat
  res$intercept <- beta0
  res$sigma.sq <- sigma.sq
  res$quantiles <- t(betaquants)
  res$quantiles.int <- beta0quants
  res$quantiles.sigma.sq <- sigma.sq.quants
  class(res) <- "BayesPLS"
  res
}