predict.BayesPLS <-
function(obj, newX, type="mean"){
  fpred <- function(x, newX){
    newX <- cbind(1, newX)
    newX%*%x
  }
  if(type=="mean"){
    drop(obj$intercept) + newX%*%obj$coef
  }else if(type=="dist"){
    betamat <- cbind(obj$trace.beta0, obj$trace.beta)
    predmat <- apply(betamat, 1, fpred, newX=newX)
    apply(predmat,1,mean)
  }
}
