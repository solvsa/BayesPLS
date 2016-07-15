predict.BayesPLS <-
function(obj, newX){
  drop(obj$intercept) + newX%*%obj$coef
}
