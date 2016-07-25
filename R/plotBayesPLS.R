plotBayesPLS <- function(obj, start=NULL, stop=NULL, thin=1){
  if(is.null(start)|is.null(stop)){
    start <- 1
    stop <- length(obj$theta$solu)
  }
  plotiter <- seq(start, stop, by=thin)
  dev.new(noRStudioGD = TRUE, height = 10, width = 15) 
  layout(matrix(1:4,2,2, byrow=TRUE))
  par(mar=c(4,4,5,1))  
  plot(obj$theta$solu[plotiter],type="l",xlab="Iteration number",
       ylab="Value",main=expression(sigma^2))
  matplot(obj$gamma$solu[plotiter,],type="l",lty=1,xlab="Iteration number",
          ylab="Value",main=expression(gamma))
  matplot(obj$nu$solu[plotiter,],type="l",lty=1,xlab="Iteration number",
          ylab="Value",main=expression(nu))
  matplot(obj$betas[plotiter,],type="l",lty=1,xlab="Iteration number",
          ylab="Value",main=expression(beta))
}
