plot.BayesPLS <- function(obj, start=NULL, stop=NULL, thin=1, which=c(1L:4L), ask=TRUE){
  if(is.null(start)|is.null(stop)){
    start <- 1
    stop <- length(obj$theta$solu)
  }
  plotiter <- seq(start, stop, by=thin)
  dev.new(noRStudioGD = TRUE, height = 10, width = 15) 
  def.par <- par(no.readonly = TRUE)
  # if(!ask){
  #   layout(matrix(1:4,2,2, byrow=TRUE))
  # }
  # par(mar=c(4,4,5,1))
  # if (ask) {
  #   oask <- devAskNewPage(TRUE)
  #   on.exit(devAskNewPage(oask))
  # }
  if (!is.numeric(which) || any(which < 1) || any(which > 4)) 
    stop("'which' must be in 1:4")
  show <- rep(FALSE, 4)
  show[which] <- TRUE
  if(!ask){
    layout(matrix(c(1:4),2,2, byrow=TRUE))
  }
  par(mar=c(5.1, 4.1, 4.1, 4.1))
  # if (ask) {
  #   oask <- devAskNewPage(TRUE)
  #   on.exit(devAskNewPage(oask))
  # }
  
  if(show[1L]){
    dev.hold()
    plot(obj$theta$solu[plotiter],type="l",xlab="Iteration number",
       ylab="Value",main=expression(sigma^2))
    dev.flush()
  }
  if(show[2L]){
    dev.hold()
    matplot(obj$gamma$solu[plotiter,],type="l",lty=1,xlab="Iteration number",
            ylab="Value",main=expression(gamma))
    dev.flush()
  }
  if(show[3L]){
    dev.hold()
    matplot(obj$nu$solu[plotiter,],type="l",lty=1,xlab="Iteration number",
            ylab="Value",main=expression(nu))
    dev.flush()
  }
  if(show[4L]){
    dev.hold()
    matplot(obj$betas[plotiter,],type="l",lty=1,xlab="Iteration number",
           ylab="Value",main=expression(beta))
  dev.flush()
  }
}
