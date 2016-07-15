DIC <-
function(obj, start=2,  stop=length(obj$theta$solu)){
      Y <- obj$Y
      X <- obj$X
      if(obj$scale){
        X <- scale(X)
        Y <- scale(Y)
      }else{
        X <- scale(X, scale=FALSE)
        Y <- scale(Y, scale=FALSE)    
      }
      pds <- .pd1(Y, X, obj, start,  stop)
      dic <- pds$Dbar + 2*pds$pd
      list(dic=dic, pd=pds$pd)
    }
