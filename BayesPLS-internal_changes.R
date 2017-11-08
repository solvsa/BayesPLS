.adjust <-
function(x,minrate=0.25,maxrate=0.5,stepfactor=1.2, doadjust=TRUE, adaptint)
    {
        x$rate <- x$accept/adaptint
        if(doadjust & (x$eps>1e-16)){
            if(x$rate<minrate){x$eps<-x$eps/stepfactor}
            if(x$rate>maxrate){x$eps<-x$eps*stepfactor}
        }
        return(list(x$eps,0,x$rate))
    }
.adjust.multi <-
function(x,minrate=0.025,maxrate=0.3,stepfactor=1.2, doadjust=TRUE, adaptint)
    {
        n.adjust <- length(x$rate)
        for(iii in 1:n.adjust){
            x$rate[iii] <- x$accept[iii]/adaptint
            if(doadjust){
                if(x$rate[iii]<minrate){x$eps[iii]<-x$eps[iii]/stepfactor}
                if(x$rate[iii]>maxrate){x$eps[iii]<-x$eps[iii]*stepfactor}
            }
        }
        return(list(eps=x$eps,accept=rep(0,n.adjust),rate=x$rate))
    }
.apost <-
function(Y, X, a, ncomp, eps)
        {
            .loglikey(Y, X, a$theta, a$gamma, a$dvek, ncomp) + 
            .loglikex(X,a$nu, a$dvek) + 
            .logminv1(a$nu, eps$nueps) + 
            .logminv2(abs(a$gamma), eps$gammaeps) +                      
            .logminv1(a$theta, eps$thetaeps)
        }
.cond.mom <-
function(a, mu, sig, subs){
        if(length(subs)==length(a$gamma)){
            mu.cond <- mu
            sig.cond <- sig
        }else{ 
            ss <- tcrossprod(sig[subs,-subs,drop=F],solve(sig[-subs,-subs,drop=F]))
            mu.cond <- mu[subs, drop=F] + ss%*%(a$gamma[-subs, drop=F]-mu[-subs])
            sig.cond <- sig[subs,subs,drop=F] - ss%*%sig[-subs, subs, drop=F]            
        }
        return(list(mu=mu.cond, sig=sig.cond))
    }
.Dbar <-
function(Y, X, obj, start=2, stop=length(obj$theta$solu)){
      ncomp <- dim(obj$gamma$solu)[2]
      niter <- stop - start + 1
      Ds <- rep(0, niter)
      for(it in 1:niter){
        i <- it + start - 1
        Ds[it] <- -2*(.loglikey(Y,X,obj$theta$solu[i], obj$gamma$solu[i,], obj$D$solu[i,,], ncomp)) 
        #- 2*(.loglikex(X, obj$nu$solu[i,,drop=FALSE],obj$D$solu[i,,]))
      }
      list(Dbar = mean(Ds), Ds=Ds)
    }
.dgamma.propos <-
function(x, a, b, eps){
    res <- rep(0,length(x))
    id <- which(x>a & x<b)
    res[id] <- eps/((-a)^eps + (b)^eps)*1/(abs(x[id])^(1-eps))    
    return(res)
}
.Dmeans <-
function(Y,X, obj, start=2, stop=length(obj$theta$solu)){
      ncomp <- dim(obj$gamma$solu)[2]
      theta <- mean(obj$theta$solu[start:stop])
      gamma <- apply(obj$gamma$solu[start:stop,,drop=FALSE],2,mean)
      nu <- apply(obj$nu$solu[start:stop,,drop=FALSE],2,mean)
      dvek <- apply(obj$D$solu[start:stop,,],c(2,3),mean)
      Dm <- -2*(.loglikey(Y, X, theta, gamma, dvek, ncomp)) 
      #- 2*(.loglikex(X, nu, dvek))
      Dm
    }
.exact.gamma <-
function(a, mus, sigs, eps, nz = 100){

        zs <- rmvnorm(nz, mus, diag(sigs))
        f <- function(x, eps){1/((abs(x)^(1-eps)))}
        const <- apply(f(zs,eps$gammaeps),2,mean)
        const <- 1/(sqrt(sigs*2*pi)*const)
        const*1/abs(a$gamma)^(1-eps$gammaeps)*exp(-1/(2*sigs)*(a$gamma-mus)^2)
    }
.fgammainv <-
function(u,a,b,eps){
    x <- -( (-a)^eps - u*((-a)^eps + b^eps) )^(1/eps)
    return(x)
}
.find.rho <-
function(Y, a, A, AAinv, eps, subsets=as.list(1:length(a$gamma))){
        if(length(a$gamma)==1){rho <- 1
        }else{
            moms <- .moments(a, A, AAinv, Y)
          
            mu.cond <- rep(0, length(subsets))
            sig.cond <- rep(0, length(subsets))
            nom <- rep(0, length(subsets))            
            #Evaluating conditional normal for gammas        
            for(i in 1:length(subsets)){
                id <- subsets[[i]]
                condmoms <- .cond.mom(a,moms$mu.gamma, a$theta*moms$sig.gamma, id)
                mu.cond[i] <- condmoms$mu
                sig.cond[i] <- condmoms$sig
                nom[i] <- dnorm(a$gamma[id], mu.cond[i], sqrt(sig.cond[i]))
            }
            denom <- .exact.gamma(a, mu.cond, sig.cond, eps)
            rho <- exp(-eps$lambda*nom/denom)
            }
        return(rho)
    }
.householder <-
function(x){
            m <- length(x)
            alpha <- sqrt(drop(crossprod(x)))
            e <- c(1,rep(0,m-1))
            u <- x - alpha*e
            v <- u/sqrt(drop(crossprod(u)))
            diag(m) - 2*v%*%t(v)
    }
.initiate <-
function(Y,X,ncomp, init.method, reorder.nu=TRUE){
        
        if(init.method == "PCR"){
            last <- .initiate.pcr(Y,X,ncomp)
        }else if(init.method == "PLS"){
            last <- .initiate.pls(Y,X,ncomp)
        }
        
        if(reorder.nu){   
            #Reordering, largest nu first
            index <- rev(order(abs(last$nu)))
            id <- which(index<=ncomp)
            last$gamma[1:ncomp] <- last$gamma[index[id]]
            last$nu <- last$nu[index]
            last$dvek <- last$dvek[,index]    
        }
        last
    }
.initiate.pcr <-
function(Y,X,ncomp){
        last.init <- list()
        p <- dim(X)[2]
        n <- dim(X)[1]
        last.init$dvek <- matrix(0,p,p)    
        S <- eigen(cov(X))
        alpha <- cor((X%*%S$vectors),Y)
        D <- sweep(S$vectors,2,sign(alpha),"*")
        Z <- X%*%D[,1:ncomp]
        beta.pcr <- D[,1:ncomp]%*%solve(t(Z)%*%(Z))%*%t(Z)%*%Y
        last.init$theta <- drop(crossprod(Y-X%*%beta.pcr))/(n-ncomp)
        last.init$gamma <- t(D[,1:ncomp])%*%beta.pcr
        last.init$nu <- S$values + 1.0e-8*S$values[1]
        last.init$dvek <- D
        last.init
    }
.initiate.pls <-
function(Y,X,ncomp){
        last.init <- list()
        p <- dim(X)[2]
        n <- dim(X)[1]
        last.init$dvek <- matrix(0,p,p) 
        S <- cov(X)
        maxcomp <- min(n-1,p-1)
        plsfit <- plsr(Y~X, ncomp=maxcomp, method="oscorespls")
        last.init$theta <- drop(crossprod(plsfit$resid[,,ncomp]))/(n-ncomp)
        if(maxcomp < (p-1)){
            R <- matrix(rnorm((p-maxcomp)*p,0,1),ncol=(p-maxcomp))
            D <- plsfit$loading.weights
            Proj <- diag(p)-D%*%solve(t(D)%*%D)%*%t(D)
            R <- Proj%*%R
            R <- orthonormalization(R)
            D <- cbind(D,R[,1:(p-maxcomp)])
        }else{
            plsfit <- plsr(Y~X, ncomp=maxcomp+1, method="oscorespls")
            D <- plsfit$loading.weights
        }
        last.init$dvek <- D
        last.init$gamma <- t(plsfit$coef[,,ncomp]%*%plsfit$loading.weights[,1:ncomp])
        nus <- diag(t(D)%*%S%*%D)
        last.init$nu <- nus
        #last.init$nu <- rep(1/sum(nus),p)
        last.init
    }
.kern <-
function(x, logSdet, Sigmainv, mu){
        -0.5*(logSdet + crossprod((x-mu),crossprod(Sigmainv,(x-mu))))
    }
.loginv <-
function(x){-log(x)}
.loglikex <-
function(X, nu, dvek){
            likex <- .logmvdnorm2(X, lambdas=nu, eigenvecs=dvek)        
            return(likex)
        }
.loglikey <-
function(Y, X, theta, gamma, dvek, ncomp)    
        {
            betavek <- tcrossprod(gamma, dvek[,1:ncomp,drop=FALSE])
            Yhat <- t(tcrossprod(betavek, X))
            likey <- .logmvdnorm(Y,mu=Yhat,sigma=sqrt(theta))
            return(likey)
        }
.logminv1 <-
function(xvec,eps){-sum(log(xvec)+eps/xvec)}
.logminv2 <-
function(xvec,eps){-(1-eps)*sum(log(xvec))}
.logmvdnorm <-
function(vars,mu=rep(0,length(vars)),sigma=1)
    {
        sigma2 <- sigma^2
        -0.5*(length(vars)*log(2*pi*sigma2)+(sigma2^(-1))*c(crossprod(vars-mu)))
    }
.logmvdnorm2 <-
function(vars, lambdas, eigenvecs, mu=rep(0,length(lambdas)))
    {
        p <- length(lambdas)
        Sigmainv <- tcrossprod(eigenvecs*rep(sqrt(1/lambdas),each=dim(eigenvecs)[1]))
        logSdet <- sum(log(lambdas))
        oo <- apply(vars,1,.kern, logSdet, Sigmainv, mu) 
        sum(oo)
  }
.ldinvgamma <- 
  function (x, shape, scale){
    n <- length(x)
    alpha <- shape
    beta <- scale
    log.density <- alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(x) - (beta/x)
  }
.logprop <-
function(Y,X,a,b=NULL,pars=c("dvek","theta","nu","gamma"), 
                        rho=NULL, delta=NULL, mom=NULL, eps, r.gamma, updates){
      n <- dim(X)[1]
      p <- dim(X)[2]
      D <- a$dvek
        m <- length(a$gamma)
        A <- tcrossprod(X,t(D[,1:m,drop=FALSE]))
        AAinv <- solve(crossprod(A))
        H1 <- tcrossprod(AAinv,A)
        H <- A%*%H1
        logprop.dvek <- logprop.theta <- logprop.nu <- logprop.gamma <- logprop.r <- 0
        if(is.null(mom)){
            rho <- .find.rho(Y, a, A, AAinv, eps)
            mom <- .moments(a, A, AAinv, Y)
            if(length(updates>0)){
                diffg <- b$gamma[updates] - mom$mu.gamma[updates]
                delta1 <- crossprod(diffg,solve(mom$sig.gamma[updates,updates]))   
                delta <- crossprod(t(delta1),diffg)
            }else{
                delta <- 0
            }        
        }
        if(c("rho")%in%pars){logprop.r <- sum(log(rho^(r.gamma)*(1-rho)^(1-r.gamma)))}

        if(c("theta")%in%pars){
            cpy <- crossprod(Y,(diag(n)-H))
            logprop.theta <- log(dinvgamma(a$theta, (n-length(updates))/2, (crossprod(t(cpy),Y) + delta)/2))
          }
        if(c("nu")%in%pars){
            invscale <- apply(D,2,.nuscale2, A=X, eps=eps$nueps)
            logprop.nu <- sum(.ldinvgamma(a$nu, n/2,invscale))
        }
        if(c("gamma")%in%pars && length(updates)>0){
            cond <- .cond.mom(a, mom$mu.gamma, a$theta*mom$sig.gamma, subs=updates)
            logprop.gamma <- dmvnorm(c(a$gamma[updates]), mean=cond$mu, sigma=cond$sig, log=TRUE)
        }
        return(logprop.dvek + logprop.theta + logprop.nu + logprop.gamma + logprop.r)
    }
.metrop <-
function(Y, X, actual, last, logprop1, logprop2, eps)
    {
        u1<-runif(1,0,1)
        logpi1<-.apost(Y, X, actual, ncomp=length(actual$gamma), eps)
        logpi2<-.apost(Y, X, last, ncomp=length(last$gamma), eps)
        logsum <- logpi1-logpi2+logprop2-logprop1       
        acc <- exp(logsum)
        accprob<-min(1,acc)
        test<-ifelse(u1<=accprob,TRUE,FALSE)
        return(test)
    }
.moments <-
function(a, A, AAinv, Y){
        sig.gamma <- AAinv
        mu.gamma <- tcrossprod(AAinv,A)%*%Y
        return(list(mu.gamma=mu.gamma, sig.gamma=sig.gamma))
    }
.nuscale <- 
  function(x,A,eps){
  x <- matrix(x,ncol=1)
  n <- dim(x)[1]
  sc <- eps + 0.5*sum((A%*%x)^2)
  cand <- rinvgamma(1, n/2, sc)
    }
.nuscale2 <- 
  function(x,A,eps){
  x <- matrix(x,ncol=1)
  sc <- eps + 0.5*sum((A%*%x)^2)
}
.pd1 <-
function(Y, X, obj, start=2, stop=length(obj$theta$solu)){
      D1 <- .Dbar(Y, X, obj, start, stop) 
      D2 <- .Dmeans(Y, X, obj, start, stop)
      pd <- (D1$Dbar - D2)
      list(pd = pd, Dbar = D1$Dbar)
    }
.QR <-
function(X){

        n <- dim(X)[1]
        if(dim(X)[2]!=n)stop("Applies only to quadratic matrices \n")
    
    
        Ident <- diag(n)
        Q <- Ident
        R <- X
        Emark <- X
        for(j in 1:(n-1)){
            Emark <- R[j:n,j:n]
            Qmark <- .householder(Emark[,1])
            Q1 <- Ident
            Q1[j:n,j:n] <- Qmark
            R <- Q1%*%R    
            Q <- Q%*%Q1
        }

        return(list(Q=Q,R=R))
    }
