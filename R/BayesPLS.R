BayesPLS <-
function(Y, X, ncomp,
              scale = TRUE,
              totiter=1000,
              start = 1,
              doinit = TRUE,
              init.method = c("PLS"),
              dotrace = TRUE,
              plotint = 10,
              thin = 10,  
              adaptint = 100,
              approx = FALSE,
              appit = 10,
              appcomp = 100,
              approp = 0.999,
              update = list(
                  update.dvek  = TRUE,    
                  update.nu    = TRUE,
                  update.theta = TRUE,
                  update.gamma = TRUE,
                  update.seq.gamma = TRUE),
              eps=list(gammaeps = 1/19,
                  nueps = 1/31,
                  thetaeps = 1/31,
                  #dvekeps = 1,
                  lambda = -log(0.001),
                  fi = 0.5),
              deps = 0,
              previousobj = NULL 
              ){
        
        #Input check of eps elements
        if(any(unlist(lapply(eps[1:3],function(x){(1/x)%%2!=1})))){stop("The inverse of gammaeps, nueps and thetaeps must be odd numbers")}

        #For plotting and updating etc.
        if(dotrace){
          dev.new(noRStudioGD = TRUE, height = 10, width = 15) 
          layout(matrix(1:4,2,2, byrow=TRUE))
          par(mar=c(4,4,2,1))
        }
        
        
        plottime<-seq(plotint,totiter,by=plotint) 
        thinind<-seq(thin,totiter-1,by=thin)
        Tt <- length(thinind)+1 
        adapttime<-seq(adaptint,totiter,by=adaptint)    
        
        n <- dim(X)[1]
        p <- dim(X)[2]
        if(p>=500){
          warning("For large p approx=TRUE is recommended\n")
        }
        
        #Centering of matrices
        X.c <- scale(X, scale=scale)
        meanX <- attr(X.c, "scaled:center")
        Y.c <- scale(Y, scale=scale)
        meanY <- attr(Y.c, "scaled:center")
        if(scale){
          sdX <- attr(X.c, "scaled:scale")
          sdY <- attr(Y.c, "scaled:scale")
        }
        covX <- cov(X.c)
        
        
        #Creating solution objects
        dvekobj    <-list(solu=array(0,dim=c(Tt,p,p))   ,eps=0 ,accept=0   ,rate=0)
        thetaobj   <-list(solu=rep(0,Tt)                ,eps=0  ,accept=0   ,rate=0)
        gammaobj   <-list(solu=matrix(0,Tt,ncomp)       ,eps=rep(0,ncomp)  ,accept=rep(0,ncomp)   ,rate=rep(0,ncomp))
        nuobj      <-list(solu=matrix(0,Tt,p)           ,eps=0  ,accept=0   ,rate=0)
        betasolu   <-matrix(0,Tt,p)
        SSE        <-rep(0,Tt)

        if(doinit){ 
          last <- .initiate(Y.c, X.c, ncomp, init.method)
          
          #Initiating storage objects
          nuobj$solu[1,] <- last$nu; 
          gammaobj$solu[1,] <- last$gamma; 
          dvekobj$solu[1,,] <- last$dvek;  
          thetaobj$solu[1] <- last$theta; 
          #dvekobj$eps  <- eps$dvekeps
          if(identical(deps%%1,0)){
            dvekobj$eps <- 2^{deps}*exp(1.2 - 0.01*n - 0.12*p - 3.2*ncomp + 0.003*n*ncomp + 0.036*p*ncomp - 0.00007*n*p)
          }else{
            dvekobj$eps <- deps
          }
          #(if(approx){dvekobj$eps <- deps*exp(0.32 - 0.02*n - 2.5*ncomp + 0.005*n*ncomp)})
        }else{
          last <- list()
          if(is.null(previousobj))stop("No previous fitted object found\n")
          #Reading initial values from saved object
          nuobj$solu[1,] <- last$nu <- previousobj$last$nu; 
          gammaobj$solu[1,] <- last$gamma <- previousobj$last$gamma; 
          dvekobj$solu[1,,] <- last$dvek <- previousobj$last$dvek;  
          thetaobj$solu[1] <- last$theta <- previousobj$last$theta; 
          dvekobj$eps  <- previousobj$D$eps
        }
        
        if(approx){
          nuprop <- cumsum(last$nu)/sum(last$nu)
          bignu <- which(nuprop<approp)
          nnu <- length(bignu)
          usenus <- min(nnu, appcomp)
        }

        k <- 1
        dosave <- FALSE
        
        for(i in 2:totiter){
          
          if(is.element(i,thinind)){
            dosave <- TRUE
            k <- k+1
          }
          
          # 1) Block-update of D, theta, nu and some gammas-------------------------------
          
          actual <- last    
          
          #----   dvek   --------
          if(update$update.dvek){   
            if(!approx){
              #Rotating the d-vectors by a random rotation 
              #(random walk on unit sphere)
              E <- matrix(rnorm(p^2,0,dvekobj$eps),p,p)+diag(p)
              Rot <- .QR(E)$Q
              neg <- which(diag(Rot)<0)
              Rot[,neg] <- (-1)*Rot[,neg]
              cand <- crossprod(Rot,actual$dvek)#*sample(c(-1,1),1)
              actual$dvek <- cand
            }else{
              if(i<appit){
                use <- p
              }else{
                use <- usenus
              }
              #Rotating the d-vectors corresponding to the largest eigenvalues 
              #by a random rotation (random walk on unit sphere)
              Rot <- diag(p)
              E <- matrix(rnorm(use^2,0,dvekobj$eps),use,use)+diag(use)
              Rot2 <- .QR(E)$Q
              neg <- which(diag(Rot2)<0)
              Rot2[,neg] <- (-1)*Rot2[,neg]
              Rot[1:use, 1:use] <- Rot2
              #browser()
              cand <- crossprod(Rot,actual$dvek[,1:use])#*sample(c(-1,1),1)
              actual$dvek[,1:use] <- cand
            }
          }
          #Derived variables
          D <- actual$dvek
          A <- tcrossprod(X.c,t(D[,1:ncomp,drop=FALSE]))
          AAinv <- solve(crossprod(A)) 
          H1 <- tcrossprod(AAinv,A)
          H <- crossprod(t(A),H1)
          
          
          
          #Individual probabilities for the gammas for joining a block update
          rho.gamma <- .find.rho(Y.c, actual, A=A, AAinv=AAinv, eps)
          u.gamma <- runif(ncomp,0,1)
          r.gamma <- as.numeric(u.gamma<rho.gamma)
          updates <- which(r.gamma==0)
          #Compute the contribution from gamma to the posterior of theta
          moms <- .moments(actual, A, AAinv, Y.c)
          if(length(updates>0)){
            diffg <- actual$gamma[updates] - moms$mu.gamma[updates]
            delta.rho1 <- crossprod(diffg,solve(moms$sig.gamma[updates,updates])) 
            delta.rho <- crossprod(t(delta.rho1),diffg)
          }else{
            delta.rho <- 0
          }
          
          
          #----   theta  -------
          if(update$update.theta){
            cp <- crossprod(Y.c,diag(n)-H)
            cand <- rinvgamma(1, (n-length(updates))/2, eps$thetaeps + (crossprod(t(cp),Y.c) + delta.rho)/2)
            actual$theta <- cand        
          }
          
          
          
          #----  nu   ---------
          if(update$update.nu){     
            cand <- actual$nu
            if(!approx){
              cand <- apply(D,2,.nuscale,A=X.c, eps=eps$nueps)
            }else{
              if(i<appit){
                use <- p
              }else{
                use <- usenus
              }
              cand[1:use] <- apply(D[,1:use],2,.nuscale,A=X.c, eps=eps$nueps)
            }
            actual$nu <- cand
          }
          

          #----   gamma   ------
          #Sample here a block update candidate for gamma_{\rho} using a conditional normal distribution
          if(update$update.gamma){
            if(length(updates)>0){
              condmoms <- .cond.mom(actual, moms$mu.gamma, actual$theta*moms$sig.gamma, subs=updates)
              cand <- mvrnorm(1, mu=condmoms$mu, Sigma = condmoms$sig)
              actual$gamma[updates] <- cand
            }
          }
          

          if(i==2){ #Accepting first candidate as startingpoint with probability = 1
            last <- actual
          }else{
            #Block acceptance of D, theta, nu and gamma_{\rho}
            logq1 <- .logprop(Y.c, X.c, a=actual, b=actual, pars=c("rho","theta","nu","gamma"), 
                             rho=rho.gamma, delta=delta.rho, mom=moms, eps=eps, r.gamma=r.gamma, updates=updates)
            logq2 <- .logprop(Y.c, X.c, a=last, b=actual, pars=c("rho","theta","nu","gamma"),
                             rho=NULL, delta=NULL, mom=NULL, eps=eps, r.gamma=r.gamma, updates=updates)
            test<-.metrop(Y.c, X.c, actual, last,logq1,logq2, eps)
            if(test){               
              last <- actual
              dvekobj$accept <- dvekobj$accept+1
              thetaobj$accept <- thetaobj$accept+1
              nuobj$accept <- nuobj$accept+1
            }
          }    

          #2) Sequential updating of gammas from mixture distribution.
          if(update$update.seq.gamma){
            actual <- last
            D <- actual$dvek
            A <- tcrossprod(X.c,t(D[,1:ncomp,drop=FALSE]))
            AAinv <- solve(crossprod(A))        
            moms <- .moments(actual, A, AAinv, Y.c)
            for(j in 1:ncomp){
              actual <- last
              mod <- runif(1,0,1)
              if(mod<eps$fi){
                cand <- rnorm(1, moms$mu.gamma[j], sqrt(moms$sig.gamma[j,j]))
              }else{
                cand <- 0
                while(cand==0){
                  u <- runif(1,0,1)
                  cand <- .fgammainv(u,-abs(moms$mu.gamma[j]), abs(moms$mu.gamma[j]), eps$gammaeps)
                }
              }
              actual$gamma[j] <- cand
              logq1 <- eps$fi*dnorm(actual$gamma[j],moms$mu.gamma[j], sqrt(moms$sig.gamma[j,j])) +
                (1-eps$fi)*.dgamma.propos(actual$gamma[j],-abs(moms$mu.gamma[j]), abs(moms$mu.gamma[j]), eps$gammaeps)
              logq2 <- eps$fi*dnorm(last$gamma[j],moms$mu.gamma[j], sqrt(moms$sig.gamma[j,j])) +
                (1-eps$fi)*.dgamma.propos(last$gamma[j],-abs(moms$mu.gamma[j]), abs(moms$mu.gamma[j]), eps$gammaeps)
              test <- .metrop(Y.c, X.c, actual,last,logq1,logq2, eps)
              if(test){
                last <- actual
                gammaobj$accept[j] <- gammaobj$accept[j]+1
              } 
              
            }
          }
          
          
          if(dosave){
            dvekobj$solu[k,,]<-last$dvek
            thetaobj$solu[k]<-last$theta
            nuobj$solu[k,]<-last$nu
            gammaobj$solu[k,]<-last$gamma
          }        
          
          
          
          #Computing beta for each iteration
          betasolu[k,] <- last$dvek[,1:ncomp,drop=FALSE]%*%last$gamma
          SSE[k] <- sum((Y.c-X.c%*%betasolu[k,])^2)

          
          #Adjustment for proposal tuning parameter(s)
          if((is.element(i,adapttime))&&(i<=totiter))
          {
            dvekobj[2:4] <- .adjust(dvekobj, minrate=0.1, maxrate=0.3, adaptint=adaptint)
            thetaobj[2:4] <- .adjust(thetaobj, doadjust=FALSE, adaptint=adaptint)
            nuobj[2:4] <- .adjust(nuobj, doadjust=FALSE, adaptint=adaptint) 
            gammaobj[2:4] <- .adjust.multi(gammaobj, doadjust=FALSE, adaptint=adaptint)            
          }
          
          #Trace plots
          if((is.element(i,plottime)) & dotrace){
            plot(thetaobj$solu[1:k],type="l",xlab="iter",ylab=expression(theta),main=paste("Acceptance rate =",thetaobj$rate))
            matplot(gammaobj$solu[1:k,],type="l",lty=1,xlab="iter",ylab=expression(gamma),main=paste("Acceptance rate =",gammaobj$rate[1]))
            matplot(nuobj$solu[1:k,],type="l",lty=1,xlab="iter",ylab=expression(nu),main=paste("Acceptance rate =",nuobj$rate))
            matplot(betasolu[1:k,],type="l",lty=1,xlab="iter",ylab=expression(beta),main="beta coefficients (for pot. scaled data)")
         }
          
          dosave <- FALSE
          cat(i,"\n")
        }
        from <- ceiling(start/thin)
        betahat <- apply(betasolu[from:k,],2,mean)
        if(scale){ 
            betahat <- betahat*sdY/sdX
        }
        beta0 <- meanY - betahat%*%meanX
        sigma.sq <- mean(thetaobj$solu[from:k])
        if(scale){
          sigma.sq <- sdY^2*sigma.sq
        }
#        df.reg <- n-1-mean(SSE[from:k])/sigma.sq
        
        res <- list(
          coefficients = betahat,
          intercept = beta0,
          sigma.sq = sigma.sq,
#          df.reg = df.reg,
          D = dvekobj,
          gamma = gammaobj,
          nu = nuobj,
          theta = thetaobj,
          betas = betasolu,
          SSEs = SSE,
          scale = scale,
          Y = Y,
          X = X,
          last = list(gamma=last$gamma,nu=last$nu,theta=last$theta,dvek=last$dvek)
        )
        class(res) <- "BayesPLS"
        res
}
