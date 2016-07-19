library(simrel)
set.seed(1430)
sim <- simrel(n=50, p=100, m=2, q=50, relpos=c(1,3), gamma=0.5, 
              R2=0.8, muY=3, muX=rnorm(100,5,1), ntest=1000)
X <- sim$X
Y <- sim$Y
ncomp <- 3

Rprof("BatesPLSeval.txt")
test <- BayesPLS(Y,X,ncomp,scale=TRUE, totiter=100, start=20, thin=5,
                 dotrace=TRUE, adaptint=50, plotint=10, approx=TRUE, 
                 approp=0.7, appit = 3)
Rprof()

summaryRprof("BatesPLSeval.txt")
