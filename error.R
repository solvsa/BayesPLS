parSetting <- structure(list(n = 50, 
                             p = 15, 
                             R2 = 0.5, 
                             relpos = c(1, 2), 
                             gamma = 0.5, 
                             q = 15,
                             ntest = 5000
), 
.Names = c("n", "p", "R2", "relpos", "gamma", "q", "ntest"))

set.seed(777)

simObj <- do.call(simrel, parSetting)

train <- data.frame(x = I(simObj$X), y = I(simObj$Y))
test <- data.frame(x = I(simObj$TESTX), y = I(simObj$TESTY)) 

bayesObj4 <- with(train, BayesPLS(y, x, ncomp = 1, scale = FALSE, init.method = "PCR",
                                 dotrace = FALSE, totiter = 50000, freeze = 0.02, 
                                 compreduce = TRUE, thin = 50, adaptint = 500, plotint = 500,
                                 approx = TRUE))