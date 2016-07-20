\name{BayesPLS}
\alias{BayesPLS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian estimation in PLS regression
}
\description{
A Bayesian estimator for the population Partial Least Squares Model 
based on the concept of relevant components and MCMC
}
\usage{
BayesPLS(Y, X, ncomp, scale = TRUE, totiter = 1000, start = 1, doinit = TRUE, init.method = c("PLS"), 
dotrace = TRUE, plotint = 10, thin = 10, adaptint = 100, approx = FALSE, appit = 100, appcomp = 50,
approp = 0.999, update = list(update.dvek = TRUE, update.nu = TRUE, 
update.theta = TRUE, update.gamma = TRUE, update.seq.gamma = TRUE), eps = list(gammaeps = 1/19, nueps = 1/31, 
thetaeps = 1/31, dvekeps = 5e-05, lambda = -log(0.001), fi = 0.5), previousobj = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
A response vector
}
  \item{X}{
A predictor matrix
}
  \item{ncomp}{
The number of latent components
}
  \item{scale}{
Should the data be scaled (standardized)? Defaults to TRUE.
}
  \item{totiter}{
Total number of MCMC iterations
}
  \item{start}{
At which result iteration the estimation of parameters should start. All previous iterations are regarded as burn-in.
Note, the value refers to iteration number before thinning.
}
  \item{doinit}{
Should the parameters be initiated using regular PLSR or PCR? Defaults to TRUE, but may be set to FALSE if a previous
run object is supplied as extra argument. In that case the final iteration of the previous run is used as initial
values for the parameters.
}
  \item{init.method}{
Method for initiating the parameters, either "PLS" or "PCR". Both use the specified ncomp as the number
of components.
}
  \item{dotrace}{
Should the traceplots be generated? Defaults to TRUE.
}
  \item{plotint}{
The frequency of trace plot updating. Defaults to every 10th generated value (after thinnning).
}
  \item{thin}{
The thinning value. Defaults to 10 in which case only every 10th value is saved.
}
  \item{adaptint}{
The step-lengths for the rotation of the eigenvectors (D) is tuned with this frequency. Defauilts to 100.
At this frequency the acceptance rate is also reported i nthe trace plots.
}
  \item{approx}{
Logical. If TRUE then only a subspace corresponding to the largest eigenvalues (nu's) are updated after \code{appit} iterations. The dimensions of the subspace is set by the minimum cumulative sum of the eigenvalues determined by \code{approp}. 
}
  \item{appit}{
  The number of initial rotations of the complete predictor space before a reduced space is used if 
\code{approx = TRUE}
  }
  \item{appcomp}{
  The maximum dimension of the reduced space if \code{approx = TRUE}.
}
  \item{approp}{
The proportion of the total variance of X (cumulative sum of eigenvalues) to be explained
by the reduced space (if \code{approx = TRUE}). Note that \code{appcomp} is the maximum dimension.
}
  \item{update}{
A logical list. Sets whether various parameters should be estimated (updated in MCMC) or not. 
Defaults to TRUE for all parameters.
}
  \item{eps}{
Some internal tuning parameters for priors and the step-length parameter for the rotation of 
the eigenvector matrix D.
}
  \item{previousobj}{
A previously fitted BayesPLS object from which starting values will be retrieved if initiate = FALSE
and the model specification and data is identical.
}
}
\details{
To be filled in.
}
\value{
\item{coefficients }{The estimated regression coefficients for the original (unscaled) predictors.}
\item{intercept }{The estimated intercept of the model.}
\item{sigma.sq }{The estimated noise variance (same as the parameter theta in the model).}
%\item{df.reg }{An estimate of the degrees of freedom of regression.}
\item{D }{The MCMC object of the eigenvectors, D. Contains a sub-element 'solu' with the sampled values.}
\item{gamma }{The MCMC object of the weights of the latent components. Contains a sub-element 'solu' with the sampled values.}
\item{nu }{The MCMC object of the eigenvalues, nu. Contains a sub-element 'solu' with the sampled values.}
\item{theta }{The MCMC object of the noise variance. Contains a sub-element 'solu' with the sampled values.}
\item{betas }{The MCMC object of the regression coefficients (for scaled predictors).}
\item{SSEs }{The MCMC object of the error sums of squares.}
\item{scale }{Logical. Were the variables standardized?}
\item{Y}{The input response data.}
\item{X}{The input predictor data.}
\item{last}{The final iteration value for all parameters.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Solve S?b?
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
library(simrel)
set.seed(1430)
sim <- simrel(n=30, p=50, q=50, relpos=c(1,3), gamma=0.5, R2=0.8, muY=3, muX=rnorm(50,5,1), ntest=1000)
X <- sim$X
Y <- sim$Y

test <- BayesPLS(Y,X,3,scale=TRUE, totiter=500, start=20, dotrace=TRUE, adaptint=50, plotint=10)
pred <- predict(test, sim$TESTX)
windows()
plot(sim$TESTY, pred, xlab="Observed", ylab="Predicted", main="Prediction plot for test data")
curve(1*x, add=TRUE, col=2)

print(DIC(test,start=20))

}
\keyword{models }% use one of  RShowDoc("KEYWORDS")
\keyword{multivariate}% __ONLY ONE__ keyword per line
\keyword{regression}
