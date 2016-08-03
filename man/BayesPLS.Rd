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
BayesPLS(Y, X, ncomp, scale = TRUE, totiter = 1000, doinit = TRUE, init.method = c("PLS"), init.ncomp = NULL, init.sort = FALSE, dotrace = TRUE, plotint = 10, thin = 10, adaptint = 100, approx = FALSE, appit = 10, appcomp = 100, approp = 0.999, update = list(update.dvek = TRUE, update.nu = TRUE, update.theta = TRUE, update.gamma = TRUE, update.seq.gamma = TRUE), eps = list(gammaeps = 1/19, nueps = 1/31, thetaeps = 1/31, lambda = -log(0.001), fi = 0.5), deps = 0, compreduce = TRUE, alpha.reduce = 0.4, freeze = 0.1, previousobj = NULL)
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
  \item{doinit}{
Should the parameters be initiated using regular PLSR or PCR? Defaults to TRUE, but may be set to FALSE if a previous
run object is supplied as extra argument. In that case the final iteration of the previous run is used as initial
values for the parameters.
}
  \item{init.method}{
Method for initiating the parameters, either "PLS" or "PCR". Both use the specified ncomp as the number
of components.
}
  \item{init.ncomp}{
The number of components used by PLS or PCR in the initiation step for the Bayes PLS model parameters.
Deafults to NULL which means that init.ncomp = ncomp. 
}
  \item{init.sort}{
Should the components found in the intiation step through PLS or PCR be sorted in a decreasing order according to the initial values of their weights? Deafults to FALSE which is most appropriate for init.method = "PLS", If init.method ="PCR" init.sort = TRUE is recommended.
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
Some internal tuning parameters for priors and proposal distributions.
}
  \item{deps}{
Sets the step-length parameter for the rotation of the eigenvector matrix D.
Defaults to 0 in which case the steps are determined as a function of n, p and ncomp.
Other integer values of deps can be used to decrease (deps > 0) or increase (deps > 0) the initial step lengths. The adjustment will then be equal to 2^{deps}. If deps
is a non-integer (if modulo deps is not zero) the value is interpreted as an actual
step length value to replace the internal default value.
}
  \item{compreduce}{
Logical. If TRUE the weights for the components (gamma-values) are tested using bayesian p-values.
If the Bayesian p.value exceeds alpha.reduce for a given component, the component is excluded from further
iterations. 
}
  \item{alpha.reduce}{
Test level for testing the significance (from zero) for the gammas. Defaults to 0.4 which means that if the proportion of sampled values above (or below) zero is in the range [0.4, 0.6], the component is regarded as irrelevant.
}
  \item{freeze}{
A number between 0 and 1 which sets the proportion of totiter to be run before any component may be removed from the model. Defaults to 0.1. The same number of iterations will be run after a removal before any other component may be removed from the model.
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
est <- estimate(test, start=20)
pred <- predict(est, sim$TESTX)
windows()
plot(sim$TESTY, pred, xlab="Observed", ylab="Predicted", main="Prediction plot for test data")
curve(1*x, add=TRUE, col=2)

print(DIC(est,start=20))

}
\keyword{models }% use one of  RShowDoc("KEYWORDS")
\keyword{multivariate}% __ONLY ONE__ keyword per line
\keyword{regression}
