\name{estimateBayesPLS}
\alias{estimateBayesPLS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation function for Bayes PLS
}
\description{
Estimates regression coefficients (with quantiles) and noise variance from a fitted Bayes PLS object.
}
\usage{
estimateBayesPLS(obj, start = 1, stop = NULL, thin = 1, probs = c(0.025, 0.975))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
A fitted Bayes PLS object as returned from a call to BayesPLS().
}
  \item{start}{
Start from this iteration number.
}
  \item{stop}{
Stop at this iteration number.
}
  \item{thin}{
Thinning interval.
}
  \item{probs}{
Quantiles to be computed from the estimated posterior distributions.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{coefficients }{The estimated regression coefficients.}
  \item{intercept }{The estimated intercept coefficient.}
  \item{sigma.sq }{The estimated noise variance.}
  \item{quantiles }{The quantiles (as set in probs) for the posterior distributions of the regression coefficients.}
  \item{quantiles.int }{The estimated quantiles for posterior distribution of the intercept coefficient.}
  \item{quantiles.sigma.sq }{The estimated quantiles for posterior distribution of the noise variance.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Solve Sæbø
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## To be filled in
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ inference }% use one of  RShowDoc("KEYWORDS")
\keyword{ model }% __ONLY ONE__ keyword per line
