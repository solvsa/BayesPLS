\name{DIC}
\alias{DIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the Deviance Information Criterion (DIC) for a fitted model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DIC(obj, start = 2, stop = length(obj$theta$solu))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
A BayesPLS object
}
  \item{start}{
Compute DIC starting with this iteration. Note, the value refers to iteration number after thinning.
}
  \item{stop}{
Compute DIC ending with this iteration. Note, the value refers to iteration number after thinning.
}
}
\details{
DIC accoring to Spiegelhalter et al. (2002)
}
\value{
  \item{dic}{The DIC value}
  \item{pd }{The estimated effective number of parameters estimated in the model.}
%% ...
}
\references{
Spiegelhalter, David J.; Best, Nicola G.; Carlin, Bradley P.; van der Linde, Angelika (2002). "Bayesian measures of model complexity and fit (with discussion)". Journal of the Royal Statistical Society, Series B 64 (4)
}
\author{
Solve Sæbø
}

\examples{
## See ?BayesPLS
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }% use one of  RShowDoc("KEYWORDS")
\keyword{ models }% __ONLY ONE__ keyword per line
