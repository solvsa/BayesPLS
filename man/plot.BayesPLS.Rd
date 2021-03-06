\name{plot.BayesPLS}
\alias{plot.BayesPLS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Trace plots for Bayes PLS object
}
\description{
This function reproduces the plots which typically are generated during the call to BayesPLS. 
}
\usage{
plot(obj, start = NULL, stop = NULL, thin = 1, which=c(1L:4L), ask=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
A fitted Bayes PLS object
}
  \item{start}{
Plot from this iteration number. Default is NULL which means from start of file.
}
  \item{stop}{
Plot to this ieration number. Deafult is NULL which means end of file.
}
  \item{thin}{
Thinning interval.
}
  \item{which}{
if a subset of the plots is required, specify a subset of the numbers 1:4.
}
  \item{ask}{
logical; if TRUE, the user is asked before each plot, see par(ask=.)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
##See ?BayesPLS
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ graphics }% use one of  RShowDoc("KEYWORDS")
\keyword{ model }% __ONLY ONE__ keyword per line
