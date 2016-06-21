\name{plotMulti}
\alias{plotMulti}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot multiple spectra}
\description{
Using this function allows plot the spectra by regions and colored by class if the user wants to
}
\usage{
plotMulti(ppm, nmrData, A, B, col)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ppm}{ppm of spectra }
  \item{nmrData}{matrix of spectra}
  \item{A}{lower limit of ppm}
  \item{B}{upper limit of ppm}
  \item{col}{vector of colors, has the size of the data spectra, if you use for example col=1 plot in a single color}
}
%\details{
%}

\author{
Julien Wist, Jessica Medina
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

See also  \code{\link{plotMetabo}}.

}
\examples{

data(coffee)

nmrData<-data$nmrData
ppm<-data$ppm

A<-1
B<-5
groups = factor(gsub('(.).', '\\1', data$param$specie))
col=c("blue","red")[groups]

plotMulti(ppm,t(nmrData), A,B,col)
plotMulti(ppm,t(nmrData), A,B,3)

}
\keyword{ plot }
\keyword{ nmr }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
