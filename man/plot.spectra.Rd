\name{plot.spectra}
\alias{plot.spectra}
\title{
methods to plot of spectra objects
}
\description{
Methods to plot spectra objects
}
\usage{
	\method{plot}{spectra}(x,...)
}
\arguments{
  \item{x}{
 a spectra object to be ploted
}
  \item{...}{
 xlimit and ylimit can be set to define ranges for the x and y axes.
}
}
\references{
www.mylims.org
}
\author{
Julien Wist
}

\seealso{
See lims.spectraCreator and plot (generic function)
}

\examples{

# example 
\dontrun{
fakeSpectra <- list(x=c(1:100),y=sin(seq(1,6.27,length.out=100)))
s <- lims.spectraCreator(fakeSpectra)
plot(s,xlimit=c(40,60))
}
}

\keyword{ lims }
\keyword{ nmr }
\keyword{ metabolomics }
\keyword{ plot }
