\name{lims.spectra.toDataframe}
\alias{lims.spectra.toDataframe}
\title{
converts spectra objects to data.frame objects
}
\description{
Converts spectra objects into data.frame objects.
}
\usage{
lims.spectra.toDataframe(spectra,full=FALSE)
}
\arguments{
  \item{spectra}{
A spectra object
}
  \item{full}{
By default full is set to FALSE, if changed to TRUE, the function will return a data.frame object by ignoring existing masks. 
}
}
\references{
www.mylims.org
}
\author{
Julien Wist
}

\examples{

# example 
fakeSpectra <- list(x=c(1:100),y=sin(seq(1,6.27,length.out=100)))
s <- lims.spectraCreator(fakeSpectra)
## will use the plot.spectra methods
plot(s) 
## while this will use the generic plot function
plot(as.numeric(lims.spectra.toDataframe(s)))

}

\keyword{ lims }
\keyword{ nmr }
\keyword{ spectra }
\keyword{ chemometrics }
\keyword{ metabolic profiling }
