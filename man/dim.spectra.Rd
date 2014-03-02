\name{dim.spectra}
\alias{dim.spectra}
\title{
method to get dimension of spectra objects
}
\description{
Method to get dimension of spectra objects
}
\usage{
	\method{dim}{spectra}(x)
}
\arguments{
  \item{x}{
 a spectra object 
}
}

\references{
www.mylims.org
}
\author{
Julien Wist
}

\seealso{
See lims.spectraCreator and dim (generic function)
}

\examples{
# example 
# fakeSpectra <- list(x=c(1:100),y=sin(seq(1,6.27,length.out=100)))
# s <- lims.spectraCreator(fakeSpectra)
# dim(s)
}

\keyword{ lims }
\keyword{ nmr }
\keyword{ metabolomics }
\keyword{ dim }
