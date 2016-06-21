\name{lims.vast}
\alias{lims.vast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Vast data scaling
}
\description{
Vast scaling factor
}
\usage{
lims.vast(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A numeric vector, array, or data frame.}
}
\value{
  Returns the Vast scaling factor on the argument, i.e. the ratio of
  the squared standard deviation to the mean.
 }
\references{
\href{http://dx.doi.org/10.1186/1471-2164-7-142}{van den Berg RA, Hoefsloot HC, Westerhuis JA, Smilde AK, van der    Werf MJ. 2006. "Centering, scaling, and transformations: improving the biological information content of metabolomics data." BMC genomics (7). 142, 2013}
}

\author{
Julien Wist, Jessica Medina
}

\seealso{
  See also \code{\link{lims.scaling}}
  }

\examples{
#normalize the intensities on a data set of NMR spectra using the Vast method
data(coffee)
vast<-lims.scaling(data$nmrData, scale=lims.vast)
plot(vast[1,], type="l") #plot the first scaled spectrum


}

\keyword{ centering }
\keyword{ scaling }
