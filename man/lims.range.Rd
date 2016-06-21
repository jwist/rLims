\name{lims.range}
\alias{lims.range}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Range interval of a data series
}
\description{
Computes the difference between the minimum and maximum of a data
series. Used for data scaling in some applications.
}
\usage{
lims.range(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A numeric vector, array, or data frame.}
}
\value{
  {x} {max(x) - min(x) over the argument provided}
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
#normalize the intensities on a data set of NMR spectra proportionally
#to intensity range
data(coffee)
range<-lims.scaling(data$nmrData, scale=lims.range)
plot(range[1,], type="l") #plot the first scaled spectrum


}

\keyword{ centering }
\keyword{ scaling }
