\name{lims.scaling}
\alias{lims.scaling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Data centering and scaling
}
\description{
Provides a more convinient interface to the core "scale" function that
performs centering and scaling of data based on functions rather than
vectors.
}
\usage{
lims.scaling(data, center, scale)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data matrix}
  \item{center}{Function used for data centering. Uses "mean" by
    default. FALSE causes data not to be centered.}
  \item{scale}{Function used for data scaling. Uses "sd" by
    default. FALSE causes data not to be scaled. "pareto", "range" and
    "vast" methods are provided in the rLims package (see "lims.pareto",
    "lims.range", "lims.vast".}
  
}
\value{
  Returns a data matrix with the original data centered and scaled,
  i.e. each observation xi of variable x is replaced by (xi - center(x)) / scale(x)
 }
\references{
\href{http://dx.doi.org/10.1186/1471-2164-7-142}{van den Berg RA, Hoefsloot HC, Westerhuis JA, Smilde AK, van der Werf MJ. 2006. "Centering, scaling, and transformations: improving the biological information content of metabolomics data." BMC genomics (7). 142, 2013}
}


\author{
Julien Wist, Jessica Medina
}

\examples{
#scale and center intensities on a data set of NMR spectra
data(coffee)

#univariance method
nmrData<-data$nmrData
univariance<-lims.scaling((nmrData), scale=sd)
plot(univariance[1,], type="l") #plot the first scaled spectrum

#vast method
vast<-lims.scaling(nmrData, scale=lims.vast)
plot(vast[1,], type="l")

#pareto method
pareto<-lims.scaling(nmrData, scale=lims.pareto)
plot(pareto[1,], type="l")

#range method
range<-lims.scaling(nmrData, scale=lims.range)
plot(range[1,], type="l")

#custom (and arguably not so clever) method
madness <- lims.scaling(nmrData, center=median, scale=function(x){
mean(x) / sd(x) })
plot(madness[1,], type="l")

}

\keyword{ centering }
\keyword{ scaling }
