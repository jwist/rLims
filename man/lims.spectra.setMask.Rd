\name{lims.spectra.setMask}
\alias{lims.spectra.setMask}
\title{
selects the regions to display
}
\description{
This function allows to select the regions to be diplayed
}
\usage{
	lims.spectra.setMask(spectraVarName,list=list())
}
\arguments{
  \item{spectraVarName}{
 the name of the spectra variable (string)
}
  \item{list}{
 a c() list object in ascending order. An empty list will select the whole spectrum
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
\dontrun{#example
fake <- list(x=c(1:100),y=sin(seq(1,6.27,length.out=100)))
s <- lims.spectraCreator(fake)
plot(s)
lims.spectra.setMask('s',list=c(10,20,40,50,80,100))
plot(s)}
}

\keyword{ lims }
\keyword{ nmr }
\keyword{ metabolomics }
\keyword{ plot }
