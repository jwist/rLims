\name{lims.getParameters}
\alias{lims.getParameters}
\title{
retrieves the parameters for a list of entries and return a data.frame
}
\description{
This function allows to retrieve the parameters of one or several entries and returns a data.frame.
}
\usage{
	lims.getParameters(entryList)
}
\arguments{
  \item{entryList}{
 an list of entries produced by lims.getJSON
}
}
\details{
Only certain entries will contains parameters. This is particularly useful for metabolic profiling, since a complete data.frame may be obtained from the lims. It is convenient to store the data along with their metadata.
}
\references{
www.mylims.org
}
\author{
Julien Wist
}

\seealso{
See package "rjson" and www.mylims.org. See also \link{lims.getJSON}, \link{lims.findSolventNames} and \link{lims.findExperimentNames}.
}

\examples{

#example:
url <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=1144728&key=QpOcIBAE2M'
entryList <- lims.getJSON(url)
lims.findExperimentNames(entryList)
lims.findSolventNames(entryList)
param <- lims.getParameters(entryList)
}

\keyword{ lims }
\keyword{ nmr }
\keyword{ metabolomics }
\keyword{ plot }
