\name{lims.findSolventNames}
\alias{lims.findSolventNames}
\title{
retrieves the names of solvents contained in one or several entries
}
\description{
This function allows to retrieve the names of the solvents for one or several entries.
}
\usage{
	lims.findSolventNames(entryList)
}
\arguments{
  \item{entryList}{
 an list of entries produced by lims.getJSON
}
}
\details{
Solvents can be used to filter and download only the experiments acquired in one particular solvent. It might be useful for profiling if several solvents are used for analysis.
}
\references{
www.mylims.org
}
\author{
Julien Wist
}

\seealso{
See package "rjson" and www.mylims.org. See also \link{lims.getJSON} and \link{lims.findExperimentNames}.
}

\examples{
#example:
url <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=513117&key=cx9QmnDpDF'
entryList <- lims.getJSON(url)
lims.findExperimentNames(entryList)
lims.findSolventNames(entryList)
}

\keyword{ lims }
\keyword{ nmr }
\keyword{ metabolomics }
\keyword{ plot }
