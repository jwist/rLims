\name{lims.findExperimentNames}
\alias{lims.findExperimentNames}
\title{
retrieves the names of experiments contained in one or several entries
}
\description{
This function allows to retrieve the names of the experiments for one or several entries.
}
\usage{
	lims.findExperimentNames(entryList)
}
\arguments{
  \item{entryList}{
 an list of entries produced by lims.getJSON
}
}
\details{
Experiment names can be used to filter and download only the type of experiment required for one analysis.
}
\references{
www.mylims.org
}
\author{
Julien Wist
}

\seealso{
See package "rjson" and www.mylims.org. See also \link{lims.getJSON} and \link{lims.findSolventNames}.
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
