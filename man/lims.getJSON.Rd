\name{lims.getJSON}
\alias{lims.getJSON}
\title{
retrieves the json entry and converts it to a list
}
\description{
This function converts the json entry into a list object
}
\usage{
	lims.getJSON(url,LOG=FALSE,...)
}
\arguments{
  \item{url}{
 the url of an entry (string) or a list of url
}
  \item{LOG}{
 a boolean to switch log ON or OFF
}
  \item{...}{
 additional parameters
}
}
\details{
All resources can be retrieved as JSON objects. This function works to retrieve JSON for entries. It is not intended (and won't work) for retrieving other kind of resources directly, such as nmr or parameters urls.
}
\references{
www.mylims.org
}
\author{
Julien Wist
}

\seealso{
See package "rjson" and www.mylims.org.
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
