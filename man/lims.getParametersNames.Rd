\name{lims.getParametersNames}
\alias{lims.getParametersNames}
\title{
retrieves the names of the parameters stored in a list of entries
}
\description{
This function allows to retrieve the names of the parameters of an entry. If several entries are present we report the union of the names.
}
\usage{
	lims.getParametersNames(jsonList)
}
\arguments{
  \item{jsonList}{
 an list of entries produced by lims.getJSON
}
}
\details{
XXX
}
\references{
www.mylims.org
}
\author{
Julien Wist
}

\seealso{
See package "rjson" and www.mylims.org. See also lims.getJSON.
}

\examples{

#example:
url1 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=8026534&key=5quA7vBcoQ'
url2 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=8405946&key=MdmsufXhXk'

url <- c(url1,url2)
entryList <- lims.getJSON(url)
names <- lims.getParametersNames(entryList)
}

\keyword{ lims }
\keyword{ nmr }
\keyword{ metabolomics }
\keyword{ plot }
