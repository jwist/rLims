\name{lims.createDataSet}
\alias{lims.createDataSet}
\title{
creates a data.frame object containing the spectra and the parameters for analysis
}
\description{
This function creates a data.frame object containing the spectra and the parameters for analysis.
}
\usage{
	lims.createDataSet(nmrList)
}
\arguments{
  \item{nmrList}{
 the result of lims.getNmrs
}
}
\references{
www.mylims.org
}
\author{
Julien Wist
}

\seealso{
See package "rjson" and www.mylims.org. See also lims.getNmrs.
}

\examples{

#example:
url1 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=8026534&key=5quA7vBcoQ'
url2 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=6528763&key=BDiayNJUS8'

url <- c(url1,url2)
entryList <- lims.getJSON(url)
names <- lims.getParametersNames(entryList)
entryData <- lims.getParameters(entryList)
lims.findExperimentNames(entryList)
lims.findSolventNames(entryList)
Filter <- list('experiment'=c('noesygpps1dcomp'),'solvent'=c('COFFEEmeoh'))
nmrList <- lims.getNmrs(entryList,entryData,Filter=Filter,OP='AND')
nmrList <- lims.getNmrs(entryList,entryData,Filter=Filter,OP='AND')
data <- lims.createDataSet(nmrList)
}

\keyword{ lims }
\keyword{ nmr }
\keyword{ metabolomics }
\keyword{ plot }
