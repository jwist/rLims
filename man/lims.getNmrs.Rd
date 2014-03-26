\name{lims.getNmrs}
\alias{lims.getNmrs}
\title{
retrieves the experiments from an list of entries and returns a list of spectra objects
}
\description{
This function retrieves the experiments from an list of entries and returns a list of spectra objects.
}
\usage{
	lims.getNmrs(entryList=entryList,...)
}
\arguments{
  \item{entryList}{
 an list of entries produced by lims.getJSON
}
  \item{...}{
 Additional parameters:
 
 -entryData: an entryData (entryData <- lims.getParameters(entryList))
 
 -Filter: a list of experiment names and solvents
 (Filter <- list('experiment'=c('zg30'),'solvent'=c('C6D6')))
 
 -LOG: a boolean, if TRUE log will be printed
 
 -dry: if true, only a dry run is performed but no real spectra is downloaded.
}
}
\references{
www.mylims.org
}
\author{
Julien Wist
}

\seealso{
See package "rjson" and www.mylims.org. See also \link{lims.getJSON}.
}

\examples{

#example:
url1 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=8026534&key=5quA7vBcoQ'
url2 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=8405946&key=MdmsufXhXk'

url <- c(url1,url2)
entryList <- lims.getJSON(url)
names <- lims.getParametersNames(entryList)
entryData <- lims.getParameters(entryList)
lims.findExperimentNames(entryList)
lims.findSolventNames(entryList)
Filter <- list('experiment'=c('noesygpps1dcomp'),'solvent'=c('COFFEEmeoh'))
nmrList <- lims.getNmrs(entryList,entryData=entryData,Filter=Filter,OP='AND')
}

\keyword{ lims }
\keyword{ nmr }
\keyword{ metabolomics }
\keyword{ plot }
