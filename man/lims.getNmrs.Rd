\name{lims.getNmrs}
\alias{lims.getNmrs}
\title{
retrieves the experiments from an list of entries and returns a list of spectra objects
}
\description{
This function retrieves the experiments from an list of entries and returns a list of spectra objects.
}
\usage{
	lims.getNmrs(entryList=entryList,factor,...)
}
\arguments{
  \item{entryList}{
 an list of entries produced by 
}
 \item{factor}{
 if there is imaginary data factor=2, if is not factor=1 
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

  \value{
  \item{info}{information about registration entry}
  \item{nmrInfo}{acquisition parameters}
  \item{param}{parameters of dataset}
  \item{spectra}{list spectra}
   }

\references{
www.mylims.org
}
\author{
Julien Wist
}

\seealso{
See package "rjson" and www.mylims.org. See also \code{\link{lims.JSON.getEntry}}.
}

\examples{

#example:
url1 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=8026534&key=5quA7vBcoQ'
url2 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=8405946&key=MdmsufXhXk'

url <- c(url1,url2)
entryList <- lims.JSON.getEntry(url)
entryData <- lims.getParameters(entryList)
Filter <- list('experiment'=c('noesygpps1dcomp'),'solvent'=c('COFFEEmeoh'))
nmrList <- lims.getNmrs(entryList,entryData=entryData,Filter=Filter,OP='AND')
}

\keyword{ lims }
\keyword{ nmr }
\keyword{ metabolomics }
\keyword{ plot }
