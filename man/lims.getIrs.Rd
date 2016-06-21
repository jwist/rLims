\name{lims.getIrs}
\alias{lims.getIrs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
lims.getIrs
}
\description{
This function 
}
\usage{
lims.getIrs(entryList = entryList, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{entryList}{
an list of entries produced by 
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
  \item{info }{information about registration entry}
  \item{nmrInfo}{acquisition parameters}
  \item{param}{parameters of dataset}
  \item{spectra}{list spectra}
    

}
\references{
 mylims.univalle.edu.co
}
\author{
Julien Wist, Jessica Medina
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{lims.getNmrs}}
}
\examples{

#example:
url1 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=7165849&key=7KuEVXcaq8'
url2 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=9319999&key=rr3Cz90H2F'

url <- c(url1,url2)
entryList <- lims.JSON.getEntry(url)
entryData <- lims.getParameters(entryList)
Filter <- list('experiment'=c('IR'))
nmrList <- lims.getIrs(entryList,entryData=entryData,Filter=Filter,OP='AND')
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ir }
\keyword{ NMR }% __ONLY ONE__ keyword per line
