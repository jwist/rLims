\name{lims}
\alias{lims}
\title{
retrieve spectra from www.mylims.org database
}
\description{
retrieve a set of spectra from www.mylims.org database, chosen by experiment names.
}
\usage{
lims(urlList, experimentList = list(), ...)
}
\arguments{
  \item{urlList}{
A list of URLs. Each URL should be the complete URL that is obtained by launching pipeline from www.mylims.org
}
  \item{experimentList}{
A list of experiment names that should be retrieved.
}
  \item{...}{
Addtional parameters to be passed to the function.
}
}
\details{
URLs can be obtained by replacing the correct information for entryID and key in "http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=6528214&key=NK2pYSFIhD"
}
\references{
www.mylims.org
}
\author{
Julien Wist
}

\examples{

## use this command to read from a file in.txt
## listOfURL<-as.vector(unlist(read.table('in.txt',sep='')))

data(urlList)
dataSet <- lims(as.vector(unlist(urlList)),c('noesygpps1dcomp'))

}

\keyword{ lims }
\keyword{ nmr }
\keyword{ metabolomics }
