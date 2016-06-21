\name{findRepetitions}
\alias{findRepetitions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
findRepetitions
}
\description{This function allows find the repetition of a sample in the nmrList
}

\usage{
findRepetitions(nmrList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nmrList}{
List of input to create matrix data
}
}
\value{
  \item{nmrList }{Description of 'comp1'}
}
\references{mylims.org}

\author{
Julien Wist, Jessica Medina
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
mylims.univalle.edu.co
}

\examples{
url1 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=8026534&key=5quA7vBcoQ'
url2 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=6528763&key=BDiayNJUS8'

url <- c(url1,url2)
entryList <- lims.JSON.getEntry(url)
entryData <- lims.getParameters(entryList)
Filter <- list('experiment'=c('noesygpps1dcomp'),'solvent'=c('COFFEEmeoh'))
nmrList <- lims.getNmrs(entryList,entryData=entryData,Filter=Filter,OP='AND')
#findRepetitions(nmrList)
}

\keyword{ nmrList }
\keyword{ entryID }

