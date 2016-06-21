\name{lims.JSON.getEntry}
\alias{lims.JSON.getEntry}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
lims.JSON.getEntry
}
\description{
Returns from an account in mylims.org a url list of each entry
}
\usage{
lims.JSON.getEntry(url, LOG = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{url}{
JSON url
}
  \item{LOG}{
log
}
  \item{\dots}{
list
}
}

\value{

 \item{json data }{return list json data}

}
\references{
mylims.org
}
\author{
Julien Wist
}

\seealso{
See package "rjson" and www.mylims.org. See also \code{\link{lims.getUserEntries}}
}
\examples{

url='http://mylims.univalle.edu.co/lims/default/sample/listJson.jsp?&queryUser=jessica.medina@correounivalle.edu.co&key=PDfEZzKo3K'
entryList <- lims.JSON.getUserEntries(url=url)
entries <- lims.JSON.getEntry(entryList)


}

\keyword{ json }
