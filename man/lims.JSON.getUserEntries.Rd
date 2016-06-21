\name{lims.JSON.getUserEntries}
\alias{lims.JSON.getUserEntries}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
lims.JSON.getUserEntries
}
\description{
Retrieves urls entries from an account in mylims.org
}
\usage{
lims.JSON.getUserEntries(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{json url}
}

\value{
 \item{urlList}{list of url}
 
}
\references{
mylims.org
}
\author{
Julien Wist
}

\seealso{
See package "rjson" and www.mylims.org. See also \code{\link{lims.JSON.getEntry}}
}
\examples{
url='http://mylims.univalle.edu.co/lims/default/sample/listJson.jsp?&queryUser=jessica.medina@correounivalle.edu.co&key=PDfEZzKo3K'
account <- list(server='mylims.univalle.edu.co',user='jessica.medina@correounivalle.edu.co',key='PDfEZzKo3K')
entryList <- lims.JSON.getUserEntries(url=url)

}

\keyword{ rjson }
\keyword{ lims }% __ONLY ONE__ keyword per line
