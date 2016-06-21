\name{lims.getUserEntries}
\alias{lims.getUserEntries}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
lims.getUserEntries
}
\description{
Using an account retrieves ulrs for each entry in mylims.org
}
\usage{
lims.getUserEntries(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{the account of some user in mylims.org}
}
\value{
  \item{urlList}{List of url}
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
 entries <- lims.getUserEntries(account=account)
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ lims }
\keyword{ json }% __ONLY ONE__ keyword per line
