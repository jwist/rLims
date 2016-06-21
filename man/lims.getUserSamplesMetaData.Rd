\name{lims.getUserSamplesMetaData}
\alias{lims.getUserSamplesMetaData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
lims.getUserSamplesMetaData
}
\description{
Contains info of the data refered as metadata

}
\usage{
lims.getUserSamplesMetaData(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{Account from an user in mylims.org }
}

\value{

\item{list of parameters}{Allows to see which experiments are in the account of the user, which nucleus, temperature and quantity of spectrums}
}
\references{
mylims.org
}
\author{
Julien Wist

}


\seealso{
See package "rjson" and www.mylims.org. See also \code{\link{lims.JSON.getEntry}},
\code{\link{lims.getUserEntries}}
}


\examples{
url='http://mylims.univalle.edu.co/lims/default/sample/listJson.jsp?&queryUser=jessica.medina@correounivalle.edu.co&key=PDfEZzKo3K'
account <- list(server='mylims.univalle.edu.co',user='jessica.medina@correounivalle.edu.co',key='PDfEZzKo3K')
userData <- lims.getUserSamplesMetaData(account=account)
}
\keyword{ lims}
\keyword{ json }
\keyword{ nmr }
