\name{sampleParam}
\alias{sampleParam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
sampleParam
}
\description{
Find the parameters with numer of entry
}
\usage{
sampleParam(name,data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{ number of the entry}
  \item{data}{data with parameters}
}

\value{
  \item{index of spectra}{position of the spectra}
 \item{param}{parameters of the selected sample}
%% ...
}
\references{
mylims.org
}
\author{
Julien Wist, Jessica Medina
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also  \code{\link{lims.getNmrs}}.
}
\examples{
data(coffee)

sampleParam("8571129",data)

}
\keyword{ entry }
\keyword{ param }% __ONLY ONE__ keyword per line
