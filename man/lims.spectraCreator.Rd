\name{lims.spectraCreator}
\alias{lims.spectraCreator}
\title{
create spectra object from x/y vector
}
\description{
Creates a object of class='spectra' from two x and y vectors.
}
\usage{
lims.spectraCreator(spec,type='nmr',mask=NULL)
}
\arguments{
  \item{spec}{
A list object of type list('x'=xVector,'y'=yVector)
}
  \item{type}{
The type is a string that can take the following values: 'nmr', 'ir'. If omited the type is 'nmr'.
}
  \item{mask}{
The mask is a boolean vector of length that matches the length of the x and y vectors. If omited the full spectra will be selected, i.e., the mask is a TRUE vector of length length(spec$x).
}
}
\references{
www.mylims.org
}
\author{
Julien Wist
}

\examples{
# first we create a x&y vector and then we cast it as a spectra object
fakeSpectra <- list(x=c(1:100),y=sin(seq(1,6.27,length.out=100)))
s <- lims.spectraCreator(fakeSpectra)
}

\keyword{ lims }
\keyword{ nmr }
\keyword{ spectra }
\keyword{ chemometrics }
\keyword{ metabolic profiling }
