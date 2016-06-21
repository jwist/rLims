\name{coffee}
\alias{data}
\docType{data}
\title{
nmr coffee data
}
\description{
NMR data of two different species of roasted coffee
}
\usage{data("coffee")}
\format{
  data frame with NMR data and parameters 
  \describe{
    \item{\code{ppm}}{ppm of NMR spectra}
    \item{\code{nmrData}}{matrix contains the intensity}
    \item{\code{param}}{parameters of each data}
  }
}
%\details{
%}

\source{
mylims.org
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(coffee)
## plot(data$ppm,data$nmrData[,1], xlim=c(10,0),type="l") ...
}
\keyword{datasets}
