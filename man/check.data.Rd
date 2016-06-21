\name{check.data}
\alias{check.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
check data
}
\description{
Check outliers and filter data
}
\usage{
check.data(lowerBound,upperBound,nmr,ppm, param)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lowerBound}{lowest ppm}
  \item{upperBound}{hightest ppm}
  \item{nmr}{nmrData}
   \item{ppm}{ppm}
   \item{param}{data parameters}
}

\details{
Delete from a list outliers using principal components
}

\value{

\item{listOfoutliers }{entryID of outliers}
 \item{nmrData}{data filtered}
  \item{param}{parameters of data filtered}
  \item{outliers}{data of outliers}


}
%\references{

%}

\author{
Julien Wist, Jessica Medina
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
www.mylims.org

}
\examples{
data(coffee)
ppm<-data$ppm
nmrData<-data$nmrData
param<-data$param
check<-check.data(1,3,t(nmrData),ppm, param)

}

\keyword{ outliers }
\keyword{ pca}% __ONLY ONE__ keyword per line

