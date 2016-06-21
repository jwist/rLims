\name{plotMetabo}
\alias{plotMetabo}

\title{
plotMetabo
}
\description{
This function makes the plots used in metabolomics: s-plot, stocsy, load plot
}
\usage{
plotMetabo(data,class,nComp,ppm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{matrix data}
  \item{class}{classification of data}
  \item{nComp}{number of components}
  \item{ppm}{ppm of spectra}
}
\details{
These functions are used to find a biomarkers in metabolic data 
}
\value{
 \item{s-plot }{s-plot where is calculated corr(t,X) vs (cov(t,X))}
 \item{loading plot}{loading1 vs loading2 are ploted}
 \item{colored loading plot}{using back scaled loadings are ploted as color code to find discriminant variables}
}

\references{
\href{http://doi.wiley.com/10.1002/cem.695}{Trygg, Johan, and Svante Wold. 2002. Orthogonal projections to latent structures (O-PLS). Journal of Chemometrics 16(3),119 128}
}

\author{
Julien Wist, Jessica Medina}


\seealso{
See package "rjson" and www.mylims.org. See also \code{\link{lims.opls.model}}


}
\examples{

data(coffee)
nmrData<-data$nmrData
species<-data$param$species
ppm<-data$ppm


train<-sample(1:dim(nmrData)[2],6)
data<-nmrData[,train]
y.train<-species[train]
test<-nmrData[,-train]
y.test<-species[-train]


plotAnalysis<-plotMetabo(t(data),y.train,2,ppm)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ loadings }
\keyword{ s-plot }% __ONLY ONE__ keyword per line
