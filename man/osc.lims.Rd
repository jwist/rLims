\name{osc.lims}
\alias{osc.lims}
\title{
Orthogonal Signal Correction 
}
\description{
Orthogonal Signal Correction method for a spectra
}
\usage{
osc.lims(predictors, responses, nComp, training, xcenter = mean, xscale = sd, ycenter = mean, yscale = sd, discriminant = TRUE, accuracy = 1e-05, maxit = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{predictors}{Data matrix of predictor variables. Observations in
    rows, variables in columns.}
  \item{responses}{Data matrix of responses variables. Observations in
    rows, variables in columns.}
  \item{nComp}{Number of components or latent variables to be included
    in the PLS model}
     \item{training}{A vector with the indices of the rows of "predictors"
    to be used as training set}
         \item{xcenter}{Function used for predictor variables centering, see
    lims.scaling. "mean" by default.}
    \item{xscale}{Function used for predictor variables scaling, see
    lims.scaling. "sd" by default.}
      \item{ycenter}{Function used for response variables centering, see
    lims.scaling. "mean" by default.}
    \item{yscale}{Function used for response variables centering, see
    lims.scaling. "sd" by default.}
      \item{discriminant}{If discriminant=TRUE the function will evaluate discriminant analysis for two classes.}
      \item{accuracy}{Maximum change between iteration cycles of the NIPALS
       algorithm. 0.00001 by default.}
       \item{maxit}{Maximum number of iteration cycles. 100 by default.}
       }
%\details{
%%  ~~ If necessary, more details than the description above ~~

\value{

 \item{X}{Corrected matrix}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
\href{http://www.sciencedirect.com/science/article/pii/S0169743998001099}{Wold, S.,Antti, H., Lindgren, F., Ohman, J. 1998. Orthogonal signal correction of near-infrared spectra. Chemometrics and Intelligent Laboratory Systems. (44),175 185}
}
\author{
Julien Wist, Jessica Medina, Andres Bernal 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See package "rjson" and www.mylims.org. See also \code{\link{lims.opls.model}}
}
\examples{
data(coffee)
nmrData<-t(data$nmrData)
composition<-matrix(as.numeric(as.character(data$param$composition)),ncol=1)
model<-osc.lims(nmrData,composition,nComp=5,training=1:11, discriminant = FALSE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ osc }
\keyword{ signal correction}
