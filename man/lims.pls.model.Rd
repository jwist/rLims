\name{lims.pls.model}
\alias{lims.pls.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Partial Least Squares Regression / Discriminant analysis
}
\description{
Method partial least squares regression / Discriminant analysis. Uses
the NIPALS algorithm as described in H. Abdi, Wiley Interdisciplinary
Review: Computational Statistics, Vol. 2, pp. 97-106
}
\usage{
lims.pls.model(predictors, responses, nComp, training, xcenter=mean, xscale=sd, ycenter=mean, yscale=sd, discriminant=TRUE, accuracy=1e-05, maxit=100)

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
   \item{xcenter}{Function used for predictor variables centering. "mean" by default.}
   \item{xscale}{Function used for predictor variables scaling. "sd" by
     default.}
   \item{ycenter}{Function used for responses variables centering. "mean" by default.}
  \item{yscale}{Function used for response variables scaling. "sd" by
    default.}
   \item{discriminant}{TRUE (default) if discriminant analysis, FALSE if
    regression}
  \item{accuracy}{Maximum change between iteration cycles of the NIPALS
    algorithm. 0.00001 by default.}
  \item{maxit}{Maximum number of iteration cycles of the NIPALS
    algorithm. 100 by default.}

}

\value{
  Returns an object of the "model" class with the following elements:
  \item{xweights}{Predictor variables weights.}
  \item{xscores}{Predictor variables scores.}
  \item{xloadings}{Predictor variables loadings.}
  \item{yweights}{Response variables loadings.}
  \item{yscores}{Response variables scores.}
  \item{coefficients}{Coefficients for PLS prediction.}
  \item{fitted.values}{Fitted values of the PLS model.}
  \item{predicted.values}{Predicted values on the validation data set
    (original matrix minus observations used for training)}
  \item{Q2}{Q2 statistic}
  \item{R2}{Correlation coefficient for the fitted values (training
    set)}
  
}

\references{
\href{http://www.sciencedirect.com/science/article/pii/S0169743901001551}{Wold, Svante, M Sj\"ostr\"om, and L Eriksson. 2001. PLS-regression: a basic tool of chemometrics. Chemometrics and intelligent laboratory (58), 109 130}


}
\author{
Julien Wist, Andres Bernal, Jessica Medina

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{lims.opls.model}}

}
\examples{
####change levels####

library(caret)
data(coffee)
nmrData<-t(data$nmrData)
composition<-matrix(as.numeric(as.character(data$param$composition)),ncol=1)
model<-lims.pls.model(nmrData,composition,nComp=5,training=1:11)
plot(composition[1:11,],model$fitted.values[,,5])
legend("topleft", legend =paste(model$R2[5]))


}

\keyword{ pls }
\keyword{ regression }% __ONLY ONE__ keyword per line
\keyword{ discriminant analysis }
