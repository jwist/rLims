\name{lims.opls.model}
\alias{lims.opls.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Orthogonal Partial to Latent Structures Model
}
\description{
Method Orthogonal Partial to Latent Structures / Discriminant analysis. Uses
the NIPALS algorithm as described in Wold et al.2002.(119).128.
}

\usage{

lims.opls.model(predictors, responses, nComp, training, xcenter =
                 mean, xscale = sd, ycenter = mean, yscale = sd,
                 discriminant = TRUE, accuracy = 1e-05, maxit = 100)

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
   \item{maxit}{Maximum number of iteration cycles of the NIPALS
    algorithm. 100 by default.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  Returns an object of the "model" class with the following elements:
  \item{xweights}{Predictor variables weights.}
  \item{xscores}{Predictor variables scores.}
  \item{xloadings}{Predictor variables loadings.}
  \item{orthoyweights}{Response variables loadings.}
  \item{orthoscores}{Response variables scores.}
   \item{ortholoadings}{Predictor variables loadings.}
  \item{yweights}{Response variables loadings.}
  \item{yscores}{Response variables scores.}
  \item{coefficients}{Coefficients for PLS prediction.}
  \item{fitted.values}{Fitted values of the PLS model.}
  \item{predicted.values}{Predicted values on the validation data set
    (original matrix minus observations used for training)}
  \item{Q2}{Q2 statistic}
  \item{R2}{Correlation coefficient for the fitted values (training
    set)}
  \item{TSS}{Total Sum of Squares}
  \item{PRESS}{Predicted Residual Sum of Squares}
 
 
}
\references{
\href{http://doi.wiley.com/10.1002/cem.695}{Trygg, Johan, and Svante Wold. 2002. Orthogonal projections to latent structures (O-PLS). Journal of Chemometrics 16(3),119 128}
}

\author{
Julien Wist, Jessica Medina,Andres Bernal 
}

\seealso{
See package "rjson" and www.mylims.org. See also \code{\link{lims.pls.model}}

}
\examples{
library(caret)
data(coffee)
nmrData<-t(data$nmrData)
composition<-matrix(as.numeric(as.character(data$param$composition)),ncol=1)
model<-lims.opls.model(nmrData,composition,nComp=5,training=1:11)

}

\keyword{ opls }
\keyword{ pls }% __ONLY ONE__ keyword per line
