\name{lims.validate}
\alias{lims.validate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
lims.validate
}
\description{
Retrieves validation models by Monte Carlo or K-fold
validation. Currently only tested for discriminant analysis with two classes.
}
\usage{
lims.validate(model, predictors, responses, method="randomp", nModels=100, nComp=1,
valsize=0.2, xcenter=mean, xscale=sd, ycenter=mean,
yscale=sd, codomain=c(-1,1), discriminant=TRUE, accuracy = 1e-05, maxit = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ A class model must be selected, this can be lims.pls.model or lims.opls.model
}
  \item{predictors}{ Matrix data containing the test and training
}
  \item{responses}{ Response matrix (y class)}
  \item{method}{ Validation method. Can be either "randomp" (default)
    for Monte Carlo validation (nModels are created with randomly
    sampled training set according to valsize) or "kfold" for K-fold
    validation (dataset is partitioned in valsize subsets and valsize
    models are built with each of these subsets as testing set) }
  \item{nModels}{ Number of models to be created, only used if
    method="randomp", 100 by default }
  \item{nComp}{ How many components are used for the model, 1 by
    default}
  \item{valsize}{ If method="randomp", the fraction of subjects to be
    used for testing. If method="kfold", the number of subsets in which
    the dataset is to be partitioned. } 
  \item{xcenter}{Function used for predictor variables centering. "mean" by default.}
   \item{xscale}{Function used for predictor variables scaling. "sd" by
     default.}
   \item{ycenter}{Function used for responses variables centering. "mean" by default.}
  \item{yscale}{Function used for response variables scaling. "sd" by
    default.}
  \item{codomain}{ Domain of the response variable. c(-1,1) by
    default. It is recommended to use the default in the current version
    of the package.}
  \item{discriminant}{TRUE (default) if discriminant analysis, FALSE if
    regression}
  \item{accuracy}{accuracy of method }
  \item{maxit}{ number of maximun iterations}
}

%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{

  \item{full.results}{list of n models, within n list are saved the model output such as: xscores, xloadings, 
  Q2, R2, fitted.values, xweights, etc}
  \item{Q2}{matrix of Q2 calculated nmodelx ncomp}
  \item{R2}{matrix of R2 calculated nmodelx ncomp}
  \item{trainig.set}{training set of the model}
%% ...
}
\references{

\href{http://www.sciencedirect.com/science/article/pii/S0169743901001551}{Wold, S., M Sjostrom, and L Eriksson. 2001. "PLS-regression: a basic tool of chemometrics". Chemometrics and intelligent laboratory (58), 109 130}

\href{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3337399/?tool=pmcentrez}{Szymanska, E.
Saccenti, E., Smilde, A. K., Westerhuis, J. A. 2012. "Double-check: validation of diagnostic statistics for PLS-DA models in metabolomics studies". Metabolomics (8), 3 16}


}
\author{
Andres Bernal
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See also \code{\link{lims.pls.model}} and \code{\link{lims.opls.model}}
}
\examples{

data(coffee)
nmrData<-t(data$nmrData)
species<-matrix(as.numeric(data$param$species),ncol=1)

nmr.validation <-lims.validate(model = lims.opls.model, predictors = nmrData, 
                                 responses = species, 
                                 nModels = 5, nComp=10, ycenter = TRUE, yscale = TRUE,
                                 xscale=FALSE,
                                 codomain = c(1,2))
  
plot(nmr.validation$full.results[[1]]$xscores, nmr.validation$full.results[[1]]$orthoscores)



}
\keyword{ validation }
\keyword{ kfold }% __ONLY ONE__ keyword per line
