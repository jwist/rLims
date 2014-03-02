\name{lims.log}
\alias{lims.log}
\title{
A log function that display and append into a file 
}
\description{
Create list object that contains "user" and "date". The date comes with a random time stamp to create unique entries.
}
\usage{
lims.log(prefix='        ',scriptName='defaultScriptName',msg=msg,file=NULL)
}
\arguments{
  \item{prefix}{
	a prefix (string) of 8 characters that can be used to highlight certain log entries. For example, log during the init part of a script may be tagged with prefix="    INIT"
}
  \item{scriptName}{
  the name of the script or the scritps that execute the log command (usefull for debugging)
}
  \item{msg}{
  a string message
}
  \item{file}{
  a filename for the log
}
}
\details{
A unique name should be given to the log file. Therefore lims.getDate() will provide the date and a random stamp.
}
\references{
www.mylims.org
}
\author{
Julien Wist
}

\examples{

lims.log(msg='dfd','test.log')
## returns the following to the standard output
## [1] "Sun Feb 23 11:30:58 2014test.logdefaultScriptName: dfd "
## and creates the file test.log with the same enty line

}

\keyword{ lims }
\keyword{ nmr }
\keyword{ metabolomics }
\keyword{ log }
