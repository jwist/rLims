\name{getJSON}
\alias{getJSON}
\title{
retrieve json object from www.mylims.org database
}
\description{
retrieve a complete entry of www.mylims.org database.
}
\usage{
getJSON(url)
}
\arguments{
  \item{url}{
URL should be the complete URL that is obtained by launching pipeline from www.mylims.org
}
}
\details{
}
\value{
}
\references{
www.mylims.org
}
\author{
Julien Wist
}
\note{
}


\seealso{
}
\examples{

## The function is currently defined as
function (url) 
{
    print(paste("<b>*** extractJSON2.R:getJSON:started with ", 
        url, "</b>"))
    json_data <- fromJSON(paste(readLines(url), collapse = ""))
    print(paste("<b>*** extractJSON2.R:getJSON:finished dim for json_data", 
        length(json_data$entry[[1]]), "</b>"))
    return(json_data)
  }
}
\keyword{ lims }
\keyword{ metabolomics }% __ONLY ONE__ keyword per line
