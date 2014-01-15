pkgname <- "rLims"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('rLims')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("lims")
### * lims

flush(stderr()); flush(stdout())

### Name: lims
### Title: retrieve spectra from www.mylims.org database
### Aliases: lims
### Keywords: lims nmr metabolomics

### ** Examples


## use this command to read from a file
## listOfURL<-as.vector(unlist(read.table('in.txt',sep='')))
#
#url1 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=8405946&key=MdmsufXhXk'
#url2 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=8026534&key=5quA7vBcoQ'
#url3 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=6528763&key=BDiayNJUS8'
#url4 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=8316869&key=P5pzErJ490'
#url5 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=6528214&key=NK2pYSFIhD'

data(urlList)
dataSet <- lims(as.vector(unlist(urlList)),c('noesygpps1dcomp'))
## dataSet <- lims(c(url1,url2,url3,url4,url5),c('noesygpps1dcomp'))




### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
