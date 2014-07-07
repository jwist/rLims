
# #Make generic function
# # This is the "constructor" function...
# # ... UseMethod should have the name of the function!
# f <- function(x,...) UseMethod("f")
# 
# #Default method
# # ... The class name should be the same as the constructor
# f.default <- function(a,b=5,c=3,...){
# 	out <- a+b+c
# 	class(out) <- "f"
# 	out # must return the object out, not the class name!
# }
# 
# # Print method
# # The "f" part of "print.f" must be the same as the class!
# print.f <- function(x,...){
# 	cat("Result for f: ")
# 	print(unclass(x)) # Must unclass to avoid infinite recursion
# 	# NextMethod(x) # Alternative, but prints the class attribute...
# }
# 
# # Summary method
# # Should return a summary object (and not print it!)
# # Need a unique class for it ("fsummary")
# summary.f <- function(object,...){
# 	res <- object
# 	class(res) <- "fsummary"
# 	res
# }
# 
# # Now need to print the summary too:
# print.fsummary <- function(x, ...) {
# 	cat("f summary!\n")
# 	# Nice summary print goes here...
# }
# 
# # Plot method
# plot.f <-function(x,p=0.3,...){ cat("PLOTTING!\n") }
# 
# # Try it out:
# 
# x <- f(3)
# x # print x
# 
# y <- summary(x) # 
# y # print summary
# 
# plot(x)

##############################
#### GENERAL FUNCTIONS ####

lims.getDate <- function() {
	today <- Sys.Date()
	today <- format(today, format="%Y_%m_%d")
	stamp <- round(runif(1,0,10000)*10000000) # this is because Sys.time() is broken by ChemoSpec
	today <- paste(today,stamp,sep="_")
	user <- try(system("whoami", intern = TRUE))
	if (exists('user')) {
		res <- list("date"=today,"user"=user)
	} else {
		res <- list("date"=today,"user"='unknown')
	}
	return(res)
}

lims.log <- function(prefix='        ',scriptName='defaultScriptName',msg=msg,file=NULL) {
	if (is.null(file)) {
		file <- paste('log_',lims.getDate()$date,'.txt',sep='')
	}
	message <- paste(date(),prefix,scriptName,': ',msg,' ',sep='')
	print(message)
	write.table(message,file=file,append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)
}

lims.loopAppend <- function(var,dest) {
	
	#check if dest exists
	if ( (a <- mget( dest, envir=.GlobalEnv, ifnotfound=list(dest='not found') )) == 'not found' ) {
		assign( dest, var,  envir=.GlobalEnv )
	} else { # if exists
		assign( dest, append( eval(parse(text=dest)), var, length( eval(parse(text=dest)) ) ), envir=.GlobalEnv)
	}
}

##############################
#### SPECTRA CLASS AND METHODS ####
##############################

##############################
## lims.spectraCreator allows to create a spectra object (S3) from a list of x and y data. 
## A mask can be attributed that will be used by plot methods and other to display
## a selected area of the spectra without replicating the data in a new variable
lims.spectraCreator <- function(spec,type='nmr',mask=NULL) {
	x_min <- range(spec[[1]])[1]
	x_max <- range(spec[[1]])[2]
	y <- spec[[2]]
	
	lines=list() # todo retrieve nmrLines if exist on lims
	
	if (is.null(mask)) {
		mask  <- rep(TRUE,length.out=length(y))
	}
	
	switch(type,
				 nmr=spectra  <- list('type'=type,'intensity'=rev(y),'xlims'=range(spec[[1]]),'nmrLines'=lines,'mask'=mask),
				 ir=spectra  <- list('type'=type,'intensity'=rev(y),'xlims'=range(spec[[1]]),'irLines'=lines,'mask'=mask),
				 gcms=,
				 'not defined yet'
	)
	
	class(spectra) <- 'spectra'
	
	return(spectra)
}

# example 1
# fakeSpectra <- list(x=c(1:100),y=sin(seq(1,6.27,length.out=100)))
# s <- lims.spectraCreator(fakeSpectra)

# example 2
# logfile <- 'lims.log'
# 
# entry <- lims.JSON.getEntry('http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=8026534&key=5quA7vBcoQ')
#
# spect <- read.table(sprintf("%s&filter=JcampToXY",entry[[1]]$nmrs[[1]]$resourceURL),sep=',',colClasses='numeric')
# spectra <- lims.spectraCreator(spec)

##############################
## this methods returns a data.frame for spectra object in such a way that data
## can be displayed using the generic plot function
lims.spectra.toDataframe <- function(spectra,full=FALSE) {
	if (full) {
		return(data.frame(t(as.numeric(spectra$intensity))))
	} else {
		return(data.frame(t(as.numeric(spectra$intensity[spectra$mask]))))
	}
}

# example 
# fakeSpectra <- list(x=c(1:100),y=sin(seq(1,6.27,length.out=100)))
# s <- lims.spectraCreator(fakeSpectra)
## will use the plot.spectra methods
# plot(s) 
## while this will use the generic plot function
# plot(as.numeric(lims.spectra.toDataframe(s)))

##############################
## this function retrieve the dimension of a spectra object
dim.spectra <- function(x) {
	return(length(x$intensity))
}

# example 
# fakeSpectra <- list(x=c(1:100),y=sin(seq(1,6.27,length.out=100)))
# s <- lims.spectraCreator(fakeSpectra)
# dim(s)

##############################
## plot.spectra is the method to plot spectra object.
plot.spectra <- function(x,...) {
	
	spectra <- x ## to maintain the same arguments than generic function
	
	y <- spectra$intensity[spectra$mask]
	x <- seq(spectra$xlims[1],spectra$xlims[2],along.with=spectra$intensity)[spectra$mask]
	
	argList<-list(...)
	
	if ( sum( names(argList) == 'xlimit' ) == 0 ) {
		xlimit <- (range(x))
	} else {
		xlimit= (range(argList$xlimit))
	}
	
	if ( sum( names(argList) == 'ylimit' ) == 0 ) {
		ylimit <- range(y)
	} else {
		ylimit=range(argList$ylimit)
	}
	
	switch(spectra$type, 
				 nmr={
				 	plot(x,(y),type='n',ylab='intensity',xlab='ppm',xlim=rev(xlimit), ylim=ylimit)
				 	lines(x,(y),xlim=rev(xlimit),ylim=ylimit,...) ## the tree dots must go in order to allow for passing generic parameters
				 }, 
				 ir={
				 	plot(x,(y),type='n',ylab='intensity',xlab='frequency',xlim=xlimit, ylim=ylimit)
				 	lines(x,(y),xlim=xlimit,ylim=ylimit,...)				 	
				 }, 
				 gcms=, 
				 'not defined yet')
	
}
# example
# plot(spectra,xlimit=c(9.4,9.7),ylimit=c(-1e4,3e5))

##############################
lims.spectra.setMask <- function(spectraVarName,list=list()) {
	
	spectra <- get(spectraVarName,envir = .GlobalEnv)
	if (class(spectra) != 'spectra') {
		stop('Spectra is not of class spectra, please check!')
	} else {
		
		ppm <- seq(spectra$xlims[1],spectra$xlims[2],along.with=spectra$intensity)
		
		if (length(list) == 0) {
			
			warning('list is empty, full region considered')
			I <- factor(rep(TRUE,length(ppm)))
			
		} else {
			
			I <- factor(findInterval(ppm,list))
			
			# evaluate how to deal with found intervals 
			# todo some cases are not treated
			if (length(levels(I)) == length(list)+1) {
				print(1)
				levels(I) <- rep(c(FALSE,TRUE),times=length(list)/2+1)[-(length(list)+2)]
				print(sum(as.logical(I)))
			} else if (length(levels(I)) == length(list)) {
				if (list[1] <= spectra$xlims[1]) {
					print(21)
					levels(I) <- rep(c(TRUE,FALSE),times=length(list)-2)
				} else if (list[length(list)] >= spectra$xlims[2]){
					print(22)
					levels(I) <- rep(c(FALSE,TRUE),times=length(list)-2)	
				}
			} else if (length(levels(I)) == length(list)-1) {
				print(3)
				levels(I) <- rep(c(TRUE,FALSE),times=length(list)-2)[-length(list)]
			}
		}
	}
	
	spectra$mask <- as.logical(I)
	assign(spectraVarName,spectra,envir = .GlobalEnv)
	
}

#example
# spectra <- lims.spectraCreator(spec[1000:2000,])
# lims.spectra.setMask('spectra',list=c(8.4,8.7,8.8,9,9.5,9.57))
# plot(spectra)

# fake <- list(x=c(1:100),y=sin(seq(1,6.27,length.out=100)))
# s <- lims.spectraCreator(fake)
# plot(s)
# lims.spectra.setMask('s',list=c(10,20,40,50,80,100))
# plot(s)


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### LIMS FUNCTIONS ####
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lims.isArg <- function(arg,argList) {
	!is.na(match(arg,names(argList)))
}

## this function get all the entries in a user account and returns a list of URL that can be retrieved with
## the function lims.JSON.getEntry
lims.JSON.getUserEntries <- function(...) {
	
	argList<-list(...)
	
	if ( length(names(argList)) == 0 ) {
		stop("Not enough arguments for function lims.JSON.getUserEntries. Provides the url for the user")
		
  } else if ( lims.isArg('url',argList) ) {
		url <- argList$url
		
	} else if ( lims.isArg('account',argList) ) {
		if ( is.na(pmatch('http://',argList$account$server)) ) {
			url=paste('http://',argList$account$server,'/lims/default/sample/listJson.jsp?&queryUser=',argList$account$user,'&key=',argList$account$key,sep='')
		} else {
			url=paste(argList$account$server,'/lims/default/sample/listJson.jsp?&queryUser=',argList$account$user,'&key=',argList$account$key,sep='')
		}
	}
		
	#url=paste('http://mylims.univalle.edu.co/lims/default/sample/listJson.jsp?&queryUser=',user,'&key=PDfEZzKo3K',sep='')
	t <- fromJSON(paste(readLines(url), collapse=""))[[1]]
	#urlList <- sapply(t,function(x) x$entryDetails)
	urlList <- lapply(t,function(x) x$entryDetails)
	return(urlList)
}

## examples
#lims.JSON.getUserEntries(url='http://mylims.univalle.edu.co/lims/default/sample/listJson.jsp?&queryUser=jessica.medina@correounivalle.edu.co&key=PDfEZzKo3K')
#lims.JSON.getUserEntries(account=list(server='http://mylims.univalle.edu.co',user='jessica.medina@correounivalle.edu.co',key='PDfEZzKo3K'))
#lims.JSON.getUserEntries(account=list(server='mylims.univalle.edu.co',user='jessica.medina@correounivalle.edu.co',key='PDfEZzKo3K'))

##%%%%%%%%%%%%

## this function gets JSON from a list of entry URL
## it should be renamed lims.JSON.getEntry
lims.JSON.getEntry<-function(url,LOG=FALSE,...) {
	
	scriptName <- 'lims.JSON.getEntry'
	argList<-list(...)
	
	if (length(url) == 0) {
		stop(paste('Please provide at least one url!'))
		
	} else if ( is.vector(url) && !is.list(url) ) {
		#json_data <- lapply( seq_along(url),function(i) c(fromJSON(paste(readLines(url[i],warn=FALSE), collapse="")) ,'url'=url[i])) # retrieve JSON
		json_data <- lapply( seq_along(url),function(i) c(fromJSON(paste(readLines(url[i],warn=FALSE), collapse=""))[[1]][[1]] ,'url'=url[i])) # retrieve JSON
		
	} else if ( is.list(url) ) {
		#stop('Please provide url as a vector!')
		json_data <- lapply( seq_along(url),function(i) c(fromJSON(paste(readLines(url[[i]],warn=FALSE), collapse=""))[[1]][[1]] ,'url'=url[i])) # retrieve JSON
	}
	
	if (LOG) {
		msg <- paste("retrieved: ",length(json_data)," json" )
		lims.log(prefix='        ',scriptName=scriptName,msg=msg,file=argList$logfile)
	}	
	
	return(json_data)	
}

#example:
# url1 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=8026534&key=5quA7vBcoQ'
# url2 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=8405946&key=MdmsufXhXk'
# url3 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=6528763&key=BDiayNJUS8'
# url4 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=8316869&key=P5pzErJ490'
# url5 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=6528214&key=NK2pYSFIhD'

# url <- c(url1,url2,url3,url4,url5)
# g1 <- lims.JSON.getEntry('http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=8026534&key=5quA7vBcoQ')
# entryList <- lims.JSON.getEntry(url1)

##%%%%%%%%%%%%

## this function retrieves the whole entries for user.
lims.getUserEntries <- function(...) {
	argList <- list(...)
	entryList <- lims.JSON.getUserEntries(account=argList$account)
	lims.JSON.getEntry(entryList)
}

##%%%%%%%%%%%%

## this function returns the list of available resources for entryList. 
## Should be in list format ex[1] and not ex[[1]]
lims.listResources <- function(entryList) {
	
	if ( length(entryList) == 0 ) {
		stop("Please provide at leaste one entry")
		
	} else if ( length(entryList) == 1 ) {
		return(names(entryList[[1]]))
		
  } else {
  	dups <- length(unique(lapply(entryList,function(x) names(x) )))
  	if ( dups == 1 ) {
  		return(names(entryList[[1]]))
  	} else {
  		a = NULL
  		for ( i in 1:dups ) {
  			a <- unique(c(a,names(entryList[[i]])))
  		}
  	}
  	return(a)
  }
}

## examples
##lims.listResources(entryList[[38]])

##%%%%%%%%%%%%

lims.getResource <- function(entryList,resources=c("batchID","catalogID")) {
	
	fun <- function(x) {ifelse(is.null(x[i][[1]]),NA,x[i][[1]])}
	
	for ( i in resources ) {
		
		tmp_res <- unlist( lapply(entryList, fun) )
		if ( sum(is.na(tmp_res)) != 0 ) { warning(paste("NA found in resource: ",i)) }
		ifelse( !exists('resource'), resource <- tmp_res, resource <- cbind(resource,tmp_res) )
	}
	
	if ( length(resources) > 1) {
		colnames(resource) <- resources
	}
	
	return(resource)
}

##%%%%%%%%%%%%

lims.getAllResource <- function(entryList) {
	
	resources <- c("entryID","catalogID","batchID","creationDate","user")
		
	fun <- function(x) {ifelse(is.null(x[i][[1]]),NA,x[i][[1]])}
	
	for ( i in resources ) {
		
		tmp_res <- unlist( lapply(entryList, fun) )
		if ( sum(is.na(tmp_res)) != 0 ) { warning(paste("NA found in resource: ",i)) }
		ifelse( !exists('resource'), resource <- tmp_res, resource <- cbind(resource,tmp_res) )
	}
	
	if ( length(resources) > 1) {
		colnames(resource) <- resources
	}
	
	return(resource)
}

##%%%%%%%%%%%%

## this function checks if parameters exist and then returns their names
lims.getParametersNames <- function(jsonList) {
	msg <- 'Check your parameters, some have empty descriptions, entryID='
	#G <- sapply(jsonList,function(x) sapply(x$parameters, function(y) y$description))
	G <- sapply(jsonList,function(x) sapply(x$parameters, function(y) 
	if (y$description=="") {warning(paste(msg,x$entryID,sep='')); y$description} else {y$description} ))
	G <- unique(unlist(G))
	return(G) # returns a list of entries
}
#example:
# parameterNames <- lims.getParametersNames(entryList)

##%%%%%%%%%%%%

lims.getParameters <- function(entryList) {
	
	param <- lapply(entryList, function(x) if (length(x$parameters) > 0)
			{data.frame("entryID"=x$entryID, sapply(x$parameters, 
			function(y) if (y$description == "") { warning(paste('empty parameter description:',x$entryID)); eval(parse(text=paste('data.frame("',"unknown",'"=y$value)',sep=''))) }
			else {eval(parse(text=paste('data.frame("',I(y$description),'"=y$value)',sep='')))} ))} 
			else {data.frame("entryID"=x$entryID)})
	param <- Reduce(function(x,y) merge(x,y,all=TRUE),param)
	
	keywords <- lapply(entryList, function(x) if (length(x$keywords) > 0) 
			{data.frame("entryID"=x$entryID, "keywords"=sapply(x$keywords, 
			function(y) data.frame(y$value)))} else {data.frame("entryID"=x$entryID)})
	keywords <- Reduce(function(x,y) merge(x,y,all=TRUE),keywords)
	
	iupacs <- lapply(entryList, function(x) if (length(x$iupacs) > 0) 
			{data.frame("entryID"=x$entryID, "iupacs"=sapply(x$iupacs, 
			function(y) data.frame(y$value)))} else {data.frame("entryID"=I(x$entryID))})
	iupacs <- Reduce(function(x,y) merge(x,y,all=TRUE),iupacs)
	
	infos <- lapply(entryList, function(x) data.frame("entryID"=x['entryID'], "catalogID"=I(x$catalogID),
																										"batchID"=I(x$batchID), "creationDate"=I(x$creationDate),
																										"userEmail"=I(x$user$email),"lastModificationDate"=I(x$lastModificationDate)))
	infos <- Reduce(function(x,y) merge(x,y,all=TRUE),infos)
	
	# do some cleanup
	if (sum(is.na(param))) {
		warning('NA found in parameters, please check')
		for (i in 1:ncol(param)) {
			levels(param[,i]) <- c(levels(param[,i]),"")
			param[,i][is.na(param[,i])] <- ""
		}
	}
	
	for (i in 1:ncol(keywords)) {
		levels(keywords[,i]) <- c(levels(keywords[,i]),"")
		keywords[,i][is.na(keywords[,i])] <- ""
	}
	
	for (i in 1:ncol(iupacs)) {
		levels(iupacs[,i]) <- c(levels(iupacs[,i]),"")
		iupacs[,i][is.na(iupacs[,i])] <- ""
	}
	# returns entryData
	return(list("infos"=infos,"params"=param,"keywords"=keywords,"iupacs"=iupacs))
}
#example:
# entryData <- lims.getParameters(entryList)

# filter out index
#lapply(seq_along(g2), function(i) if (!is.na(match(i,index))) {g2[[i]]})

##%%%%%%%%%%%%

## this function organize all the data in a user account
## filtering the samples should be done here
lims.getUserSamplesMetaData <- function(...) {
	argList <- list(...)
	entries <- lims.getUserEntries(account=argList$account)
	entryData <- lims.getParameters(entries)
	
	lims.findExperimentNames(entries)
	lims.findSolventNames(entries)
	lims.findNucleusType(entries)
	lims.findByTemperature(entries)
	lims.findNumberOfSpectra(entries)
	
	## we have to sort the row because Reduce/Merge pair function in getParameters 
	## shuffles the rows
	F <- match(entryData$params$entryID,lims.getResource(entries,c("entryID")))
	entries <- entries[F]
	return(res=list("samples"=entries,"metaData"=entryData))
}

##%%%%%%%%%%%%

lims.selectSamples <- function(userData,index,exclude=FALSE) {
	userData$samples <- userData$samples[index]
	
	if ( exclude == TRUE ) { index <- -index }
	
	for ( i in c("infos","params","iupacs","keywords") ) {
		if (ncol(userData$metaData[i][[1]]) == 1) {
			userData$metaData[i][[1]] <- userData$metaData[i][[1]][[1]][index]
		} else {
			userData$metaData[i][[1]] <- userData$metaData[i][[1]][index,]
		}
	}
	
	return(userData)
}

##%%%%%%%%%%%%

lims.findExperimentNames <- function(entryList) {
	G <- lapply(entryList,function(x) lapply(x$nmrs, function(y) y$experiment))
	print(paste("number of spectra ordered by experiment type:"))
	print(table(unlist(G)))
	G <- unique(unlist(G))
	return(G)
}

##%%%%%%%%%%%%

lims.findSolventNames <- function(entryList) {
	G <- lapply(entryList,function(x) lapply(x$nmrs, function(y) y$solvent))
	print(paste("number of spectra ordered by solvents:"))
	print(table(unlist(G)))
	G <- unique(unlist(G))
	return(G)
}

##%%%%%%%%%%%%

lims.findNucleusType <- function(entryList) {
	G <- lapply(entryList,function(x) lapply(x$nmrs, function(y) y$nucleus))
	print(paste("number of spectra ordered by nucleus:"))
	print(table(unlist(G)))
	G <- unique(unlist(G))
	return(G)
}

##%%%%%%%%%%%%

lims.findByTemperature <- function(entryList) {
	G <- lapply(entryList,function(x) lapply(x$nmrs, function(y) y$temperature))
	print(paste("number of spectra ordered by temperature:"))
	print(table(unlist(G)))
	G <- unique(unlist(G))
	return(G)
}

##%%%%%%%%%%%%

lims.findNumberOfSpectra <- function(entryList) {
	G <- lapply(entryList,function(x) lapply(x$nmrs, function(y) y$resourceURL))
	G <- length(unique(unlist(G)))
	print(paste("total number of spectra: ",G))
	return(G)
}

##%%%%%%%%%%%%

lims.getNmrs <- function(entryList=entryList,...) {
	
	functionName <- 'lims.getNmrs'
	argList<-list(...)
	
	if (!is.na(match('LOG',names(argList)))) {
		LOG <- argList$LOG
	} else {
		LOG <- FALSE
	}
	
	if (!is.na(match('dry',names(argList)))) {
		dry <- argList$dry
	} else {
		dry <- FALSE
	}
	
	if (!is.na(match('OP',names(argList)))) {
		OP <- argList$OP
	} else {
		OP <- 'AND'
	}
	
	if (!is.na(match('Filter',names(argList)))) {
		Filter <- argList$Filter
	} else {
		Filter <- list()
	}
	
	t <- list()
	s <- list()
	p <- list()
	info <- list()
	
	if (sum(is.na(match(names(Filter),names(entryList[[1]]$nmrs[[1]])))) != 0) {
		warning('lims.getNmrs: some names in parameter filter are not correct, please check')
	}
	
	if (length(entryList) > 0) {
		for (i in 1:length(entryList)) {
			x <- entryList[[i]]
			if (length(x$nmrs) > 0) {
				for (j in 1:length(x$nmrs)) {
					y <- x$nmrs[[j]]
					
					z <- NULL
					for (i in 1:length(Filter)) {
						zt <- paste('!is.na(match(y[names(Filter)[',i,']],Filter[[',i,']]))',sep='')
						if (is.null(z)) {
							z <- zt
						} else {
							if (OP == 'AND') {
								z  <- paste(z,zt,sep=' & ')
							} else if (OP == 'OR') {
								z  <- paste(z,zt,sep=' | ')
							}
						}
					}
					
					if (eval(parse(text=z))) { # test for filter
						
						# retrieve spectra and store them as a list to avoid dimension problems
						if (!dry) {
							spect <- read.table(sprintf("%s&filter=JcampToXY",y$resourceURL),sep=',',colClasses='numeric')
							s <- rbind(s, list('spectra'=lims.spectraCreator(spect[1:(dim(spect)[1]/2),] ))) # taking only real part of spectrum
							if (LOG) {
								msg <- paste('spectra: ',x$entryID,' / ', y$resourceURL,sep='')
								lims.log(prefix='        ',scriptName=functionName,msg=msg,file=argList$logfile)
							}
							
						}
						# merge infos
						if (!is.na(match('entryData',names(argList)))) {
							entryData = argList$entryData
							p <- rbind(p, entryData$params[entryData$params$entryID == x$entryID,])
							info  <- rbind(info, entryData$info[entryData$info$entryID == x$entryID,])
							t <- rbind(t,list(
								'entryID'=x$entryID,
								'experiment'=y$experiment,
								'temperature'=y$temperature,
								'nucleus'=y$nucleus,
								'solvent'=y$solvent,	
								'resourceURL'=y$resourceURL
							))
						}
					}
				}} #nmrs
		}} #entry
	
	# returns nmrData
	if (!dry) {
		return(list('spectra'=s,'nmrInfo'=data.frame(t),'param'=data.frame(p),'info'=data.frame(info)))
	} else {
		return(list('nmrInfo'=data.frame(t),'param'=data.frame(p),'info'=data.frame(info)))
	}
}
#example
# Filter <- list('experiment'=c('zg30','noesygpps1dcomp'),'solvent'=c('C6D6','COFFEEmeoh'))
# Filter <- list('experiment'=c('noesygpps1dcomp'),'solvent'=c('COFFEEmeoh'))
# nmrList <- lims.getNmrs(entryList,entryData=entryData,Filter=Filter,OP='AND')

##%%%%%%%%%%%%

lims.getSpectra <- function(metaData=metaData,...) {
	
	entryList <- metaData$samples
	entryData <- metaData$metaData
	
	functionName <- 'lims.getSpectra'
	argList<-list(...)
	
	if (!is.na(match('LOG',names(argList)))) {
		LOG <- argList$LOG
	} else {
		LOG <- FALSE
	}
	
	if (!is.na(match('dry',names(argList)))) {
		dry <- argList$dry
	} else {
		dry <- FALSE
	}
	
	if (!is.na(match('OP',names(argList)))) {
		OP <- argList$OP
	} else {
		OP <- 'AND'
	}
	
	if (!is.na(match('Filter',names(argList)))) {
		Filter <- argList$Filter
	} else {
		Filter <- list()
	}
	
	t <- list()
	s <- list()
	p <- list()
	info <- list()
	
	if (sum(is.na(match(names(Filter),names(entryList[[1]]$nmrs[[1]])))) != 0) {
		warning('lims.getSpectra: some names in parameter filter are not correct, please check')
	}
	
	if (length(entryList) > 0) {
		for (i in 1:length(entryList)) {
			x <- entryList[[i]]
			if (length(x$nmrs) > 0) {
				for (j in 1:length(x$nmrs)) {
					y <- x$nmrs[[j]]
					
					z <- NULL
					for (i in 1:length(Filter)) {
						zt <- paste('!is.na(match(y[names(Filter)[',i,']],Filter[[',i,']]))',sep='')
						if (is.null(z)) {
							z <- zt
						} else {
							if (OP == 'AND') {
								z  <- paste(z,zt,sep=' & ')
							} else if (OP == 'OR') {
								z  <- paste(z,zt,sep=' | ')
							}
						}
					}
					
					if (eval(parse(text=z))) { # test for filter
						
						# retrieve spectra and store them as a list to avoid dimension problems
						if (!dry) {
							spect <- read.table(sprintf("%s&filter=JcampToXY",y$resourceURL),sep=',',colClasses='numeric')
							s <- rbind(s, list('spectra'=lims.spectraCreator(spect[1:(dim(spect)[1]/2),] ))) # taking only real part of spectrum
							if (LOG) {
								msg <- paste('spectra: ',x$entryID,' / ', y$resourceURL,sep='')
								lims.log(prefix='        ',scriptName=functionName,msg=msg,file=argList$logfile)
							}
							
						}
						# merge infos

						p <- rbind(p, entryData$params[entryData$params$entryID == x$entryID,])
						info  <- rbind(info, entryData$infos[entryData$infos$entryID == x$entryID,])
						t <- rbind(t,list(
							'entryID'=x$entryID,
							'experiment'=y$experiment,
							'temperature'=y$temperature,
							'nucleus'=y$nucleus,
							'solvent'=y$solvent,	
							'resourceURL'=y$resourceURL
						))
						
					}
				}} #nmrs
		}} #entry
	
	# returns nmrData
	if (!dry) {
		return(list('spectra'=s,'nmrInfo'=data.frame(t),'param'=data.frame(p),'info'=data.frame(info)))
	} else {
		return(list('nmrInfo'=data.frame(t),'param'=data.frame(p),'info'=data.frame(info)))
	}
}
#example

##%%%%%%%%%%%%

lims.findSpectraSize <- function(data) {
	return(sapply(data$spectra,function(x) length(x$intensity)))
}

##%%%%%%%%%%%%

# this function creates data ready for analysis
# 1) flatten spectra into spectra objects
# 2) add relevant information
lims.createDataSet <- function(nmrList) {
	# check for dimension problem
	t <- unlist(lapply(nmrList$spectra,function(x) dim(x)))
	if (min(t) != max(t)) {
		stop('Cannot create data set, check size of spectra')
	} else {
		
		s <- nmrList[['spectra']]
		
		for (i in 1:nrow(s)) {
			ifelse(!exists('nmrData'),nmrData <- as.numeric(lims.spectra.toDataframe(s[[i]])),
						 nmrData <- rbind(nmrData, as.numeric(lims.spectra.toDataframe(s[[i]]))))
		}
		
		ppm <- seq(s[[1]]$xlims[1],s[[1]]$xlims[2],along.with=s[[1]]$intensity)[s[[1]]$mask]
	}
	return(list('nmrData'=nmrData,'ppm'=ppm,'param'=nmrList$param,'nmrInfo'=nmrList$nmrInfo))
}
# example
# data <- lims.createDataSet(nmrList)

##%%%%%%%%%%%%

lims <- function(urlList,experimentList=list(),...) {
	
	# retrieve optional arguments
	argList<-list(...)
	
	today <- lims.getDate()$date
	logfile <- paste(today,'_lims.log',sep='')
	
	entryList <- lims.JSON.getEntry(urlList,LOG=FALSE)
	
	# this will produce a warning if empty parameter description (not necessary)
	#lims.getParametersNames(entryList)
	
	# new list of nmr url	
	resourceURLList <- sapply(entryList,function(x) x$nmrs[[1]]$resourceURL)
	
	if (!is.na(match('old',names(argList)))) {
		F <- match(resourceURLList,argList$old$nmrInfo$resourceURL)
	}
	
	entryData <- lims.getParameters(entryList)
	
	#spect <- read.table(sprintf("%s&filter=JcampToXY",entryList[[1]]$nmrs[[1]]$resourceURL),sep=',',colClasses='numeric')
	#Filter <- list('experiment'=c('zg30','noesygpps1dcomp'),'solvent'=c('C6D6','COFFEEmeoh'))
	Filter <- list('experiment'=experimentList)
	nmrList <- lims.getNmrs(entryList,entryData=entryData,Filter=Filter,OP='AND',dry=FALSE,LOG=TRUE,logfile=logfile)
	
	data <- lims.createDataSet(nmrList)
	
	return(data)
}

#install.packages('devtools')
# library(devtools)
# install_github("jwist/rLims")
# library(rLims)
# url='http://mylims.univalle.edu.co/lims/default/sample/listJson.jsp?&queryUser=jessica.medina@correounivalle.edu.co&key=PDfEZzKo3K'
# account <- list(server='mylims.univalle.edu.co',user='jessica.medina@correounivalle.edu.co',key='PDfEZzKo3K')
# entryList <- lims.JSON.getUserEntries(url=url)
# entries <- lims.JSON.getEntry(entryList)
# entries <- lims.getUserEntries(account=account)
# entryData <- lims.getParameters(entries) ## parameters are not ordered as sample
# ## the following two lines are for ordering entries according to entryData
# F <- match(entryData$params$entryID,lims.getResource(entries,c("entryID")))
# entries <- entries[F]
# userData <- list("samples"=entries,"metaData"=entryData)
# userData <- lims.getUserSamplesMetaData(account=account)
# Filter <- list('experiment'=c('noesygpps1dcomp'))
# nmrList <- lims.getNmrs(entries,entryData=entryData,Filter=Filter,OP='AND',dry=FALSE)
# 
# 
# url='http://mylims.univalle.edu.co/lims/default/sample/listJson.jsp?&queryUser=jessica.medina@correounivalle.edu.co&key=PDfEZzKo3K'
# metaData <- lims.getUserSamplesMetaData(url=url)
# 
# 
# account <- list(server='mylims.univalle.edu.co',user='jessica.medina@correounivalle.edu.co',key='PDfEZzKo3K')
# metaData <- lims.getUserSamplesMetaData(account=account)
# metaData <- lims.selectSamples(metaData,c(6:8))
# lims.findExperimentNames(metaData$samples)
# Filter <- list('experiment'=c('noesygpps1dcomp'))
# data <- lims.getSpectra(metaData, Filter=Filter,OP='AND',dry=FALSE)
# lims.findSpectraSize(data) ## check size of spectra before creating the dataset
# dataSet <- lims.createDataSet(data)


#### DISPLAY FUNCTIONS ####


diffplot <- function(index=ppm,data=nmrData,png=NULL,...) {
	
	# get arguments 
	args <- as.list( sys.call() ) # get names
	argList<-list(...) # get args
	print(nrow(data))
	ifelse(is.null(args$main), title <- "" , title <- argList$main)
	ifelse(is.null(args$sub), caption <- "" , caption <- argList$sub)
	if (is.null(args$xlab)) { args$xlab <- "Chemical Shift [ppm]" }
	if (is.null(args$ylab)) { args$ylab <- "Relative Intensity" }
	if (is.null(args$col)) {
		if (is.null(nrow(data))) { 
			group <- 1 
		} else {
			group <- (c(1:nrow(data)) %% 9 ) + 1
		}
	} else {
		group = argList$col
	}
	
	# print
	ifelse(is.null(args$xlim), limitsx <- range(ppm) , limitsx <- argList$xlim)
	
	#  png(paste("./images/noise_",log,"_",param$name[i],".png",sep=""), bg="transparent", width=700, height=300)
	F <- index > min(limitsx) & index < max(limitsx)
	
	if (is.null(nrow(data))) {
		ifelse(is.null(args$ylim), limitsy <- range(c(min(data[F]),max(data[F])*1.3)) , limitsy <- argList$ylim)
		matplot(index[F],data[F], type="l",ylim=limitsy,xlab=args$xlab,ylab=args$ylab,col=group,xlim=rev(range(index[F])),main=title,sub=caption,cex.sub=0.7,cex.main=0.8)
	} else {
		ifelse(is.null(args$ylim), limitsy <- range(c(min(data[,F]),max(data[,F])*1.3)) , limitsy <- argList$ylim)
		matplot(index[F],t(data[,F]), type="l",ylim=limitsy,xlab=args$xlab,ylab=args$ylab,col=group,xlim=rev(range(index[F])),main=title,sub=caption,cex.sub=0.7,cex.main=0.8)
	}
	return(res=list('xlim'=limitsx,'ylim'=limitsy))
}

saveDevice <- function (fn) {
	dev.copy(device= pdf, file=paste(fn, ".pdf", sep=""),
					 width=7, height =5)
	dev.off()
	dev.copy(device= png, bg="transparent", file=paste(fn, ".png", sep=""),
					 width=1000, height = 400)
	dev.off()
}

# Examples
# region <- c(0.65,1)
# samples <- c(1,2,6)
# color <- c('red','green','black')
# main <- "Dominican Republic Study"
# caption <- ""
# par(mar=c(4.5,4.5,3,1))
# 
# res<-diffplot(limits,index=ppm,data=nmrData[samples,],xlim=region,col=color, sub=caption, main=main)
# 
# legend_txt  <- paste(color,": ",param$name[samples],param$name[samples],' ')
# legend(max(res$xlim),max(res$ylim),legend_txt,cex=0.5)
# 
# saveDevice("./images/a1")



##############################
# ## new function to do it simpler (not in use)
# lims.get <- function(entryList=entryList) {
# 	
# 	n <- list()
# 	for (i in 1:length(entryList)) {
# 		n <- c(n, unlist(setdiff(names(entryList[[i]]),n)))
# 	}
# 	
# 	for (i in 1:length(n)) {
# 		for (j in 1:length(entryList)) {
# 			if (length(entryList[[j]][[i]]) == 1) {
# 				if (length(entryList[[j]][[i]][[1]]) == 1) {
# 					lims.loopAppend(unlist(entryList[[j]][[i]]),'res')
# 				} else {
# 					lims.loopAppend("list",'res')
# 				}
# 			} else {
# 				lims.loopAppend(NA,'res')
# 			}
# 		}
# 	}
# 	
# 	t <- as.matrix(res)
# 	dim(t) <- c(length(entryList),length(n))
# 	t <- data.frame(t)
# 	colnames(t) <- n
# 	return(t)
# }
# #remove(res)
# #lims.get(entryList)




