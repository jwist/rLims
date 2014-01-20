
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
## lims.spectra allows to create a spectra object (S3) from a list of x and y data. 
## A mask can be attributed that will be used by plot methods and other to display
## a selected area of the spectra without replicating the data in a new variable
lims.spectra <- function(spec,type='nmr',mask=NULL) {
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
# s <- lims.spectra(fakeSpectra)

# example 2
# logfile <- 'lims.log'
# 
# entry <- lims.getJSON('http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=8026534&key=5quA7vBcoQ')
#
# spect <- read.table(sprintf("%s&filter=JcampToXY",entry[[1]]$nmrs[[1]]$resourceURL),sep=',',colClasses='numeric')
# spectra <- lims.spectra(spec)

##############################
## this methods returns a data.frame for spectra object in such a way that data
## can be displayed using the generic plot function
data.frame.spectra <- function(spectra,full=FALSE) {
	if (full) {
		return(data.frame(t(as.numeric(spectra$intensity))))
	} else {
		return(data.frame(t(as.numeric(spectra$intensity[spectra$mask]))))
	}
}

# example 
# fakeSpectra <- list(x=c(1:100),y=sin(seq(1,6.27,length.out=100)))
# s <- lims.spectra(fakeSpectra)
## will use the plot.spectra methods
# plot(s) 
## while this will use the generic plot function
# plot(as.numeric(data.frame.spectra(s)))

##############################
## this function retrieve the dimension of a spectra object
dim.spectra <- function(spectra) {
	return(length(spectra$intensity))
}

# example 
# fakeSpectra <- list(x=c(1:100),y=sin(seq(1,6.27,length.out=100)))
# s <- lims.spectra(fakeSpectra)
# dim(s)

##############################
## plot.spectra is the method to plot spectra object.
plot.spectra <- function(spectra,xlimit=NULL,ylimit=NULL) {
	y <- spectra$intensity[spectra$mask]
	x <- seq(spectra$xlims[1],spectra$xlims[2],along.with=spectra$intensity)[spectra$mask]
	
	if (is.null(xlimit)) {
		xlimit <- (range(x))
	} else {
		xlimit= (range(xlimit))
	}
	
	if (is.null(ylimit)) {
		ylimit <- range(y)
	} else {
		ylimit=range(ylimit)
	}
	
	switch(spectra$type, 
				 nmr={
				 	plot(x,(y),type='n',ylab='intensity',xlab='ppm',xlim=rev(xlimit), ylim=ylimit)
				 	lines(x,(y),col=1,xlim=rev(xlimit),ylim=ylimit)
				 }, 
				 ir={
				 	plot(x,(y),type='n',ylab='intensity',xlab='frequency',xlim=xlimit, ylim=ylimit)
				 	lines(x,(y),col=1,xlim=xlimit,ylim=ylimit)				 	
				 }, 
				 gcms=, 
				 'not defined yet')
	
}
# example
# plot(spectra,xlimit=c(9.4,9.7),ylimit=c(-1e4,3e5))

##############################
lims.spectra.setMask <- function(spectra,list=list()) {

	if (class(spectra) != 'spectra') {
		stop('Spectra is not of class spectra, please check!')
	} else {
		
		if (length(list) == 0) {
			warning('list is empty, no mask created')
		} else {
			
			ppm <- seq(spectra$xlims[1],spectra$xlims[2],along.with=spectra$intensity)
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
  #unclass(spectra)
	#assign("spectra[['mask']]",as.logical(I),envir = .GlobalEnv)
	#spectra$mask <<- as.logical(I)
	
	#return(as.logical(I))
	spectra$mask <- as.logical(I)
	return(spectra)
}

#example
#spectra <- lims.spectra(spec[1000:2000,])
# lims.spectra.setMask(spectra,list=c(8.4,8.7,8.8,9,9.5,9.57))
# plot(spectra)

# fake <- list(x=c(1:100),y=sin(seq(1,6.27,length.out=100)))
# s <- lims.spectra(fake)
# plot(s)
# lims.spectra.setMask(s,list=c(10,20,40,50,60,100))
# plot(s)

##############################
#### LIMS FUNCTIONS ####
##############################

##############################
lims.getJSON<-function(url,LOG=FALSE,...) {
	
	scriptName <- 'getJSON'
	argList<-list(...)
	
	if ( is.vector(url) & !is.list(url) ){
		if (length(url) == 0) {
			stop(paste('Please provide at least one url!'))
		} else {
			json_data <- lapply( seq_along(url),function(i) c(fromJSON(paste(readLines(url[i]), collapse=""))[[1]][[1]] ,'url'=url[i])) # retrieve JSON
			if (LOG) {
				msg <- paste("retrieved: ",length(json_data)," json" )
				lims.log(prefix='        ',scriptName=scriptName,msg=msg,file=argList$logfile)
			}
			return(json_data)	
		}
	} else {stop('Please provide url as a vector!')}

}

#example:
# url1 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=8026534&key=5quA7vBcoQ'
# url2 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=8405946&key=MdmsufXhXk'
# url3 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=6528763&key=BDiayNJUS8'
# url4 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=8316869&key=P5pzErJ490'
# url5 <- 'http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=6528214&key=NK2pYSFIhD'

# url <- c(url1,url2,url3,url4,url5)
# g1 <- lims.getJSON('http://mylims.univalle.edu.co/lims/Lims?action=GetJSON&table=entry&id=8026534&key=5quA7vBcoQ')
# entryList <- lims.getJSON(url)

##############################
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

##############################
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


##############################
lims.getNmrs <- function(entryList=entryList,entryData=entryData,Filter=list(),OP='AND',dry=FALSE,LOG=FALSE,...) {
	
	functionName <- 'lims.getNmrs'
	argList<-list(...)
	
	t <- list()
	s <- list()
	p <- list()
	info <- list()
	
	if (sum(is.na(match(names(Filter),names(entryList[[1]]$nmrs[[1]])))) != 0) {
		warning('lims.getNmrs: some names in parameter filter are not correct, please check')
	}
	
	for (i in 1:length(entryList)) {
		x <- entryList[[i]]
		
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
					s <- rbind(s, list('spectra'=lims.spectra(spect[1:(dim(spect)[1]/2),] ))) # taking only real part of spectrum
					if (LOG) {
						msg <- paste('spectra: ',x$entryID,' / ', y$resourceURL,sep='')
						lims.log(prefix='        ',scriptName=functionName,msg=msg,file=argList$logfile)
					}
					
				}
				# merge infos
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
	}
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
# nmrList <- lims.getNmrs(entryList,entryData,Filter=Filter,OP='AND')


##############################
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
			ifelse(!exists('nmrData'),nmrData <- as.numeric(data.frame.spectra(s[[i]])),
						 nmrData <- rbind(nmrData, as.numeric(data.frame.spectra(s[[i]]))))
		}
		
		ppm <- seq(s[[1]]$xlims[1],s[[1]]$xlims[2],along.with=s[[1]]$intensity)[s[[1]]$mask]
	}
	return(list('nmrData'=nmrData,'ppm'=ppm,'param'=nmrList$param,'nmrInfo'=nmrList$nmrInfo))
}
# example
# data <- lims.createDataSet(nmrList)



##############################


lims <- function(urlList,experimentList=list(),...) {

	# retrieve optional arguments
	argList<-list(...)
	
	today <- lims.getDate()$date
	logfile <- paste(today,'_lims.log',sep='')

	entryList <- lims.getJSON(urlList,LOG=FALSE)
	
	# new list of nmr url	
	resourceURLList <- sapply(entryList,function(x) x$nmrs[[1]]$resourceURL)
	
	if (!is.na(match('old',names(argList)))) {
		F <- match(resourceURLList,argList$old$nmrInfo$resourceURL)
	}
	
	entryData <- lims.getParameters(entryList)

	#spect <- read.table(sprintf("%s&filter=JcampToXY",entryList[[1]]$nmrs[[1]]$resourceURL),sep=',',colClasses='numeric')
	#Filter <- list('experiment'=c('zg30','noesygpps1dcomp'),'solvent'=c('C6D6','COFFEEmeoh'))
	Filter <- list('experiment'=experimentList)
	nmrList <- lims.getNmrs(entryList,entryData,Filter=Filter,OP='AND',dry=FALSE,LOG=TRUE,logfile=logfile)

	data <- lims.createDataSet(nmrList)

	return(data)
}

# lims.getNmrs <- function(entryList,filterList) {
# 	#"spectra"=read.table(sprintf("%s&filter=JcampToXY",y$resourceURL),sep=',',colClasses='numeric')
# 	t <- (lapply(entryList,function(x) unlist(lapply(x$nmrs,
# 				function(y) if (match(y$experiment,filterList[[1]]) != 0) {list("entryID"=x$entryID,
# 				"url"=I(y$resourceURL), "temperature"=y$temperature,
# 				"nucleus"=I(y$nucleus), "solvent"=I(y$solvent), "experiment"=I(y$experiment),
# 				#"spectra"=read.table(sprintf("%s&filter=JcampToXY",y$resourceURL),sep=',',colClasses='numeric'),
# 				"parameters"=ggg$params[ggg$params$entryID == x$entryID,]  )} ),recursive=FALSE) ))
# 	return(t)
# }
# #example
# t <- lims.getNmrs(g2,filterList=list(list('zg30',"noesygpps1dcomp"),list('C6D6','COFFEEmeoh')))

# 
# lims.getList <- function(y) {
# 	print(paste('ddd:',y[[1]],length(y[[1]])))
# 	if (length(y) != 0) {
# 		eval(parse(text=paste('data.frame("',y$description,'"=y$value)',sep='')))
# 	} else {
# 		data.frame("NA"=NA)
# 	}	
# }
# 
# lims.getListRec <- function() {
# 	if (!is.list(value)[[1]] == 1) {
# 		value <- paste(strsplit(g2[[1]]$keywords[[4]]$value,'\\.')[[1]],collapse='_')
# 		print(value)
# 		eval(parse(text=paste('data.frame("','d','"="',unlist(value[[1]][1]),'")',sep='')))
# 	} else if (length(value) == 0) {
# 		data.frame("NA"=NA)
# 	} else if (length(value) > 1) {
# 		value <- paste(strsplit(g2[[1]]$keywords[[4]]$value,'\\.')[[1]],collapse='_')
# 		#data.frame( sapply(seq_along(value),lims.getListRec) )
# 		data.frame("NA"=NA)
# 	}
# }


#sapply(seq_along(entryList), function(i) sapply(seq_along(entryList[[i]]), function(i) if(!is.list(entryList[[1]][i][[1]][1])) {eval(parse(text=paste('data.frame("',names(entryList[[1]][i]),'"="',I(entryList[[1]][i][[1]][1]),'")',sep='')))} ))


