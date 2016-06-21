###this is taken from chemospec !! ########


isWholeNo <- function(x, tol = .Machine$double.eps^0.5)  {
  
  # Taken from the help to is.integer()
  # Bryan Hanson, DePauw Univ, Nov 2009
  # Used by binSpec in ChemoSpec
  
  abs(x - round(x)) < tol
  
}


check4Gaps<-function (x, y = NULL, tol = 0.01, plot = FALSE, silent = FALSE, 
                      ...) 
{
  len.x <- length(x)
  p <- x[2] - x[1]
  d1 <- c()
  d2 <- c()
  d1i <- c()
  d2i <- c()
  dend <- data.frame(d2 = NA, d2i = NA)
  dbeg <- data.frame(d1 = NA, d1i = NA)
  for (n in 1:(len.x - 1)) {
    t <- x[n + 1] - x[n]
    if (!isTRUE(all.equal(t, p, tolerance = tol))) {
      dend <- rbind(dend, data.frame(d2 = x[n], d2i = n))
      dbeg <- rbind(dbeg, data.frame(d1 = x[n + 1], d1i = n + 
                                       1))
    }
  }
  chk <- dim(dend)[1]
  if (chk > 2) {
    dend <- dend[-1, ]
    dbeg <- dbeg[-1, ]
    tmp <- data.frame(beg.freq = NA, end.freq = NA, size = NA, 
                      beg.indx = NA, end.indx = NA)
    first <- data.frame(beg.freq = x[1], end.freq = dend[1, 
                                                         1], size = NA, beg.indx = 1, end.indx = dend[1, 2])
    last <- data.frame(beg.freq = dbeg[nrow(dbeg), 1], end.freq = x[length(x)], 
                       size = NA, beg.indx = dbeg[nrow(dbeg), 2], end.indx = length(x))
    for (n in 1:(nrow(dbeg) - 1)) {
      tmp <- rbind(tmp, data.frame(beg.freq = dbeg[n, 1], 
                                   end.freq = dend[n + 1, 1], size = NA, beg.indx = dbeg[n, 
                                                                                         2], end.indx = dend[n + 1, 2]))
    }
    tmp <- tmp[-1, ]
    df <- rbind(first, tmp, last)
    df[, 3] <- abs(df[, 2] - df[, 1])
  }
  if (chk == 2) {
    dend <- dend[-1, ]
    dbeg <- dbeg[-1, ]
    first <- data.frame(beg.freq = x[1], end.freq = dend[1, 
                                                         1], size = NA, beg.indx = 1, end.indx = dend[1, 2])
    last <- data.frame(beg.freq = dbeg[nrow(dbeg), 1], end.freq = x[length(x)], 
                       size = NA, beg.indx = dbeg[nrow(dbeg), 2], end.indx = length(x))
    df <- rbind(first, last)
    df[, 3] <- abs(df[, 2] - df[, 1])
  }
  if ((chk == 1) && (!silent)) 
    cat("No gaps were found by check4Gaps\nNo plot will be made\n")
  if (chk == 1) 
    df <- FALSE
  if ((chk > 1) && (plot)) {
    if (missing(y)) 
      stop("No y values provided; cannot plot!")
    plot(x, y, type = "l", col = "black", main = "Gaps in Frequency Data", 
         ylab = "", xlab = "marked regions are skipped in data set", 
         ...)
    ybottom <- min(y) - 0.1 * diff(range(y))
    ytop <- max(y) + 0.1 * diff(range(y))
    for (n in 1:(nrow(df) - 1)) {
      lines(x = c(df[n, 2], df[n + 1, 1]), y = c(y[df[n, 
                                                      5]], y[df[n + 1, 4]]), lty = 2, col = "white")
      rect(xleft = df[n, 2] + p, ybottom, xright = df[n + 
                                                        1, 1] - p, ytop, density = 15, col = "pink")
    }
  }
  df
}



binData<-function (x = NULL, y = NULL, bin.ratio = 2) 
{
  if (bin.ratio <= 1) 
    stop("bin.ratio must > 1")
  if (!isWholeNo(bin.ratio)) 
    stop("bin.ratio must be an integer > 1")
  if (!is.null(y) && !is.null(x)) {
    if (!identical(length(x), length(y))) 
      stop("x and y vectors in binData have different lengths")
  }
  chk <- check4Gaps(x, silent = TRUE)
  if (length(chk) > 1) 
    stop("The data being binned has gaps and cannot be binned accurately")
  br <- bin.ratio
  if (!is.null(x)) 
    len <- length(x)
  if (!is.null(y)) 
    len <- length(y)
  no.bins <- len/br
  if (!isWholeNo(no.bins)) {
    chop <- NULL
    n <- 0
    while (n < br) {
      n <- n + 1
      l <- len - n
      no.b <- l/br
      if (isWholeNo(no.b)) {
        chop <- n
        break
      }
    }
    rem <- c(1:chop)
    if (!is.null(x)) 
      x <- x[-rem]
    if (!is.null(y)) 
      y <- y[-rem]
    if (!is.null(x)) 
      len <- length(x)
    if (!is.null(y)) 
      len <- length(y)
    no.bins <- len/br
  }
  b.x <- c(rep(NA, no.bins))
  b.y <- c(rep(NA, no.bins))
  cnt <- seq(1, len, br)
  inc <- br - 1
  if (!is.null(x)) {
    for (n in 1:no.bins) {
      r <- c(cnt[n]:(cnt[n] + inc))
      b.x[n] <- mean(x[r])
    }
  }
  if (!is.null(y)) {
    for (n in 1:no.bins) {
      r <- c(cnt[n]:(cnt[n] + inc))
      b.y[n] <- sum(y[r])
    }
  }
  if (!is.null(y)) 
    res <- data.frame(sum.y = b.y)
  if (!is.null(x)) 
    res <- data.frame(mean.x = b.x)
  if (!is.null(y) && !is.null(x)) 
    res <- data.frame(mean.x = b.x, sum.y = b.y)
  res
}


####################################################################################################################

#### extract intensity ###########################

nmr.bin<-function(ppm,nmrData,Nbin){
  
  bin = Nbin
  for (K in 1:dim(nmrData)[2]){
    z <- binData(seq(1:dim(nmrData)[1]),nmrData[,K],bin)
    zshift <- seq(to=ppm[length(ppm)], from=ppm[1], length.out=dim(z)[1])
    b <- data.frame(I(zshift))
    b <- cbind(b,z$sum.y)
    colnames(b) <- c("shift","intensity")
    ifelse(K == 1, binned.data<-data.frame(I(b)),binned.data<-cbind(binned.data,I(b)))
  }
  res <- list("binned.data"=binned.data)
  return(res)
}

extract.nmr.bin<-function(binned.data){
  
  dataBin<-matrix(ncol=dim(binned.data)[2], nrow=dim(binned.data)[1])
  for (i in 1:dim(binned.data)[2]){
    dataBin[,i]<- binned.data[,i]$intensity
  }
  return(dataBin)
}


######cut data #########
nmr.crop <- function(lowerBound,upperBound,nmrData,ppm){
     nmr.cropped <- nmrData[,lowerBound:upperBound]
     ppm.cropped<- ppm[lowerBound:upperBound]
    return(list("nmr"=nmr.cropped,"ppm"=ppm.cropped))
}

spect.slim <- function(lowerBound,upperBound,nmr,ppm){
  F <- ppm< lowerBound | ppm > upperBound
  ppm <- ppm[F]
  nmrData<- nmr[,F]
  spect <- list('nmrData'=nmrData,'ppm'=ppm)
  return(spect)
}


check.data <- function(lowerBound,upperBound,nmr,ppm, param){
  F <- ppm> lowerBound & ppm < upperBound
  ppm <- ppm[F]
  nmrData<- nmr[,F]
  #spect <- list('nmrData'=nmrData,'ppm'=ppm)
  pca<-prcomp(nmrData, center=TRUE,scale=TRUE)
  
  #std <- lapply(components,function(comp){
  #comppca <- princomp(cbind(x[comp], y[comp]));
  #commppca<-princomp(t(pca$Scores), scale=TRUE, center=TRUE)   
  #std<- list(
  #      "w" = comppca$sdev[[1]] * 2, #shouldn't be 4? ("within 2 stds of the mean")
  #     "h"= comppca$sdev[[2]] * 2,
  #    "x" = comppca$center[[1]],
  #   "y" = comppca$center[[2]],
  #  "a" = (atan(comppca$loading[[1,2]]/comppca$loading[[1,1]])*(180/pi))
  #  * -(comppca$loading[[1,1]] * comppca$loading[[2,2]] - comppca$loading[[2,1]] * comppca$loading[[1,2]])
  # );
  
  sc<-apply(pca$x[,1:2],1,sd)
  ms<-mean(apply(pca$x[,1:2],1,sd))
  o<-sc<2*ms
  outliers<-list("entryID"=param$entryID[!o])
  nmrData<-nmr[o,]
  
  
  # nmr.noiseStdev<-apply(nmrData,2,sd) 
  return(list("listOutliers"=outliers,"outliers"=nmr[!o,], "nmrData"=nmrData, "param"=param[o,]))
}



###################################################################################################################


##check reference #######


#for (K in c(1:10)){
#	#I = K*2-1
#	#F <- data[,I] > A & data[,I] < B
#	F <- data[,K]$shift > A & data[,K]$shift < B
#	#chunk <- list(shift=data[F,I],intensity=data[F,I+1])
#	chunk <- list(shift=data[,K]$shift[F],intensity=data[,K+1]$intensity[F])
#	peak <- pp(chunk,3,noiseLevel)
#	peak$integral <- sum(chunk$intensity)/length(chunk$intensity)
	##plot(range(chunk$shift),range(chunk$intensity),type="n")
	##lines(chunk$shift,chunk$intensit, col=K)
	##points(peak$shift,peak$intensity,col=K)
	#ifelse(K == 1, trigoneline <- peak, trigoneline <- cbind(trigoneline, peak))
	#ifelse(K == 1, Int <- max(chunk$intensity), ifelse(Int < max(chunk$intensity),Int <- max(chunk$intensity), Int))
	# deriv <- Deriv1(chunk$shift,chunk$intensity)
	# peakd <- pp(deriv,3,3e9)
	#deriv$intensity <- -deriv$intensity
	 #valleyd <- pp(deriv,3,3e9)
	 #d <- mean(peakd$shift)-mean(valleyd$shift)
	 #center <- (mean(peakd$shift)+mean(valleyd$shift))/2
	 #fit1 <- nls(intensity~(a/b)*exp(-(shift-c)^2/(2*b^2)),data=chunk, start=list(a=(1/sqrt(2*pi)) / d, b=d, c=center),trace=TRUE,control=nls.contr(maxiter = 50, tol = 		1e-05, minFactor = 1/4096,printEval = FALSE, warnOnly = FALSE))
	 #plot(chunk$shift,chunk$intensity)
	 #lines(chunk$shift,(coef(fit1)[1]/coef(fit1)[2])*exp(-(chunk$shift-coef(fit1)[3])^2/(2*coef(fit1)[2]^2)))
	 #lines(chunk$shift,chunk$intensity,col=2)
	 #T <- abs(2*sqrt(2*log(2))*coef(fit1)[2])*400
	 #ifelse(K == 1, FWHM <- T,FWHM <- cbind(FWHM, T))
	 #T <- coef(fit1)[3]
	 #ifelse(K == 1, truePeakPosition <- T,truePeakPosition <- cbind(truePeakPosition,T))
	 #T <- (coef(fit1)[1]/coef(fit1)[2])*exp(0)
	 #ifelse(K == 1, truePeakHeight <- T,truePeakHeight <- cbind(truePeakHeight,T))
	 #T <- sum((coef(fit1)[1]/coef(fit1)[2])*exp(-(chunk$shift-coef(fit1)[3])^2/(2*coef(fit1)[2]^2)))/length(chunk$shift)
	 #ifelse(K == 1, truePeakIntensity <- T,truePeakIntensity <- cbind(truePeakIntensity,truePeakIntensity=T))
	 
#}
#colnames(trigoneline) <- seq(1,K)
#hist(FWHM[truePeakIntensity > 4.5e5],6,main='Full Width at Half Maximum for Trigoneline Peak', xlab='Hz', ylab='number of spectra')


#####################################################################################################

#checkData<-function(nmr,param,date){
#	stamp <- round(runif(1,0,100)*1000)
#	postscript(paste(path,'/R/preprocess/', date$date,'-',stamp,'_prepData2.checkData.ps',sep=''),"ps")
#	loopCounter <- 0
#	repeat{
#		loopCounter <- loopCounter + 1
		
#		nmr.noise.sd.left <- nmr.noiseStdev(nmr,1, 10500)
#		nmr.noise.sd.right <- nmr.noiseStdev(nmr,115000,dim(nmr)[1])
		
#		distribution.radius <- 2*mean(apply(nmr.noise.sd.left,2,sd),apply(nmr.noise.sd.right,2,sd))
#		distribution.x <- apply(nmr.noise.sd.left,2,mean) # calculate center of circle (x)
#		distribution.y <- apply(nmr.noise.sd.right,2,mean) # calculate center of circle (x)
		
		## calculate distance to center of circle
#		nmr.noise.sd.distance <- apply(cbind(nmr.noise.sd.left,nmr.noise.sd.right),1,function(v) sqrt((v[1]-distribution.x)^2+(v[2]-distribution.y)^2))
		
	#	F <- nmr.noise.sd.distance < distribution.radius
 	#	plot(nmr.noise.sd.left,nmr.noise.sd.right)
	#	points(nmr.noise.sd.left[F],nmr.noise.sd.right[F],pch=20)
 #		textxy(nmr.noise.sd.left[!F],nmr.noise.sd.right[!F],param$catalogID[!F], dcol=2)
	#	textxy(nmr.noise.sd.left[F],nmr.noise.sd.right[F],param$catalogID[F], dcol=1)
	#	draw.circle(distribution.x,distribution.y,distribution.radius,nv=100,border=NULL,col=NA,lty=1,lwd=1)

		## delete bad spectra
	#	print(paste("/", param$catalogID[-F]," /"))
	#	print(paste('number of deleted spectra:', sum(!F)))
	#	print(paste('loop counter:', loopCounter))
	#	print(paste('distribution radius', distribution.radius))
	#	F <- nmr.noise.sd.distance < distribution.radius
	#	nmr <- nmr[,F]
	#	param <- param[F,]
	#	print("loop fertig")

	#	if(distribution.radius < 150){break}
	
#	} ## end of repeat
#	dev.off()
#	res <- list("nmr"=nmr,"param"=param)
#	return(res)

# }






















