

nmr.bin<-function(ppm,nmrData,Nbin){
  gaps <- findGaps(ppm)
  print(gaps)
  if (gaps == "no gaps"){
    print('gaps = 0')
    for (K in 1:dim(nmrData)[2]){
      print(paste("K", K))
      z <- binData(seq(1:dim(nmrData)[1]), nmrData[, K], Nbin)
      ifelse(K == 1, dat <- z$sum.y, dat <- rbind(dat, z$sum.y))
      if (K == 1) {shift <- seq(to = ppm[length(ppm)], from = ppm[1], length.out = dim(z)[1])}
    }
    res <- list("binned.data" = dat, "binned.ppm" = shift)
  } else {
    print("gaps > 0")
    for (K in 1:dim(nmrData)[2]){
      print(paste("K", K))
      for (L in 1:(dim(gaps)[1] + 1)){
        ifelse(L == 1, begInt <- 1, begInt <- gaps$indOfJumps[L - 1] + 1)
        ifelse(L == dim(gaps)[1] + 1, endInt <- length(ppm), endInt <- gaps$indOfJumps[L])
        z <- binData(seq(from = begInt, to = endInt), nmrData[begInt:endInt, K], Nbin)
        ifelse(L == 1, dataChk <- z$sum.y, dataChk <- c(dataChk, z$sum.y))
        if (K == 1){
          zshift <- seq(to = ppm[endInt], from = ppm[begInt], length.out=dim(z)[1])
          ifelse(L == 1, ppmChk <- zshift, ppmChk <- c(ppmChk, zshift))
        }
      }
      if (K == 1) {shift <- ppmChk}
      ifelse(K == 1, dat <- dataChk, dat <- rbind(dat, dataChk))
    }
    res <- list("binned.data" = I(dat), "binned.ppm" = shift)
  }
  return(res)
}
