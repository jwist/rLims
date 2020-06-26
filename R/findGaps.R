

findGaps <- function(ppm){
  dw <- c()
  for (i in 1:(length(ppm) - 1)){
    dw[i] <- abs(ppm[i + 1] - ppm[i])
  }
  nIntervals <- sum(dw > (2 * min(dw)))
  indOfJumps <- which(dw > (2 * min(dw)))

  if (nIntervals > 0){
    endInt <- c()
    begInt <- c()
    for (j in 1:length(indOfJumps)){
      begInt[j] <- ppm[indOfJumps[j]]
      endInt[j] <- ppm[indOfJumps[j] + 1]
    }
    res <- data.frame(nIntervals, indOfJumps, begInt, endInt)
  } else {
    res <- "no gaps"
  }
  return(res)
}
