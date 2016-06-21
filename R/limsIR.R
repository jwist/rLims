lims.getIrs <- function(entryList=entryList,...) {
  
  functionName <- 'lims.getIrs'
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
  
  if (sum(is.na(match(names(Filter),names(entryList[[1]]$irs[[1]])))) != 0) {
    warning('lims.getIrs: some names in parameter filter are not correct, please check')
  }
  
  if (length(entryList) > 0) {
    for (i in 1:length(entryList)) {
      x <- entryList[[i]]
      if (length(x$irs) > 0) {
        for (j in 1:length(x$irs)) {
          y <- x$irs[[j]]
          
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
              #  spect <- read.table(sprintf("%s&filter=JcampToXY",entry[[1]]$irs[[1]]$resourceURL),sep=',',colClasses='numeric')
              #  print(spect)
              # print(y$resourceURL)
              # print(sum(is.na(spect[1:(dim(spect)[1]/2),][[1]])))
              #   print(range(spect[1:(dim(spect)[1]/2),][[1]]))
              #   plot(lims.spectraCreator(spect[1:(dim(spect)[1]/2),] ))
              
              s <- rbind(s, list('spectra'=lims.spectraCreator(data.frame(spect[1:(dim(spect)[1]),] ), type="ir"))) # taking only real part of spectrum
              #  print(s)
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
        }} #irss
    }} #entry
  
  # returns nmrData
  if (!dry) {
    return(list('spectra'=s,'nmrInfo'=data.frame(t),'param'=data.frame(p),'info'=data.frame(info)))
  } else {
    return(list('nmrInfo'=data.frame(t),'param'=data.frame(p),'info'=data.frame(info)))
  }
}


msc.lims<-function(x){
  xc<-apply(x,2,function(x) x-mean(x))
  r<-apply(x,1,function(x) mean(x))
  rc<-r-mean(r)
  b<- solve(t(rc)%*%rc)%*%t(rc)%*%xc
  xco<-matrix(nrow=dim(xc)[1], ncol=dim(xc)[2])
  for (i in 1:dim(xc)[2]){
    xco[,i]<-xc[,i]/b[,i] +mean(r)
  }
  xco
}




osc.lims <- function(predictors, responses, nComp, training, xcenter = mean, xscale = sd, ycenter = mean,
                            yscale = sd, discriminant = TRUE, accuracy = 1e-5, maxit = 100){
  
  ###Center and Scale####
  #X <- lims.scaling(as.matrix(predictors[training,]), center = xcenter, scale = xscale);
  #Y <- lims.scaling(as.matrix(responses[training,]), center = ycenter, scale = yscale);
  
  X<-as.matrix(predictors)
  Y<-as.matrix(responses)
  ###Create output arrays
  xscores <- c();
  yscores <- c();
  xweights <- c();
  xloadings <- c();
  yweights <- c();
  B <- c();
  coefficients <- array(0, c(dim(X)[2], dim(Y)[2], nComp));
  
  ###Variables for loop control
  iterations <- c();
  deltas <- c();
  
  if (dim(Y)[2] == 1){
    maxit <- 1;
  }
  
  ###loop on components (aka latent variables)
  for (j in 1:nComp) {
    u <- matrix(Y[,1],ncol=1) #/ sqrt(sum(Y * Y)); #initial guess for Y score #/ sqrt(sum(Y * Y)) is apparently not necessary, double check
    
    delta = 1000000;
    lastt = NULL;
    iteration = 0;
    
    while (delta > accuracy && iteration < maxit){
      iteration = iteration + 1
      t=prcomp(X,scale=FALSE,center=FALSE)$x
      tnew = 1-  (Y%*%solve(t(Y)%*%(Y))%*%t(Y))%*%t  ##affects in iteration
      w = ginv(X) %*% tnew # / c(crossprod(u)); # X weight is X transposed projected on y score# denominator is unnecessary since we normalize next
      t = (X %*% w) #/ c(crossprod(w));
      
      p = t(t(t)%*%X) %*% (ginv(t%*%tnew))
      #b = (crossprod(t,u));
      if (!is.null(lastt)){
        delta = sqrt(c(crossprod(t - lastt)));
      }
      lastt = t;
    }
    
    iterations = c(iterations, iteration);
    deltas = c(deltas, delta);
    
    ##Save current component results        
    xweights = cbind(xweights, w);
    xscores = cbind(xscores, t);
    xloadings = cbind(xloadings, p);
    #B = c(B,b);
    #coefficients[,,j] = ginv(t(xloadings)) %*% (B * t(yweights));
    
    ##Deflate
    E = X- t%*%t(p)
  }
  ###quality assess
  #X <- lims.scaling(as.matrix(predictors), center = xcenter, scale = xscale);
  #Y <- lims.scaling(as.matrix(responses), center = ycenter, scale = yscale);
  #X.validation <- as.matrix(X[!(1:dim(X)[1] %in% training),]);
  #Y.validation <- as.matrix(Y[!(1:dim(Y)[1] %in% training),]);
  #X <- as.matrix(X[training,]);
  #Y <- as.matrix(Y[training,]);
  
  #fitted.values <- array(0,c(dim(Y),nComp))
  #predicted.values <- array(0,c(dim(Y.validation),nComp))
  #R2 = c();
  #Q2 = c();
  #PRESS=c();  
  #RESS=c();
  
  #fitted.filter = TRUE;
  #predicted.filter = TRUE;
  
  # for (i in 1:nComp){
  #   fitted.values[,,i] = X %*% coefficients[,,i];
  #   predicted.values[,,i] = X.validation %*% coefficients[,,i];
  #   
  #   if (discriminant){
  #     fitted.filter = !(abs(fitted.values[,,i]) > 1 & fitted.values[,,i] * Y > 0);
  #     predicted.filter = !(abs(predicted.values[,,i]) > 1 & predicted.values[,,i] * Y.validation > 0);
  #   }
  #   
  #   RESS= c(RESS,sum(((Y - fitted.values[,,i]) * fitted.filter)^2))
  #   PRESS= c(PRESS,sum(((Y.validation - predicted.values[,,i]) * predicted.filter)^2))
  #   # R2 = c(R2, 1 - (sum((Y - fitted.values[,,i])^2) / TSS));
  #   # Q2= Q2 = c(Q2, 1 - (sum((Y.validation - predicted.values[,,i])^2) / TSS)); #according to Szymánska et al, Metabolomics, (2012). May be inadequate for classification, in that case graphic analysis is better (prediction vs. observation plot)
  #   R2 = c(R2, 1 - RESS[i] / sum((Y - mean(Y))^2));
  #   Q2 = c(Q2, 1 - PRESS[i] / sum((Y.validation - mean(Y.validation))^2));
  #   

  ###return results
  result <- list("xweights" = xweights,"xscores" = xscores,
                 "xloadings" = xloadings, "X"=E,"iterations" = iterations, "deltas" = deltas);
  
  ### results is of a class for multivariate regression models (testing)
  class(result) = "model";
  result;
}
