###statistics###


##Moore-Penrose generalized inverse of a matriz
ginv<-function(X, tol = sqrt(.Machine$double.eps))
{
  ## Generalized Inverse of a Matrix
  dnx <- dimnames(X)
  if(is.null(dnx)) dnx <- vector("list", 2)
  s <- svd(X)
  nz <- s$d > tol * s$d[1]
  structure(
    if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
    dimnames = dnx[2:1])
}

##Standard scale functions

lims.range <- function(x) { max(x) - min(x)};
lims.pareto <- function(x) { sqrt(sd(x)) };
lims.vast <- function(x) { sd(x)^2 / mean(x) };

#Scale and center data
lims.scaling <-function(data, center=mean, scale=sd){
  
  if(is.logical(center)){
    if (center == FALSE){
      center <- function(x){0};
    }
    else{
      if (center == TRUE){
        center <- mean;
      }
    }
  }
  if(is.logical(scale)){
    if (scale == FALSE){
      scale <- function(x){1};
    }
    else{
      if (scale ==TRUE){
        scale <- sd;
      }
    }
  }
  
  #     data <- t(data);
  #     
  #     scale(data, center = apply(data,2,mean), scale = apply(data,2,scale)));
  
  apply(data,2,function(x) (x - center(x)) / scale(x));
}



lims.pls.model <- function(predictors, responses, nComp, training, xcenter = mean, xscale=sd, ycenter = mean, yscale=sd, discriminant = TRUE, accuracy = 1e-5, maxit = 100){
  
  
  ###Center and Scale####
  # X<-predictors
  #Y<-responses
  X <- lims.scaling(as.matrix(predictors[training,]), center = xcenter, scale = xscale);
  Y <- lims.scaling(as.matrix(responses[training,]), center = ycenter, scale = yscale);
  
  ###Create output arrays
  xscores <- c();
  yscores <- c();
  xweights <- c();
  xloadings <- c();
  yweights <- c();
  F <-c();
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
      w = crossprod(X, u)# / c(crossprod(u)); # X weight is X transposed projected on y score# denominator is unnecessary since we normalize next
      
      w = w / sqrt(c(crossprod(w))); #normalize X weight
      t = (X %*% w)# / c(crossprod(w));  # X score is X projected on X weight #denominator is most likely unnecessary as we already normalized w, double check
      t = t / sqrt(c(crossprod(t)))
      q = crossprod(Y, t)# / c(crossprod(t)); #Y weight is Y transposed projected on X score
      q = q / sqrt(c(crossprod(q)))
      u = (Y %*% q)# / c(crossprod(q)); # new Y score
      if (!is.null(lastt)){
        delta = sqrt(c(crossprod(t - lastt)));
      }
      lastt = t;
    }
    
    iterations = c(iterations, iteration);
    deltas = c(deltas, delta);
    p = crossprod(X, t)# / c(crossprod(t)); #compute X loading
    b = crossprod(t, u);
    
    
    ##Deflate
    X = X - tcrossprod(t, p);
    Y = Y - tcrossprod(t, q); #Y  #Deflation on Y is allegedly not necesary, see Westerhuis, 1998 and/or Dayal and MacGregor 1997
    f =  tcrossprod(t, q);
    
    
    ##Save current component results        
    xweights = cbind(xweights, w);
    yweights = cbind(yweights, q);
    xscores = cbind(xscores, t);
    yscores = cbind(yscores, u);
    xloadings = cbind(xloadings, p);
    F= cbind(F,f)
    B = c(B,b);
    coefficients[,,j] = ginv(t(xloadings)) %*% (B * t(yweights));
    
  }
  
  ############ VIP ####################
  #SSX
  #Sum of squares of the X block. For component number A, 
  #it is the X residual Sum of Squares after component A.
  
  #SSY
  #Sum of squares of the Y block. For component number A, 
  #it is the Y residual Sum of Squares after component A.
  
  
  VIP<-matrix(nrow=ncol(X), ncol=nComp)
  vip<-c()
  sx= X[1:nComp,] #xloadings%*%t(xloadings) 
  
  ssx<-apply(sx,2, function(x) sum((x - mean(x))^2))
  ssy<-(apply(F,2, function(x) sum((x - mean(x))^2)))
  
  for (i in c(1:nComp)){
    VIP[,i]<-  ((xweights[,i])^2)*ssy[i]
  }
  vip<- (length(xweights[,1])/sum(ssy)*rowSums(VIP))^0.5 #(apply(VIP,1,sum)) 
  
  #     #VIP <- function(object) {
  #       
  SS <- c(yweights)^2 * colSums(xscores^2)
  Wnorm2 <- colSums(xweights^2)
  SSW <- sweep(xweights^2, 2, SS / Wnorm2, "*")
  #vip<- sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
  #     #}
  #     
  #     
  #     ## VIPjh returns the VIP of variable j with h components
  #  #   VIPjh <- function(object, j, h) {
  #       for (i in c(1:nComp)){
  #       b <- c(yweights)
  #       T <- xscores
  #       SS <- b^2 * colSums(T^2)
  #       W <- xweights
  #       Wnorm2 <- colSums(W^2)
  #       sqrt(nrow(W) * sum(SS * W[j,]^2 / Wnorm2) / sum(SS))
  #     }
  #     
  #     
  #     
  #     
  
  
  
  ###quality assess
  X <- lims.scaling(as.matrix(predictors), center = xcenter, scale = xscale);
  Y <- lims.scaling(as.matrix(responses), center = ycenter, scale = yscale);
  X.validation <- as.matrix(X[!(1:dim(X)[1] %in% training),]);
  Y.validation <- as.matrix(Y[!(1:dim(Y)[1] %in% training),]);
  X <- as.matrix(X[training,]);
  Y <- as.matrix(Y[training,]);
  
  fitted.values <- array(0,c(dim(Y),nComp))
  predicted.values <- array(0,c(dim(Y.validation),nComp))
  R2 = c();
  Q2 = c();
  PRESS=c();  
  RESS=c();
  
  fitted.filter = TRUE;
  predicted.filter = TRUE;
  
  for (i in 1:nComp){
    fitted.values[,,i] = X %*% coefficients[,,i];
    predicted.values[,,i] = X.validation %*% coefficients[,,i];
    
    if (discriminant){
      fitted.filter = !(abs(fitted.values[,,i]) > 1 & fitted.values[,,i] * Y > 0);
      predicted.filter = !(abs(predicted.values[,,i]) > 1 & predicted.values[,,i] * Y.validation > 0);
    }
    
    RESS= c(RESS,sum(((Y - fitted.values[,,i]) * fitted.filter)^2))
    PRESS= c(PRESS,sum(((Y.validation - predicted.values[,,i]) * predicted.filter)^2))
    # R2 = c(R2, 1 - (sum((Y - fitted.values[,,i])^2) / TSS));
    # Q2= Q2 = c(Q2, 1 - (sum((Y.validation - predicted.values[,,i])^2) / TSS)); #according to Szymánska et al, Metabolomics, (2012). May be inadequate for classification, in that case graphic analysis is better (prediction vs. observation plot)
    R2 = c(R2, 1 - RESS[i] / sum((Y - mean(Y))^2));
    Q2 = c(Q2, 1 - PRESS[i] / sum((Y.validation - mean(Y.validation))^2));
    
    
  }
  
  ## max1=max(responses[training,])
  ## min1=min(responses[training,])
  ## max= max(fitted.values)
  ## min= min(fitted.values)
  ## scal=function(x)(max1-min1)/(max-min)*(x-min)+ min1
  ## fitted.values<-(scal(fitted.values))
  
  ##  max1=max(responses[training,])
  ##  min1=min(responses[training,])
  ##  max= max(fitted.values)
  ##  min= min(fitted.values)
  ##  scal=function(x)(max1-min1)/(max-min)*(x-min)+ min1
  ##  fitted.values<-(scal(fitted.values))
  
  
  
  ###return results
  result <- list("xweights" = xweights, "xscores" = xscores, "xloadings" = xloadings,
                 "yweights" = yweights, "yscores" = yscores, "iterations" = iterations, "deltas" = deltas,
                 "coefficients" = coefficients, "fitted.values" = fitted.values,"ssx"=ssx,
                 "predicted.values" = predicted.values, "Q2" = Q2, "R2" = R2, "PRESS"=PRESS, 
                 "RESS" = RESS,"F"=F, "vip"=vip);
  
  ### results is of a class for multivariate regression models (testing)
  class(result) = "model";
  result;
}

#prediction with multiariate regression model (testing)
prediction <- function(model,predictor,nComp,...) UseMethod("prediction")

prediction.model <- function(model,predictor,nComp,...){
  predictor %*% model$coefficients[,,nComp]
}

#Cross-validation
#If codomain (possible values of response variable) is given, it is used to homogenously sample from each value
#e.g if codomain=(-1,1) (default), training consists of 80% of positive (1) and 20% of negative (-1) responses
#WARNING I assume that codomain, if given, has at least two elements (because doing this makes no sense otherwise)
#method: randomp builds nModels with randomly selected valsize * sample points for validation
#        lpo (leave p out) builts a model with each subest of valsizez elements as validation set NOT IMPLEMENTED
#        kfold partitions the set radomly in valsize equal sized subsets and builts a model with each of these subsets as validation set
lims.validate <- function(model, predictors, responses, method="randomp", nModels=100, nComp=1, valsize=0.2, xcenter = mean, xscale=sd,
                          ycenter = mean, yscale=sd, codomain = c(-1,1), discriminant = TRUE,
                          accuracy = 1e-5, maxit = 100){
  
  results <- list();
  if (method=="randomp"){
    results$Q2 <- matrix(nrow = nModels, ncol = nComp);
    results$R2 <- matrix(nrow = nModels, ncol = nComp);
    results$training.set <- c();
    results$full.results <- c();
    if(is.vector(codomain) == FALSE){
      responses.by.value = list(seq_along(responses));
      #results$training.set <- matrix(nrow = nModels, ncol = round(nrow(predictors) * (1-valsize)));
    }
    else{
      responses.by.value = lapply(codomain, function(x) which(responses == x));
     # results$training.set <- matrix(nrow = nModels,
                             #        ncol = Reduce(function(x,y) round(length(x)*(1-valsize) + round(length(y)*(1-valsize))),
                                #                   responses.by.value));
    }
    for (i in 1:nModels){
      training = Reduce(c, lapply(responses.by.value, function(y) sample(y, round(length(y) * (1-valsize)))));
      current.model <- model(predictors = predictors, responses = responses, nComp = nComp,
                             training = training, xcenter = xcenter, ycenter = ycenter, 
                             xscale=xscale, yscale=yscale, accuracy = accuracy, maxit = maxit,
                             discriminant = discriminant);
      results$training.set[[i]]=training;
      results$Q2[i,] = current.model$Q2;
      results$R2[i,] = current.model$R2;
      results$full.results = c(results$full.results, list(current.model));
    }
  }
  else if(method=="kfold"){
    results$Q2 <- matrix(nrow = valsize, ncol = nComp);
    results$R2 <- matrix(nrow = valsize, ncol = nComp);
    results$full.results <- c();
    results$training.set <- c(); #matrix(nrow = valsize, ncol = nrow(predictors) - floor(nrow(predictors)/valsize));
    
    if(is.vector(codomain) == FALSE){
      responses.by.value = list(seq_along(responses));
      #results$training.set <- matrix(nrow = nModels, ncol = round(nrow(predictors) * (1-valsize)));
    }
    else{
      responses.by.value = lapply(codomain, function(x) which(responses == x));
      # results$training.set <- matrix(nrow = nModels,
      #        ncol = Reduce(function(x,y) round(length(x)*(1-valsize) + round(length(y)*(1-valsize))),
      #                   responses.by.value));
    }
    
    partition1 <- split(sample(responses.by.value[[1]], length(responses.by.value[[1]]), replace=FALSE),1:valsize);
    partition2 <- split(sample(responses.by.value[[2]], length(responses.by.value[[2]]), replace=FALSE),1:valsize);
    partition <- list();
    for (i in 1:valsize){
      partition[[i]] = c(partition1[[i]], partition2[[i]]);
    }
#    partition <- Reduce(rbind,lapply(responses.by.value, function(y) split(sample(y,length(y), replace=FALSE),1:valsize)))
 #     c(split(sample(nrow(predictors),nrow(predictors),replace=FALSE),1:valsize);
    
    for (i in 1:valsize){
      training <- c(1:nrow(predictors))[-partition[[i]]];

      current.model <- model(predictors = predictors, responses = responses, nComp = nComp,
                             training = training, xcenter = xcenter, ycenter = ycenter, 
                             xscale=xscale, yscale=yscale, accuracy = accuracy, maxit = maxit,
                             discriminant = discriminant);
      
      #if(length(partition[[i]]) == floor(nrow(predictors)/valsize)){
        results$training.set[[i]] = training;
     # }
     # else{
      #  results$training.set[i,] = training; c(training,0);
     # }
      
      results$Q2[i,] = current.model$Q2;
      results$R2[i,] = current.model$R2;
      results$full.results = c(results$full.results, list(current.model));
    }
  }
  results;
}
