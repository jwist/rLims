lims.opls.model <- function(predictors, responses, nComp, training, xcenter = mean, xscale = sd, ycenter = mean,
                            yscale = sd, discriminant = TRUE, accuracy = 1e-5, maxit = 100){
  
  ###Center and Scale####
  X <- lims.scaling(as.matrix(predictors[training,]), center = xcenter, scale = xscale);
  Y <- lims.scaling(as.matrix(responses[training,]), center = ycenter, scale = yscale);
  
  ###Create output arrays
  VIPx<-matrix(nrow=ncol(X), ncol=nComp)
  VIPy<-matrix(nrow=ncol(X), ncol=nComp)
  VIPxortho<-matrix(nrow=ncol(X), ncol=nComp)
  VIPyortho<-matrix(nrow=ncol(X), ncol=nComp)
  xscores <- c();
  yscores <- c();
  xweights <- c();
  xloadings <- c();
  yweights <- c();
  orthoscores<-c();
  orthoweights<-c();
  ortholoadings<-c();
  orthoyweights<-c();
  F<-c();
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
      t = (X %*% w) / c(crossprod(w));  # X score is X projected on X weight #denominator is most likely unnecessary as we already normalized w, double check
      t = t / sqrt(c(crossprod(t)))  ##affects in iteration
      q = crossprod(Y, t)# / c(crossprod(t)); #Y weight is Y transposed projected on X score
      q = q / sqrt(c(crossprod(q)))
      u = (Y %*% q) / c(crossprod(q)); # new Y score
      p = crossprod(X, t) / c(crossprod(t)); #compute X loading
      wortho= p - t(((t(w)%*%p)%*% solve(t(w)%*%w))%*% t(w)) ## solve
      wortho = wortho / sqrt(c(crossprod(wortho))); #normalize X weight
      #wortho= wortho%*%solve(norm(wortho))
      tortho = (X %*% wortho) %*% solve(t(wortho) %*% wortho)
      #tortho = X %*% wortho %*% solve(t(wortho) %*% wortho)
      #cortho = (Y)*(tortho)  # %*% solve(t(tortho) %*% tortho)
      #cortho = (Y* tortho)# / c(crossprod(q)); # new Y score ..
      b = (crossprod(t,u));
      if (!is.null(lastt)){
        delta = sqrt(c(crossprod(t - lastt)));
      }
      lastt = t;
    }
    
    iterations = c(iterations, iteration);
    deltas = c(deltas, delta);
    portho = t(X) %*% tortho %*% solve(t(tortho) %*% tortho) ## affects in normalization of weights fitting the model
    
    
        
    ##Deflate
    X = X - tcrossprod(tortho, portho);
   # Y = Y - tcrossprod(t, q);   #Deflation on Y is allegedly not necesary, see Westerhuis, 1998 and/or Dayal and MacGregor 1997
    
    ##Save current component results        
    xweights = cbind(xweights, w);
    yweights = cbind(yweights, q);
    xscores = cbind(xscores, t);
    yscores = cbind(yscores, u);
    orthoscores<-cbind(orthoscores,tortho)
    orthoweights<-cbind(orthoweights,wortho)
    ortholoadings<-cbind(ortholoadings,portho)
    #orthoyweights<-cbind(orthoyweights,cortho)
    xloadings = cbind(xloadings, p);
    B = c(B,b);
    coefficients[,,j] = ginv(t(xloadings)) %*% (B * t(yweights));
  #  F<-cbind(F,Y)   
        
  }
  
  #### VIP predictive #######################

  # sx= X
  # 
  # ssx<-apply(sx,2, function(x)sum( (x - mean(x) )^2))
  # #print((ssx))
  # ssy<-(apply(F,2, function(x)sum( (x - mean(x) )^2)))
  # 
  # for (i in c(1:nComp)){
  #   VIPx[,i]<-  ((xloadings[,i])^2)*ssx[i]
  #   VIPy[,i]<-  ((xloadings[,i])^2)*ssy[i]
  #   VIPxortho[,i]<-  ((ortholoadings[,i])^2)*ssx[i]
  #   VIPyortho[,i]<-  ((ortholoadings[,i])^2)*ssy[i]
  #   
  # }
  # 
  # vip<- (length(xweights[,1])*((apply(VIPx,1,sum)/sum(ssx)) + (apply(VIPy,1,sum)/sum(ssy))))^0.5 
  # 
  # vipT<- (length(xweights[,1])/2*((apply(VIPx,1,sum)/sum(ssx)) + (apply(VIPxortho,1,sum)/sum(ssx)) 
  #         + (apply(VIPyortho,1,sum)/sum(ssy)) +(apply(VIPy,1,sum)/sum(ssy))))^0.5 
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
  
  ## max1=max(Y)
  ##min1=min(Y)
  ##max= max(fitted.values)
  ##min= min(fitted.values)
  ##scal=function(x)(max1-min1)/(max-min)*(x-min)+ min1
  ##fitted.values<-(scal(fitted.values))
  # this is from world ####
  # Porto<-ortholoadings[,1]
  # Weights<-orthoweights[,1]
  # Xnew<-X.validation
  #  tnew<-as.matrix(Xnew)%*%Weights%*%solve((t(Weights)%*%Weights))
  #  Xnew<- Xnew-(tnew%*%t(Porto))
  #  input.opls<-prcomp(Xnew,scale=TRUE,center=TRUE)
  #  pred<-input.opls$x[,1]
  #pred<-(input.opls$x[,1]*-1< mean(input.opls$x[,1]))
  #input.opls<-prcomp(t,scale=TRUE,center=TRUE)
  #  
  #  test<-as.matrix(testData)
  # Y.test1<-(testy)
  #  levels(Y.test1)<-c("1","-1")
  #  t<-as.matrix(Xnew)%*%(Porto)
  #
  
  
  
  ###return results
  result <- list("xweights" = xweights,"xscores" = xscores, "orthoscores"=orthoscores,
                 "xloadings" = xloadings, "yweights" = yweights,"ortholoadings"=ortholoadings, 
                 "yscores" = yscores, "iterations" = iterations, "deltas" = deltas, "coefficients" = coefficients,
                 "fitted.values" = fitted.values, "predicted.values" = predicted.values, "Q2" = Q2, "R2" = R2,
                 "RESS"=RESS, "PRESS"=PRESS);
  
  ### results is of a class for multivariate regression models (testing)
  class(result) = "model";
  result;
}
