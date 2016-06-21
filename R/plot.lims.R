



#opls <- function(data,...) UseMethod("opls")

#nComp<-5  ##R
opls.default <- function(data,class,nComp){   ## quitar trainID          
  class<-class
  y<-class
  X<-as.matrix(data)  ## as.matrix(data)[trainID,]
  nClass <- length(levels(y))
  levels(y)<-c(1,-1,0)
  Y<-as.numeric(as.character(y))
  nClass <- length(levels(y))
  T = c()
  P = c()
  C = c()
  W = c()
  U = c()
  Tortho = c()
  Portho = c()
  Wortho = c()
  Cortho = c()
  
  for (j in 1:nComp) {
    w = (t(X) %*% Y) %*% solve(t(Y) %*% Y)
    w=w%*%solve(norm(w))
    t = (X %*% w) %*% solve(t(w) %*% w)
    c = t(Y) %*% t %*% solve(t(t) %*% t)
    u = Y %*% c %*% solve(t(c) %*% c)
    p = (t(X) %*% t) %*% solve(t(t) %*% t)
    wortho = p - w
    wortho= p - t(((t(w)%*%p)%*% solve(t(w)%*%w))%*% t(w))
    wortho= wortho%*%solve(norm(wortho))
    tortho = X %*% wortho %*% solve(t(wortho) %*% wortho)
    portho = t(X) %*% tortho %*% solve(t(tortho) %*% tortho)
    cortho = t(Y) %*% tortho %*% solve(t(tortho) %*% tortho)
    X = X - tortho %*% t(portho)
    T = matrix(c(T, t))
    P = matrix(c(P, p))
    C = matrix(c(C, c))
    W = matrix(c(W, w))
    U = matrix(c(U,u))
    Tortho = matrix(c(Tortho, tortho))
    Portho = matrix(c(Portho, portho))
    Wortho = matrix(c(Wortho, wortho))
    Cortho = matrix(c(Cortho, cortho))
  }
  T = matrix(T, ncol = nComp)
  T = scale(T, scale = FALSE, center = TRUE)
  P = matrix(P, ncol = nComp)
  C = matrix(C, ncol = nComp)
  W = matrix(W, ncol = nComp)
  U = matrix(U, ncol = nComp)
  Tortho = matrix(Tortho, ncol = nComp)
  Portho = matrix(Portho, ncol = nComp)
  Wortho = matrix(Wortho, ncol = nComp)
  Xortho = Tortho %*% t(Portho)
  results.opls <- list("T"=T[, nComp], "To"= Tortho[, 1], "Xvar"= paste(X),"C"=C, "loadings"= Portho,"wortho"= Wortho, "P"= P,"weights"=W )
  return(results.opls)
    }

#model<-opls.default((data),species,5)

#A<-7.6
#B<-7.8

selectPeak<-function(ppm, data,A,B,zoom){

F<-A< ppm & ppm< B
r<-as.matrix(data[,F])
g<-ppm[F]
p<-g[which(r==max(r), TRUE)]
areaPlot<-matplot(ppm, t(data),xlim=c(B,A),ylim=c(0,mean(r[1,]*zoom)), type="l")
peak<-match(p,ppm)
return(list("VariPeak"=peak,"peaks"=p,"area"=areaPlot))
}

#f<-selectPeak(ppm.binned,binned.data,3,3.5,100)

plotMetabo<-function(data,class,nComp,ppm){
  model1<-  opls.default((data),class,nComp) 
  data.Scale<-scale(data, scale=TRUE, center=TRUE)
  model<-  opls.default((data.Scale),class,nComp)
  sdData<-apply(data,2,sd)  
  load<-model$P[,1]*-1
  Porto1<-model$P[,1]
  Wortho1<-model$weights[,1]
  Portho<-model1$loadings
  print(dim(Portho))
#dataR<-t(scale(binned.data, center=TRUE, scale=TRUE))
#corData<-as.matrix(dataR)%*%(as.matrix(t(dataR)))
#Wortho<-model.Scale$weights
colorSd<-sdData*model$P[,1]
##### este si va 
##cor(model$C[1,])
###sd*weights vs color weights ##Loadings(sd*weight vs weights) way1 paper
###loading<-data.frame("x"=(ppm.binned),"y"=(colorSd), "Coefficients"=abs(Wortho1))  ##abs wortho paper
###d <- ggplot(x=x,y=y, data=loading)
#d + geom_line(data=loading, aes(x=x,y=y ,colour=r), size=0.8) + scale_colour_gradient( low="blue",high="red")
##rgb.palette <- colorRampPalette(c("blue", "cyan","green","orange", "red"), space = "rgb")

##plot.loading<-d + ggtitle("")+ geom_line(data=loading, aes(x=x,y=y ,colour=Coefficients), size=0.8)+ scale_colour_gradientn(colours=rgb.palette(dim(Portho)[1]), breaks=dim(Portho)[2]) + coord_cartesian(xlim = c(max(ppm.binned),min(ppm.binned))) + theme_bw() + ylab("Coeficientes PLS (au)") + xlab("ppm") + scale_x_reverse()


#+ theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))



###matrix correlation

#corData = cor((data), method = c("pearson"))
#x <- 1:nrow(CorCol)
#y <- 1:ncol(CorCol)
#filled.contour(x, y, CorCol, color = (rainbow))  
#cor1<-data.frame("x"=ppm.binned,"y"=colorSd, "r"=as.numeric(signal))
#d <- ggplot(x=x,y=y, data=cor1)
#d + geom_line(data=loading, aes(x=x,y=y ,colour=r), size=0.8) + scale_colour_gradient( low="blue",high="red")
#rgb.palette <- colorRampPalette(c("blue", "cyan","green","orange", "red"), space = "rgb")
#stocsyM<-d + ggtitle("STOCSY para 2.1ppm way2 paper") + geom_line(data=cor1, aes(x=x,y=y ,colour=r), size=0.8)+ scale_colour_gradientn(colours=rgb.palette(dim(Portho)[1]), breaks=dim(Portho)[2])+ coord_cartesian(xlim = c(10,-1)) 



##loadings vs color weight  ## parecido simca
#loading1<-data.frame("x"=seq(1:dim(Portho)[1]),"y"=load, "r"=(Wortho1))
#d <- ggplot(x=x,y=y, data=loading)
#d + geom_line(data=loading, aes(x=x,y=y ,colour=r), size=0.8) + scale_colour_gradient( low="blue",high="red")
#rgb.palette <- colorRampPalette(c("blue", "cyan","green","orange", "red"), space = "rgb")
#plot.loading1<-d + ggtitle("Loadings(load vs weight)SIMCA")+ geom_line(data=loading1, aes(x=x,y=y ,colour=r), size=0.8)+ scale_colour_gradientn(colours=rgb.palette(dim(Portho)[1]), breaks=dim(Portho)[2])

##da lo mismo con Wortho1
p1<-Portho[,1]
p2<-Portho[,2]
plot(p1,p2,pch=16, main="Loading plot", xlab="LV1", ylab="LV2")
#text(p1[which(p1<mean(p1-0.4*sd(p1)),T)],p2[which(p1<mean(p1-0.4*sd(p1)),T)],round(ppm[which(p1<mean(p1-0.4*sd(p1)),T)],2), cex=0.6, pos=3,col=4)
#text(p1[which(p1>mean(p1+sd(p1)),T)],p2[which(p1>mean(p1+sd(p1)),T)],round(ppm[which(p1>mean(p1+sd(p1)),T)],2), cex=0.6, pos=3,col=4)


#text(Portho[,1],Portho[,2],ppm.binned, cex=0.5, pos=4)

#stocsy  ##wiklund no va
spec<-cov(as.numeric(model$T),as.matrix(data))  ##revisar
r<-(as.numeric(cor(model$T,as.matrix(data))))
#cor<-data.frame("x"=seq(1:length(spec)),"y"=(as.numeric(spec)), "r"=(r))
#d <- ggplot(x=x,y=y, data=cor)
#d + geom_line(data=loading, aes(x=x,y=y ,colour=r), size=0.8) + scale_colour_gradient( low="blue",high="red")
#rgb.palette <- colorRampPalette(c("blue", "cyan","green","orange", "red"), space = "rgb")
#stocsy<-d + ggtitle("STOCSY Wiklund") + geom_line(data=cor, aes(x=x,y=y ,colour=r), size=0.8)+ scale_colour_gradientn(colours=rgb.palette(dim(Portho)[1]), breaks=dim(Portho)[2])

#splot

plot(spec,r,pch=16, main="S plot",xlab="cov(T,X)", ylab="cor(T,X)")
text(spec[which(spec<mean(spec-sd(spec)),T)[,2]],r[which(spec<mean(spec-sd(spec)),T)[,2]],round(ppm[which(spec<mean(spec-sd(spec)),T)[,2]],2), cex=0.6, pos=3,col=4)
#text(spec,r,ppm.binned, cex=0.5, pos=4)
text(spec[which(spec>mean(spec+sd(spec)),T)[,2]],r[which(spec>mean(spec+sd(spec)),T)[,2]],round(ppm[which(spec>mean(spec+sd(spec)),T)[,2]],2), cex=0.6, pos=3,col=4)

###"loadingsSD"=plot.loading
#plot.results <- list("splot"=Splot, "Lplot"=Lplot)
return(list("load1"=p1,"load2"=p2))
}

#plotsAnalysis<-plotMetabo((binned.data),param$country,5,ppm.binned,corData[,1098])
#http://www.r-bloggers.com/fast-computation-of-cross-validation-in-linear-models/
#http://stackoverflow.com/questions/7402313/generate-sets-for-cross-validation-in-r
#B<- (solve((X)%*%t(X)))%*% ((X)%*%t(Y))

