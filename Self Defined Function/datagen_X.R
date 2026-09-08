
####data generating with simulated X,Z,W


datagen_X<-function(n,ui,vi,sigmae,tau,t.jump.w,alpha0,alphaz,
                    alphax,alphat,lambda.m,t.jump.m,betax,betaz,                   betaw,deltam,lambda.y,t.jump.y,etax,etaz,etaw,etam,delta1,delta2,cen){
  #Generate X and trt 
  # Generate the corrlation matrix
  k=NROW(betax)
  mu<- rep(0,k)
  corMat <- rcorrmatrix(k)
  
  # Generate data from normal distribution
  X <- mvrnorm(n, mu=mu, Sigma=corMat, empirical = TRUE)
  # Generate binary trt
  Z <- rbinom(n,1,0.5)

  # Generate W0 and W1
  SS=length(wT)
  W0_tmp<-matrix(0,n,SS)
  W1_tmp<-matrix(0,n,SS)
  
  alpha.all=c(alpha0,alphaz,alphax,alphat,1,1)
  for (s in 1:SS){
   ei=rnorm(n)*sigmae
   W0_tmp[,s]=c(cbind(1,0,X,wT[s],ui,ei)%*%alpha.all)
   W1_tmp[,s]=c(cbind(1,1,X,wT[s],ui,ei)%*%alpha.all)
   
  }
  
  ###Generate M
  M0=NULL
  M1=NULL
  p=ncol(X)
  ####Prepare the risk0/1 for M
  ID=1:n
  beta.all=c(betaz,betax,1,deltam)
  risk0=c(exp(cbind(0,X,vi,ui)%*%beta.all))
  risk1=c(exp(cbind(1,X,vi,ui)%*%beta.all))
  for (i in 1:n){  
  riskw0=exp(W0_tmp[i,]*betaw)
  riskw1=exp(W1_tmp[i,]*betaw)
  lambda.risk0<-myprod(t.jump.m,lambda.m,t.jump.w,riskw0)[[2]]
  t.jump.risk0<-myprod(t.jump.m,lambda.m,t.jump.w,riskw0)[[1]]
  lambda.risk1<-myprod(t.jump.m,lambda.m,t.jump.w,riskw1)[[2]]
  t.jump.risk1<-myprod(t.jump.m,lambda.m,t.jump.w,riskw1)[[1]]

  M0_tmp=mysimRec(t.jump.risk0,lambda.risk0,tau,risk0,recurrent=TRUE)
  M1_tmp=mysimRec(t.jump.risk1,lambda.risk1,tau,risk1,recurrent=TRUE)
  
  M0_tmp=M0_tmp[which(M0_tmp$ID==i),]
  M1_tmp=M1_tmp[which(M1_tmp$ID==i),]
  
  M0=rbind(M0,M0_tmp)
  M1=rbind(M1,M1_tmp)
  }
 
  M0=M0[order(M0$ID,M0$stoptime),]
  M1=M1[order(M1$ID,M1$stoptime),]

  ####Simulate T^{zM^z'W^z''} one at a time####

  T000=T010=T001=T011=T100=T110=T101=T111=rep(0,n)
  risk0.y=c(exp(cbind(0,X,vi,ui)%*%c(etaz,etax,delta1,delta2)))
  risk1.y=c(exp(cbind(1,X,vi,ui)%*%c(etaz,etax,delta1,delta2)))
  for (i in 1:n){
    tmpM0=M0[which(M0$ID==i),]
    tmpM1=M1[which(M1$ID==i),]
    tmpW0=W0_tmp[i,]
    tmpW1=W1_tmp[i,]
    lambda.risk0=exp((0:length(tmpM0$stoptime))*etam)
    t.jump.risk0=tmpM0$stoptime
    lambda.risk1=exp((0:length(tmpM1$stoptime))*etam)
    t.jump.risk1=tmpM1$stoptime
    
    lambda.riskw0=exp(tmpW0*etaw)
    lambda.riskw1=exp(tmpW1*etaw)
    lambda.risk00w=myprod(t.jump.w,lambda.riskw0,t.jump.risk0,lambda.risk0)[[2]]
    t.jump.risk00w=myprod(t.jump.w,lambda.riskw0,t.jump.risk0,lambda.risk0)[[1]]
    lambda.risk10w=myprod(t.jump.w,lambda.riskw1,t.jump.risk0,lambda.risk0)[[2]]
    t.jump.risk10w=myprod(t.jump.w,lambda.riskw1,t.jump.risk0,lambda.risk0)[[1]]
    lambda.risk01w=myprod(t.jump.w,lambda.riskw0,t.jump.risk1,lambda.risk1)[[2]]
    t.jump.risk01w=myprod(t.jump.w,lambda.riskw0,t.jump.risk1,lambda.risk1)[[1]]
    lambda.risk11w=myprod(t.jump.w,lambda.riskw1,t.jump.risk1,lambda.risk1)[[2]]
    t.jump.risk11w=myprod(t.jump.w,lambda.riskw1,t.jump.risk1,lambda.risk1)[[1]]
    
    my00=myprod(t.jump.y,lambda.y,t.jump.risk00w,lambda.risk00w)
    my10=myprod(t.jump.y,lambda.y,t.jump.risk10w,lambda.risk10w)	
    my01=myprod(t.jump.y,lambda.y,t.jump.risk01w,lambda.risk01w)
    my11=myprod(t.jump.y,lambda.y,t.jump.risk11w,lambda.risk11w)	
    
    T000[i]=mysimRec(my00[[1]],my00[[2]],tau,risk0.y[i],recurrent=FALSE)$stoptime
    T100[i]=mysimRec(my00[[1]],my00[[2]],tau,risk1.y[i],recurrent=FALSE)$stoptime
    T010[i]=mysimRec(my10[[1]],my10[[2]],tau,risk0.y[i],recurrent=FALSE)$stoptime
    T110[i]=mysimRec(my10[[1]],my10[[2]],tau,risk1.y[i],recurrent=FALSE)$stoptime
    
    T001[i]=mysimRec(my01[[1]],my01[[2]],tau,risk0.y[i],recurrent=FALSE)$stoptime
    T101[i]=mysimRec(my01[[1]],my01[[2]],tau,risk1.y[i],recurrent=FALSE)$stoptime
    T011[i]=mysimRec(my11[[1]],my11[[2]],tau,risk0.y[i],recurrent=FALSE)$stoptime
    T111[i]=mysimRec(my11[[1]],my11[[2]],tau,risk1.y[i],recurrent=FALSE)$stoptime
  }

  ####Simulate Z####
  ID1=ID[which(Z==1)]
  ID0=ID[which(Z==0)]
  Mdat=rbind(M0[which(M0$ID%in%ID0),],M1[which(M1$ID%in%ID1),]) 
  Mdat=Mdat[which(Mdat$stoptime<=tau),]
  Mdat=Mdat[order(Mdat$ID,Mdat$stoptime),]
  
  W0=data.frame(cbind(ID,W0_tmp))
  colnames(W0) <- c('ID','W1','W2','W3','W4','W5')
  W0_long <- gather(W0, stoptime, W,W1:W5, factor_key=TRUE)
 
  levels(W0_long$stoptime)[levels(W0_long$stoptime)=="W1"] <- 0
  levels(W0_long$stoptime)[levels(W0_long$stoptime)=="W2"] <- 2
  levels(W0_long$stoptime)[levels(W0_long$stoptime)=="W3"] <- 4
  levels(W0_long$stoptime)[levels(W0_long$stoptime)=="W4"] <- 6
  levels(W0_long$stoptime)[levels(W0_long$stoptime)=="W5"] <- 8
  
  W0_long$event=3
  W1=data.frame(cbind(ID,W1_tmp))
  colnames(W1) <- c('ID','W1','W2','W3','W4','W5')
  W1_long <- gather(W1, stoptime, W, W1:W5, factor_key=TRUE)
  levels(W1_long$stoptime)[levels(W1_long$stoptime)=="W1"] <- 0
  levels(W1_long$stoptime)[levels(W1_long$stoptime)=="W2"] <- 2
  levels(W1_long$stoptime)[levels(W1_long$stoptime)=="W3"] <- 4
  levels(W1_long$stoptime)[levels(W1_long$stoptime)=="W4"] <- 6
  levels(W1_long$stoptime)[levels(W1_long$stoptime)=="W5"] <- 8
  W1_long$event=3
  
  Wdat=rbind(W0_long[which(W0_long$ID%in%ID0),],W1_long[which(W1_long$ID%in%ID1),]) 
  
  T=T000*(1-Z)+T111*Z

  if (cen){
    C=runif(n,tau/2,tau)
  }
  Y=pmin(C,T)
  D=as.numeric(T<=C)
  
  ydata=data.frame(ID,Y,D,Z,X)
  names(ydata)=c("ID","Y","D",paste("X",as.character(1:(p+1))))
  Mdat$W=''
  mdata=rbind(Mdat,Wdat)
  #remove duplicated tau
  mdata<-mdata[-which(mdata$stoptime==10),]
  mdata=mdata[order(mdata$ID,mdata$stoptime),]
  
  adata=merge(ydata,mdata,by="ID")
  adata=adata[which(adata$stoptime<adata$Y),]
  adata$D=adata$event
  adata$Y=adata$stoptime
  adata=adata[,c("ID","Y","D",paste("X",as.character(1:(p+1))),"W")]
  ydata$D=ydata$D*2
  ydata$W=''
  alldata=rbind(adata,ydata)
  alldata=alldata[order(alldata$ID,alldata$Y),]
  ###first is Z, all others are X
  names(alldata)=c("ID","stoptime","event",paste("X",as.character(1:(p+1)),sep=""),"W")
  ydata$vi=vi
  ydata$ui=ui
  data=list(alldata=alldata,edat=data.frame(ID,T000,T010,T001,T011,T100,T110,T101,T111),mdat=merge(data.frame(Mdat[,-4]),ydata[,-c(2,3)],by="ID"), wdat=Wdat)
}
