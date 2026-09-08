####Function used to generate recurrent events M_Q1 and M_Q3
datagen_M_CD4<-function(X,T0,T1,W0_tmp, W1_tmp,t.jump.w, ui,vi,tau,t.jump.m,lambda.m,betax,betaz,
                             betaw,deltam){
  
  ###Generate M
  M00=NULL
  M11=NULL

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
    
    M00_tmp=mysimRec(t.jump.risk0,lambda.risk0,tau,risk0,recurrent=TRUE)
    M11_tmp=mysimRec(t.jump.risk1,lambda.risk1,tau,risk1,recurrent=TRUE)
    
    M00_tmp=M00_tmp[which(M00_tmp$ID==i),]
    M11_tmp=M11_tmp[which(M11_tmp$ID==i),]
    
    M00=rbind(M00,M00_tmp)
    M11=rbind(M11,M11_tmp)
  }
  
  
  M00=M00[order(M00$ID,M00$stoptime),]
  M11=M11[order(M11$ID,M11$stoptime),]

  M00$M=0
  M11$M=1
  
  return(list(M00,M11))
}

