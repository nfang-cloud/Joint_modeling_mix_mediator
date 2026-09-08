

library(MASS)
library(survival)
library(dplyr)
library(reda)
library(clusterGeneration)
###load in self-defined functions#####
source("Lambdainv.R")
source("mysimRec.R")
source("myprod.R")
source("datagenT01_M_multi.R")
source("S.R")
source("datagen_X.R")

#######################compute true effects using large data#####################################
#################################################################################################
##### TRUE parameter ######
n=10000
k=1
tau=10
t.jump.m1=c(2,4,6)
t.jump.m2=c(2,4,6)
t.jump.y=c(2,4,6)
lambda.m1=c(0.65,0.5,0.45,0.35)
lambda.m2=c(0.5,0.35,0.3,0.25)
lambda.y=c(0.033,0.084,0.1,0.6)

betax1=rep(0.2,k)
betax2=rep(0.15,k)
betaz1=0.35
betaz2=0.4
etaz=0.35
etax=rep(0.15,k)
etam1=0.18
etam2=0.15
delta1=0.7
deltax1=0.7
deltax2=0.84
delta2_1=0
delta2_2=0

shape=1.5
scale=1

est=cbind(log1_r1=log(lambda.m1[1]),log1_r2=log(lambda.m1[2]),log1_r3=log(lambda.m1[3]),log1_r4=log(lambda.m1[4]),
          log2_r1=log(lambda.m2[1]),log2_r2=log(lambda.m2[2]),log2_r3=log(lambda.m2[3]),log2_r4=log(lambda.m2[4]),
          log_h1=log(lambda.y[1]),log_h2=log(lambda.y[2]),log_h3=log(lambda.y[3]),log_h4=log(lambda.y[4]),
          betax1=betax1,betax2=betax2,
          betaz1=betaz1,betaz2=betaz2,
          deltax1=deltax1,deltax2=deltax2,
          etax=etax,etaz=etaz,delta1=delta1,etam1=etam1,etam2=etam2)
est=as.data.frame(est)

#####Generate a large dataset#######


set.seed(8888)
v_gamma=rgamma(n,shape = 1.5, scale = 1)
vi=log(v_gamma)
data=datagen_X(n,vi,tau,betaz1,betaz2,betax1,betax2,lambda.m1,lambda.m2,t.jump.m1,t.jump.m2,
               etaz,etax,etam1,etam2,delta1,deltax1,deltax2,delta2_1,delta2_2,
               lambda.y,t.jump.y,cen=TRUE)
X2 <-data %>% group_by(ID) %>% filter(row_number(X2)==1) %>% select(ID,X2)
X=X2$X2
tseq=seq(1,8,by=1)

#####Get the true value####
mymed<-function(X,n,est,tseq,B){
  t.jump.m1=c(2,4,6)
  t.jump.m2=c(2,4,6)
  t.jump.y=c(2,4,6)
  tau=10
  T0=0
  T1=1
  betaz1=est$betaz1
  betaz2=est$betaz2
  betax1=est$betax1
  betax2=est$betax2
  lambda.m1=exp(c(est$log1_r1,est$log1_r2,est$log1_r3,est$log1_r4))
  lambda.m2=exp(c(est$log2_r1,est$log2_r2,est$log2_r3,est$log2_r4))
  lambda.y=exp(c(est$log_h1,est$log_h2,est$log_h3,est$log_h4))
  delta_X1=est$deltax1
  delta_X2=est$deltax2
  etaz=est$etaz
  etax=est$etax
  etam1=est$etam1
  etam2=est$etam2
  delta1=est$delta1
  delta2_1=0
  delta2_2=0
  
  
  #Generate the M_00/M_11
  avgres=matrix(data=0,nrow=n,ncol=8*length(tseq))  
  avgressq= matrix(data=0,nrow=n,ncol=8*length(tseq))
  
  for (b in 1:B){
    
    set.seed(1111+b)
    v_gamma=rgamma(n,shape = 1.5, scale = 1)
    vi=log(v_gamma)
    M_T=datagenT01_M_multi(X,T0,T1,vi,
                           tau,betaz1,betaz2,betax1,betax2,lambda.m1,lambda.m2,
                           t.jump.m1,t.jump.m2,delta_X1,delta_X2)
    M_00=M_T[[1]]
    M_01=M_T[[2]]
    M_10=M_T[[3]]
    M_11=M_T[[4]]
    R_00=M_00[which(M_00$event==1),][,-c(3,5)]
    R_11=M_11[which(M_11$event==1),][,-c(3,5)]
    R_01=M_01[which(M_01$event==1),][,-c(3,5)]
    R_10=M_10[which(M_10$event==1),][,-c(3,5)]
    
    
    tmpres=avgres
    etam1=c(etam1,delta2_1)
    etam2=c(etam2,delta2_2)
    S000b=S100b=S011b=S111b=S001b=S101b=S010b=S110b=matrix(ncol=length(tseq),nrow=n)
    
    for (i in 1:n){
      Xi=X[i]
      R_00i=R_00[which(R_00$ID==i),]
      R_11i=R_11[which(R_11$ID==i),]
      R_01i=R_01[which(R_01$ID==i),]
      R_10i=R_10[which(R_10$ID==i),]
      
      R_00_1i=R_00i[which(R_00i$event_Tp==1),]
      R_00_2i=R_00i[which(R_00i$event_Tp==2),]
      R_11_1i=R_11i[which(R_11i$event_Tp==1),]
      R_11_2i=R_11i[which(R_11i$event_Tp==2),]
      
      R_01_1i=R_01i[which(R_01i$event_Tp==1),]
      R_01_2i=R_01i[which(R_01i$event_Tp==2),]
      R_10_1i=R_10i[which(R_10i$event_Tp==1),]
      R_10_2i=R_10i[which(R_10i$event_Tp==2),]
      
      
      lambda.risk00_1=exp((0:length(R_00_1i$stoptime))*(etam1[1]+etam1[2]*vi[i]))
      t.jump.risk00_1=R_00_1i$stoptime
      
      lambda.risk00_2=exp((0:length(R_00_2i$stoptime))*(etam2[1]+etam2[2]*vi[i]))
      t.jump.risk00_2=R_00_2i$stoptime
      
      lambda.risk11_1=exp((0:length(R_11_1i$stoptime))*(etam1[1]+etam1[2]*vi[i]))
      t.jump.risk11_1=R_11_1i$stoptime
      
      lambda.risk11_2=exp((0:length(R_11_2i$stoptime))*(etam2[1]+etam2[2]*vi[i]))
      t.jump.risk11_2=R_11_2i$stoptime
      
      lambda.risk01_1=exp((0:length(R_01_1i$stoptime))*(etam1[1]+etam1[2]*vi[i]))
      t.jump.risk01_1=R_01_1i$stoptime
      
      lambda.risk01_2=exp((0:length(R_01_2i$stoptime))*(etam2[1]+etam2[2]*vi[i]))
      t.jump.risk01_2=R_01_2i$stoptime
      
      lambda.risk10_1=exp((0:length(R_10_1i$stoptime))*(etam1[1]+etam1[2]*vi[i]))
      t.jump.risk10_1=R_10_1i$stoptime
      
      lambda.risk10_2=exp((0:length(R_10_2i$stoptime))*(etam2[1]+etam2[2]*vi[i]))
      t.jump.risk10_2=R_10_2i$stoptime
      
      
      
      lambda.risk00<-myprod(t.jump.risk00_1,lambda.risk00_1,t.jump.risk00_2,lambda.risk00_2)[[2]]
      t.jump.risk00<-myprod(t.jump.risk00_1,lambda.risk00_1,t.jump.risk00_2,lambda.risk00_2)[[1]]
      
      lambda.risk11<-myprod(t.jump.risk11_1,lambda.risk11_1,t.jump.risk11_2,lambda.risk11_2)[[2]]
      t.jump.risk11<-myprod(t.jump.risk11_1,lambda.risk11_1,t.jump.risk11_2,lambda.risk11_2)[[1]]
      
      
      lambda.risk01<-myprod(t.jump.risk01_1,lambda.risk01_1,t.jump.risk01_2,lambda.risk01_2)[[2]]
      t.jump.risk01<-myprod(t.jump.risk01_1,lambda.risk01_1,t.jump.risk01_2,lambda.risk01_2)[[1]]
      
      lambda.risk10<-myprod(t.jump.risk10_1,lambda.risk10_1,t.jump.risk10_2,lambda.risk10_2)[[2]]
      t.jump.risk10<-myprod(t.jump.risk10_1,lambda.risk10_1,t.jump.risk10_2,lambda.risk10_2)[[1]]
      
      my00=myprod(t.jump.y,lambda.y,t.jump.risk00,lambda.risk00)
      my11=myprod(t.jump.y,lambda.y,t.jump.risk11,lambda.risk11)
      my01=myprod(t.jump.y,lambda.y,t.jump.risk01,lambda.risk01)
      my10=myprod(t.jump.y,lambda.y,t.jump.risk10,lambda.risk10)
      
      
      S000b=S(Xi,T0,vi[i],my00, etaz,etax,delta1,tseq)#modified
      S100b=S(Xi,T1,vi[i],my00, etaz,etax,delta1,tseq)#modified
      S011b=S(Xi,T0,vi[i],my11, etaz,etax,delta1,tseq)#modified
      S111b=S(Xi,T1,vi[i],my11, etaz,etax,delta1,tseq)#modified
      
      S001b=S(Xi,T0,vi[i],my01, etaz,etax,delta1,tseq)#modified
      S101b=S(Xi,T1,vi[i],my01, etaz,etax,delta1,tseq)#modified
      S010b=S(Xi,T0,vi[i],my10, etaz,etax,delta1,tseq)#modified
      S110b=S(Xi,T1,vi[i],my10, etaz,etax,delta1,tseq)#modified
      
      tmpres[i,]= c(S000b,S011b,S100b,S111b,S001b,S101b,S010b,S110b)
      
    }
    
    tmpressq=tmpres^2
    
    avgres=avgres+(tmpres-avgres)/b
    avgressq=avgressq+(tmpressq-avgressq)/b
    
  }
  avgsd=sqrt(avgressq-avgres^2)
  
  avg_S=apply(avgres,2,mean,na.rm=TRUE)
  sd_S=apply(avgres,2,sd,na.rm=TRUE)
  
  S000_n=avg_S[1:length(tseq)]
  S011_n=avg_S[(length(tseq)+1):(length(tseq)*2)]
  S100_n=avg_S[(length(tseq)*2+1):(length(tseq)*3)]
  S111_n=avg_S[(length(tseq)*3+1):(length(tseq)*4)]
  S001_n=avg_S[(length(tseq)*4+1):(length(tseq)*5)]
  S101_n=avg_S[(length(tseq)*5+1):(length(tseq)*6)]
  S010_n=avg_S[(length(tseq)*6+1):(length(tseq)*7)]
  S110_n=avg_S[(length(tseq)*7+1):(length(tseq)*8)]
  
  
  S000_sd=sd_S[1:length(tseq)]
  S011_sd=sd_S[(length(tseq)+1):(length(tseq)*2)]
  S100_sd=sd_S[(length(tseq)*2+1):(length(tseq)*3)]
  S111_sd=sd_S[(length(tseq)*3+1):(length(tseq)*4)]
  S001_sd=sd_S[(length(tseq)*4+1):(length(tseq)*5)]
  S101_sd=sd_S[(length(tseq)*5+1):(length(tseq)*6)]
  S010_sd=sd_S[(length(tseq)*6+1):(length(tseq)*7)]
  S110_sd=sd_S[(length(tseq)*7+1):(length(tseq)*8)]
  
  NDE01=S100_n-S000_n
  NDE10=S011_n-S111_n
  NIE01=S011_n-S000_n
  NIE10=S100_n-S111_n
  
  NDE=(NDE01-NDE10)/2
  NIE=(NIE01-NIE10)/2
  TE=NDE+NIE
  
  NIE01_1_1=S011_n-S001_n
  NIE01_1_0=S010_n-S000_n
  NIE10_1_1=S101_n-S111_n
  NIE10_1_0=S100_n-S110_n
  
  NIE01_2_1=S011_n-S010_n
  NIE01_2_0=S001_n-S000_n
  NIE10_2_1=S110_n-S111_n
  NIE10_2_0=S100_n-S101_n
  
  NIE01_1=(NIE01_1_1+NIE01_1_0)/2
  NIE01_2=(NIE01_2_1+NIE01_2_0)/2
  NIE10_1=(NIE10_1_1+NIE10_1_0)/2
  NIE10_2=(NIE10_2_1+NIE10_2_0)/2
  
  NIE1=(NIE01_1-NIE10_1)/2
  NIE2=(NIE01_2-NIE10_2)/2
  
  
  return(as.data.frame(list(NDE=NDE,NIE=NIE,TE=TE,NIE1=NIE1,NIE2=NIE2,
                            S000=S000_n,S011=S011_n,S100=S100_n,S111=S111_n,
                            S001=S001_n,S010=S010_n,S101=S101_n,S110=S110_n,
                            S000_sd=S000_sd,S011_sd=S011_sd,S100_sd=S100_sd,S111_sd=S111_sd,
                            S001_sd=S001_sd,S010_sd=S010_sd,S101_sd=S101_sd,S110_sd=S110_sd
                            
  )))
  
}


truevalue=mymed(X,n,est,tseq,B=1000)
truevalue_in=c(truevalue$NDE,truevalue$NIE,truevalue$NIE1,truevalue$NIE2)
truevalue3=truevalue_in$x
###handle each simulated dataset

allest=NULL
allsd=NULL
myest=NULL

for (i in 1:200){
  myBest=NULL
  myest1=read.csv(paste("est2_multi",as.character(i),".csv",sep=""))
	myest=c(myest1$NDE,myest1$NIE,myest1$NIE1, myest1$NIE2)
	for (batch in 1:5){
	  skip_to_next <- FALSE
	  tryCatch({myBest_temp=read.csv(paste("Boot2_multilres",as.character(i),"batch=",as.character(batch),".csv",sep=""))}, error = function(e) { skip_to_next <<- TRUE})
	  
	  if(skip_to_next) { next }
	  if (!is.na(sum(myBest_temp))){
	  myBest=rbind(myBest,myBest_temp)
	  }
	  
	}
	
	if (!is.null(myBest)){
	myBest[abs(myBest)>1]<-NA
	myest.sd=apply(myBest,2,sd,na.rm=TRUE)
	allest=rbind(allest,myest)
	allsd=rbind(allsd,myest.sd)
	}
}

bias=apply(allest,2,mean,na.rm=TRUE)-truevalue3
ese=apply(allsd*is.finite(allsd),2,mean,na.rm=TRUE)
meSE=apply(allsd,2,median, na.rm=TRUE)
sd=apply(allest,2,sd,na.rm=TRUE)
cr=apply((abs(allest-rep(1,nrow(allest))%o%truevalue3)/allsd)<1.96,2,mean,na.rm=TRUE)

###might change file name to make it clear the parameter setting
write.csv(rbind(truevalue3,bias,sd,meSE,ese,cr),"sim2_multi_NIENDE_sum.csv")



