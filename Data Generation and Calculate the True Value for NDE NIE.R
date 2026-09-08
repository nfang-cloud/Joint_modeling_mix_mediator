library(MASS)
library(survival)
library(dplyr)
library(reda)
library(clusterGeneration)
library(tidyr)

###load in self-defined functions#####
source("Lambdainv.R")
source("mysimRec.R")
source("myprod.R")
source("datagen_X.R")


###################### Setting 1 Generate the dataset##################
sigmae<-1
n=200
k=1
tau=10
t.jump.m=c(2,4,6)
t.jump.y=c(2,4,6)
wT<-seq(0,8,by=2)
t.jump.w<-wT[-1]
lambda.m=c(0.2,0.5,0.6,0.45)
lambda.y=c(0.033,0.084,0.136,0.2)
alpha0=0
alphaz=0.6
alphax=0.2
alphat=0.2
betax=rep(0.2,k)
betaz=0.35
betaw=0.1
deltam=0.5
etaz=0.35
etax=rep(0.15,k)
etam=0.25
etaw=0.2
delta1=1
delta2=-0.5
sigmav=0.7
sigmau=1

##Generate 1000 datasets for potential T ##
for (iter in 1:1000){
set.seed(8888+iter)
  
  vi=rnorm(n)*sigmav
  ui=rnorm(n)*sigmau
 
  data=datagen_X(n,ui,vi,sigmae,tau,t.jump.w,alpha0,alphaz,
                  alphax,alphat,lambda.m,t.jump.m,betax,betaz,             betaw,deltam,lambda.y,t.jump.y,etax,etaz,etaw,etam,delta1,delta2,cen=TRUE)  

 write.csv(data$edat,paste("part3_sim_edat",as.character(iter),".csv",sep=""),row.names=F)
}


###Compute the real value of NDE NIE NIE_m NIE_w####
tseq=seq(2,8,by=2)
NIE_m=NIE_w=NDE=TE=matrix(data=NA,nrow=1000,ncol=length(tseq))
avgres=matrix(data=0,nrow=length(tseq),ncol=8) 
for (i in 1:1000){
edat=read.csv(paste("part3_sim_edat",as.character(i),".csv",sep=""))[,-1]
s=1
for (t in tseq){
ReturnT<-ifelse(edat>t,1,0)
S_T<-apply(ReturnT,2,mean,na.rm=TRUE)
avgres[s,]=S_T
s=s+1
  }

  S000_n=avgres[,1]
  S010_n=avgres[,2]
  S001_n=avgres[,3]
  S011_n=avgres[,4]
  S100_n=avgres[,5]
  S110_n=avgres[,6]
  S101_n=avgres[,7]
  S111_n=avgres[,8]
  
  NDE01=S100_n-S000_n
  NDE10=S011_n-S111_n
  
  NDE[i,]=(NDE01-NDE10)/2

  NIE01_m_1=S011_n-S001_n
  NIE01_m_0=S010_n-S000_n
  NIE10_m_1=S101_n-S111_n
  NIE10_m_0=S100_n-S110_n
  
  NIE01_w_1=S011_n-S010_n
  NIE01_w_0=S001_n-S000_n
  NIE10_w_1=S110_n-S111_n
  NIE10_w_0=S100_n-S101_n
  
  NIE01_m=(NIE01_m_1+NIE01_m_0)/2
  NIE01_w=(NIE01_w_1+NIE01_w_0)/2
  NIE10_m=(NIE10_m_1+NIE10_m_0)/2
  NIE10_w=(NIE10_w_1+NIE10_w_0)/2
  
  NIE_m[i,]=(NIE01_m-NIE10_m)/2
  NIE_w[i,]=(NIE01_w-NIE10_w)/2
  TE[i,]=NDE[i,]+NIE_m[i,]+NIE_w[i,]
}

NDE_SE=apply(NDE,2,sd,na.rm=TRUE)/sqrt(1000) 
NIE_m_SE=apply(NIE_m,2,sd,na.rm=TRUE)/sqrt(1000)
NIE_w_SE=apply(NIE_w,2,sd,na.rm=TRUE)/sqrt(1000)
NDE_Mean=apply(NDE,2,mean,na.rm=TRUE)
NIE_m_Mean=apply(NIE_m,2,mean,na.rm=TRUE)
NIE_w_Mean=apply(NIE_w,2,mean,na.rm=TRUE)

  truevalue1=c(NDE_Mean,NIE_m_Mean,NIE_w_Mean)
  write.csv(truevalue1,"true model1 NIE NDE_no SE.csv")
  truevalue1_SE=c(NDE_SE,NIE_m_SE,NIE_w_SE)
  write.csv(truevalue1_SE,"true model1 NIE NDE SE.csv")







###################### Setting 2 Generate the dataset##################
sigmae<-1
n=200
k=1
tau=10
t.jump.m=c(2,4,6)
t.jump.y=c(2,4,6)
wT<-seq(0,8,by=2)
t.jump.w<-wT[-1]
lambda.m=c(0.2,0.5,0.6,0.45)
lambda.y=c(0.033,0.084,0.136,0.2)
alpha0=0
alphaz=0.6
alphax=0.2
alphat=0.2
betax=rep(0.2,k)
betaz=0.35
betaw=0.1
deltam=0.5
etaz=0.35
etax=rep(0.15,k)
etam=0.25
etaw=0.2
delta1=1
delta2=-0.5
sigmav=0.7
sigmau=1

##Generate 1000 datasets for potential T ##
for (iter in 1:1000){
set.seed(8888+iter)
  
 v_gamma=rgamma(n,shape = 2, scale = 1.8)
 vi=log(v_gamma)
 u_gamma=rgamma(n,shape = 1.5, scale = 1)
 ui=log(u_gamma)
 
  data=datagen_X(n,ui,vi,sigmae,tau,t.jump.w,alpha0,alphaz,
                  alphax,alphat,lambda.m,t.jump.m,betax,betaz,             betaw,deltam,lambda.y,t.jump.y,etax,etaz,etaw,etam,delta1,delta2,cen=TRUE)  

 write.csv(data$edat,paste("part3_sim2_edat",as.character(iter),".csv",sep=""),row.names=F)
}


###Compute the real value of NDE NIE NIE_m NIE_w####
tseq=seq(2,8,by=2)
NIE_m=NIE_w=NDE=TE=matrix(data=NA,nrow=1000,ncol=length(tseq))
avgres=matrix(data=0,nrow=length(tseq),ncol=8) 
for (i in 1:1000){
edat=read.csv(paste("part3_sim2_edat",as.character(i),".csv",sep=""))[,-1]
s=1
for (t in tseq){
ReturnT<-ifelse(edat>t,1,0)
S_T<-apply(ReturnT,2,mean,na.rm=TRUE)
avgres[s,]=S_T
s=s+1
  }

  S000_n=avgres[,1]
  S010_n=avgres[,2]
  S001_n=avgres[,3]
  S011_n=avgres[,4]
  S100_n=avgres[,5]
  S110_n=avgres[,6]
  S101_n=avgres[,7]
  S111_n=avgres[,8]
  
  NDE01=S100_n-S000_n
  NDE10=S011_n-S111_n
  
  NDE[i,]=(NDE01-NDE10)/2

  NIE01_m_1=S011_n-S001_n
  NIE01_m_0=S010_n-S000_n
  NIE10_m_1=S101_n-S111_n
  NIE10_m_0=S100_n-S110_n
  
  NIE01_w_1=S011_n-S010_n
  NIE01_w_0=S001_n-S000_n
  NIE10_w_1=S110_n-S111_n
  NIE10_w_0=S100_n-S101_n
  
  NIE01_m=(NIE01_m_1+NIE01_m_0)/2
  NIE01_w=(NIE01_w_1+NIE01_w_0)/2
  NIE10_m=(NIE10_m_1+NIE10_m_0)/2
  NIE10_w=(NIE10_w_1+NIE10_w_0)/2
  
  NIE_m[i,]=(NIE01_m-NIE10_m)/2
  NIE_w[i,]=(NIE01_w-NIE10_w)/2
  TE[i,]=NDE[i,]+NIE_m[i,]+NIE_w[i,]
}

NDE_SE=apply(NDE,2,sd,na.rm=TRUE)/sqrt(1000) 
NIE_m_SE=apply(NIE_m,2,sd,na.rm=TRUE)/sqrt(1000)
NIE_w_SE=apply(NIE_w,2,sd,na.rm=TRUE)/sqrt(1000)
NDE_Mean=apply(NDE,2,mean,na.rm=TRUE)
NIE_m_Mean=apply(NIE_m,2,mean,na.rm=TRUE)
NIE_w_Mean=apply(NIE_w,2,mean,na.rm=TRUE)

  truevalue2=c(NDE_Mean,NIE_m_Mean,NIE_w_Mean)
  write.csv(truevalue2,"true model2 NIE NDE_no SE.csv")
  truevalue2_SE=c(NDE_SE,NIE_m_SE,NIE_w_SE)
  write.csv(truevalue2_SE,"true model2 NIE NDE SE.csv")

