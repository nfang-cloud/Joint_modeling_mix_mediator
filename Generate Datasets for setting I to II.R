library(MASS)
library(survival)
library(dplyr)
library(reda)
library(clusterGeneration)
library(tidyr)

#####load in self-defined functions#####
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

##Generate 200 datasets##
for (iter in 1:200){
set.seed(8888+iter)
  vi=rnorm(n)*sigmav
  ui=rnorm(n)*sigmau
  data=datagen_X(n,ui,vi,sigmae,tau,t.jump.w,alpha0,alphaz,
                  alphax,alphat,lambda.m,t.jump.m,betax,betaz,
   betaw,deltam,lambda.y,t.jump.y,etax,etaz,etaw,etam,delta1,delta2,cen=TRUE)  

 write.csv(data$alldata,paste("part3_sim1_alldata",as.character(iter),".csv",sep=""),row.names=F)
}



#################################################################################################################Setting II######################################

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
  
  #####Generate 200 datasets######

  for (iter in 1:200){
    set.seed(8888+iter)
    v_gamma=rgamma(n,shape = 2, scale = 1.8)
    vi=log(v_gamma)
    u_gamma=rgamma(n,shape = 1.5, scale = 1)
    ui=log(u_gamma)
   
    data=datagen_X(n,ui,vi,sigmae,tau,t.jump.w,alpha0,alphaz,
                   alphax,alphat,lambda.m,t.jump.m,betax,betaz,
                   betaw,deltam,lambda.y,t.jump.y,etax,etaz,etaw,etam,delta1,delta2,cen=TRUE)  
    
    write.csv(data$alldata,paste("part3_sim2_alldata",as.character(iter),".csv",sep=""),row.names=F)
  }
  
  