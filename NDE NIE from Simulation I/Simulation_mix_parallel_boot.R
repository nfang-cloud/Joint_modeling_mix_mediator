
library(MASS)
library(survival)
library(dplyr)

###load in self-defined functions###
source("Lambdainv.R")
source("mysimRec.R")
source("myprod.R")
source("datagen_M_CD4.R")
source("S.R")
source("mymed.R")



###setting running parameters###


for (simseed in 0:999){

bootbatch=simseed%%5+1
simseed=simseed%/%5

####read in simulated datasets by simseed####

data=read.csv(paste("part3_sim_alldata",as.character(simseed+1),".csv",sep=""),header=TRUE)
X2 <-data %>% group_by(ID) %>% filter(row_number(X2)==1) %>% select(ID,X2)
X=X2$X2
n=length(X)


####read in estimated parameter by simseed####

est=read.table(paste("estI",as.character(simseed+1),".txt",sep=""),header=FALSE,row.names=1)
est = as.data.frame(t(est))
est_m=c(est$log_lam_m1,est$log_lam_m2,est$log_lam_m3,est$log_lam_m4,est$log_lam_y1,
est$log_lam_y2,est$log_lam_y3,est$log_lam_y4,est$betax,est$betaz,
est$betaw,est$etaz,est$etax,est$etam,est$etaw,est$alpha0,est$alphaz,
est$alphax,est$alphat,est$delta1,est$delta2,est$deltam,est$log_varv,est$log_vare,est$log_varu) 


####read in estimated variance by simseed####

estvar=read.table(paste("covI",as.character(simseed+1),".txt",sep=""),header=FALSE,row.names=1,na.strings = ".")


###define parameters

tseq=seq(2,8,by=2)
###number of sampling 
B=10000
###number of bootstrap 
M1=2
M2=0
###get point estimate (save in a vector form)
#myest=mymed(X,n,est,tseq,B)


myBest=matrix(data=NA,nrow=M1-M2,ncol=length(tseq)*3)


###compute variance via Bootstrap
if (!is.na(sum(estvar))){
  for (m in 1:(M1-M2)){
    set.seed(8888+(bootbatch-1)*20+m+M2)
    estb=mvrnorm(1,mu=est_m,Sigma=estvar)
    est=as.data.frame(t(estb))
    mymed_m=mymed(X,n,est,tseq,B)
    myBest[m,]=c(mymed_m$NDE,mymed_m$NIE_m,mymed_m$NIE_w)
  }
}
###save results


write.csv(myBest,paste("Boot1_mix",as.character(simseed+1),"batch=", as.character(bootbatch),"M=", as.character(M1),".csv",sep=""),row.names=FALSE)
}
