
library(MASS)
library(survival)
library(dplyr)
###load in self-defined functions###
source("Lambdainv.R")
source("mysimRec.R")
source("myprod.R")
source("datagenT01_M_multi.R")
source("S.R")
source("mymed.R")


###setting running parameters###


for (simseed in 0:999){
  
bootbatch=simseed%%5+1
simseed=simseed%/%5

####read in simulated datasets by simseed####

data=read.csv(paste("sim2_multi",as.character(simseed+1),".csv",sep=""),header=TRUE)
X2 <-data %>% group_by(ID) %>% filter(row_number(X2)==1) %>% select(ID,X2)
X=X2$X2
n=length(X)

####read in estimated  parameters by simseed####

est=read.table(paste("estIImulti",as.character(simseed+1),".txt",sep=""),header=FALSE,row.names=1)
est = as.data.frame(t(est))

est_m=c(est$log1_r1,est$log1_r2,est$log1_r3,est$log1_r4,
        est$log2_r1,est$log2_r2,est$log2_r3,est$log2_r4,
        est$log_h1,est$log_h2,est$log_h3,est$log_h4,
        est$betax1,est$betax2,est$betaz1,est$betaz2,est$deltax1,est$deltax2,
        est$etax,est$etaz,est$delta1,est$etam1,est$etam2)

####read in estimated variance by simseed####

estvar=read.table(paste("covIImulti",as.character(simseed+1),".txt",sep=""),header=FALSE,row.names=1,na.strings = ".")


###define parameters
tseq=seq(1,8,by=1)
###number of sampling 
B=10000
###number of bootstrap 
M=20

###get point estimate (save in a vector form)
myest=mymed(X,n,est,tseq,B)

myBest=matrix(data=NA,nrow=M,ncol=length(myest$NDE)*4)


###compute variance via Bootstrap
if (!is.na(sum(estvar))){
for (m in 1:M){
  
	set.seed(8888+(bootbatch-1)*M+m)
	estb=mvrnorm(1,mu=est_m,Sigma=estvar)
	est=as.data.frame(t(estb))
	mymed_m=mymed(X,n,est,tseq,B)
	myBest[m,]=c(mymed_m$NDE,mymed_m$NIE,mymed_m$NIE1,mymed_m$NIE2)
	
}
}
###save results

write.csv(myBest,paste("Boot2_multilres",as.character(simseed+1),"batch=", as.character(bootbatch),".csv",sep=""),row.names=FALSE)

}

