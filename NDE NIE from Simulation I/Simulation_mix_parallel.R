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



####read in simulated datasets by simseed####
for (simseed in 0:199){

data=read.csv(paste("part3_sim_alldata",as.character(simseed+1),".csv",sep=""),header=TRUE)
X2 <-data %>% group_by(ID) %>% filter(row_number(X2)==1) %>% select(ID,X2)
X=X2$X2


####read in estimated parameter  by simseed####

est=read.table(paste("estI",as.character(simseed+1),".txt",sep=""),header=FALSE,row.names=1)
est = as.data.frame(t(est))


###define parameters
tseq=seq(2,8,by=2)  
###number of sampling 
B=10000
n=length(X)


###get point estimate (save in a vector form)
myest=mymed(X,n,est,tseq,B)

###save results

write.csv(myest,paste("est1_mix",as.character(simseed+1),".csv",sep=""),row.names=FALSE)
}