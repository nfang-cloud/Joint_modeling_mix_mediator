
library(MASS)
library(survival)
library(dplyr)
###load in self-defined functions###
source("Lambdainv.R")
source("mysimRec.R")
source("datagenT01_M_CD4.R")
source("myprod.R")
source("S.R")
source("mymed_T01.R")

####check real data###########
####read in data##############

covdata=read.csv("covdata.csv")

covdata$trt=covdata$RANDGRP-1
covdata$gender=covdata$GENDER-1
covdata$hemobl=covdata$HEMOBL-12
covdata$stratum=covdata$STRATUM-1
covdata$prevoi=covdata$PREVOI

Z=covdata[,"trt"]
T0=0
T1=1
tseq=seq(1,15,by=1)
X=as.matrix(covdata[,c("prevoi","gender","stratum","hemobl")])
tau=15
n=nrow(X)

############################Import the Estimation######################
 
 est=read.table(paste("partIII_estI_real_cd4100_quantile",".txt",sep=""),header=FALSE,row.names=1)
 est = as.data.frame(t(est))
 
t.jump.m=c(2.767,6.333, 10.2)
t.jump.y=c(4.583, 8.6, 11.617)


#######################est1 from the real Estimates###############
 
 
est1=c(est$log_lam_m1,est$log_lam_m2,est$log_lam_m3,est$log_lam_m4,
        est$log_lam_y1,est$log_lam_y2,est$log_lam_y3,est$log_lam_y4,
        est$betax1,est$betax2,est$betaz,est$betax3,est$betax4,est$betaw,
        est$etax1,est$etax2,est$etaz,est$etax3,est$etax4,est$etam,est$etaw,
        est$alpha0,est$alphax1,est$alphax2,est$alphaz,est$alphax3,est$alphax4,
        est$alphat,est$delta1,est$delta2,est$deltam,est$log_varv,est$log_vare,est$log_varu)
 
 

##Mediation Analysis##
 B=100000/20
 for (batch in 1:20){
Bres1_new=NULL
 Best=est1
 Bsigmav=sqrt(exp(Best[32]))
 Bsigmae=sqrt(exp(Best[33]))
 Bsigmau=sqrt(exp(Best[34]))
 
 Blambda.m=exp(Best[1:4])
 Blambda.y=exp(Best[5:8])
 
 Bbetaz=Best[11]
 Bbetax=Best[c(9:10,12:13)]
 Bbetaw=Best[14]
 
 Betaz=Best[17] 
 Betax=Best[c(15:16,18:19)]
 Betam=Best[20]
 Betaw=Best[21]
 
 Balpha0=Best[22]
 Balphaz=Best[25]
 Balphax=Best[c(23:24,26:27)]
 Balphat=Best[28]
 
 Bdelta1=Best[29]
 Bdelta2=Best[30]
 Bdeltam=Best[31]
 
 ###from the mymed function, need to add these besides the mymed function###
 Bres1_new=mymed_T01(X,n,Z,T0,T1,Bsigmav,Bsigmae,Bsigmau,tau,
                     Bbetaz,Bbetax,Bbetaw,Blambda.m,t.jump.m,
                     Betaz,Betax,Betam,Betaw,
                     Balpha0,Balphaz,Balphax,Balphat,wT,
                     Bdelta1,Bdelta2,Bdeltam,
                     Blambda.y,t.jump.y,tseq,B,seed=1111+batch*B)
 
 
 write.csv(Bres1_new,paste("Res_Mix_Trt_new",as.character(1111),"batch=",as.character(batch),".csv",sep=""),row.names=FALSE)

}
 