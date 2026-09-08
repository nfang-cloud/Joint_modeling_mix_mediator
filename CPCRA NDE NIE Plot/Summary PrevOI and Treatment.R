
library(sas7bdat)
library(dplyr)
library(survival)

MM=20
tseq=seq(0,15*30,by=30)

#########################################################
##########First summary the real estimates of NIE/NDE####
#########################################################
resPrevOI_NDE=matrix(data=NA,nrow=MM,ncol=length(tseq))

resPrevOI_NIE_m=matrix(data=NA,nrow=MM,ncol=length(tseq))
resPrevOI_NIE_w=matrix(data=NA,nrow=MM,ncol=length(tseq))
resPrevOI_TE=matrix(data=NA,nrow=MM,ncol=length(tseq))
res_PrevOI=NULL
for (mm in 1:MM){
  res_PrevOI=read.csv(paste("Res_Mix_PrevOI_new",as.character(1111),"batch=",as.character(mm-1),".csv",sep=""))
  resPrevOI_NDE[mm,]=c(0,res_PrevOI$NDE)
  resPrevOI_NIE_m[mm,]=c(0,res_PrevOI$NIE_m)
  resPrevOI_NIE_w[mm,]=c(0,res_PrevOI$NIE_w)
  resPrevOI_TE[mm,]=c(0,res_PrevOI$TE)
  
}

RPrevOI_NDE=apply(resPrevOI_NDE,2,mean,na.rm=TRUE)

RPrevOI_NIE_m=apply(resPrevOI_NIE_m,2,mean,na.rm=TRUE)
RPrevOI_NIE_w=apply(resPrevOI_NIE_w,2,mean,na.rm=TRUE)
RPrevOI_TE=apply(resPrevOI_TE,2,mean,na.rm=TRUE)


################################################
##########Summary the bootstrap of NIE/NDE####
################################################
######average the batch####
Ba_tmp=20
simseed=1000
tmpresPrevOI_NDE=matrix(data=NA,nrow=Ba_tmp,ncol=length(tseq))

tmpresPrevOI_NIE_m=matrix(data=NA,nrow=Ba_tmp,ncol=length(tseq))
tmpresPrevOI_NIE_w=matrix(data=NA,nrow=Ba_tmp,ncol=length(tseq))
tmpresPrevOI_TE=matrix(data=NA,nrow=Ba_tmp,ncol=length(tseq))
AresPrevOI_NDE=matrix(data=NA,nrow=simseed,ncol=length(tseq))

AresPrevOI_NIE_m=matrix(data=NA,nrow=simseed,ncol=length(tseq))
AresPrevOI_NIE_w=matrix(data=NA,nrow=simseed,ncol=length(tseq))
AresPrevOI_TE=matrix(data=NA,nrow=simseed,ncol=length(tseq))
SEresPrevOI_NDE=matrix(data=NA,nrow=simseed,ncol=length(tseq))

SEresPrevOI_NIE_m=matrix(data=NA,nrow=simseed,ncol=length(tseq))
SEresPrevOI_NIE_w=matrix(data=NA,nrow=simseed,ncol=length(tseq))
SEresPrevOI_TE=matrix(data=NA,nrow=simseed,ncol=length(tseq))


BresPrevOI_tmp=NULL
for (sim in 1:simseed){
  for (mm in 1:Ba_tmp){
    skip_to_next <- FALSE
    tryCatch({ BresPrevOI_tmp=read.csv(paste("Bres_Mix_PrevOI_new",as.character(1111+sim-1),"batch=",as.character(mm),".csv",sep=""))}, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }
    if (!is.na(sum(BresPrevOI_tmp))){
    tmpresPrevOI_NDE[mm,]=c(0,BresPrevOI_tmp$NDE)
    tmpresPrevOI_NIE_m[mm,]=c(0,BresPrevOI_tmp$NIE_m)
    tmpresPrevOI_NIE_w[mm,]=c(0,BresPrevOI_tmp$NIE_w)
    tmpresPrevOI_TE[mm,]=c(0,BresPrevOI_tmp$TE)
    }
  }
  AresPrevOI_NDE[sim,]=apply(tmpresPrevOI_NDE,2,mean,na.rm=TRUE)
  
  AresPrevOI_NIE_m[sim,]=apply(tmpresPrevOI_NIE_m,2,mean,na.rm=TRUE)
  AresPrevOI_NIE_w[sim,]=apply(tmpresPrevOI_NIE_w,2,mean,na.rm=TRUE)
  AresPrevOI_TE[sim,]=apply(tmpresPrevOI_TE,2,mean,na.rm=TRUE)
  SEresPrevOI_NDE[sim,]=apply(tmpresPrevOI_NDE,2,sd,na.rm=TRUE)
  
  SEresPrevOI_NIE_m[sim,]=apply(tmpresPrevOI_NIE_m,2,sd,na.rm=TRUE)
  SEresPrevOI_NIE_w[sim,]=apply(tmpresPrevOI_NIE_w,2,sd,na.rm=TRUE)
  SEresPrevOI_TE[sim,]=apply(tmpresPrevOI_TE,2,sd,na.rm=TRUE)
  
}


SDPrevOI_NDE=apply(AresPrevOI_NDE,2,sd,na.rm=TRUE)
SDPrevOI_NIE_m=apply(AresPrevOI_NIE_m,2,sd,na.rm=TRUE)
SDPrevOI_NIE_w=apply(AresPrevOI_NIE_w,2,sd,na.rm=TRUE)
SDPrevOI_TE=apply(AresPrevOI_TE,2,sd,na.rm=TRUE)

SENDE_PrevOI=apply(SEresPrevOI_NDE,2,mean,na.rm=TRUE)
SENIE_m_PrevOI=apply(SEresPrevOI_NIE_m,2,mean,na.rm=TRUE)
SENIE_w_PrevOI=apply(SEresPrevOI_NIE_w,2,mean,na.rm=TRUE)
SETE_PrevOI=apply(SEresPrevOI_TE,2,mean,na.rm=TRUE)

QLPrevOI_NDE=apply(AresPrevOI_NDE,2,quantile, probs=0.025,na.rm=TRUE)
QUPrevOI_NDE=apply(AresPrevOI_NDE,2,quantile, probs=0.975,na.rm=TRUE)

QLPrevOI_NIE_m=apply(AresPrevOI_NIE_m,2,quantile, probs=0.025,na.rm=TRUE)
QUPrevOI_NIE_m=apply(AresPrevOI_NIE_m,2,quantile, probs=0.975,na.rm=TRUE)

QLPrevOI_NIE_w=apply(AresPrevOI_NIE_w,2,quantile, probs=0.025,na.rm=TRUE)
QUPrevOI_NIE_w=apply(AresPrevOI_NIE_w,2,quantile, probs=0.975,na.rm=TRUE)

QLPrevOI_TE=apply(AresPrevOI_TE,2,quantile, probs=0.025,na.rm=TRUE)
QUPrevOI_TE=apply(AresPrevOI_TE,2,quantile, probs=0.975,na.rm=TRUE)

##################################################
#####################Plot########################

##################Alternative methods#################
covdata=read.sas7bdat("covar.sas7bdat")
covdata$trt=covdata$RANDGRP-1
covdata$gender=covdata$GENDER-1
covdata$hemobl=covdata$HEMOBL-12
covdata$stratum=covdata$STRATUM-1
covdata$prevoi=covdata$PREVOI
Z=covdata[,"prevoi"]
surv_PrevOI=NULL
surv_PrevOI=read.sas7bdat("surv.sas7bdat")
surv_PrevOI$event=surv_PrevOI$event*(1/2)
surv_PrevOI$PrevOI=Z

############Cox#############
surv_PrevOI=cbind(surv_PrevOI,covdata)
cox_fit_PrevOI <- coxph(Surv(fup*30,event)~PrevOI+trt+gender+hemobl+stratum, data=surv_PrevOI)
# Create the new data  
data0_PrevOI=data1_PrevOI=surv_PrevOI
data0_PrevOI$PrevOI=0
data1_PrevOI$PrevOI=1

fit2PrevOI_1<-survfit(cox_fit_PrevOI,newdata=data1_PrevOI)
fit2PrevOI_0<-survfit(cox_fit_PrevOI,newdata=data0_PrevOI)
CoxPrevOI_S1<-summary(fit2PrevOI_1, time=tseq)
CoxPrevOI_S0<-summary(fit2PrevOI_0, time=tseq)
fit2_PrevOI<-apply(CoxPrevOI_S1[["surv"]]-CoxPrevOI_S0[["surv"]],1,mean,na.rm=TRUE)

######Treatment#########################

MM=20
tseq=seq(0,15*30,by=30)

#########################################################
##########First summary the real estimates of NIE/NDE####
#########################################################
res_NDE=matrix(data=NA,nrow=MM,ncol=length(tseq))

res_NIE_m=matrix(data=NA,nrow=MM,ncol=length(tseq))
res_NIE_w=matrix(data=NA,nrow=MM,ncol=length(tseq))
res_TE=matrix(data=NA,nrow=MM,ncol=length(tseq))
res=NULL

for (mm in 1:MM){
  res=read.csv(paste("Res_Mix_trt_new",as.character(1111),"batch=",as.character(mm-1),".csv",sep=""))
  res_NDE[mm,]=c(0,res$NDE)

  res_NIE_m[mm,]=c(0,res$NIE_m)
  res_NIE_w[mm,]=c(0,res$NIE_w)
  res_TE[mm,]=c(0,res$TE)
  
}

R_NDE=apply(res_NDE,2,mean,na.rm=TRUE)

R_NIE_m=apply(res_NIE_m,2,mean,na.rm=TRUE)
R_NIE_w=apply(res_NIE_w,2,mean,na.rm=TRUE)
R_TE=apply(res_TE,2,mean,na.rm=TRUE)


################################################
##########Summary the bootstrap of NIE/NDE####
################################################
######average the batch####
Ba_tmp=10
simseed=1000
tmpres_NDE=matrix(data=NA,nrow=Ba_tmp,ncol=length(tseq))

tmpres_NIE_m=matrix(data=NA,nrow=Ba_tmp,ncol=length(tseq))
tmpres_NIE_w=matrix(data=NA,nrow=Ba_tmp,ncol=length(tseq))
tmpres_TE=matrix(data=NA,nrow=Ba_tmp,ncol=length(tseq))
Ares_NDE=matrix(data=NA,nrow=simseed,ncol=length(tseq))

Ares_NIE_m=matrix(data=NA,nrow=simseed,ncol=length(tseq))
Ares_NIE_w=matrix(data=NA,nrow=simseed,ncol=length(tseq))
Ares_TE=matrix(data=NA,nrow=simseed,ncol=length(tseq))
SEres_NDE=matrix(data=NA,nrow=simseed,ncol=length(tseq))

SEres_NIE_m=matrix(data=NA,nrow=simseed,ncol=length(tseq))
SEres_NIE_w=matrix(data=NA,nrow=simseed,ncol=length(tseq))
SEres_TE=matrix(data=NA,nrow=simseed,ncol=length(tseq))


Bres_tmp=NULL
for (sim in 1:simseed){
  for (mm in 1:Ba_tmp){
    skip_to_next <- FALSE
    tryCatch({ Bres_tmp=read.csv(paste("Bres_Mix_Trt_new",as.character(1111+sim-1),"batch=",as.character(mm),".csv",sep=""))}, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { next }
    if (!is.na(sum(Bres_tmp))){
    tmpres_NDE[mm,]=c(0,Bres_tmp$NDE)
    tmpres_NIE_m[mm,]=c(0,Bres_tmp$NIE_m)
    tmpres_NIE_w[mm,]=c(0,Bres_tmp$NIE_w)
    tmpres_TE[mm,]=c(0,Bres_tmp$TE)
    }
  }
  Ares_NDE[sim,]=apply(tmpres_NDE,2,mean,na.rm=TRUE)
  Ares_NIE_m[sim,]=apply(tmpres_NIE_m,2,mean,na.rm=TRUE)
  Ares_NIE_w[sim,]=apply(tmpres_NIE_w,2,mean,na.rm=TRUE)
  Ares_TE[sim,]=apply(tmpres_TE,2,mean,na.rm=TRUE)
  
  SEres_NDE[sim,]=apply(tmpres_NDE,2,sd,na.rm=TRUE)
  SEres_NIE_m[sim,]=apply(tmpres_NIE_m,2,sd,na.rm=TRUE)
  SEres_NIE_w[sim,]=apply(tmpres_NIE_w,2,sd,na.rm=TRUE)
  SEres_TE[sim,]=apply(tmpres_TE,2,sd,na.rm=TRUE)
  
}
SD_NDE=apply(Ares_NDE,2,sd,na.rm=TRUE)
SD_NIE_m=apply(Ares_NIE_m,2,sd,na.rm=TRUE)
SD_NIE_w=apply(Ares_NIE_w,2,sd,na.rm=TRUE)
SD_TE=apply(Ares_TE,2,sd,na.rm=TRUE)

SENDE=apply(SEres_NDE,2,mean,na.rm=TRUE)
SENIE_m=apply(SEres_NIE_m,2,mean,na.rm=TRUE)
SENIE_w=apply(SEres_NIE_w,2,mean,na.rm=TRUE)
SETE=apply(SEres_TE,2,mean,na.rm=TRUE)

QL_NDE=apply(Ares_NDE,2,quantile, probs=0.025,na.rm=TRUE)
QU_NDE=apply(Ares_NDE,2,quantile, probs=0.975,na.rm=TRUE)

QL_NIE_m=apply(Ares_NIE_m,2,quantile, probs=0.025,na.rm=TRUE)
QU_NIE_m=apply(Ares_NIE_m,2,quantile, probs=0.975,na.rm=TRUE)

QL_NIE_w=apply(Ares_NIE_w,2,quantile, probs=0.025,na.rm=TRUE)
QU_NIE_w=apply(Ares_NIE_w,2,quantile, probs=0.975,na.rm=TRUE)

QL_TE=apply(Ares_TE,2,quantile, probs=0.025,na.rm=TRUE)
QU_TE=apply(Ares_TE,2,quantile, probs=0.975,na.rm=TRUE)

##################################################
#####################Plot########################

##################Alternative methods#################
covdata=read.sas7bdat("covar.sas7bdat")
covdata$trt=covdata$RANDGRP-1
covdata$gender=covdata$GENDER-1
covdata$hemobl=covdata$HEMOBL-12
covdata$stratum=covdata$STRATUM
covdata$cd4bl=covdata$CD4BL/100
covdata$prevoi=covdata$PREVOI
Z=covdata[,"trt"]
surv=NULL
surv=read.sas7bdat("surv.sas7bdat")
surv$event=surv$event*(1/2)
surv$trt=Z



############Cox#############
surv=cbind(surv,covdata)
cox_fit <- coxph(Surv(fup*30,event)~trt+gender+hemobl+stratum+prevoi, data=surv)
# Create the new data  
data0=data1=surv
data0$trt=0
data1$trt=1

fit2_1<-survfit(cox_fit,newdata=data1)
fit2_0<-survfit(cox_fit,newdata=data0)
Cox_S1<-summary(fit2_1, time=tseq)
Cox_S0<-summary(fit2_0, time=tseq)
fit2<-apply(Cox_S1[["surv"]]-Cox_S0[["surv"]],1,mean,na.rm=TRUE)

###Quantile CI##

pdf("Effects_Mix.pdf")
par( mfrow= c(2,2) )
plot(c(R_NDE,R_TE)~c(tseq,tseq),type="n",ylim=c(-0.10,0.15)
     ,ylab="Effect (Difference in Survival Probability)",xlab="Time (days)",pch=8,cex.lab=1)
lines(R_NDE~tseq,col="red",lty=1,lwd=2)
lines(I(QL_NDE)~tseq,col="red",lty=2,lwd=2)
lines(I(QU_NDE)~tseq,col="red",lty=2,lwd=2)

lines(R_TE~tseq,col="black",lty=1,lwd=2)
lines(I(QL_TE)~tseq,col="black",lty=2,lwd=2)
lines(I(QU_TE)~tseq,col="black",lty=2,lwd=2)

lines(fit2~tseq,col="green",lty=1,lwd=2)
title(main = "A. Trt (ddC vs ddI)", 
      cex.main = 1.2,   font.main= 2)
legend("bottomleft",c("NDE","TE","Cox"),col=c("red","black","green"),lty=1,lwd=2,cex=0.8)
abline(h=0)

plot(c(R_NIE_m,R_NIE_w)~c(tseq,tseq),type="n",ylim=c(-0.03,0.03)
     ,ylab="Effect (Difference in Survival Probability)",xlab="Time (days)",pch=8,cex.lab=1)

lines(R_NIE_m~tseq,col="black",lty=1,lwd=2)
lines(I(QL_NIE_m)~tseq,col="black",lty=2,lwd=2)
lines(I(QU_NIE_m)~tseq,col="black",lty=2,lwd=2)

lines(R_NIE_w~tseq,col="purple",lty=1,lwd=2)
lines(I(QL_NIE_w)~tseq,col="purple",lty=2,lwd=2)
lines(I(QU_NIE_w)~tseq,col="purple",lty=2,lwd=2)
title(main = "B. Subgroup NIEs Trt (ddC vs ddI)", 
      cex.main = 1.2,   font.main= 2)

legend("bottomleft",c("NIE_PrevOI","NIE_CD4"),col=c("black","purple"),lty=1,lwd=2,,cex=0.8)
abline(h=0)



plot(c(RPrevOI_NDE,RPrevOI_TE)~c(tseq,tseq),type="n",ylim=c(-0.35,0.15)
     ,ylab="Effect (Difference in Survival Probability)",xlab="Time (days)",pch=8,cex.lab=1)
lines(RPrevOI_NDE~tseq,col="red",lty=1,lwd=2)
lines(I(QLPrevOI_NDE)~tseq,col="red",lty=2,lwd=2)
lines(I(QUPrevOI_NDE)~tseq,col="red",lty=2,lwd=2)

#lines(RPrevOI_NIE~tseq,col="blue",lty=1,lwd=2)
#lines(I(QLPrevOI_NIE)~tseq,col="blue",lty=2,lwd=2)
#lines(I(QUPrevOI_NIE)~tseq,col="blue",lty=2,lwd=2)

lines(RPrevOI_TE~tseq,col="black",lty=1,lwd=2)
lines(I(QLPrevOI_TE)~tseq,col="black",lty=2,lwd=2)
lines(I(QUPrevOI_TE)~tseq,col="black",lty=2,lwd=2)

lines(fit2_PrevOI~tseq,col="green",lty=1,lwd=2)
title(main = "C. PrevOI", 
      cex.main = 1.2,   font.main= 2)
legend("bottomleft",c("NDE","TE","Cox"),col=c("red","black","green"),lty=1,lwd=2,,cex=0.8)
abline(h=0)


plot(c(RPrevOI_NIE_m,RPrevOI_NIE_w)~c(tseq,tseq),type="n",ylim=c(-0.15,0.15)
     ,ylab="Effect (Difference in Survival Probability)",xlab="Time (days)",pch=8,cex.lab=1)

lines(RPrevOI_NIE_m~tseq,col="black",lty=1,lwd=2)
lines(I(QLPrevOI_NIE_m)~tseq,col="black",lty=2,lwd=2)
lines(I(QUPrevOI_NIE_m)~tseq,col="black",lty=2,lwd=2)


lines(RPrevOI_NIE_w~tseq,col="purple",lty=1,lwd=2)
lines(I(QLPrevOI_NIE_w)~tseq,col="purple",lty=2,lwd=2)
lines(I(QUPrevOI_NIE_w)~tseq,col="purple",lty=2,lwd=2)
title(main = "D. Subgroup NIEs PrevOI", 
      cex.main = 1.2,   font.main= 2)
legend("bottomleft",c("NIE_PrevOI","NIE_CD4"),col=c("black","purple"),lty=1,lwd=2,,cex=0.8)
abline(h=0)
dev.off()
