
library(MASS)
library(survival)
library(dplyr)
library(reda)
library(clusterGeneration)


###Import True Values######   
truevalue_in=read.csv(paste("true model1 NIE NDE_no SE",".csv",sep=""))
truevalue1=truevalue_in$x
truevalue_SE_in=read.csv(paste("true model1 NIE NDE SE",".csv",sep=""))
truevalue1_SE=truevalue_SE_in$x

###Import the NDE NIE estimates#####

allest=NULL
allsd=NULL
myest=NULL

for (i in 1:200){
  myBest=NULL
  myest1=read.csv(paste("est1_mix",as.character(i),".csv",sep=""))[c(2,4,6,8),]
  myest=c(myest1$NDE,myest1$NIE_m,myest1$NIE_w)
  for (batch in 1:5){
    for (m in seq(2,20,by=2)){
      skip_to_next <- FALSE
      tryCatch({myBest_temp=read.csv(paste("Boot1_mix",as.character(i),"batch=",as.character(batch),"M=",as.character(m),".csv",sep=""))}, error = function(e) { skip_to_next <<- TRUE})
      
      if(skip_to_next) { next }
      if (!is.na(sum(myBest_temp))){
        myBest=rbind(myBest,myBest_temp)
      }
    }
  }
  
  
  if (!is.null(myBest)){
    myest.sd=apply(myBest,2,sd,na.rm=TRUE)
    allest=rbind(allest,myest)
    allsd=rbind(allsd,myest.sd)
  }
}

bias=apply(allest,2,mean,na.rm=TRUE)-truevalue1
ese=apply(allsd*is.finite(allsd),2,mean,na.rm=TRUE)
meSE=apply(allsd,2,median, na.rm=TRUE)
sd=apply(allest,2,sd,na.rm=TRUE)
cr=apply((abs(allest-rep(1,nrow(allest))%o%truevalue1)/allsd)<1.96,2,mean,na.rm=TRUE)

###might change file name to make it clear the parameter setting
write.csv(rbind(truevalue1,truevalue1_SE,bias,sd,meSE,ese,cr),"sim1_mix_NIENDE_sum.csv")




