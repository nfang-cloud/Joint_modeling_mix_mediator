
####################################################################
##########################True parameters###########################
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

################Computing log() for lambda.m, lambda.y and sigmav, sigmau, sigmae###############

log_lam_m1=log(lambda.m[1])
log_lam_m2=log(lambda.m[2])
log_lam_m3=log(lambda.m[3])
log_lam_m4=log(lambda.m[4])

log_lam_y1=log(lambda.y[1])
log_lam_y2=log(lambda.y[2])
log_lam_y3=log(lambda.y[3])
log_lam_y4=log(lambda.y[4])

log_varv=log(sigmav^2)
log_vare=log(sigmae^2)
log_varu=log(sigmau^2)



##############################################################################
###################################Setting I###################################
True_est=cbind(log_lam_m1,log_lam_m2,log_lam_m3,log_lam_m4,
          log_lam_y1,log_lam_y2,log_lam_y3,log_lam_y4,
          betax,betaz,betaw,etaz,etax,etam,etaw,alpha0,
          alphaz,alphax, alphat, 
          delta1,delta2,deltam,log_varv,log_vare,log_varu)



#####Import the estimation and covariance######


est=NULL
var=NULL
for (iter in 1:200){
  est.next=NULL
  var.next=NULL
  cov=NULL
  est.next=read.table(paste("estI",as.character(iter),".txt",sep=""))
  est.next=as.data.frame(t(est.next))
  cov=read.table(paste("covI",as.character(iter),".txt",sep=""),na.strings=c("."))
  var.next=t(diag(as.matrix(cov)))
  est=rbind(est,est.next[2,])
  var=rbind(var,var.next)
}
names(est)<-c("log_lam_m1","log_lam_m2","log_lam_m3","log_lam_m4",
              "log_lam_y1","log_lam_y2","log_lam_y3","log_lam_y4",
              "betax","betaz","betaw","etaz",
              "etax","etam","etaw","alpha0","alphaz",
              "alphax", "alphat","delta1","delta2",
              "deltam","log_varv","log_vare","log_varu") 
names(var)<-c("log_lam_m1","log_lam_m2","log_lam_m3","log_lam_m4",
              "log_lam_y1","log_lam_y2","log_lam_y3","log_lam_y4",
              "betax","betaz","betaw","etaz",
              "etax","etam","etaw","alpha0","alphaz",
              "alphax", "alphat","delta1","delta2",
              "deltam","log_varv","log_vare","log_varu") 
    
est=as.data.frame(sapply(est, as.numeric))

################Computing bias, empirical SD, median SE, CR###############
est_mean <- colMeans(est)
est_mean=apply(est,2,mean)
est_bias=apply(est-rep(1,nrow(est))%o%c(True_est),2,mean)
est_SD=apply(est,2,sd)
est_CR=apply((abs(est-rep(1,nrow(est))%o%c(True_est))/sqrt(var))<1.96,2,mean,na.rm=TRUE)
est_meSE=apply(sqrt(var),2,median, na.rm=TRUE)
est_ESE=apply(sqrt(var),2,mean,na.rm=TRUE)
my_eva<-data.frame(Mean=est_mean,Bias=est_bias,SD=est_SD,meSE=est_meSE, ESE=est_ESE,CR=est_CR)
write.csv(my_eva,paste("part3_sim1_estbialog",".csv",sep=""),row.names = TRUE)


##############################################################################
###################################Setting II###################################
True_est=cbind(log_lam_m1,log_lam_m2,log_lam_m3,log_lam_m4,
            log_lam_y1,log_lam_y2,log_lam_y3,log_lam_y4,
            betax,betaz,betaw,etaz,etax,etam,etaw,alpha0,
            alphaz,alphax, alphat, 
            delta1,delta2,deltam,log_varv,log_vare,log_varu)

#####Import the estimation and covariance######

est=NULL
var=NULL
for (iter in 1:200){
    est.next=NULL
    var.next=NULL
    cov=NULL
    est.next=read.table(paste("estII",as.character(iter),".txt",sep=""))
    est.next=as.data.frame(t(est.next))
    cov=read.table(paste("covII",as.character(iter),".txt",sep=""),na.strings=c("."))
    var.next=t(diag(as.matrix(cov)))
    est=rbind(est,est.next[2,])
    var=rbind(var,var.next)
  }
  names(est)<-c("log_lam_m1","log_lam_m2","log_lam_m3","log_lam_m4",
                "log_lam_y1","log_lam_y2","log_lam_y3","log_lam_y4",
                "betax","betaz","betaw","etaz",
                "etax","etam","etaw","alpha0","alphaz",
                "alphax", "alphat","delta1","delta2",
                "deltam","log_varv","log_vare","log_varu") 
  names(var)<-c("log_lam_m1","log_lam_m2","log_lam_m3","log_lam_m4",
                "log_lam_y1","log_lam_y2","log_lam_y3","log_lam_y4",
                "betax","betaz","betaw","etaz",
                "etax","etam","etaw","alpha0","alphaz",
                "alphax", "alphat","delta1","delta2",
                "deltam","log_varv","log_vare","log_varu") 
  
  est=as.data.frame(sapply(est, as.numeric))

################Computing bias, empirical SD, median SE, CR###############
est_mean <- colMeans(est)
est_mean=apply(est,2,mean)
est_bias=apply(est-rep(1,nrow(est))%o%c(True_est),2,mean)
est_SD=apply(est,2,sd)
est_CR=apply((abs(est-rep(1,nrow(est))%o%c(True_est))/sqrt(var))<1.96,2,mean,na.rm=TRUE)
est_meSE=apply(sqrt(var),2,median, na.rm=TRUE)
est_ESE=apply(sqrt(var),2,mean,na.rm=TRUE)
my_eva<-data.frame(Mean=est_mean,Bias=est_bias,SD=est_SD,meSE=est_meSE, ESE=est_ESE,CR=est_CR)
write.csv(my_eva,paste("part3_sim2_estbialog",".csv",sep=""),row.names = TRUE)


