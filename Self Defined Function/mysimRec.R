
###length of lambda is of one longer than of jump.t

###using Cinlar's inversion Method to generate Non-homogeneous Poisson process
mysimRec<-function(t.jump,lambda,tau,risk,recurrent=TRUE){
  n=length(risk)
  ID=1:n
  s=rep(0,n)
  dat=NULL
  T=rep(0,n)
  
  sss=1:n
  while(length(sss)>0){
    u=runif(n)
    s=s-log(u)
    ID=sss
    T=Lambdainv(s[sss],t.jump,lambda,tau,risk[sss])
    sss=sss[which(T<tau)]
    if (!recurrent){sss=NULL}
    event=as.numeric(T<tau)
    tmp=data.frame(ID=ID,stoptime=T,event=event)
    dat=rbind(dat,tmp)
  }
  dat[order(dat$ID,dat$stoptime),]
  return(dat[order(dat$ID,dat$stoptime),])
}

