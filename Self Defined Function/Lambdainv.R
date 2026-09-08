###Find vector T such that Lambda(T)=s for a vector input s
Lambdainv<-function(s,t.jump,lambda,tau,risk){
  ###compute cumulative 
  n=length(risk)
  dup=0
  if (n==1){
    dup=1
    s=c(s,s)
    n=2*n
    risk=c(risk,risk)
  }
  Lambda=cumsum((c(t.jump,tau)-c(0,t.jump))*lambda)
  Lambda=risk%o%c(0,Lambda,Inf)-s%o%rep(1,length(lambda)+2)
  t.jump.exp=c(0,t.jump,tau)
  
  index=max.col(Lambda>=0,"first")
  T=rep(tau,n)
  ss=which(index<=length(t.jump.exp))
  if (length(ss)>1){
    upper=0-diag(Lambda[ss,index[ss]-1])
    lower=diag(Lambda[ss,index[ss]]-Lambda[ss,index[ss]-1])
    T[ss]=t.jump.exp[index[ss]-1]+(t.jump.exp[index[ss]]-t.jump.exp[index[ss]-1])*upper/lower
  }
  if (length(ss)==1){
    upper=0-Lambda[ss,index[ss]-1]
    lower=Lambda[ss,index[ss]]-Lambda[ss,index[ss]-1]
    T[ss]=t.jump.exp[index[ss]-1]+(t.jump.exp[index[ss]]-t.jump.exp[index[ss]-1])*upper/lower
    
  }
  if (dup==1){T=T[1]}
  T
}
