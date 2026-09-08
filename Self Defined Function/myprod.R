###product of two step function
myprod<-function(t.jump1,lambda1,t.jump2,lambda2){
  f1=stepfun(y=c(0,lambda1),x=c(0,t.jump1))
  f2=stepfun(y=c(0,lambda2),x=c(0,t.jump2))
  t.jump=sort(union(t.jump1,t.jump2))
  lambda=f1(c(0,t.jump))*f2(c(0,t.jump))
  return(list(t.jump,lambda))
}