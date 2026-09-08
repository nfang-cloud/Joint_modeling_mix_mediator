####Function used to generate recurrent events M_Q1 and M_Q3
datagenT01_M_multi<-function(X,T0,T1,vi,
                           tau,betaz1,betaz2,betax1,betax2,lambda.m1,lambda.m2,
                           t.jump.m1,t.jump.m2,delta_X1, delta_X2){
  
  ####Simulate M^z(t)
  ID=1:n
  beta.all1=c(betaz1,betax1,delta_X1)
  risk1_0=c(exp(cbind(T0,X,vi)%*%beta.all1))
  risk1_1=c(exp(cbind(T1,X,vi)%*%beta.all1))
  
  beta.all2=c(betaz2,betax2,delta_X2)
  risk2_0=c(exp(cbind(T0,X,vi)%*%beta.all2))
  risk2_1=c(exp(cbind(T1,X,vi)%*%beta.all2))
  
  ####generate survival time until all individual has follow-up at lest time tau
  M1_0=mysimRec(t.jump.m1,lambda.m1,tau,risk1_0,recurrent=TRUE)
  M1_1=mysimRec(t.jump.m1,lambda.m1,tau,risk1_1,recurrent=TRUE)
  
  M2_0=mysimRec(t.jump.m2,lambda.m2,tau,risk2_0,recurrent=TRUE)
  M2_1=mysimRec(t.jump.m2,lambda.m2,tau,risk2_1,recurrent=TRUE)
  
  M1_0$event_Tp=1
  M1_1$event_Tp=1
  M2_0$event_Tp=2
  M2_1$event_Tp=2
  M1_0$M=0
  M1_1$M=1
  M2_0$M=0
  M2_1$M=1
  M1_0[which(M1_0$event==0),]$event_Tp=0
  M1_1[which(M1_1$event==0),]$event_Tp=0
  M2_0[which(M2_0$event==0),]$event_Tp=0
  M2_1[which(M2_1$event==0),]$event_Tp=0
  
  M_00=rbind(M1_0,M2_0)
  M_01=rbind(M1_0,M2_1)
  M_10=rbind(M1_1,M2_0)
  M_11=rbind(M1_1,M2_1)
  
  M_00=M_00[order(M_00$ID,M_00$stoptime),]
  M_01=M_01[order(M_01$ID,M_01$stoptime),]
  M_10=M_10[order(M_10$ID,M_10$stoptime),]
  M_11=M_11[order(M_11$ID,M_11$stoptime),]
  
  M_00=M_00[!duplicated( M_00[ , c("ID","stoptime","event","event_Tp")]),]
  M_01=M_01[!duplicated( M_01[ , c("ID","stoptime","event","event_Tp")]),]
  M_10=M_10[!duplicated( M_10[ , c("ID","stoptime","event","event_Tp")]),]
  M_11=M_11[!duplicated( M_11[ , c("ID","stoptime","event","event_Tp")]),]

  return(list(M_00,M_01,M_10,M_11))
}

