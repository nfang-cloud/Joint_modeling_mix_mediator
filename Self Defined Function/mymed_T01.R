####Mediation analysis
mymed_T01<-function(X,n,Z,T0,T1,
                    sigmav,sigmae,sigmau,tau,
                    betaz,betax,betaw,lambda.m,t.jump.m,
                    etaz,etax,etam,etaw,
                    alpha0,alphaz,alphax,alphat,wT,
                    delta1,delta2,deltam,
                    lambda.y,t.jump.y,tseq,B,seed){
  wT=seq(0,20,by=2)
  t.jump.w<-wT[-1]
  
  avgres=matrix(data=0,nrow=n,ncol=8*length(tseq))  
  avgressq=matrix(data=0,nrow=n,ncol=8*length(tseq))
  
  for (b in 1:B){
    
    set.seed(seed+b)
    vi=rnorm(n,0,sigmav)
    ui=rnorm(n,0,sigmau)
    
    #### Simulate W
    SS=length(wT)
    W0_tmp<-matrix(0,n,SS)
    W1_tmp<-matrix(0,n,SS)
    alpha.all=c(alpha0,alphaz,alphax,alphat,1,1)
    for (s in 1:SS){
      ei=rnorm(n)*sigmae
      W0_tmp[,s]=c(cbind(1,0,X,wT[s],ui,ei)%*%alpha.all)
      W1_tmp[,s]=c(cbind(1,1,X,wT[s],ui,ei)%*%alpha.all)
    }
    
    #Generate M

    M_T=datagenT01_M_CD4(X,T0,T1,ui,vi,tau,
                      W0_tmp,W1_tmp,t.jump.w,t.jump.m,
                      lambda.m,betax,betaz,betaw,deltam)
    R_00=M_T[[1]]
    R_11=M_T[[2]]
    
    tmpres=avgres

    S000b=S100b=S011b=S111b=S001b=S101b=S010b=S110b=matrix(ncol=length(tseq),nrow=n)
    
    for (i in 1:n){
      Xi=X[i,]
      R_00i=R_00[which(R_00$ID==i),]
      R_11i=R_11[which(R_11$ID==i),]
      W_0i=W0_tmp[i,]
      W_1i=W1_tmp[i,]
      
      lambda.risk0=exp((0:length(R_00i$stoptime))*etam)
      t.jump.risk0=R_00i$stoptime
      lambda.risk1=exp((0:length(R_11i$stoptime))*etam)
      t.jump.risk1=R_11i$stoptime
      
      lambda.riskw0=exp(W_0i*etaw)
      lambda.riskw1=exp(W_1i*etaw)
      lambda.risk00w=myprod(t.jump.w,lambda.riskw0,t.jump.risk0,lambda.risk0)[[2]]
      t.jump.risk00w=myprod(t.jump.w,lambda.riskw0,t.jump.risk0,lambda.risk0)[[1]]
      lambda.risk10w=myprod(t.jump.w,lambda.riskw1,t.jump.risk0,lambda.risk0)[[2]]
      t.jump.risk10w=myprod(t.jump.w,lambda.riskw1,t.jump.risk0,lambda.risk0)[[1]]
      lambda.risk01w=myprod(t.jump.w,lambda.riskw0,t.jump.risk1,lambda.risk1)[[2]]
      t.jump.risk01w=myprod(t.jump.w,lambda.riskw0,t.jump.risk1,lambda.risk1)[[1]]
      lambda.risk11w=myprod(t.jump.w,lambda.riskw1,t.jump.risk1,lambda.risk1)[[2]]
      t.jump.risk11w=myprod(t.jump.w,lambda.riskw1,t.jump.risk1,lambda.risk1)[[1]]
      
      my00=myprod(t.jump.y,lambda.y,t.jump.risk00w,lambda.risk00w)
      my10=myprod(t.jump.y,lambda.y,t.jump.risk10w,lambda.risk10w)	
      my01=myprod(t.jump.y,lambda.y,t.jump.risk01w,lambda.risk01w)
      my11=myprod(t.jump.y,lambda.y,t.jump.risk11w,lambda.risk11w)
      
      S000b=S(Xi,T0,vi[i],ui[i],my00, etaz,etax,delta1,delta2,tseq)#modified
      S100b=S(Xi,T1,vi[i],ui[i],my00, etaz,etax,delta1,delta2,tseq)#modified
      S011b=S(Xi,T0,vi[i],ui[i],my11, etaz,etax,delta1,delta2,tseq)#modified
      S111b=S(Xi,T1,vi[i],ui[i],my11, etaz,etax,delta1,delta2,tseq)#modified
      
      S001b=S(Xi,T0,vi[i],ui[i],my01, etaz,etax,delta1,delta2,tseq)#modified
      S101b=S(Xi,T1,vi[i],ui[i],my01, etaz,etax,delta1,delta2,tseq)#modified
      S010b=S(Xi,T0,vi[i],ui[i],my10, etaz,etax,delta1,delta2,tseq)#modified
      S110b=S(Xi,T1,vi[i],ui[i],my10, etaz,etax,delta1,delta2,tseq)#modified
      
      tmpres[i,]= c(S000b,S011b,S100b,S111b,S001b,S101b,S010b,S110b)
      
    }
    
    tmpressq=tmpres^2
    
    avgres=avgres+(tmpres-avgres)/b
    avgressq=avgressq+(tmpressq-avgressq)/b
    
  }
  avgsd=sqrt(avgressq-avgres^2)
  
  avg_S=apply(avgres,2,mean,na.rm=TRUE)
  sd_S=apply(avgres,2,sd,na.rm=TRUE)
  
  S000_n=avg_S[1:length(tseq)]
  S011_n=avg_S[(length(tseq)+1):(length(tseq)*2)]
  S100_n=avg_S[(length(tseq)*2+1):(length(tseq)*3)]
  S111_n=avg_S[(length(tseq)*3+1):(length(tseq)*4)]
  S001_n=avg_S[(length(tseq)*4+1):(length(tseq)*5)]
  S101_n=avg_S[(length(tseq)*5+1):(length(tseq)*6)]
  S010_n=avg_S[(length(tseq)*6+1):(length(tseq)*7)]
  S110_n=avg_S[(length(tseq)*7+1):(length(tseq)*8)]
  
  
  S000_sd=sd_S[1:length(tseq)]
  S011_sd=sd_S[(length(tseq)+1):(length(tseq)*2)]
  S100_sd=sd_S[(length(tseq)*2+1):(length(tseq)*3)]
  S111_sd=sd_S[(length(tseq)*3+1):(length(tseq)*4)]
  S001_sd=sd_S[(length(tseq)*4+1):(length(tseq)*5)]
  S101_sd=sd_S[(length(tseq)*5+1):(length(tseq)*6)]
  S010_sd=sd_S[(length(tseq)*6+1):(length(tseq)*7)]
  S110_sd=sd_S[(length(tseq)*7+1):(length(tseq)*8)]
  
  NDE01=S100_n-S000_n
  NDE10=S011_n-S111_n
  #NIE01=S011_n-S000_n
  #NIE10=S100_n-S111_n
  
  NDE=(NDE01-NDE10)/2
  #NIE=(NIE01-NIE10)/2
  
  
  NIE01_m_1=S011_n-S001_n
  NIE01_m_0=S010_n-S000_n
  NIE10_m_1=S101_n-S111_n
  NIE10_m_0=S100_n-S110_n
  
  NIE01_w_1=S011_n-S010_n
  NIE01_w_0=S001_n-S000_n
  NIE10_w_1=S110_n-S111_n
  NIE10_w_0=S100_n-S101_n
  
  NIE01_m=(NIE01_m_1+NIE01_m_0)/2
  NIE01_w=(NIE01_w_1+NIE01_w_0)/2
  NIE10_m=(NIE10_m_1+NIE10_m_0)/2
  NIE10_w=(NIE10_w_1+NIE10_w_0)/2
  
  NIE_m=(NIE01_m-NIE10_m)/2
  NIE_w=(NIE01_w-NIE10_w)/2
  TE=NDE+NIE_m+NIE_w
  
  return(as.data.frame(list(NDE=NDE,TE=TE,NIE_m=NIE_m,NIE_w=NIE_w,
                            S000=S000_n,S011=S011_n,S100=S100_n,S111=S111_n,
                            S001=S001_n,S010=S010_n,S101=S101_n,S110=S110_n,
                            S000_sd=S000_sd,S011_sd=S011_sd,S100_sd=S100_sd,S111_sd=S111_sd,
                            S001_sd=S001_sd,S010_sd=S010_sd,S101_sd=S101_sd,S110_sd=S110_sd
                            
  )))
  
}
