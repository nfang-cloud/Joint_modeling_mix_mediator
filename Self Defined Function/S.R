S<-function(Xi,Q,vi,ui,myprod,etaz,etax,delta1,delta2,tseq){
  
  myjump=myprod[[1]]
  myjump_gap=c(myjump[1],myjump[2:length(myjump)]-myjump[1:(length(myjump)-1)])
  myrisk=myprod[[2]]
  sss=length(tseq)
  S_seq=rep(0,sss)
  
  for (seq in 1:sss){
    
    j_prime=length(myjump)
    sum_j=matrix(0,j_prime+1,1)
    
    if (j_prime!=0){
      for (j_p in 1:j_prime){
        T=tseq[seq]-myjump[j_p]
        Cbj=myrisk[j_p]
        if (tseq[seq]<myjump[j_p]) {
          sum_j[j_p + 1,1]=Cbj*(tseq[seq]-max(0,myjump[j_p-1]))
          break
        }
        else if (tseq[seq]>=myjump[j_p]){sum_j[j_p + 1,1]=Cbj*myjump_gap[j_p]}
        
      }
    }
    
    
    sum_2=sum(sum_j)
    S=exp(-exp(etaz*Q + etax%*%Xi+delta1*vi+delta2*ui)*sum_2)
    
    S_seq[seq]=S
  }
  
  return(S_seq)
}
