### Implementation of DtACI method from https://arxiv.org/abs/2208.08401

vecZeroMax <- Vectorize(function(x){max(x,0)})
vecZeroMin <- Vectorize(function(x){min(x,0)})

pinball <- function(u,alpha){
  alpha*u - vecZeroMin(u)
}

dtaci <- function(betas,alpha,gammas,sigma=1/1000,eta=2.72,alphaInit=alpha,etaAdapt=FALSE,etaLookback=500){
  T <- length(betas)
  k <- length(gammas)
  meanAlphaSeq <- rep(0,T)
  meanErrSeq <- rep(0,T)
  meanGammas <- rep(0,T)
  lossSeq <- rep(0,T)
  
  expertAlphas <- rep(alphaInit,k)
  expertWs <- rep(1,k)
  expertProbs <- rep(1/k,k)
  for(t in 1:T){
    if(t > etaLookback & etaAdapt){
      eta <- sqrt((log(2*k*etaLookback) + 1)/sum(lossSeq[(t-etaLookback):(t-1)]^2))
    }
    meanAlphaSeq[t] <- sum(expertProbs*expertAlphas)
    meanErrSeq[t] <- as.numeric(meanAlphaSeq[t]>betas[t])
    meanGammas[t] <- sum(expertProbs*gammas)
    
    expertLosses <- pinball(betas[t] - expertAlphas,alpha)
    lossSeq[t] <- sum(expertLosses*expertProbs)
    
    expertAlphas <- expertAlphas + gammas*(alpha-as.numeric(expertAlphas>betas[t]))
    
    if(eta < Inf){
      expertBarWs <- expertWs*exp(-eta*expertLosses)
      expertNextWs <- (1-sigma)*expertBarWs/sum(expertBarWs) + sigma/k
      
      expertProbs <- expertNextWs/sum(expertNextWs)
      expertWs <- expertNextWs
    }else{
    }
    
  }
  return(list(meanAlphaSeq,meanErrSeq,meanGammas))
}

