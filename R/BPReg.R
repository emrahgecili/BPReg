

BPReg<- function(y,x,time,id,n.samples,r.init,delta.init,method=c("BEN","BAL","BR","BL")){

  library(MCMCpack)
  library(statmod)
  library(MASS)
  library(glmnet)
  library(matrixcalc)
  library(corpcor)
  library(msm)
  library(robustbase)
  library(mvtnorm)
  library(brms)
  library(MCMCvis)
  library(MCMCvis)




  if (method=="BEN")

  {
    BEN.MCMC(y,x,time,id,n.samples,r.init,delta.init)
  }

  else if (method=="BAL")

  {

    BAL.MCMC(y,x,time,id,n.samples,r.init,delta.init)

  }

  else if (method=="BL")

  {

    BL.MCMC(y,x,time,id,n.samples,r.init,delta.init)

  }
  else if (method=="BR")

  {

    BR.MCMC(y,x,time,id,n.samples,r.init,delta.init)

  }
}
