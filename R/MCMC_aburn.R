MCMC_aburn=function(burnin){
  Aft_burnin=cbind(MCMCsamp$alpha[(burnin + 1) : dim(MCMCsamp[[1]])[1], ],MCMCsamp$beta[(burnin + 1) :dim(MCMCsamp[[1]])[1], ],MCMCsamp$tau2[(burnin + 1) : dim(MCMCsamp[[1]])[1]],
        MCMCsamp$phi1[(burnin + 1) : dim(MCMCsamp[[1]])[1]],MCMCsamp$phi2[(burnin + 1) : dim(MCMCsamp[[1]])[1]])

  colnames(Aft_burnin)<-c("Int","age",colnames(x),"tau2","phi1","phi2")
  return(Aft_burnin)
  }
