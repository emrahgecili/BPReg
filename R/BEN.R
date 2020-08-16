
# BEN MCMC with lambda as a parameter
BEN.MCMC <- function(y, x,time,id,n.samples,r.init,delta.init){


  new_solve<-function(A){return(chol2inv(chol(A)))}

  func_invGau<-function(theta,chi){ chisq1<- rnorm(1, 0, 1)^2;
  y1    <- theta + 0.5*theta/chi * ( theta*chisq1 - sqrt(4*theta*chi*chisq1 + theta^2*chisq1^2) );
  y2    <- theta^2/y1;
  out_1 <- runif(1) < theta/(theta+y1);
  value <- out_1*y1+(1-out_1)*y2;
  return(value);}

  nsubj  <- length(unique(id))   # number of subjects
  ntotal <- length(y)            # total number of observations
  Time <- tapply(time, id, function(x) x)
  data2 <- data.frame(cbind(id, x))
  DM <- split(data2[, -1], data2$id)
  YM <- tapply(y, id, function(x) x)

  nobs <- tapply(id, id, function(x) length(x))
  n <- nrow(x)
  p <- ncol(x)
  #set initial values from
  beta <- solve(crossprod(x) + diag(p), crossprod(x, y))
  tau2 <- sum((y - x %*% beta)^2) / (n - p - 1)
  a2 <- rep(1, p)
  lambda12 <- 1
  lambda2 <- 1     #comes into account for EN
  alpha<- c(mean(y),1)

  ##now assume that r1=r2=1, del1=del2=1
  r <- r.init
  delta <- delta.init
  phi1<- rtnorm(1, 1,1,0,Inf)
  phi2<- rtnorm(1, 1,1,0,Inf)

  can.sd <- rep(0.1,2)
  acc <- rep(0,2)
  att <- 0

  keep <- matrix(0, n.samples, 2 * p + 4+1+2) #+1 for lambda2 for EN paramter

  for (s in 1 : n.samples){

    #########   updating alpha   ###############################################
    tmpa1 <- 0;
    tmpa2 <- 0;
    Vii<-NULL
    Viii<-NULL
    for (i in 1:nsubj)
    {
      # calculate mean and variance
      Timei <- Time[[i]]
      ni    <- nobs[i]
      Wi<- cbind(rep(1,ni),Timei)

      Ji <- matrix(1, ni, ni)
      Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))
      Ii <- diag(ni)
      Vii[[i]]<- phi1*Ji+phi2*Ri+Ii      #vi in derivation is tau2*Vii

      Xi <- DM[[i]];
      Yi <-  YM[[i]];
      alpha_ma<- tau2*Vii[[i]]
      tmpa  <- new_solve(alpha_ma) %*% as.matrix(Wi);
      tmpa1 <- tmpa1 + t(Yi-(as.matrix(Xi) %*% as.vector(beta))) %*% tmpa;
      tmpa2 <- tmpa2 + t(Wi) %*% tmpa;
      #emr<-emr+((new_solve(t(Wi)%*%new_solve(alpha_ma)%*% Wi)%*% t(Wi))%*%as.numeric((t(Yi-(as.matrix(Xi) %*% as.vector(beta))))%*%
      #         new_solve(alpha_ma)%*%(Yi-(as.matrix(Xi)%*% as.vector(beta)))))
    }

    Var_alp <- new_solve(tmpa2);
    #Var_j<- make.positive.definite(Var_j, tol=1e-3)
    Mu_alp <- Var_alp %*% t(tmpa1);
    alpha <- rmvnorm(1, Mu_alp,Var_alp);

    #########   updating Beta   ###############################################
    tmp2 <- 0;
    tmp3 <- 0;
    Vii<-NULL
    for (i in 1:nsubj)
    {
      # calculate mean and variance
      Timei <- Time[[i]]
      ni    <- nobs[i]
      Wi<- cbind(rep(1,ni),Timei)

      Ji <- matrix(1, ni, ni)
      Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))
      Ii <- diag(ni)
      Vii[[i]]<- phi1*Ji+phi2*Ri+Ii      #vi in derivation is tau2*Vii

      Xi <- DM[[i]];
      Yi <-  YM[[i]];
      tauma<- tau2*Vii[[i]]
      tmp  <- new_solve(tauma) %*% as.matrix(Xi);
      tmp3 <- tmp3 + t(Yi-(as.matrix(Wi) %*% as.matrix(t(alpha)))) %*% tmp;
      tmp2 <- tmp2 + t(Xi) %*% tmp;
    }

    Var_j <- new_solve(new_solve(tau2*diag(((1/a2)+lambda2)^-1))+tmp2);
    #Var_j<- make.positive.definite(Var_j, tol=1e-3)
    Mu_j <- Var_j %*% t(tmp3);
    beta <- rmvnorm(1, Mu_j,Var_j);

    ############# updating a2 #############
    for (j in 1:p)
    {
      Inva2_1 <- sqrt(lambda12 * tau2)/(abs(beta[j]));
      a2[j]   <- 1/func_invGau( Inva2_1, lambda12);
    }

    # update lambda1^2
    lambda12 <-rgamma(1, ncol(x)+r, sum(a2)/2 + delta)

    # update lambda2
    lambda2 <-rgamma(1, ncol(x)/2+r, sum(beta^2)/(2*tau2) + delta)


    ######### updating tau2 ########################################################################
    tau2_scale <- 0;
    for(i in 1:nsubj)
    {
      Yi<- YM[[i]];  #can be taken out no need these again and again
      Xi<- DM[[i]];
      Timei <- Time[[i]]
      ni    <- nobs[i]
      Wi<- cbind(rep(1,ni),Timei)

      Ji <- matrix(1, ni, ni)
      Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))
      Ii <- diag(ni)
      Vii[[i]]<- phi1*Ji+phi2*Ri+Ii      #vi in derivation is tau2*Vii
      #calculate scale parameter
      tmp_tau <- Yi-(as.matrix(Xi)) %*% (as.vector(beta))-(as.matrix(Wi) %*% as.matrix(t(alpha)))
      tau2_scale <- tau2_scale + t(tmp_tau)%*% new_solve(Vii[[i]]) %*%(tmp_tau);
    }
    tau2_sc <- tau2_scale+(t(as.vector(beta)) %*% diag(1/a2+lambda2,p)) %*% as.vector(beta)
    tau2_sc <- tau2_sc/(sum(nobs)+p);
    tau2 <- 1/rchisq(1, sum(nobs)+p)*(sum(nobs)+p)*tau2_sc;
    tau2 <- tau2[1,1];

    #####################
    # Update phi1
    log.post <- function(Yi, Xi,Wi, beta,alpha, tau2,phi1,phi2){
      tmp7 <- 0;
      tmp8 <- 0;

      for (i in 1:nsubj)
      {
        Yi<- YM[[i]];
        Xi<- DM[[i]];
        Timei <- Time[[i]]
        ni    <- nobs[i]
        Wi<- cbind(rep(1,ni),Timei)

        Ji <- matrix(1, ni, ni)
        Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))
        Ii <- diag(ni)
        Viii[[i]]<- tau2*(phi1*Ji+phi2*Ri+Ii)
        # calculate mean and variance
        tmp7<- tmp7 + as.numeric(determinant(Viii[[i]],logarithm = TRUE)$'modulus')
        tmp_tau <- Yi-(as.matrix(Xi)) %*% (as.vector(beta))-(as.matrix(Wi) %*% as.matrix(t(alpha)))
        tmp8 <- tmp8 + t(tmp_tau)%*% new_solve(Viii[[i]]) %*%(tmp_tau);
      }

      like<- -0.5*(tmp7+tmp8)
      return(like)
    }

    log.cur.post <- log.post(Yi,Xi,Wi,beta,alpha,tau2,phi1,phi2)
    canphi1 <- rtnorm(1, phi1,1,0,Inf)
    log.can.post <- log.post(Yi, Xi, Wi,beta,alpha, tau2,canphi1,phi2)
    R <- exp(log.can.post - log.cur.post)
    if (runif(1) < R){
      phi1 <- canphi1
      acc[1] <- acc[1] + 1}


    ############################################
    # Update phi2
    log.post2 <- function(Yi, Xi, Wi,beta,alpha, tau2,phi1,phi2){
      tmp9 <- 0;
      tmp10 <- 0;

      for (i in 1:nsubj)
      {
        Yi<- YM[[i]];
        Xi<- DM[[i]];
        Timei <- Time[[i]]
        ni    <- nobs[i]
        Wi<- cbind(rep(1,ni),Timei)

        Ji <- matrix(1, ni, ni)
        Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))
        Ii <- diag(ni)
        Viii[[i]]<- tau2*(phi1*Ji+phi2*Ri+Ii)
        # calculate mean and variance
        tmp9<- tmp9 + as.numeric(determinant(Viii[[i]],logarithm = TRUE)$'modulus')
        tmp_tau <- Yi-((as.matrix(Xi)) %*% (as.vector(beta)))-(as.matrix(Wi) %*% as.matrix(t(alpha)))
        tmp10 <- tmp10 + t(tmp_tau)%*% new_solve(Viii[[i]]) %*%(tmp_tau);
      }

      like2<- -0.5*(tmp9+tmp10)
      return(like2)
    }

    log.cur.post2 <- log.post2(Yi,Xi,Wi,beta,alpha,tau2,phi1,phi2)
    canphi2 <- exp(rnorm(1, log(phi2),can.sd[2]))
    log.can.post2 <- log.post2(Yi, Xi,Wi, beta,alpha, tau2,phi1,canphi2)
    R <- exp(log.can.post2 - log.cur.post2)
    if (runif(1) < R){
      phi2 <- canphi2
      acc[2] <- acc[2] + 1
    }

    if (i < 3000 & att > 50){
      can.sd[(acc / att) < 0.3] <- can.sd[(acc / att) < 0.3] * 0.8
      can.sd[(acc / att) > 0.6] <- can.sd[(acc / att) > 0.6] * 1.2

      acc <- rep(0, p + 1)
      att <- 0
    }

    keep[s,] <- c(beta, tau2, lambda12,a2,phi1,phi2,lambda2,alpha)

  }

  list(beta = keep[ , 1 : p], tau2 = keep[ ,p+1], lambda12 = keep[ ,p+2],
       a2 = keep[ , (p+3): (2 * p + 2)],phi1 = keep[ , (2*p+3)],phi2 = keep[ , (2*p+4)],lambda2 = keep[ , (2*p+4+1)],alpha=keep[,(2*p+4+2):(2*p+4+3)])

}
