
Inf_criteria<-function(M){
post_sum<-posterior_summary(M, probs = c(0.025, 0.975),robust =FALSE)

new_solve<-function(A){return(chol2inv(chol(A)))}
beta_hat<- post_sum[3:(dim(x)[2]+2),1]
alpha_hat<- t(post_sum[1:2,1])
tau2_hat<- post_sum[(dim(x)[2]+3),1]
phi1_hat<- post_sum[(dim(x)[2]+4),1]
phi2_hat<- post_sum[(dim(x)[2]+5),1]


N=length(unique(id))
nsubj  <- length(unique(id))  # number of subjects
ntotal <- length(y)            # total number of observations
Time <- tapply(time, id, function(x) x)
data2 <- data.frame(cbind(id, x))
DM <- split(data2[, -1], data2$id)
YM <- tapply(y, id, function(x) x)
nobs <- tapply(id, id, function(x) length(x))
ntot<-dim(x)[1]
npar <- sum(beta_hat!=0)+sum(alpha_hat!=0) + 3   #3 for tau2,phi1, and phi2

tmp7 <- 0
tmp8 <- 0
Viii<- NULL

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
  Viii[[i]]<- tau2_hat*(phi1_hat*Ji+phi2_hat*Ri+Ii)
  # calculate mean and variance
  tmp7<- tmp7 + as.numeric(determinant(Viii[[i]],logarithm = TRUE)$'modulus')
  tmp_tau <- Yi-(as.matrix(Xi)) %*% (as.vector(beta_hat))-(as.matrix(Wi) %*% as.matrix(t(alpha_hat)))
  tmp8<- tmp8 + t(tmp_tau)%*% new_solve(Viii[[i]]) %*%(tmp_tau);
}

like<- -0.5*(tmp7+tmp8+dim(x)[1]*log(2*pi))
aic <- -2* like + 2*npar
bic <- -2* like + log(ntot)*npar
inf_cri<- cbind(aic,bic,like)
colnames(inf_cri)<-c("AIC","BIC","loglike")
return(inf_cri)
}


