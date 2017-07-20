#####
swank.true.Ts<- function(n.trials, beta, mu1, tau1, rho1, tau2, nclusters = 1, cluster_spacing = 1, model = "linear"){
  
  ## nonlinearity
  ## threshold effects
  ## clustered T1s
  
 
  library(mvtnorm)
  
  cov1 <- rho1*tau1[1]*tau1[2]
  
  T1 <- NULL
  sizes <- rep(n.trials %/% nclusters, nclusters)
  sizes[nclusters] <- sizes[nclusters] + n.trials %% nclusters
  
  for(j in 1:length(sizes)){
    
    T1 <- rbind(T1, rmvnorm(sizes[j], mean = (mu1 + cluster_spacing * j)/nclusters, 
                            sigma = matrix(c(tau1[1], cov1, cov1, tau1[2])/nclusters, nrow = 2)))
    
  }
  
  if(model == "cosine") {
    
    T2 <- rnorm(n.trials, beta[1] + beta[2]*T1[,1] - 2 * cos(T1[,1]^2), tau2)
    
  } else if(model == "cubic-threshold") {
    
    n1 <- beta[1] + beta[2] * quantile(T1[,1], .25) + quantile(T1[,1], .25)^2 + quantile(T1[,1], .25)^3
    n2 <- quantile(T1[,1], .75)
    
    T2 <- rnorm(n.trials, ifelse(T1[,1] < quantile(T1[,1], .25), beta[1] + beta[2] * T1[,1] + T1[,1]^2 + T1[,1]^3,
                                 ifelse(T1[,1] < quantile(T1[,1], .75), n1,
                                        n1 - beta[2] * (T1[,1] - n2) - (T1[,1] - n2)^2 - (T1[,1]-n2)^3)), tau2)/5 - 2
    
  } else if(model == "step") {
    
    T2 <- rnorm(n.trials, ifelse(T1[,1] < quantile(T1[,1], .25), beta[1],
                                 ifelse(T1[,1] < quantile(T1[,1], .75), beta[2] * 2, beta[1])), tau2)
    
    
  } else if(model == "cubic") {
    
    T2 <- rnorm(n.trials,  beta[1] + beta[2] * T1[,1] + 4.25 * (T1[,1] - quantile(T1[,1], .75))^2 - 2 * (T1[,1] - median(T1[,1]))^3, tau2)
    
  }
  
  else  {
    
    T2 <- rnorm(n.trials, beta[1] + beta[2]*T1[,1] + beta[3]*T1[,2], tau2)
    
  }
  
  ## rescale by factor to make sd(T2) ~= 1
  
  T2 <- T2/sd(T2)
  
  return(data.frame(T1.a = T1[,1], T1.b = T1[,2], T2 = T2))
  
}