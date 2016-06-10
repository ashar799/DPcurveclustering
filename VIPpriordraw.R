priordraw = function(r, si,N,D, sig2.dat) {

  
  tau2t = matrix(data = NA, nrow = 1, ncol = D)
  betahatt = matrix(data = NA, nrow = 1, ncol = D)
  sigma2t <- 0
  beta0t <- 0
  lambda2t <- 0
  
  

  
  
  lambda2t <- rgamma(1,shape = r, rate = si)
  for ( i in 1:D)  {
  tau2t[1, i] <- rgamma(1, shape = 1, rate = lambda2t)
  } 
  ## Approximating the Jeffery's prior by a Gamma distribution 
 
  sigma2t <- mean(rinvgamma(100, shape = 1, scale = 1))
  
  
  beta0t <- mean(rnorm(100, 0, sd = sig2.dat))
  ## 
  scaleof <- 0
  
  scaleof <- sqrt(abs(sigma2t/lambda2t))
  
  for ( i in 1 :D) {
    
    betahatt[1, i] <- urlaplace(1, location = 0, scale = scaleof) 
  }
  
  list('beta0' = beta0t,'betahat' = betahatt, 'lambda2'= lambda2t, 'tau2'= tau2t, 'sigma2' = sigma2t)
  

}