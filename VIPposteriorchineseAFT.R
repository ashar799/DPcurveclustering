posteriorchineseAFT = function(c,Y,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K,D, r, si, Time,N, sig2.dat) {
  
  
  Ytemp <- matrix(NA, nrow = N, ncol = D)
  ctemp <- c
  
  ## This can't be parallelized !!!!!
  for(l in 1:N)  {
    
    temp <- ctemp[l]
    cminus <- ctemp
    cminus[l] <- NA
    ## The table function helps converting the data point specific indicator variables to class specific indicator variables
    g <- table(factor(cminus, levels = 1:K))
    active <- which(g!=0)
    
    
    
    kminus <- length(active)
    ## Two Auxilary Variables
    ## The name of the auxilary variables are taken to be one and two more than the maximum value in the already active cluster set
    active <- append(active, max(active)+1)
    active <- append(active, max(active)+1)
    
    
    
    ## If the observation was singelton (i.e no other point was associated with it then we assign to kminus +1 parameter)
    if(length(which(cminus==temp))==0 || length(which(cminus==temp))==1 )  
    {
      ## The kminus+1 parameter gets the value of the temporary variable
      ctemp[l] <- active[kminus+1]
      beta0[active[kminus+1]] <- beta0[temp]
      betahat[active[kminus+1], 1:D] <- betahat[temp, 1:D]
      sigma2[active[kminus+1]] <- sigma2[temp]
      lambda2[active[kminus+1]] <- lambda2[temp]
      tau2[active[kminus+1], 1:D] <- tau2[temp, 1:D]
      
      ## Also the second auxilary variable should be drawn from the prior distribution  
      
      source('VIPpriordraw.R')
      priorone <- NA
      priorone <- priordraw(r, si,N,D, sig2.dat)
      beta0[active[kminus+2]] <- priorone$beta0 
      sigma2[active[kminus+2]] <- priorone$sigma2
      betahat[active[kminus+2],1:D] <- priorone$betahat 
      lambda2[active[kminus+2]] <- priorone$lambda2 
      tau2[active[kminus+2], 1:D] <- priorone$tau2
      
      ## As we have to deal with centred matrices and if this point is alone in its cluster then
      for ( k in 1:D){
        Ytemp[l,k] <- 0
      }
      
    }else{
      ## We have to deal with centred matrices
      clust <- which(ctemp == temp)
      tempmatrix <- Y[clust,1:D]
      sd.tempmatrix <- apply(tempmatrix, 2, function(x) sd(x))
      mean.tempmatrix <- apply(tempmatrix, 2, mean)
      
      for ( k in 1:D){
       if (sd.tempmatrix[k] == 0){
         sd.tempmatrix[k] = 1
       }
      }
      
      for ( k in 1:D){
        Ytemp[l,k] <- (Y[l,k] - mean.tempmatrix[k])/(sd.tempmatrix[k])
      }
      
      source('VIPpriordraw.R')
      priortwo <- NA
      priortwo <- priordraw(r, si,N,D, sig2.dat)
      beta0[active[kminus+1]] <- priortwo$beta0 
      sigma2[active[kminus+1]] <- priortwo$sigma2
      betahat[active[kminus+1],1:D] <- priortwo$betahat 
      lambda2[active[kminus+1]] <- priortwo$lambda2 
      tau2[active[kminus+1], 1:D] <- priortwo$tau2
      
      
      source('VIPpriordraw.R')
      priorthree <- NA
      priorthree <- priordraw(r, si,N,D, sig2.dat)
      beta0[active[kminus+2]] <- priorthree$beta0 
      sigma2[active[kminus+2]] <- priorthree$sigma2
      betahat[active[kminus+2],1:D] <- priorthree$betahat 
      lambda2[active[kminus+2]] <- priorthree$lambda2 
      tau2[active[kminus+2], 1:D] <- priorthree$tau2
    }
    
    
    
    #######################################################
    
    
    posterior <- matrix(NA, nrow = length(active), ncol = 1)
    
    
    
    
    ## Calculating the probabalities for drawing the value of c_i from the active classes
    for (j in 1:kminus) {
      
        posterior[j] <- g[active[j]] /(N-1+alpha)  *  dnorm(x = That[l], mean = beta0[active[j]] + betahat[active[j],1:D] %*% as.vector(t(Ytemp[l,1:D])), sd = sqrt(sigma2[active[j]]) )
    }
    
      posterior[kminus+1] <- (0.5 * alpha) /(N-1+alpha) *  dnorm(x = That[l], mean = beta0[active[kminus+1]] + betahat[active[kminus+1],1:D] %*% as.vector(t(Ytemp[l,1:D])), sd = sqrt(sigma2[active[kminus+1]]) )
      posterior[kminus+2] <- (0.5 * alpha) /(N-1+alpha) *  dnorm(x = That[l], mean = beta0[active[kminus+2]] + betahat[active[kminus+2],1:D] %*% as.vector(t(Ytemp[l,1:D])), sd = sqrt(sigma2[active[kminus+2]]) )
    
    
    ## Calculating the normalization constant for probabilities
    normalization <- sum(posterior) 
    
    if (normalization < 1e-200 || normalization ==Inf){
      ctemp[l] <- sample(active, 1, prob = rep(1,length(active)), replace = TRUE)
    } else {  
      ctemp[l] <- sample(active, 1, prob= posterior, replace = TRUE)
    }
    
  }
  
  c <- ctemp
  ## Delete those observations that are not associcated with no data point
  g <- table(factor(c, levels = 1:K))
  inactive <- which(g==0)
  
  for ( i in 1:length(inactive)){
    beta0[inactive[i]] <- NA 
    sigma2[inactive[i]] <- NA
    betahat[inactive[i],1:D] <- NA 
    lambda2[inactive[i]] <- NA
    tau2[inactive[i], 1:D] <- NA
  }
  
  
  
  list('indicator' = c, 'beta0' = beta0,'betahat2'= betahat, 'sigma2'=sigma2, 'lambda2'= lambda2, 'tau2' = tau2)
  
}
