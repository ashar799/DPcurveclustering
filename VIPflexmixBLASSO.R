flexmixBlasso = function(Y,That, F,K, r, si, N, D, sig2.dat, c, beta0, betahat, sigma2, lambda2, tau2,c.true ) {
  
  source('VIPpriordraw.R')
  G <- F
  
  
 ### For the time data points we use flexmix to assign the clusters
    
 ## The cross validation folds for choosing lambda
 data <- data.frame(time, Y)
 fo <- sample(rep(seq(10), length = nrow(data)))
 gr.flx <- flexmix(time ~ ., data = data, k = F, cluster =c.true, model = FLXMRglmnet(foldid = fo, adaptive= FALSE), control = list(iter.max = 500))
 c <- clusters(gr.flx)
 
  
  
  prior.numclust <- table(factor(c, levels = 1:K))
  prior.activeclass<- which(prior.numclust!=0)
  
 
 
 
 
 

 
 
  ### The means are set using the k-means
  for ( i in 1:length(prior.activeclass)){
    
#    lclust <- which(c == prior.activeclass[i])
#     
#     reg.blas <- 0
#     
#     sum <- c(0)
#     
#     coeff <- 0
#     
#     Ytemp <-  matrix(NA, nrow = length(lclust), ncol = D)
#     
#     Ytemp <- scale(Y[lclust,1:D], center = TRUE, scale = TRUE)
#     
#     
#     ### Part where I use the MONOMVN PACKAGE
#     
#     Ttemp <- as.vector(That[lclust])
#     
#     ntemp <- length(lclust)
#     
#     reg.blas <- blasso(Ytemp, Ttemp, T = 300,thin =  50, RJ = TRUE, mprior = 0.0 ,normalize = TRUE, verb = 0)
#     
#     sum <- summary(reg.blas, burnin= 100)
#     
#     ## Selecting those features which are relevant
#     
#     coeff <- unlist(lapply(strsplit(sum$coef[3,], split = ":"), function(x) as.numeric(unlist(x)[2])))
#     
#     beta0[prior.activeclass[i]] <- coeff[1]
#     
#     indexplusone <- D+1
#     
#     ind <- 2:indexplusone
#     
#     betahat[prior.activeclass[i], ] <- coeff[ind]
#     
#     ta <- unlist(lapply(strsplit(sum$tau2i[3,], split = ":"), function(x) as.numeric(unlist(x)[2])))
#     
#     tau2[prior.activeclass[i],] <- ta
#     
#     sigma2[prior.activeclass[i]] <- sum$s2[3]
#     
#     lambda2[prior.activeclass[i]] <- sum$lambda2[3]
    
    source('VIPpriordraw.R')
    
    priorone <- priordraw(r, si,N,D, sig2.dat)
    beta0[prior.activeclass[i]] <- priorone$beta0 
    sigma2[prior.activeclass[i]] <- priorone$sigma2
    betahat[prior.activeclass[i],1:D] <- priorone$betahat 
    lambda2[prior.activeclass[i]] <- priorone$lambda2 
    tau2[prior.activeclass[i], 1:D] <- priorone$tau2
    
    
  }
  
  ## Deleting those values which are no longer relevant
  g <- table(factor(c, levels = 1:K))
  inactive <- which(g==0)
  
  for ( i in 1:length(inactive)){
   beta0[inactive[i]] <- NA 
    sigma2[inactive[i]] <- NA
    betahat[inactive[i],1:D] <- NA 
    lambda2[inactive[i]] <- NA
    tau2[inactive[i], 1:D] <- NA
  }
  
  
  
  list('c' = c, 'beta0'=beta0, 'betahat'= betahat, 'sigma2' =sigma2, 'lambda2' = lambda2, 'tau2'= tau2)  
  
}