## Testing the DPplusAFT model with Gibb's Sampling
### Script to test if the model works with higher dimensions

rm(list = ls())
setwd('/home/bit/ashar/Dropbox/Code/DPmixturemodel/DPplusAFTVIP')


#################################### SIMULATED DATA PROPERTIES ####################################################

############### Patient Data ##################################################
Np = 120

#### Division of per patient sample
sample.patient <- rep(1,Np) + rpois(Np,2)

#### Number of molecular data per patient ####################################
N = sum(sample.patient)

## Number of Patient Clusters
Z =2

## Order of Molecular data, this is used for simulation
order <- sample(N)

## Number of Clusters on Molecular Level
F = 2

## Distribution of the points within three clusters

p.dist = c(sum(sample.patient[1:60]), sum(sample.patient[61:120]))

## Total Number of features D
D = 40

## Total Percentage of irrelevant feature

prob.noise.feature = 0.50

## Total Percentage of censoring

prob.censoring = 0.05


## Overlap between Cluster of molecular Data of the relevant features

prob.overlap = 0.05

## Percentage of Noise/Overlap in Time Data

prob.noise = 0.05

## Actual Number of Components and dimension  which are relevant
rel.D = 20
## Actual Number of Irrelevant Componenets
irrel.D = 20


#### Use Patient Level Pairing information 
X.info <- matrix(NA,nrow = N, ncol = 1)
count = 1
for( i in 1:Np){
  for ( j in 1:sample.patient[i]){
    X.info[count] <- i
    count <- count +1
  }
}


########################################Simulated Data##############################################
A <- MixSim(MaxOmega = prob.overlap ,K = F, p = rel.D, int =c(0,10), lim = 1e08)
data.mu = array(data = NA, dim =c(F,D))
data.S = array(data = NA, dim =c(F,D,D))

for( i in 1:F){
  data.mu[i,1:rel.D] <- A$Mu[i,1:rel.D]
  data.S[i,1:rel.D,1:rel.D] <- A$S[1:rel.D,1:rel.D,i]
}

## The relevant data is genereated first
Y.rel.list <- list(0)
for ( i in 1:F){
  Y.rel.list[[i]] <- mvrnorm(n = p.dist[i], mu = data.mu[i,1:rel.D], Sigma = data.S[i,1:rel.D,1:rel.D])
}


## Irrelevant features
Y.irrel.list <- list(0)
for ( i in 1:F){
  mean <- runif(irrel.D,0,10)
  Y.irrel.list[[i]] <- mvrnorm(n =p.dist[i], mu = mean, Sigma = diag(x =1, nrow = irrel.D, ncol = irrel.D))
}


### Combining the data with relevant and irrelevant columns
data.old <- list(0) 
for (i in 1:F){
  data.old[[i]] <-  cbind(Y.rel.list[[i]], Y.irrel.list[[i]]) 
}

############################################### MAKING X from the clusters data #####################
X.old <- c(0)
for (i in 1:F){
  X.old<- rbind(X.old, data.old[[i]])
} 
X.old<- X.old[-1,]

X.old <- as.matrix(X.old)
c.true.x <- c(c(rep(1,p.dist[1])),c(rep(2,p.dist[2])))

#### Shuffle the points such that same patients fall into different clusters ########################
shuffle <- sample(N)

X <- X.old[shuffle,]
c.true.X <- c.true.x[shuffle]



### A little Visualization of the Y Data ##############
pc <- prcomp(X)
pc.pred <- predict(pc,newdata = X)
plot(pc.pred[,1], pc.pred[,2], pch = 19,col = c.true.X)





##################################################Making Y by using means of the data set ###############################
Y.list <- matrix(0,nrow = Np, ncol = D)
for ( j in 1:Np ){
  ind <- which(X.info==j)
  if (length(ind) ==1){
    Y.list[j,1:D] <- as.vector(X[ind,1:D])
  } else {
    Y.list[j,1:D] <- as.vector(apply(X[ind,1:D],2,mean))
  }
}

## True Labels for the points
c.true <- c(c(rep(1,60)),c(rep(2,60)))




### Make the Features Uncorrelated
data.old <-list(0)
data.new <- list(0)


for ( i in 1:F){
ind <- which(c.true ==i)
data.old[[i]] <- Y.list[ind,1:D]
}

######################################################################################### Making the irrelevant features independent from the dependent features #############

for (f in 1:F){
  
  V <- data.old[[f]]
  
  rel.V <- as.matrix(V[,1:rel.D])
  
  obj.qr <- qr(V)
  
  rk <- obj.qr$rank
  
  alpha <- qr.Q(obj.qr)[,1:rel.D]
  
  gamma <- qr.Q(obj.qr)[,(1+rel.D):rk]
  
  matT <- matrix(runif(n = rel.D*(rk -rel.D), min = -0.0005, max= 0.0005), nrow = rel.D, ncol = (rk -rel.D))
  
  matP <- t(matT) %*% matT
  
  max.eig <- eigen(matP)$values[1]
  
  max.corr <- sqrt(max.eig)/sqrt(1 + max.eig)
  
  linear.space <- gamma + alpha %*% matT
  
  irrel.V <- matrix(NA, nrow = nrow(V), ncol = irrel.D)
  
  for ( i in 1: irrel.D){
    
    matTemp <- matrix(runif(n = (rk -rel.D), min = -1.5, max= 1.5),  nrow = (rk-rel.D), ncol =1)
    irrel.V[,i] <- as.vector(linear.space %*% matTemp)
    
  }
  
  ## Checking if the covariance is indeed small
  
  cov.mat <- cov(rel.V,irrel.V)
  
  boxplot(cov.mat)
  
  ## Building the full data matrix
  
  V.full <- cbind(rel.V, irrel.V)
  
  
  levelplot(cov(V.full))
  
  data.new[[f]] <- V.full
  
}

#######################################################



############################################### MAKING Y from the clusters data #####################3
Y <- c(0)
for (i in 1:F){
  Y <- rbind(Y, data.new[[i]])
} 
Y <- Y[-1,]
### A little Visualization of the Y Data ##############
pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
plot(pc.pred[,1], pc.pred[,2], pch = 19,col = c.true)



## The Co-efficients have to be obtained from uniform distribution between [1,10]
beta.list <- list(0)
for ( i in 1:F){
  beta.list[[i]] <- runif(rel.D, min = -5, max = 5)
}

Y.rel.sc.list <- list(0)
for ( i in 1:F){
  ind <- which(c.true ==i)
  Y.rel.sc.list[[i]] <- scale(Y[ind,1:rel.D],center = TRUE, scale = FALSE)
}

## The pure time is generated
time.pur.list <- list(0)
for ( i in 1:F){
  time.pur.list[[i]] <- t(beta.list[[i]]) %*% t(Y.rel.sc.list[[i]])
}


time.noise.list <- list(0)
for ( i in 1:F){
  time.noise.list[[i]] <- rnorm(60, mean = 50*i^2, sd = i)
}

time.list <- list(0)
for ( i in 1:F){
  time.list[[i]] <- time.pur.list[[i]] + time.noise.list[[i]]
}


#######################################MAKING TIME from cluster data ########################################################
time <- c(0)
for (i in 1:F){
  time <- cbind(time, time.list[[i]])
} 
time <- time[,-1]
time <- as.vector(time)


####################################### Adding CENSORING INFORMATION  ################################################################
## Adding the censoring information
censoring <- rep(1,Np)



############################# PARAMETERS for GIBB's SAMPLING ######################################
iter = 100
iter.burnin = 100
iter.thin  =5

################################# GIBBS SAMPLING  ###################################################

Time <- cbind(time, censoring) 
D = NCOL(Y)
N = NROW(Y)
K = 10*as.integer(N)

## HYPER PRIORS
## Hyper parameters of the DP
shape.alpha <- 2
rate.alpha <- 1



source('rchinese.R')
alpha  = rgamma(1, shape = shape.alpha, rate = rate.alpha )
c <-  rchinese(N,alpha)
f <- table(factor(c, levels = 1:max(c)))


#Sparsity controlling hyperparameter of the BAYESIAN LASSO MODEL
r =1
si = 1.78

lambda2 <- numeric(K)
tau2 = matrix(data = NA, nrow = K, ncol = D)
betahat = matrix(data = NA, nrow = K, ncol = D)
sigma2 <- rep(NA, K)
beta0 <- rep(NA, K)
That <-  time

## Fitting a linear model to the whole model
Ysc <- scale(Y[1:N,1:D], center = TRUE, scale =TRUE)
lm.data <- lm(time ~ Ysc)
sig2.dat <-  var(lm.data$residuals)


## Set Some Initial Values for the Cluster Parameters


disclass <- table(factor(c, levels = 1:K))
activeclass <- which(disclass!=0)


for ( j in 1:length(activeclass)){
  source('VIPpriordraw.R')
  priorone <- priordraw(r, si, N, D, sig2.dat)  
  beta0[activeclass[j]] <- priorone$beta0 
  sigma2[activeclass[j]] <- priorone$sigma2
  betahat[activeclass[j],1:D] <- priorone$betahat 
  lambda2[activeclass[j]] <- priorone$lambda2 
  tau2[activeclass[j], 1:D] <- priorone$tau2
}


## Initialization part for the parmaters with FlexMix
source('VIPflexmixBLASSO.R')
t.c <- list(0)
t.sigma2 <- list(0)
t.betahat <- list(0)
t.beta0 <- list(0)
t.lambda2 <- list(0)
t.tau2 <- list(0)

rand <- c(0)
for ( q in 1:50){
  km <- flexmixBlasso(Y,That, F,K, r, si, N, D, sig2.dat, c, beta0, betahat, sigma2, lambda2, tau2)
  
  t.c[[q]] <- km$c
  rand[q] <- adjustedRandIndex(t.c[[q]], c.true)
  t.sigma2[[q]] <- km$sigma2
  t.betahat[[q]] <- km$betahat
  t.beta0[[q]] <- km$beta0
  t.lambda2[[q]] <- km$lambda2
  t.tau2[[q]] <- km$tau2
  
}

c <- t.c[[which.max(rand)]]
sigma2 <- t.sigma2[[which.max(rand)]] 
betahat <- t.betahat[[which.max(rand)]] 
beta0 <- t.beta0[[which.max(rand)]] 
lambda2 <- t.lambda2[[which.max(rand)]] 
tau2 <- t.tau2[[which.max(rand)]] 


## Adjusted Rand INDEX measure
randindex <- adjustedRandIndex(c.true,as.factor(c))




## BURNIN AND GIBBs sampling 
source('VIPpriordraw.R')
source('VIPposteriorchineseAFT.R')
source('VIPposteriortimeparameters.R')


cognate <- NA
param <- NA
paramtime <- NA
loglike<- rep(0, iter)  
timeparam <- NA
time.predicted <- c(0)
cindex <- c(0)
likli <- c(0)

source('VIPlikelihood.R')
print(loglikelihood(c,Y,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, D, r, si, Time,N, sig2.dat) )
      
      
#################### BURNIN PHASE ###################################################
print("BURNIN...PHASE")
for (o in 1:iter.burnin) {
        
        
  ################## PARAMETERS OF THE DP Mixture Model ######################################################
  ## Updating the parameters based on the observations 
  
  paramtime <- posteriortimeparameters(c, That, lambda2,tau2,sigma2,beta0, betahat, Y, K, r, si, Time, N, sig2.data)
  beta0 <- paramtime$beta0
  betahat <- paramtime$betahat
  sigma2 <- paramtime$sigma2
  lambda2 <- paramtime$lambda2
  tau2 <- paramtime$tau2
  
  
  ################# INDICATOR VARIABLE ##################################################################
  ## Updating the indicator variables and the parameters
  source('VIPposteriorchineseAFT.R')
  cognate <- posteriorchineseAFT(c,Y,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K,D, r, si, Time,N, sig2.dat)
  c <- cognate$indicator
  beta0 <- cognate$beta0
  betahat <- cognate$betahat
  sigma2 <- cognate$sigma2
  lambda2 <- cognate$lambda2
  tau2 <- cognate$tau2
  
  ########################### The Concentration Parameter #################################################################
  source('posterioralpha.R') 
  # Updating the concentration parameter
  alpha <- posterioralpha(c, N, alpha, shape.alpha, rate.alpha)
  
  ########################### Print RandIndex and Log-likelihood ####################################################
  print(adjustedRandIndex(c.true,c))
  print(loglikelihood(c,Y,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K, D, r, si, Time,N, sig2.dat))     
  print(o/iter.burnin)
      
  } 
      

############## GIBBS SAMPLING WITH THINNING ######################################################
beta0.list <- list(0)
betahat.list <- list(0) 
sigma2.list <- list(0)
lambda2.list <- list(0)
tau2.list <- list(0)
c.list <- list(0)


print("GIBB'S SAMPLING")
count = 1
for (o in 1:iter) {
  
  
  ################## PARAMETERS OF THE DP Mixture Model ######################################################
  ## Updating the parameters based on the observations 
  
  paramtime <- posteriortimeparameters(c, That, lambda2,tau2,sigma2,beta0, betahat, Y, K, r, si, Time, N, sig2.data)
  beta0 <- paramtime$beta0
  betahat <- paramtime$betahat
  sigma2 <- paramtime$sigma2
  lambda2 <- paramtime$lambda2
  tau2 <- paramtime$tau2
  
  
  ################# INDICATOR VARIABLE ##################################################################
  ## Updating the indicator variables and the parameters
  source('VIPposteriorchineseAFT.R')
  cognate <- posteriorchineseAFT(c,Y,alpha,That, beta0, betahat, sigma2, lambda2, tau2, K,D, r, si, Time,N, sig2.dat)
  c <- cognate$indicator
  beta0 <- cognate$beta0
  betahat <- cognate$betahat
  sigma2 <- cognate$sigma2
  lambda2 <- cognate$lambda2
  tau2 <- cognate$tau2
  
  ########################### The Concentration Parameter #################################################################
  source('posterioralpha.R') 
  # Updating the concentration parameter
  alpha <- posterioralpha(c, N, alpha, shape.alpha, rate.alpha)
  
  
  if(o%% iter.thin == 0 ){
    beta0.list[[count]] <- beta0
    betahat.list[[count]] <- betahat  
    sigma2.list[[count]] <- sigma2
    lambda2.list[[count]] <- lambda2
    tau2.list[[count]] <- tau2
    c.list[[count]] <- c
    count <- count +1
  }
  
  
  
  
  print(o/iter) 

} 

########## ANLAYSING THE OUTPUT #######################################################

