###########################################################################
##### This file Runs the Actual VIP data set ##############################
###########################################################################

rm(list =ls())
setwd("~/Dropbox/Code/DPmixturemodel/DPplusAFTVIP")
load("/home/bit/ashar/ExpressionSets/ONE_VIEW/VIPdataset/FINAL_LAST/DataVIPOWNsignature.RData")
X <- relev$Y.train
time.old <- log(30*(relev$time.train)+1)
censoring.old <- relev$censoring.train
Y.info <- relev$Y.info
Np <- length(table(Y.info$Case))


#################################################################
####### Actual Data Matrix #######################################
##################################################################
D <-  ncol(X)


Y.list <- matrix(0,nrow = Np, ncol = D)
for ( j in 1:Np ){
  ind <- which(levels(as.factor(Y.info$Case))[j] == Y.info$Case)
  if (length(ind) ==1){
    Y.list[j,1:D] <- as.vector(X[ind,1:D])
  } else {
    Y.list[j,1:D] <- as.vector(apply(X[ind,1:D],2,mean))
  }
}

time.old <- c(0)
for ( j in 1:Np ){
  ind <- which(levels(as.factor(Y.info$Case))[j] == Y.info$Case)
  time.old[j] <- mean(relev$time.train[ind])
}
time <- log(30*(time.old)+1)
censoring <- rep(1,Np)
That <- time


############################################### MAKING Y from the clusters data #####################3
Y <- Y.list
N <- nrow(Y)
c.true <- ((time > mean(time))+1)

### A little Visualization of the Y Data ##############
pc <- prcomp(Y)
pc.pred <- predict(pc,newdata = Y)
plot(pc.pred[,1], pc.pred[,2], pch = 19,col = c.true)


############################# PARAMETERS for GIBB's SAMPLING ######################################
iter = 100
iter.burnin = 100
iter.thin  =5

################################# GIBBS SAMPLING  ###################################################
K = 10*as.integer(N)

## HYPER PRIORS
## Hyper parameters of the DP
shape.alpha <- 2
rate.alpha <- 1

k =2
F =2


source('rchinese.R')
alpha  = rgamma(1, shape = shape.alpha, rate = rate.alpha )
c <- c.true

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



## Initialization part for the parmaters with FlexMix
source('VIPflexmixBLASSO.R')
km <- flexmixBlasso(Y,That, F,K, r, si, N, D, sig2.dat, c, beta0, betahat, sigma2, lambda2, tau2,c.true)
c <- km$c
sigma2 <- km$sigma2
betahat<- km$betahat
beta0<- km$beta0
lambda2 <- km$lambda2
tau2 <- km$tau2
  




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

Time <- cbind(time,censoring)

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

count <- count -1


########## ANLAYSING THE OUTPUT #######################################################
####### Testing the LogRank Statistic #################################################

######## Compute Summary Statistics about Genes #######################################

##### Class Assignments ########################
c.matrix <- matrix(NA, nrow = N, ncol = count)
for ( i in 1:count){
  temp <- as.factor(c.list[[i]])
  h <- length(levels(temp))
  levels(temp) <- c(1:h)
  c.matrix[,i] <- temp
}
c.final <-  ((apply(c.matrix,1,mean) > 1.01)+1)


############ Time Covariate Slopes FOR Relevant Clusters ############
list.betahat <- list(0)

for ( i in 1:count){
  active <- as.numeric(levels(as.factor(c.list[[i]])))
  list.betahat[[i]] <- (betahat.list[[i]][active,] != 0) +0
}

Q <- 2
matrix.betahat <- array(data = NA, dim =c(Q,count,D))

for ( z in 1:Q){
  for ( x  in 1:count){
    matrix.betahat[z,x,1:D] <- list.betahat[[x]][z,1:D]
  }
}


final.betahat <- matrix(NA, nrow =Q, ncol = D)
for ( z in 1:Q){
  sum <- c(0) 
  for ( x in 1:count){
    sum <- sum + matrix.betahat[z,x,1:D]
  }
 sum <- sum/count
 final.betahat[z,1:D] <- sum 
}

final.betahat[2,1:D] <- rep(0,D)
colnames(final.betahat) <- colnames(X)
### Probability of betahat of genes FOR ONE SIMULATION
##colnames(final.betahat) =  c(rep("relevant",rel.D),rep("irrelevant",irrel.D))
rownames(final.betahat) <- c("Good Prognosis","Bad prognosis")


surv.days <- Surv(exp(time),censoring)
surv.days.fit <- survfit(surv.days ~ c.final)
logrank <- survdiff(surv.days ~ c.final)


############ Some GGplots ###################################
p1 <- ggsurv(surv.days.fit, main = " m-DPMM \n Kaplan Meier for PFS (45) Patients \n 60 Gene DPMM signature \n Mean Good Prognosis (30) 385d \n Mean Poor prognosis(15) 76d", surv.col = c("Red","Black"), cens.col ="Blue")  
pdf('time-DPMM.pdf')
p1
hmcols<-colorRampPalette(c("white","black"))(128)
heatmap.2(t(final.betahat) ,dendrogram="none", margins=c(6,10),col = hmcols, main = "Feature Selection \n Peripheral Cells ", cexCol = 0.85, cexRow = 0.7, Rowv = FALSE, Colv = FALSE) 
dev.off()


