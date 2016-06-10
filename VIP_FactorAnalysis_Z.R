##### This file takes the VIP signature and then calculates Average Expression Vector for each patient #############

rm(list =ls())
setwd("~/Dropbox/Code/DPmixturemodel/DPplusAFTVIP")
load("/home/bit/ashar/ExpressionSets/ONE_VIEW/VIPdataset/DPMM-FactorAnalysis/DataVIPOWNsignatureZ.RData")
X <- relev$Y.train
Y.info <- relev$Y.info
Np <- length(table(Y.info$Case))


table.cases <-  table(Y.info$Case)
Y.joint <- matrix(NA, nrow = length(table.cases), ncol = ncol(X))
rownames(Y.joint) <- names(table.cases)
colnames(Y.joint) <- colnames(X)
survival.joint <- c(0)
censoring.joint <- rep(1,)
for ( i in 1:length(table.cases)){
  case <- names(table.cases)[i]
  
  
  if(table.cases[case] ==1){
    Y.joint[i,] <- X[Y.info$Case == case,]
    survival.joint[i] <- log(30*as.numeric(Y.info$PFS[Y.info$Case ==case])+1)
  } else {
    temp.matrix <- apply(t(X[Y.info$Case == case,]),1,mean)
    Y.joint[i,] <-t(temp.matrix)
    survival.joint[i] <- log(30*mean(as.numeric(Y.info$PFS[Y.info$Case ==case]))+1)
  }
  
}

pc <- prcomp(Y.joint)
pc.pred <- predict(pc,newdata = Y.joint)
plot(pc.pred[,1], pc.pred[,2], pch = 19)

relev <- list('Y.train' =Y.joint, 'time.train' = survival.joint )

save(relev, file = 'DataVIPOWNsignaturePLUSMeanSamplesZ.RData')
