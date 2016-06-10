###########This file plays around with the VIP data to check for good and bad prognosis ###############
#######################################################################################################

### The File takes 218 samples ###############################################################################
#### This file does Prefiltering also based on those genes which discriminate PvsZ and ARE releated to survival
##### This script CLUSTERS Periphery Cells ON THE PvsZ signature and PFS survival  ###########################
library('affy')
library('xlsx')
library('limma')
library('survival')
rm(list =ls())

## load data ##
load("/home/abidata/Dinis/VIP/Main.Data/ExpressionConsole.normalised.RData")


## The phenoData
tab <- read.xlsx(file  = '/home/abidata/Dinis/VIP/Sample.Info/15-07-14 VIP Sample Collection PvsZ DEG.xlsx', sheetIndex =1)
list.patients <- tab[,2]


### Only those patients who have a NON NA is PFS
pheno  <- pData(eset.ec.norm) 
pheno.ready <- pheno[!is.na(pheno$PFS),]
list.patients.final <- list.patients[(list.patients %in% pheno.ready[,3])]


## Include ONLY those features that have correspoding annotation ~ 27148 out of 70000
## Include ONLY those Samples which HAVE PFS and PvsZ Annotation
exprs <- eset.ec.norm[featureNames(eset.ec.norm) %in% rownames(anno.GENENAME),sampleNames(eset.ec.norm) %in% list.patients.final]
v <- t(exprs(exprs))
labs.v <- (pData(exprs)$PFS > median(as.numeric(pData(exprs)$PFS))) +0 
pc <- prcomp(v)
pc.pred <- predict(pc,newdata = v)
plot(pc.pred[,1], pc.pred[,2], pch = 19,col = c("red","green")[as.factor(labs.v)])
plot(pc.pred[,1], pc.pred[,2], pch = 19,col = c("red","green")[as.factor(pData(exprs)$Sex)])


###### What is the REASON FOR THIS BATCH EFFECT #########################
######### Let's Look at the Features which distinguish P and Z########
########################################################################


setwd("~/Dropbox/Code/DPmixturemodel/DPplusAFTVIP")
load("/home/bit/ashar/ExpressionSets/ONE_VIEW/VIPdataset/DPMM-FactorAnalysis/DataVIPOWNsignature.RData")
X <- relev$Y.train
Y.info <- relev$Y.info
Np <- length(table(Y.info$Case))

labs.x <-  Y.info$Case
pc <- prcomp(X)
pc.pred <- predict(pc,newdata = X)
plot(pc.pred[,1], pc.pred[,2], pch = 19)

