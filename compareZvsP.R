###########################################################################################################
###### Comparing DPMM model from both Z and P cells ########################################################
############################################################################################################

rm(list =ls())
load("/home/bit/ashar/ExpressionSets/ONE_VIEW/VIPdataset/DPMM-FactorAnalysis/DPMM-meansamples.RData")
c.periphery <- c.final
samples.periphery <- rownames(Y)

load("/home/bit/ashar/ExpressionSets/ONE_VIEW/VIPdataset/DPMM-FactorAnalysis/ZVIPDPMM.RData")
c.center <- c.final
samples.center <- rownames(Y)

###### Common samples

samples.periphery.common <- samples.periphery[samples.periphery %in% samples.center]
c.periphery.common <- c.periphery[samples.periphery %in% samples.center]

#### Adjusted Rand Index

