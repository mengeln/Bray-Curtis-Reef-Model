library(randomForest)
library(reshape2)
library(plyr)

BC <- function(x, y){
  total <- sum(x + y)
  numerator <- sum(abs(x - y))
  numerator/total
}

abundanceBC <- function(dat, predictors, mod, 
                        occurence, preds) {
  
  #Get predictors
  group_probs <- predict(mod, newdata = predictors[,preds], type = 'prob')
  predicted_abund <- group_probs %*% occurence
  row.names(predicted_abund) <- predictors$SampleID
  
  #Merge with observations
  predflat <- melt(predicted_abund)
  names(predflat) <- c("SampleID", "Species", "predicted_abundance")
  predflat$Species <- as.character(predflat$Species)
  comb <- join(dat, predflat, by=c("SampleID", "Species"))
  
  #Calculate BC
  ddply(comb, .(SampleID), summarise,
        BCdis = BC(Abundance, predicted_abundance))
}




