library(randomForest)
library(reshape2)
library(plyr)

## Make a model to determine if Reference status can be classified directly by fish abundance

#Read in data
species <- read.csv("Bio_all_species_level.csv",
                    stringsAsFactors=FALSE)
habitat <- read.csv("habitat_fish_trim.csv",
                    stringsAsFactors=FALSE)

#Format
fishdat <- dcast(species[species$Assemblage =="Fish", ],
                 SampleID ~ Species,
                 fill=0,
                 fun.aggregate=mean, na.rm=TRUE,
                 value.var="Abundance")
fishdat$Status <-  as.factor(
  habitat$ref_1980[match(fishdat$SampleID, habitat$SampleID)])


#Bootstrap data to make Ref/Non-Ref balanced
fishdat <- ddply(fishdat, .(Status), function(x){
  x[sample(1:nrow(x), 250, replace=TRUE), ]
})

#Make a classification model
fishmod <- randomForest(data=na.omit(fishdat[, -1]), 
                        Status ~ ., 
                        ntree=5000,
                        importance=TRUE)

#Extract the most important species from the model
impfish <- importance(fishmod)
goodfish <- names(which(impfish[order(impfish[, 3], decreasing=TRUE), 3] > 60))
