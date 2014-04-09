library(randomForest)
library(reshape2)
library(plyr)
library(ggmap)
source("abundanceBC.r")

### Build group membership predictive model ####

#Read in data
species <- read.csv("Bio_all_species_level.csv",
                    stringsAsFactors=FALSE)
habitat <- read.csv("habitat_fish_trim.csv",
                    stringsAsFactors=FALSE)

habitat <- ddply(habitat, .(ref_1980), function(x){
  set.seed(12345)
  x$Cal <- runif(nrow(x)) > 0.3
  x
})
#  Filter out species by strong predictors
species <- species[species$Species %in% goodfish, ]
species$Ref <- habitat$ref_1980[match(species$SampleID, habitat$SampleID)]
species$Cal <- habitat$Cal[match(species$SampleID, habitat$SampleID)]

#Get sample by taxon matrix
fishCal <- species[which(species$Assemblage =="Fish" &
                           species$Ref == "Yes" &
                           species$Cal), ]
speciesM <- acast(fishCal,
                  SampleID ~ Species,
                  fill=0,
                  fun.aggregate=mean, na.rm=TRUE,
                  value.var="Abundance")

#Cluster analysis set apriori to n groups
groups <- kmeans(speciesM, 10)$cluster
habitat$SiteCode <- sapply(strsplit(habitat$SampleID, "_"), "[", 1)
habitat$group <- as.factor(groups[match(habitat$SampleID, names(groups))])

#create model
rfdat <- habitat[habitat$ref_1980 == "Yes" & habitat$Cal,
                 !grepl("ref", names(habitat)) &
                   !grepl("calval", names(habitat)) &
                   !names(habitat) %in% c("SampleID",
                                          "Year",
                                          "SiteCode",
                                          "Cal")]
rfdat$Island_Mainland <- as.factor(rfdat$Island_Mainland)

rfmod <- randomForest(group ~ ., 
                      data=rfdat,
                      ntree=5000, 
                      proximity=T)

## Other stuff
preds <- names(rfdat)[names(rfdat) != "group"]
calibration <- apply(speciesM, 2, function(x)tapply(x, rfdat$group, mean))

val_abund <- species[which(species$Assemblage =="Fish" #&
                           #                              ((!species$Cal & species$Ref == "Yes")|
                           #                                 species$Ref == "No"))
                           ), ]
val_preds <- habitat[
#   (!habitat$Cal & habitat$ref_1980 == "Yes")|
#                        habitat$ref_1980 == "No"
  ,
                     !grepl("ref", names(habitat)) &
                       !grepl("calval", names(habitat)) &
                       !names(habitat) %in% c("Year",
                                              "SiteCode",
                                              "Cal", "group")]
val_preds$Island_Mainland <- as.factor(val_preds$Island_Mainland)


### test ####
test <- abundanceBC(val_abund, val_preds, 
                    mod = rfmod,
                    occurence = calibration,
                    preds = preds)

validated <- join(test, habitat)

t.test(validated$BCdis[validated$ref_1980 == "Yes" & !validated$Cal],
       validated$BCdis[validated$ref_1980 == "No"])

base <- get_map("Catalina Island", zoom=8)
ggmap(base) + 
  geom_point(data=validated, aes(Longitude, Latitude, colour=BCdis), size=4) +
  scale_color_continuous(high="red", low="green")

ggplot(validated, aes(ref_1980, BCdis, fill=Cal)) + geom_boxplot()
