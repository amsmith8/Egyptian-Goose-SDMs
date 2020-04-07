## Chukar modeling - site level factors 
## Author: Austin M. Smith

###################################################
### code chunk number 1: set up environment 
###################################################

rm(list=ls(all=TRUE))
setwd("/Users/austinsmith/Desktop/A.chukar_dismo")

######## ########  Uploaded Packages ######## ########

# data processing
library("dismo")
library("maptools")
library("raster")
library("rgdal")
library("rgeos")
#library(sf)
library(sp)

# Monitor timing 
library("tictoc")

###-----------------------------------------------------------------------------------------------------------------------

###################################################
### code chunk number 4: Collect pts and data
################################################### 

# ***** The following is needed if running script as a stand-alone run.
super.stack <- readRDS("super.stack.rds")
Model_CRS <- crs(super.stack)

#Chukar range polygon file 
Ac.poly <- readShapePoly("/Users/austinsmith/Desktop/A.chukar_dismo/Alectoris_chukar/Alectoris_chukar.shp") # via Bird Life
crs(Ac.poly) <- Model_CRS
plot(Ac.poly)

### Seperate the polygon into native and naturalized regions
Native <- subset(Ac.poly, Ac.poly$OBJECTID == 36 )  # historical  native range for A. chukar. Similar to Christensen 1970
Naturalized <- subset(Ac.poly, Ac.poly$OBJECTID != 36 ) # recognized regions of naturalized populations

library(rnaturalearth)
world.land <- ne_countries(returnclass = c("sp", "sf"))
crs(world.land) <- Model_CRS
par(1,1,1,1)
plot(world.land, main = "Global Chukar Distribution")

# Check on map
plot(Native,add = T , col = "blue")
plot(Naturalized, add = T, col = "red")

# Remove Antarctica
ws  <-  crop(world.land, c(-180, 180, -62, 83.57027))
crs(ws) <- Model_CRS
plot(ws)

# Remove Greenland
greenland0 <- readShapePoly("/Users/austinsmith/Downloads/GRL_adm/GRL_adm0.shp" )
crs(greenland0) <- Model_CRS
ws <- ws - greenland0
#plot(ws, col="green")

#### Calculations to determine how many points
sum(area(Ac.poly)) / sum(area(ws)) # determine ratio of background.pts : Native.pts

###################
set.seed(5421)
chukar.present <- spsample(Native, 1000, type = 'random')
chukar.present.coords <- data.frame(chukar.present@coords) # data frame of lon lat coords 
names(chukar.present.coords) <- c("lon", "lat")

k <- 5
k.means <- kmeans(chukar.present.coords, k)
occ.grp <- k.means$cluster

chukar.present.covariates <- extract(super.stack, chukar.present.coords) # extract values from raster stack 
chukar.present <- cbind(chukar.present.coords, chukar.present.covariates) # combine coords and covars into one df 
plot(Native, col ="gray", main = "Spatial bins of  native Chukar occurrences " )
points(chukar.present.coords, pch=21, bg=occ.grp)

ws.background.pts <- spsample(ws, 11000, type = 'random')
#points(ws.background.pts)
ws.background.coords <- data.frame(ws.background.pts@coords)
names(ws.background.coords) <- c("lon", "lat")
ws.background.covariates <- extract(super.stack, ws.background.coords)
ws.background <- cbind(ws.background.coords,ws.background.covariates)

#combine data 
pts.id <-  c(rep(1, nrow(chukar.present)), 
             rep(0, nrow(ws.background))) 

sdm.data <- data.frame(rbind(chukar.present, ws.background )) 
sdm.data <- data.frame(cbind(sdm.data, pts.id))
sdm.data <- na.omit(sdm.data) # remove N/A's


###-----------------------------------------------------------------------------------------------------------------------


###################################################
### code chunk number 5: Run algorithms 
################################################### 

######## ########  Uploaded Packages ######## ########


# Problem with files on computer.  this uploads rJava for the session on current computer 
dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_241.jdk/Contents/Home/jre/lib/server/libjvm.dylib')

#Sys.setenv(NOAWT=TRUE)
# Maxent -dependency 
library("rJava")

# Monitor timing 

library("tictoc")

# Modeling  
library("caret")
library("caretEnsemble")
library("dismo")
library("gbm")
library("kernlab")
library("neuralnet")
library("nnet")
library("randomForest")

###_________________________________________________


# Create data.frame to hold AUC scores 
auc.scores <- data.frame(matrix(NA,nrow=5,ncol=6))
auc.scores[,1]<-c("SVM", "GBM", "RF", "ANN", "MaxEnt")
colnames(auc.scores)<-c("model","Fold.1","Fold.2","Fold.3","Fold.4","Fold.5")

###### 

sensSpec.score <- data.frame(matrix(NA,nrow=5,ncol=6))
sensSpec.score[,1]<-c("SVM", "GBM", "RF", "ANN", "MaxEnt")
colnames(sensSpec.score)<-c("model","Fold.1","Fold.2","Fold.3","Fold.4","Fold.5")

########

# Model stacks 
model.stack.svm <- stack()
model.stack.gbm <- stack()
model.stack.rf <- stack()
model.stack.ann <- stack()
model.stack.maxent <- stack()


###-----------------------------------------------------------------------------------------------------------------------


##### Run Models 
######## ######### #########

for (i in 1:k) {
  # set up partition for iterations
  train <- sdm.data[occ.grp != i,] # all folds excluding the i-th 
  test <- sdm.data[occ.grp == i,] # only the i-th
  

  # SVM
  ######## ######## #######
  tic("svm")
  svm <- ksvm(pts.id ~ ., data=train[3:28] )
  #evaluate
  eval.svm <-  evaluate(p=test[test$pts.id == 1,], test[test$pts.id == 0,], svm)
  auc.scores[1,1+i] <- eval.svm@auc
  sensSpec.score[1,1+i] <- threshold(eval.svm, stat='spec_sens')
  pred.svm <- predict(super.stack, svm)
  model.stack.svm <- stack(model.stack.svm, pred.svm) # add rasters
  toc()
  
  # GBM
  ######## ######## #######
  tic("gbm")
  #gbm <- gbm(pts.id ~ ., data=train , n.trees = 100)
  gbm <- gbm(pts.id ~ ., data=train[3:28] , n.trees = 100)
  #evaluate
  eval.gbm <-  evaluate(p=test[test$pts.id == 1,], test[test$pts.id == 0,], gbm,n.trees = 100)
  auc.scores[2,1+i] <- eval.gbm@auc
  sensSpec.score[2,1+i] <- threshold(eval.gbm, stat='spec_sens')
  pred.gbm <- predict(super.stack, gbm, n.trees = 100)
  model.stack.gbm <- stack(model.stack.gbm, pred.gbm) # add rasters
  toc()
  
  # RF
  ######## ######## #######  
  tic("rf")
  rf <- randomForest(pts.id ~ ., data = train[3:28] )
  #rf <- randomForest(pts.id ~ ., data = train )
  #evaluate
  eval.rf <-  evaluate(p=test[test$pts.id == 1,], test[test$pts.id == 0,], rf)
  auc.scores[3,1+i] <- eval.rf@auc
  sensSpec.score[3,1+i] <- threshold(eval.rf, stat='spec_sens')
  pred.rf <- predict(super.stack, rf)
  model.stack.rf <- stack(model.stack.rf, pred.rf) # add rasters
  toc()
  
  # ANN
  ######## ######## #######  
  tic("ann")
  ann <- neuralnet(pts.id ~ ., data=train[3:28] , stepmax = 1e+07)
  #evaluate
  eval.ann <-  evaluate(p=test[test$pts.id == 1,], test[test$pts.id == 0,], ann)
  auc.scores[4,1+i] <- eval.ann@auc
  sensSpec.score[4,1+i] <- threshold(eval.ann, stat='spec_sens')
  pred.ann <- predict(super.stack, ann)
  model.stack.ann <- stack(model.stack.ann, pred.ann) # add rasters
  toc()
  
  # MaxEnt
  ######## ######## #######  
  tic("maxent")
  train.occ <- subset(train, pts.id == 1)
  train.bg <- subset(train, pts.id == 0)
  maxent <- maxent(x = super.stack, p = train.occ[1:2], a = train.bg[1:2] ) # background,  regression
  test.occ <- subset(test, pts.id == 1)
  test.bg <- subset(test, pts.id == 0)
  #evaluate
  eval.maxent <-  evaluate(p=test.occ[1:2], a <-test.bg[1:2], maxent, x = super.stack)
  auc.scores[5,1+i] <- eval.maxent@auc
  sensSpec.score[5,1+i] <- threshold(eval.maxent, stat='spec_sens')
  pred.maxent<- predict(super.stack, maxent)
  model.stack.maxent <- stack(model.stack.ann, pred.maxent ) # add rasters
  toc()
}

###-----------------------------------------------------------------------------------------------------------------------


###################################################
### code chunk number 7: save-output
###################################################
saveRDS(auc.scores,"AUC.Scores.Eurasia.rds")
saveRDS(sensSpec.score,"Threshholds.Eurasia.rds")


###-----------------------------------------------------------------------------------------------------------------------

# Important objects

# The points used for model building
saveRDS(chukar.present, "chukar.present.pts.rds")
saveRDS(ws.background.pts, "ws.background.pts.rds"  )

#Saster stacks of predicted models
# saveRDS(model.stack.svm, "stack.svm.Eurasia.rds" )
# saveRDS(model.stack.gbm,"stack.gbm.Eurasia.rds"  )
# saveRDS(model.stack.rf, "stack.rf.Eurasia.rds" )
# saveRDS(model.stack.ann,"stack.ann.Eurasia.rds" )
# saveRDS(model.stack.maxent, "stack.maxent.Eurasia.rds" ) 


names(model.stack.svm) <- c("svmFold1", "svmFold2", "svmFold3", "svmFold4", "svmFold5")
names(model.stack.gbm) <- c("gbmFold1", "gbmFold2", "gbmFold3", "gbmFold4", "gbmFold5")
names(model.stack.rf) <- c("rfFold1", "rfFold2", "rfFold3", "rfFold4", "rfFold5")
names(model.stack.ann) <- c("annFold1", "annFold2", "annFold3", "annFold4", "annFold5")
names(model.stack.maxent) <- c("maxentFold1", "maxentFold2", "maxentFold3", "maxentFold4", "maxentFold5", "?")

plot(model.stack.maxent)

writeRaster(model.stack.svm, filename=names(model.stack.svm), bylayer=TRUE,format="GTiff")
writeRaster(model.stack.gbm, filename=names(model.stack.gbm), bylayer=TRUE,format="GTiff")
writeRaster(model.stack.rf, filename=names(model.stack.rf), bylayer=TRUE,format="GTiff")
writeRaster(model.stack.ann, filename=names(model.stack.ann), bylayer=TRUE,format="GTiff")
writeRaster(model.stack.maxent, filename=names(model.stack.maxent), bylayer=TRUE,format="GTiff")

write.csv(chukar.present, 'Chukar.present.data.csv'  )
write.csv(ws.background.pts, 'Chukar.background.data.csv' )


