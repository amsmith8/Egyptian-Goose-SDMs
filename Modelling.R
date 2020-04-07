

# Step1: Set up working environemnt
####### ######## ######## ####### ######## ########

setwd("/Users/austinsmith/Desktop/Model_40" )

library(raster)
library(RColorBrewer)
library("dismo")
library("maptools")
library("raster")
library("rgdal")
library("rgeos")
#library(sf)
library(sp)


# Step 2 : Prepare data and functions
####### ######## ######## ####### ######## ########

# Stored folders of tiff files
svm.tiffs <- list.files(path = "./model.svm.folds", full.names= TRUE, pattern = ".tif")
gbm.tiffs <- list.files(path = "./model.gbm.folds", full.names= TRUE, pattern = ".tif")
rf.tiffs <- list.files(path = "./model.rf.folds", full.names= TRUE, pattern = ".tif")
ann.tiffs <- list.files(path = "./model.ann.folds", full.names= TRUE, pattern = ".tif")
maxent.tiffs <- list.files(path = "./model.maxent.folds", full.names= TRUE, pattern = ".tif")


# Function 1
#--------------------------------------------------
# function to stack raster layers
listToStack <- function(tiff.list){
  model.stack <- stack()
  for (i in 1:length(tiff.list)){
  r <- raster(tiff.list[i])
  model.stack <- stack(model.stack, r)}
  model.stack
}


#compile tiffs to stacks  
model.stack.svm <- listToStack(svm.tiffs)
model.stack.gbm <- listToStack(gbm.tiffs)
model.stack.rf <- listToStack(rf.tiffs)
model.stack.ann <- listToStack(ann.tiffs)
model.stack.maxent <- listToStack(maxent.tiffs)


# rasterize Ac.poly
Ac.poly <- readShapePoly("/Users/austinsmith/Desktop/A.chukar_dismo/Alectoris_chukar/Alectoris_chukar.shp") # via Bird Life
crs(Ac.poly) <- crs(model.stack.rf$rfFold1)
#r <- rasterize(Ac.poly, Final.model, field=1)
#plot(Ac.poly, add = T)


# Plot colors  
palette<-colorRampPalette(c("navy","goldenrod1","red"))
palette2 <- colorRampPalette(c("skyblue2","firebrick"))


# Step 3 : Accuracy - function,s images, plots
####### ######## ######## ####### ######## ########


# Function 2
#--------------------------------------------------
# Function to take stack, and compute the accuracy of each prediction given theshhold data 
stack.accuracy <- function(model.stack, polygn, thresholds, main.title = "Title"){
  list <- c()
  plot.stack <- stack()
  stack <- model.stack
  thresholds1 <- thresholds
 
  for (i in 1:5){
    Final.model <- (stack[[i]] > thresholds1[1,i] )
    #plot(Final.model, col=palette2(2), legend = FALSE)
    #plot(polygn, add = T)
    r <- rasterize(Ac.poly, Final.model, field=1)
    a1 <- Final.model
    a1[a1== 1 ] <- 0
    r2 <- merge(r,a1)
    # Divide the number of error
    error <- cellStats(r2 != Final.model , sum)/8667720 # This is the number of cells (land) used in model
    acc <- 1-error # accuracy
    list <- c(list, acc)
    plot.stack <- stack(plot.stack,Final.model)
  }
  names(plot.stack) <- c("Fold 1","Fold 2","Fold 3","Fold 4","Fold 5")
  plot(plot.stack,
       col=palette2(2), 
       ylim = c(-65,90), 
       legend = FALSE, 
       box.lty=0 )
  legend("bottomright", 
         title = main.title,
         legend=c("Suitable","Not Suitable"),
         fill = c("firebrick","skyblue2"), 
         cex=1.5)
  list
}

# Calculate 
valid.stack.svm <- stack.accuracy(model.stack = model.stack.svm , Ac.poly, Threshholds.Eurasia[1,2:6], main.title = "SVM Folds")
valid.stack.gbm <- stack.accuracy(model.stack =  model.stack.gbm, Ac.poly, Threshholds.Eurasia[2,2:6], main.title = "GBM Folds")
valid.stack.rf <- stack.accuracy(model.stack =  model.stack.rf, Ac.poly, Threshholds.Eurasia[3,2:6], main.title = "RF Folds")
valid.stack.ann<- stack.accuracy(model.stack =  model.stack.ann, Ac.poly, Threshholds.Eurasia[4,2:6], main.title = "ANN Folds")
valid.stack.maxent <- stack.accuracy(model.stack =  model.stack.maxent, Ac.poly, Threshholds.Eurasia[5,2:6], main.title = "MaxEnt Folds")


# valid.stack.svm
# valid.stack.gbm
# valid.stack.rf
# valid.stack.ann
# valid.stack.maxent 

# list statistics 
means <- c( mean(valid.stack.svm), mean(valid.stack.gbm),mean(valid.stack.rf),
            mean(valid.stack.ann),mean(valid.stack.maxent))

sdev <- c( sd(valid.stack.svm), sd(valid.stack.gbm),sd(valid.stack.rf),
           sd(valid.stack.ann),sd(valid.stack.maxent))


data.frame.valids <- as.data.frame(rbind(valid.stack.svm,valid.stack.gbm,valid.stack.rf,
                           valid.stack.ann, valid.stack.maxent))
row.names(data.frame.valids) <- c("SVM", "GBM", "RF", "ANN", "MaxEnt")
colnames(data.frame.valids) <- c("Fold 1","Fold 2","Fold 3","Fold 4","Fold 5")
data.frame.valids


AUC.table <- as.data.frame( as.table(as.matrix(AUC.Scores.Eurasia[2:6])) )
AUC.table$Var1 <- c("SVM", "GBM", "RF", "ANN", "MaxEnt")
#convert data frame to table for strip charts
Valid.table<- as.data.frame( as.table(as.matrix(data.frame.valids)) )
Valid.table <- Valid.table[,-2]


####
# Basic bar plots of means +/- se with jittered points
plot1 <- ggbarplot(AUC.table[,-2], x = "Var1", y = "Freq", ylim = c(0.875,1),position = position_dodge(0.8),
                 #title = "a)",
                 xlab = "",
                 ylab = "AUC",
                 fill = "Var1",
                 palette  =  c("orange","forestgreen","blue","purple", "red"),
                 legend = "none",
                 add = c("mean_sd", "jitter"))
#ggpar(qw1,legend = "none")

plot2<- ggbarplot(Valid.table, x = "Var1", y = "Freq", ylim = c(0.875,1),position = position_dodge(0.8),
                 #title = "b)",
                 xlab = "",
                 ylab = "Accuracy",
                 fill = "Var1",
                 palette  =  c("orange","forestgreen","blue","purple", "red"),
                 legend = "none",
                 add = c("mean_sd", "jitter"))
#ggpar(plot2,legend = "none")

#stack the plots to same image
ggarrange(plot1, plot2, legend = NULL,
          labels = c("Model testing results", " Model predictions"),
          ncol = 1, nrow = 2)


###-------------------------------------------------------------------------------------------------

# Step 4 : Normalize - function,s images, plots
####### ######## ######## ####### ######## ########


# Function 3
#--------------------------------------------------
# # function to normalize all layers so each layer is weighted the same during model building
normalize <- function(x) {
  min <- raster::minValue(x)
  max <- raster::maxValue(x)
  return((x - min) / (max - min))
}


# Function 4
#--------------------------------------------------
# function to tranform thresholds normalized rasters)
normalizeThreshhold <- function(r,x) {
  max <- cellStats(r, max)
  min <- cellStats(r, min)
  return((x - min) / (max - min))
}


# Function 5
#--------------------------------------------------
# function to normailize thresholds in raster stack
normalize.stack.threshholds <- function(norm.stack, thrshlds){
  norm.thresh <- c()
  for (i in 1:5){
    norm.thresh[i] <- normalizeThreshhold(norm.stack[[i]], thrshlds[[i]])
  }
  norm.thresh
}


## Have all models on same scale
svm.scale <- normalize(model.stack.svm )
gbm.scale <- normalize(model.stack.gbm )
rf.scale <- normalize(model.stack.rf )
ann.scale <- normalize(model.stack.ann )
maxent.scale <- normalize(model.stack.maxent )


#create new thresholds
norm.thresh.svm <-  normalize.stack.threshholds(model.stack.svm, Threshholds.Eurasia[1, 2:6])
norm.thresh.gbm <-  normalize.stack.threshholds(model.stack.gbm, Threshholds.Eurasia[2, 2:6])
norm.thresh.rf <-  normalize.stack.threshholds(model.stack.rf, Threshholds.Eurasia[3, 2:6])
norm.thresh.ann <-  normalize.stack.threshholds(model.stack.ann, Threshholds.Eurasia[4, 2:6])
norm.thresh.maxent <-  normalize.stack.threshholds(model.stack.maxent, Threshholds.Eurasia[5, 2:6])



plot(svm.scale, col = palette(20), main = "SVM Folds Suitbaility Results")
plot(gbm.scale, col = palette(20), main = "GBM Folds Suitbaility Results")
plot(rf.scale, col = palette(20), main = "RF Folds Suitbaility Results")
plot(ann.scale, col = palette(20), main = "ANN Folds Suitbaility Results")
plot(maxent.scale, col = palette(20), main = "MaxEnt Folds Suitbaility Results")

###-------------------------------------------------------------------------------------------------

# Step 5 : Means - calculate, plot
####### ######## ######## ####### ######## ########


# Calculate the mean of normalized stacks 
mean.svm.scale <- mean(svm.scale)
mean.gbm.scale <- mean(gbm.scale)
mean.rf.scale <- mean(rf.scale)
mean.ann.scale <- mean(ann.scale)
mean.maxent.scale <- mean(maxent.scale)



plot(mean.svm.scale, col = palette(20), main = "SVM Mean Suitbaility Result")
plot(mean.gbm.scale, col = palette(20), main = "GBM Mean Suitbaility Result")
plot(mean.rf.scale, col = palette(20), main = "RF Mean Suitbaility Result")
plot(mean.ann.scale, col = palette(20), main = "ANN Mean Suitbaility Result")
plot(mean.maxent.scale, col = palette(20), main = "MaxEnt Mean Suitbaility Result")



#calculate mean thresholds 
mean.norm.thresh.svm  <- mean(norm.thresh.svm )
mean.norm.thresh.gbm  <- mean(norm.thresh.gbm )
mean.norm.thresh.rf  <- mean(norm.thresh.rf )
mean.norm.thresh.ann  <- mean(norm.thresh.ann )
mean.norm.thresh.maxent  <- mean(norm.thresh.maxent )


mean.stack <- stack(mean.svm.scale,mean.gbm.scale,mean.rf.scale,mean.ann.scale,mean.maxent.scale)
norm.mean.thresh <- data.frame(mean.norm.thresh.svm, mean.norm.thresh.gbm, mean.norm.thresh.rf, mean.norm.thresh.ann, mean.norm.thresh.maxent) 

#Calculated means 
means.acc.scores <- stack.accuracy(mean.stack,thresholds = norm.mean.thresh)
means.acc.scores


# Mean Plots
plot(mean.svm.scale > mean.norm.thresh.svm, main = "SVM Averaged Prediction", col=palette2(2), legend = FALSE)
legend("bottomleft", 
       legend=c("Suitable","Not Suitable"),
       fill = c("firebrick","skyblue2"), 
       cex=1.0)
plot(mean.gbm.scale > mean.norm.thresh.gbm, main = "GBM Averaged Prediction", col=palette2(2), legend = FALSE)
legend("bottomleft", 
       legend=c("Suitable","Not Suitable"),
       fill = c("firebrick","skyblue2"), 
       cex=1.0)
plot(mean.rf.scale > mean.norm.thresh.rf, main = "RF Averaged Prediction", col=palette2(2), legend = FALSE)
legend("bottomleft", 
       legend=c("Suitable","Not Suitable"),
       fill = c("firebrick","skyblue2"), 
       cex=1.0)
plot(mean.ann.scale > mean.norm.thresh.ann, main = "ANN Averaged Prediction", col=palette2(2), legend = FALSE)
legend("bottomleft",
       legend=c("Suitable","Not Suitable"),
       fill = c("firebrick","skyblue2"), 
       cex=1.0)
plot(mean.maxent.scale > mean.norm.thresh.maxent, main = "MaxEnt Averaged Prediction", col=palette2(2), legend = FALSE)
legend("bottomleft", 
       legend=c("Suitable","Not Suitable"),
       fill = c("firebrick","skyblue2"), 
       cex=1.0)


###-------------------------------------------------------------------------------------------------


# Step 6 :  Final Model - Ensemble 
####### ######## ######## ####### ######## ########

ensemble.threshhold <- rowMeans(norm.mean.thresh )

ensemble <- mean(mean.stack) 

plot(ensemble, col = palette(20), main = "Ensemble Suitbaility Result")
plot(ensemble >ensemble.threshhold, col = palette2(2), main = "Ensemble Prediction", legend = FALSE)
legend("bottomleft", 
       legend=c("Suitable","Not Suitable"),
       fill = c("firebrick","skyblue2"), 
       cex=1.0)


Final.model <- (ensemble > ensemble.threshhold)
r <- rasterize(Ac.poly, Final.model, field=1)
a1 <- Final.model
a1[a1== 1 ] <- 0
r2 <- merge(r,a1)
# Divide the number of error
error <- cellStats(r2 != Final.model , sum)/8667720 # This is the number of cells (land) used in model
acc <- 1-error # accuracy
acc

