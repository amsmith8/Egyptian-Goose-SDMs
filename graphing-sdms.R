# # GLM 
# for(i in 1:k){
#   
#   plot(crop( model.stack.glm[[i]] > sensSpec.score[1,i+1], c(-165, -68, 19.5, 68)), main='GLM N. America prediction')
#   plot(Naturalized, add = T)
# }
stack.means <-rowMeans(sensSpec.score[,-1])


dev.off()
plot(ws, main = "Global Chukar Distribution")
# Check on map 
plot(Native, add = T, col = "blue")
plot(Naturalized, add = T, col = "red")
legend("bottom", legend=c("Native", "Naturalized"),
       fill=c("blue", "red"),cex=0.8)




### Calculate mean vlues for threshold, store in vector  
stack.means <-rowMeans(sensSpec.score[,-1])

### Mean value for model stack
#mean.glm <-mean(model.stack.glm)
mean.svm <-mean(model.stack.svm)
mean.gbm <-mean(model.stack.gbm)
mean.rf <- mean(model.stack.rf)
mean.ann <-mean(model.stack.ann)
mean.maxent <-mean(model.stack.maxent)

mean.stack <- stack( mean.svm,mean.gbm,
                 mean.rf,mean.ann,mean.maxent)

plot(model.stack.rf$layer.1.1 > sensSpec.score[3,2])

plot(mean.stack)


# function to tranform thresholds normalized rasters)
normalizeThreshhold <- function(r,x) {
  max <- cellStats(r, max)
  min <- cellStats(r, min)
  return((x - min) / (max - min))
}

#store and transform threshholds so all P/A are same scale 
norm.thresh.means <- rep(NA, 5)
for (i in 1:6){
  norm.thresh.means[i] <- normalizeThreshhold(mean.stack[[i]], stack.means[i])
}
norm.thresh.means 



## Have all models on same scale 
#glm.scale <- normalize(model.stack.glm ) 
svm.scale <- normalize(model.stack.svm )
gbm.scale <- normalize(model.stack.gbm )
rf.scale <- normalize(model.stack.rf )
ann.scale <- normalize(model.stack.ann )
maxent.scale <- normalize(model.stack.maxent ) 

# Calculate the mean of normalized stacks 
mean.glm.scale <- mean(glm.scale)
mean.svm.scale <- mean(svm.scale)
mean.gbm.scale <- mean(gbm.scale)
mean.rf.scale <- mean(rf.scale)
mean.ann.scale <- mean(ann.scale)
mean.maxent.scale <- mean(maxent.scale)

mean.glm.scale[1]

#par(mfrow=c(,1))
#plot(mean.glm.scale > norm.thresh.means[1], main = "GLM")
plot(mean.svm.scale > norm.thresh.means[1], main = "SVM")
plot(mean.gbm.scale > norm.thresh.means[2], main = "GBM")
plot(mean.rf.scale > norm.thresh.means[3], main = "RF")
plot(mean.ann.scale > norm.thresh.means[4], main = "ANN")
plot(mean.maxent.scale > norm.thresh.means[5], main = "MaxEnt")


ensemble <- mean( mean.svm.scale,mean.gbm.scale,
                 mean.rf.scale,mean.ann.scale,mean.maxent.scale)


plot(ensemble,axes=FALSE, box=FALSE, main = "Final Model Raw")
plot(ensemble > mean(norm.thresh.means), main = "Mean P/A prediction for Chukar")

plot(crop(ensemble, c(-165, -68, 19.5, 68)), main='Ensemble N. America prediction') 
plot(Naturalized, add = T)

Final.model <- (ensemble > mean(norm.thresh.means))

Final.model <- (ensemble > 0.25)

plot(Final.model, axes=FALSE, box=FALSE, legend=FALSE, main = "Final Model P/A")
plot(Native, add = T)
plot(Naturalized, add = T)

unique(values(Final.model))

r <- rasterize(Ac.poly, Final.model, field=1)
plot(r)






#us <- rasterize(spdf, Final.model, field=1)
#plot(us)

freq(a1)
a1<- Final.model 
plot(a1)
a1[a1== 1 ] <- 0
plot(a1)

r2 <- merge(r,a1)
plot(r2)

unique(values(r2))

plot(r2 != Final.model, legend=FALSE, main = "Model Error" )

a.diff <-r2 != Final.model

plot(a.diff)

plot(error.plot, main = "Model Error")





#Divide the number of error 
cellStats(r2 != Final.model, sum)/8667720 # This is the number of cells (land) used in model 








#######

#calculations

native.crop <- crop(land.cover, Native)
native.mask <- mask(land.cover, Native)
plot(native.mask, col = LCcolors,)
frequencies <- freq(native.mask, progress = 'text', merge =T)
f <- c(frequencies[1:17,2])


native.freq <- data.frame(LCTypes, frequencies[1:17,2] )

pie(f, labels = native.freq$LCTypes, col = LCcolors, main="Landcover distributiion- Chukar native range")


naturalized.crop <- crop(land.cover, Naturalized)
naturalized.mask <- mask(land.cover, Naturalized)
plot(naturalized.mask, col = LCcolors,)
nat.frequencies <- freq(naturalized.mask, progress = 'text', merge =T)
nat.f <- c(nat.frequencies[1:14,2])

piepercent<- round(100*nat.f/sum(nat.f), 1)
class(LCTypes)

pie(nat.f, labels = piepercent,radius = 1, cex = .95,
    # labels = c("WAT","ENF", "EBF", "MIXF",
    #                   "CSHB", "OSHB", "WSAV", "SAV", "GRA","WET",
    #                   "CROP",  "CNVM", "SNOW", "BAR"), 
    col = c("blue","seagreen3","seagreen4","forestgreen",
            "darkolivegreen3","darkolivegreen3","darkkhaki","khaki3","orange1","mediumaquamarine",
            "goldenrod1","goldenrod2","white","gray"), 
    main="Landcover distributiion- Chukar naturalized range")
 




###
library(corrplot)

A <- cor(sdm.data[3:28])
corrplot(A)
