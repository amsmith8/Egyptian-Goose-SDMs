###################################################
### code chunk number 1: set up environment 
###################################################


#rm(list=ls(all=TRUE))


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

######## ########  Set working directory ######## ########

# set wd so files can be pulled 
# Dependent on which computer to use 

#work.dir<- setwd("/Users/austinsmith/Desktop/A.chukar_dismo")

getwd() # check to see if correct

###-----------------------------------------------------------------------------------------------------------------------

###################################################
### code chunk number 2: rasters & raw data
###################################################

###  BioClim
########  ########  ########

tic("BioClim data")
# Pull from the web
bioclim.data <- getData(name = "worldclim",var = "bio",res = 2.5 )
plot(bioclim.data$bio7)
#names(bioclim.data)
#set crs for project. Refer to/ set this for all layers
Model_CRS <- crs( bioclim.data$bio1 )
#Model_CRS # check
toc()

######## ########     Landcover Raster File     ######## ########
######## ########  ##########################  ######## ########
tic("landcover - stacking")
land.cover <- raster("/Users/austinsmith/Desktop/MCD12C1.A2018001.006.2019200161458_MOD12C1.tif")
crs(land.cover) <- Model_CRS
#plot( land.cover )

# Layer currently in gardient form
#change colors to match MODIS
LCcolors <- c("blue","seagreen3","seagreen4","springgreen2","springgreen1","forestgreen",
              "darkolivegreen3","darkolivegreen3","darkkhaki","khaki3","orange1","mediumaquamarine",
              "goldenrod1","red","goldenrod2","white","gray")
# Create plot to check
par(mar=c(1,1,1,1))
plot(land.cover , col = LCcolors,
     #xlim=c(17, 122.9893), ylim=c(21.56769, 60), # Chukar native/ historic range
     axes = FALSE,
     box = FALSE,
     legend = FALSE )

### Create raster stack
LCTypes <- c("WAT","ENF", "EBF", "DNF", "DBF", "MIXF",
             "CSHB", "OSHB", "WSAV", "SAV", "GRA","WET",
             "CROP", "URB", "CNVM", "SNOW", "BAR")

#### Function
raster_values_stack <- function(r){
  raster.stack <- stack() # make empty stack
  max <- cellStats(r, max) # identifies max value
  for (i in 0:max){
    m <- mask(r, r == i, maskvalue = FALSE) # isolates pixels with value i
    raster.stack <- stack(raster.stack, m)
  }
  raster.stack # produces final stack
}

#######
landcover.data <- raster_values_stack(land.cover)
names(landcover.data) <- LCTypes

Binary.Split <- function( layer ) {
  max <- cellStats(layer, max) # measures the max value in the grid, should be the landcover identifier.
  layer[layer== max ] <- 1 # convert identifier to 1
  layer[is.na(layer[])] <- 0 # masked cells (cells with measure NA) are now 0s (not ID)
  layer # creates raster
}

# Set each layer to binary identity.
# The following was unable to run ewithout error in loop. Needed to complete individually. May be
# a package problem.

ENF <- Binary.Split( landcover.data$ENF )
EBF <- Binary.Split( landcover.data$EBF )
DNF <- Binary.Split( landcover.data$DNF )
DBF <- Binary.Split( landcover.data$DBF )
MIXF <- Binary.Split( landcover.data$MIXF )
CSHB <- Binary.Split( landcover.data$CSHB )
OSHB <- Binary.Split( landcover.data$OSHB )
WSAV <- Binary.Split( landcover.data$WSAV )
SAV <- Binary.Split( landcover.data$SAV )
GRA <- Binary.Split( landcover.data$GRA )
WET <- Binary.Split( landcover.data$WET )
CROP <- Binary.Split( landcover.data$CROP )
URB <- Binary.Split( landcover.data$URB )
CNVM <- Binary.Split( landcover.data$CNVM )
SNOW <- Binary.Split( landcover.data$SNOW )
BAR <- Binary.Split( landcover.data$BAR)
WAT <- Binary.Split( landcover.data$WAT )

##### Stack of binary rasters
landcover.data<- stack(  ENF, EBF, DNF, DBF, MIXF,
                         CSHB, OSHB, WSAV, SAV, GRA,
                         WET, CROP, URB, CNVM, SNOW,
                         BAR, WAT )

# To free memory, remove individual files, keeps the stack
rm(list = c("WAT","ENF", "EBF", "DNF", "DBF", "MIXF",
            "CSHB", "OSHB", "WSAV", "SAV", "GRA","WET",
            "CROP", "URB", "CNVM", "SNOW", "BAR"))
toc()

######## ########     Terrain Raster File     ######## ########
######## ########  ##########################  ######## ########

### Elevation
########  ########  ########
elevation <- raster("/Users/austinsmith/Downloads/ETOPO1_Bed_g_geotiff.tif")
crs(elevation) <- Model_CRS

# Remove  ocean data
elevation <- mask(elevation, elevation >= 0, maskvalue = FALSE)
# Remove Antarctica
elevation <- crop(elevation, c(-180, 180, -62, 83.57027))
#plot(elevation)

###  Slope
########  ########  ########
tic("slope")
slope <- terrain(elevation, opt ="slope", unit = "degrees", neighbors = 8)
toc()
###  Aspect
########  ########  ########
tic("aspect")
aspect <- terrain(elevation, opt ="aspect", unit = "degrees", neighbors = 8)
toc()
###  Terrain Roughness Index
########  ########  ########
tic("TRI")
TRI <- terrain(elevation, opt ="TRI", unit = "degrees", neighbors = 8) #Terrain Roughness Index
toc()

### code chunk number 3:  create stack
###################################################
terrain.data <- stack(elevation, slope, aspect, TRI)
names(terrain.data) <- c("elevation", "slope", "aspect", "TRI")

#resampled.terr <- projectRaster( terrain.data, bioclim.data, method = 'ngb' )

rm(list = c("elevation", "slope", "aspect", "TRI"))

###-----------------------------------------------------------------------------------------------------------------------

###################################################
### code chunk number 3: data cleaning & combined
###################################################

###  Rescample & align
########  ########  ########
tic("resample 1")
resampled.lc <- projectRaster( landcover.data, bioclim.data, method = 'ngb' ) # this aligns the land cover stack to the BioClim stack
toc()

rm(list = c("landcover.data"))

tic("resample 2")
resampled.terr <- projectRaster( terrain.data , bioclim.data, method = 'ngb' )
toc()

rm(list = c("terrain.data"))

###  Stack all sub-stacks
########  ########  ########
tic("normalize and finish already")
temp.stack <- stack(bioclim.data, resampled.lc,resampled.terr)
temp.stack # check to see if done
crs(temp.stack) <- Model_CRS

# Saves files to disk so we can call on it in other files/script/etc
#  Files were manually placed in a folder called "super.brick"

#writeRaster(temp.stack, filename=names(temp.stack), bylayer=TRUE,format="GTiff") # so we can call on it in other files/script/etc

super.stack <- stack()
stack.wd <- setwd("/Users/austinsmith/Desktop/A.chukar_dismo/super.brick/")
getwd()
#tif.files <- list.files(stack.wd , pattern = ".tif")

files <- list.files( full.names= TRUE, pattern = ".tif")

super.stack <- stack()
for (i in 1:length(files)){
  r <- raster(files[i])
  super.stack <- stack(super.stack, r)
}

###  Normalize layers 
########  ########  ########

# # function to normalize all layers so each layer is weighted the same during model building
normalize <- function(x) {
  min <- raster::minValue(x)
  max <- raster::maxValue(x)
  return((x - min) / (max - min))
}
# 
# # normalize data
super.stack <- normalize(super.stack)

setwd("/Users/austinsmith/Desktop/A.chukar_dismo")
saveRDS(super.stack,"super.stack.rds")

# ###### PCA
# 
# library(RStoolbox)
# 
# tic("PCA")
# pca.stack <- rasterPCA(bioclim.data, nSamples = 100000,  spca = TRUE)
# toc()
# tic("rasters")
# writeRaster(pca.stack$map, filename=names(pca.stack$map), bylayer=TRUE,format="GTiff")
# saveRDS(pca.stack,"pca.stack.rds")
# toc()
# 
# a <- raster("./PC1.tif")
# plot(a)

###-----------------------------------------------------------------------------------------------------------------------
