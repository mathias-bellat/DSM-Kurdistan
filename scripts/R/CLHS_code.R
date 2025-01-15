####################################################################
# Script for Cluster Latin hypercube sampling                      #
#                                                                  #                                                  
# Author: Mathias Bellat                                           #
# Affiliation : Tubingen University                                #
# Creation date : 01/08/2024                                       #
# E-mail: mathias.bellat@uni-tubingen.de                           #
####################################################################

# 00 Preparation ###############################################################
# 0.1 Prepare environment ======================================================

# Folder check
getwd()

# Set folder direction
setwd()

# Clean up workspace
rm(list = ls(all.names = TRUE))

# 0.2 Install packages =========================================================

# Load packages
install.packages("pacman")
library(pacman) #Easier way of loading packages
pacman::p_load(sf, mapview,raster, clhs) # Specify required packages and download it if needed

# 0.3 Show session infos =======================================================

sessionInfo()

# 1 Import and prepare the data  -----------------------------------------------
# 1.1 Select the the files #####################################################
Area <- st_read("./data/Sampling_area_2022.gpkg", layer = "Sampling_area_2022")  # Change the area for 2022 and 2023 campaigns
DEM <- raster("./data/DEM.tif")
Geomorpho <- raster("./data/Geomorphology.tif")
Soil_erosion <- raster("./data/RUSLE_map.tif")
Soil_prediction <- raster("./data/Soil_classification_prediction.tif")
Slope <- raster("./data/Slope.tif")
Wetness <- raster("./data/TWI.tif")

# 1.2 Clean and prepare the DEM ################################################

# Crop the DEM to the study area
DEM <- crop(DEM, Area, inverse=FALSE, updatevalue=NA, updateNA=TRUE)

# Mask the DEM to the study area
DEM <- mask(DEM, Area, inverse=FALSE, updatevalue=NA, updateNA=TRUE)

# 1.3 Clean and prepare the covariates #########################################

# Create a list of the different variables
list <- list(Geomorpho, Soil_erosion, Soil_prediction, Slope, Wetness) 
names(list) <- c("Geomorpho","Soil_erosion","Soil_prediction","Slope", "Wetness")

# Set the CRS WGS84 UTM38N to all raster
ZoneUTM <- c("+init=epsg:32638")
for (i in names(list)) { 
  r <- list[[i]]
  list[[i]] <- projectRaster(r, crs = ZoneUTM)
}

# Resample according to the DEM 
# With NGB only for the categorical raster
list[[1]]  <- resample(list[[1]], DEM, "ngb")
list[[3]]  <- resample(list[[3]], DEM, "ngb")

# With bilinear only for the numerical raster
list[[2]]  <- resample(list[[2]], DEM, "bilinear")
list[[4]]  <- resample(list[[4]], DEM, "bilinear")
list[[5]]  <- resample(list[[5]], DEM, "bilinear")

# 1.4 Stack all the predictors #################################################
stack <- stack(DEM, list$Geomorpho, list$Soil_erosion,list$Soil_prediction,
              list$Slope, list$Wetness)  # stack all the data in a stacked raster

plot(stack$DEM)
# 1.5 Save and export the predictors ###########################################
writeRaster(stack ,filename = "./export/predictors.tif", format="GTiff",overwrite=T)

rm(list = ls())

# 2 Create the Conditioned Latin Hypercube Sampling  ---------------------------
# 2.1 Import and prepare the files #############################################
predictors <- raster::stack("./export/predictors.tif")
names(predictor@layers) <- c("DEM","Geomorpho","Soil_erosion","Soil_prediction","Slope", "Wetness" )

# Convert the file into Data Frame
PredForMap  <- as(predictors, 'SpatialPixelsDataFrame')
PredForMap <- data.frame(PredForMap)
stack_frame <- PredForMap[complete.cases(PredForMap),] # Remove NA values

names(stack_frame) <- c("DEM","Geomorpho","Soil_erosion","Soil_prediction","Slope", "Wetness", "x", "y")
stack_frame$Geomorpho <- round(stack_frame$Geomorpho, digits=0) # If ever the categorical classification did not work
stack_frame$Soil_prediction <- round(stack_frame$Soil_prediction, digits=0) # If ever the categorical classification did not work
stack_frame <- na.omit(stack_frame) 

save.image("./export/Pre-process.RData") # Save the image

# 2.2 Conditioned Latin Hypercube sampling parameters ##########################

# Select the predictors for the Conditioned Latin Hypercube
preds <- c("DEM","Geomorpho","Soil_erosion","Soil_prediction","Slope","Wetness") 

# Set the size of the sampling 
c = 110  

# Set a seed
set.seed(1070) 

# Run the sampling
res <- clhs(stack_frame[, preds], size = c,  tdecrease = 0.95, iter = 50000, progress = FALSE, simple = FALSE) 

# 2.3 Export the results #######################################################

CLHS_sampled_res <- stack_frame[res$index_samples, ] #fit the results with the line in the table
CLHS_sampled_df <- SpatialPoints(coords =  CLHS_sampled_res[, c("x", "y")], proj4string=CRS("+init=epsg:32638")) #convert to spatial dataframe point

plot(CLHI_sampled_df)
write.table(CLHI_sampled_df, "./export/CLHS_2022.csv", col.names = TRUE, sep = ";", row.names = FALSE, 
          fileEncoding = "UTF-8")
