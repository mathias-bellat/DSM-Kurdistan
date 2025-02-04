####################################################################
# This script is  the soil depth prediction for the                #
# Northern Kurdistan region (Iraq) soils                           #                       
#                                                                  #                                                   
# Author: Mathias Bellat                                           #
# Affiliation : Tubingen University                                #
# Creation date : 17/06/2024                                       #
# E-mail: mathias.bellat@uni-tuebingen.de                          #
####################################################################


# 0 Environment setup ##########################################################

# 0.1 Prepare environment ======================================================

# Folder check
getwd()

# Set folder direction
setwd()

# Clean up workspace
rm(list = ls(all.names = TRUE))

# 0.2 Install packages =========================================================

install.packages("pacman")        #Install and load the "pacman" package (allow easier download of packages)
library(pacman)

pacman::p_load(ggplot2, terra, mapview, sf, raster, dplyr, corrplot, viridis, rsample, caret, parsnip, tmap, 
               DescTools, patchwork, quantregForest, ncdf4, wesanderson, cblindplot, grid, gridExtra, ggspatial, cowplot)

# 0.3 Show session infos =======================================================

sessionInfo()

# 01 Import data sets ##########################################################

# 01.1 Import soils depths =====================================================

# From the data accessible at the https://doi.org/10.1594/PANGAEA.973714
# We added 25 NULL results from bare rock point spoted in remote sensing images
depth <- read.csv ("./data/Soil_depth.csv", sep=";")
depth$Depth <- as.numeric(depth$Depth)

# 01.2 Plot the depth information ==============================================

# Short overview of the values
head(depth)

# Histogram of the values
ggplot(depth, aes(x = Depth)) + 
  geom_histogram(breaks = seq(0, 100, by = 10), alpha = 0.5, fill = 'steelblue', color = 'black') +
  labs(title ="Number of sample for each soil depth in cm", x = "Depth in cm", y = "Number of samples") +
  xlim(0, 100) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) 

# 01.3 Transform and plot the data =============================================

depth$sqrt <- sqrt(depth$Depth)

# Histogram of the sqrt values
ggplot(depth, aes(x = sqrt)) + 
  geom_histogram(breaks = seq(0, 10, by = 0.5), alpha = 0.5, fill = 'steelblue', color = 'black') +
  labs(title ="Number of sample for each soil depth in cm", x = "Sqrt-transformed soil depth in cm", y = "Number of samples") +
  xlim(0, 11) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) 

# 01.4 Plot the depth spatial information ======================================

# Create a spatial dataframe with WGS84 UTM 38 N cordinates
sp_df <- st_as_sf(depth, coords = c("X", "Y"), crs = 32638)

# Plot the values of depth regarding the elevation
mapview(sp_df, zcol = "Depth", legend =TRUE, map.types ="Esri.WorldShadedRelief") 

# 01.5 Import covariates raster ================================================

Landsat <- raster::stack(list.files("./data/Landsat/", full.names = TRUE))
names(Landsat)

Terrain <- raster::stack(list.files("./data/Terrain/", full.names = TRUE))
names(Terrain)

Others <- raster::stack(list.files("./data/Others/", full.names = TRUE))
names(Others)

LST <- raster::stack(list.files("./data/LST/", full.names = TRUE))
names(LST)

DSM_cov_name <- read.table("./data/Covariates_names_DSM.txt")


df_names <- data.frame()
for (i in 1:length(names(Terrain))) {
  col_name <- names(Terrain)[i]
  if (col_name %in% DSM_cov_name$V2) {
    transformed_name <- DSM_cov_name$V1[DSM_cov_name$V2 == col_name]
  } else {
    transformed_name <- paste0("TE.26")
  }
  df_names[i,1] <- transformed_name
  df_names[i,2] <- names(Terrain)[i]
}

t <- nrow(df_names)
for (i in 1:length(names(Landsat))) {
  c <- paste0("LA5.",i)
  df_names[i+t,1] <- c
  df_names[i+t,2] <- names(Landsat)[i]
}

t <- nrow(df_names)
for (i in 1:length(names(LST))) {
  c <- paste0("LST.",i)
  df_names[i+t,1] <- c
  df_names[i+t,2] <- names(LST)[i]
}

t <- nrow(df_names)
for (i in 1:length(names(Others))) {
  col_name <- names(Others)[i]
  if (col_name %in% DSM_cov_name$V2) {
    transformed_name <- DSM_cov_name$V1[DSM_cov_name$V2 == col_name]
  } else {
    transformed_name <- "OT.10"  #Only one new variable
  }
  df_names[i+t,1] <- transformed_name
  df_names[i+t,2] <- names(Others)[i]
}

write.table(df_names,"./data/Covariates_names_soil_depth.txt")
all_cov_name <- rbind(DSM_cov_name, df_names)
all_cov_name <- distinct(all_cov_name)
x <- raster::stack(Terrain, Landsat, LST, Others)
names(x) <- df_names[,1]
x <- rast(x)

terra::writeRaster(x, "./data/Stack_layers_soil_depth.tif", overwrite = TRUE)
covariates <- stack("./data/Stack_layers_soil_depth.tif")

# 01.6 Plot the covariates maps ================================================

reduce <- aggregate(covariates, fact=10, fun=mean)

plot(reduce)

# 01.7 Extract the values ======================================================

# Extract the values of each band for the sampling location
df_cov <- raster::extract(covariates, sp_df, method='simple')

# 01.8 Save the data ===========================================================

# Create a csv out of it
write.csv(data.frame(df_cov), "./data/df_cov_soil_depth.csv")

# 02 Check the data ############################################################

# 02.1 Import the data and merge ===============================================
Covdfgis    <- read.csv("./data/df_cov_soil_depth.csv")
ID    <- 1:nrow(Covdfgis)
SoilCov <- cbind(ID, depth, Covdfgis[,-c(1)])

#Check the data
names(SoilCov)
str(SoilCov)

# 02.2 Check and prepare the data ==============================================

# Prepare data for ML modelling
SoilCovML <- SoilCov[,-c(1:7)]

# 02.3 Plot and export the correlation matrix ==================================
pdf("./export/Correlation.pdf",    # File name
    width = 20, height = 20,  # Width and height in inches
    bg = "white",          # Background color
    colormodel = "cmyk")   # Color model 


# Correlation of the data
corrplot(cor(SoilCovML),  method = "color", col = viridis(200), 
         type = "upper", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 45, # Text label color and rotation
         number.cex = 0.7, # Size of the text labels
         cl.cex = 0.7, # Size of the color legend text
         cl.lim = c(-1, 1)) # Color legend limits

dev.off()

# 02.4 Export and save data ====================================================
save(covariates, SoilCov, file = "./export/save/Pre_process_soil_depth.RData")
rm(list = ls())

# 03 Tune the model ############################################################
load(file = "./export/save/Pre_process_soil_depth.RData")

# 03.1 Create the train and test set ===========================================

#R Remove the site name and x, y columns
df_soilCov <- SoilCov[,c(1,5:length(SoilCov))]

# Set a random seed
set.seed(1070)

# preprocess the layers
preproc <- preProcess(df_soilCov[,4:length(df_soilCov)], method=c("range"))
df_soilCovTrans <- predict(preproc, df_soilCov[,4:length(df_soilCov)])
df_soilCovTrans <- cbind(df_soilCovTrans, df_soilCov[,1:3])

# Create a 80% split
set.seed(1070)
df_soil <- initial_split(df_soilCovTrans, prop = 0.8, strata = ID)
trainData <- training(df_soil)
testData <- testing(df_soil)

# Implement the index
row.names(testData) <- testData$ID
row.names(trainData) <- trainData$ID

# Create the splits
X_train <- trainData[,-c(27:29)]
y_train <- trainData[,c(28)]
y_train_sqrt <- trainData[,c(29)]
X_test <- testData[,-c(27:29)]
y_test <- testData[,c(28)]
y_test_sqrt <- testData[,c(29)]

# 03.2 Create the model train controls with cross-validation  ==================
TrainControl <- trainControl(method="repeatedcv", 10, 3, allowParallel = TRUE, savePredictions=TRUE, verboseIter = TRUE)

# 03.3 Set the grid for mtry and nodesize ======================================

tuneGrid <- expand.grid(mtry = c(1, 3, 5, 7, 10, 12, 15))
nodesize_values <- c(3, 4, 5, 6, 7, 8, 9, 10, 11, 12)

all_results <- list()
set.seed(1070)

# 03.4 Run the model tunning ===================================================
# Loop over different nodesize values and different mtry
for (nodesize in nodesize_values) {
  start_time <- Sys.time()
  qrf_model <- train(x = X_train, y = y_train_sqrt,
                  method = "qrf",
                  trControl = TrainControl,
                  tuneGrid = tuneGrid,
                  ntree = 500,     # Number of trees
                  metric = "RMSE",      # Specify the metric you want to optimize
                  nodesize = nodesize     # Minimum node size
)
  end_time <- Sys.time()
  print(end_time - start_time)
  
# Store the results with the corresponding nodesize
results <- qrf_model$results
results$nodesize <- nodesize  # Add the nodesize information
all_results[[as.character(nodesize)]] <- results
}

final_results <- do.call(rbind, all_results)

# View the combined results for all mtry and nodesize combinations
print(final_results)

# 03.5 Plot the different optimization =========================================

gg1 <- ggplot(final_results, aes(x = mtry, y = Rsquared, color = as.factor(nodesize))) +
  geom_line() +
  geom_point() +
  labs(x = "mtry", y = "R²", color = "nodesize") +
  theme_minimal()

gg2 <- ggplot(final_results, aes(x = mtry, y = RMSE, color = as.factor(nodesize))) +
  geom_line() +
  geom_point() +
  labs(x = "mtry", y = expression("RMSE (" *cm^0.5* ")"), color = "nodesize") +
  theme_minimal()


# Combine R2 and RMSE in one plot
combined_plot <- wrap_plots(gg1, gg2, ncol = 2)
plot(combined_plot)

# Save and export the plot
ggsave("./export/Optimisation of the model.pdf", combined_plot, width = 20, height = 10)
ggsave("./export/Optimisation of the model.png", combined_plot, width = 20, height = 10)

# 04.1 Produce the final model #################################################

nodesize <- 6
mtry  <- expand.grid(mtry = 1)

set.seed(1070)
qrf_final <- train(x = X_train, y = y_train_sqrt,
                   method = "qrf",
                   trControl = TrainControl,
                   tuneGrid = mtry,
                   ntree = 500,     # Number of trees
                   metric = "RMSE",      # Specify the metric you want to optimize
                   nodesize = nodesize )


results_qrf <- qrf_final$results

# View R² and RMSE values for final model
print(results_qrf[, c("mtry", "RMSE", "Rsquared")])

final_model <- qrf_final$finalModel

# 04.2 Show covariates importance ==============================================

importance_df <- as.data.frame(varImp(final_model))

# Convert row names (variables) to a column for ggplot
importance_df$scale <- (importance_df$Overall/sum(importance_df$Overall)*100)
importance_df$Variable <- rownames(importance_df)

# Plot using ggplot2
gg1 <- ggplot(importance_df, aes(x = reorder(Variable, scale), y = scale)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  coord_flip() +
  xlab("Covariates") +
  ylab("Importance scaled (%)") +
  ggtitle("Variable Importance from QRF Model")

gg1

# Save and export the plot
ggsave("./export/Variables importance.png", gg1, width = 10, height = 8)
ggsave("./export/Variables importance.pdf", gg1, width = 10, height = 8)


# 04.3 Predict the test set ====================================================
predictions <- predict(final_model,  X_test)

# 04.4 Get the metrics  ========================================================
metrics <- postResample(predictions[,2], y_test_sqrt)
ccc <- CCC(y_test_sqrt, predictions[,2])
PICP <- (y_test_sqrt >= predictions[,1]) & (y_test_sqrt <= predictions[,3])
PICP <- mean(PICP)*100

Final_stats <- as.data.frame(t(c(metrics, ccc$rho.c$est, PICP)))
colnames(Final_stats) <- c("RMSE", "R²", "MAE", "CCC", "PICP")

# Print the metrics
print(Final_stats)

# Write it
write.table(Final_stats, file="./export/save/Final_stats.txt", row.names = FALSE)

# 04.5 Save the model and train and test sets ==================================
save(trainData, testData, final_model, qrf_model, all_results, qrf_final, file = "./export/save/Model_soil_depth.RData")
rm(list = ls())

# 05 Predict the map ###########################################################
# 05.1 Import the previous data ================================================
load(file = "./export/save/Model_soil_depth.RData")
stack_raster <- stack("./data/Stack_layers_soil_depth.tif")
cov <- read.csv("./data/df_cov_soil_depth.csv")
cov <- cov[,-1]
cov[] <- lapply(cov , as.numeric)
  
# 05.2 Normalise the values from the rasters ===================================
process_layer <- function(layer) {
  # Convert in numeric 
  layer <- as.numeric(layer)
  # Replace the NAs by median
  median_value <- median(layer, na.rm = TRUE)
  layer[is.na(layer)] <- median_value
  return(layer)
}

scaling_params <- lapply(1:ncol(cov), function(i) {
  list(min = min(cov[[i]], na.rm = TRUE), max = max(cov[[i]], na.rm = TRUE))
})

scale_layer <- function(layer, min_val, max_val) {
  (layer - min_val) / (max_val - min_val)
}

stack_scaled <- stack() 

# Apply transformation
for (i in 1:nlayers(raster_stack)) {
  min_val <- scaling_params[[i]]$min
  max_val <- scaling_params[[i]]$max
  
  # Normalised the layer
  layer_processed <- calc(raster_stack[[i]], process_layer)
  layer_scaled <- calc(layer_processed, function(x) scale_layer(x, min_val, max_val))
  stack_scaled <- addLayer(stack_scaled, layer_scaled)
  cat(round((i/nlayers(raster_stack))*100, 1),"% \n")
}

stack_scaled <- rast(stack_scaled)
sum(is.na(values(stack_scaled)))
names(stack_scaled) <- names(raster_stack)
terra::writeRaster(stack_scaled, "./export/Soil_depth/Stack_raster_normalised_soil_depth.tif", overwrite = TRUE)

# 05.3 Prediction map for soil depth ===========================================
rm(list = ls(all.names = TRUE))
load("./export/save/Model_soil_depth.RData")
raster_stack_normalised <- stack("./export/Soil_depth/Stack_raster_normalised_soil_depth.tif")

# Divide in block for allowing the computer to run it
block_info <- blockSize(raster_stack_normalised )

# Create an empty raster for storing predicted values
predicted_raster <- raster_stack_normalised [[1]]  # Use the first layer as a template
predicted_raster <- writeStart(predicted_raster, "./export/Prediction_depth_sqrt.tif", overwrite = TRUE)

# Loop through each block of the raster stack

for (i in 1:block_info$n) {
  # Read block of raster data
  start_time <- proc.time()
  block <- getValuesBlock(raster_stack_normalised , row = block_info$row[i], nrows = block_info$nrows[i])
  block <- as.data.frame(block)
  
  # Process the block 
  predicted_block <- predict(final_model, block,  what = c(0.05, 0.5, 0.95))  
  
  # Write the predicted block to the output raster
  predicted_raster <- writeValues(predicted_raster, predicted_block[,2], block_info$row[i])
  end_time <- proc.time()
  print(end_time - start_time)
  print(i)
}

# Close the creation of the raster
predicted_raster <- writeStop(predicted_raster)

# Rescale the raster and export it
predicted_raster_sqrt <- (predicted_raster)^2
writeRaster(predicted_raster_sqrt,"./export/Prediction_depth.tif", format = "GTiff",overwrite=T)

# 05.4 Uncertainty map of soil depth ===========================================

# Create an empty raster for storing predicted values
uncertainty_raster <- raster_stack_normalised [[1]]  # Use the first layer as a template
uncertainty_raster <- writeStart(uncertainty_raster, "./export/Uncertainty_depth_sqrt.tif", overwrite = TRUE)

# Loop through each block of the raster stack

for (i in 1:block_info$n) {
  # Read block of raster data
  start_time <- proc.time()
  block <- getValuesBlock(raster_stack_normalised, row = block_info$row[i], nrows = block_info$nrows[i])
  block <- as.data.frame(block)
  
  # Process the block 
  predicted_block <- predict(final_model, block,  what = c(0.05, 0.5, 0.95))  
  
  # Write the predicted block to the output raster
  values <- (predicted_block[,3] - predicted_block[,1]) /predicted_block[,2]
  uncertainty_raster <- writeValues(uncertainty_raster, values, block_info$row[i])
  end_time <- proc.time()
  print(end_time - start_time)
  print(i)
}

# Close the creation of the raster 
uncertainty_raster <- writeStop(uncertainty_raster)
uncertainty_raster_sqrt <- (uncertainty_raster)^2
writeRaster(uncertainty_raster_sqrt ,"./export/Uncertainty_depth.tif", format = "GTiff",overwrite=T)


# 06 Visualisation of the prediction ###########################################
# 06.1 Plot the depth maps =====================================================
rm(list = ls())
load(file = "./export/save/Pre_process_soil_depth.RData")
predicted_raster <- raster("./export/Prediction_depth.tif")
uncertainty_raster <- raster("./export/Uncertainty_depth.tif")
survey <- st_read("./data/Survey_Area.gpkg", layer = "Survey_Area")
sp_df <- st_as_sf(SoilCov, coords = c("X", "Y"), crs = 32638)

predicted_raster_resize <- aggregate(predicted_raster, fact=5, fun=mean)
uncertainty_resize <- aggregate(uncertainty_raster, fact=5, fun=mean)

# Plot continuous values of the map
mapview(predicted_raster_resize, at = seq(0,100,20), legend = TRUE, layer.name = "Soil depth prediction (cm)",col.regions = wes_palettes$Zissou1Continuous) +
  mapview(sp_df, zcol = "Depth", at = seq(0,100,20),  legend = TRUE, layer.name = "Soil samples depth (cm)", col.regions = wes_palettes$Zissou1Continuous)

mapview(uncertainty_resize, legend = TRUE)


# 06.3 Visualise for colorblind ================================================

wes_palette_hcl <- function(palette_name, n = 7) {
  wes_colors <- wes_palette(palette_name, n = n, type = "continuous")
  hex_colors <- as.character(wes_colors)  
  return(hex_colors)
}

palette_wes <- rev(wes_palette_hcl("Zissou1", n = 7))
palette_blue <- colorRampPalette(c("lightblue", "blue", "darkblue"))(7)

# Replace the cvd palette with palette_wes for depht and palette blue for uncertainty
gg1 <- cblind.plot(uncertainty_resize, cvd = palette_blue)
gg2 <- cblind.plot(uncertainty_resize, cvd = "deuteranopia")
gg3 <- cblind.plot(uncertainty_resize, cvd = "tritanopia")
gg4 <- cblind.plot(uncertainty_resize, cvd = "protanopia")

grid <- grid.arrange(
  arrangeGrob(gg1, bottom = textGrob("Original", gp = gpar(fontsize = 14))),
  arrangeGrob(gg2, bottom = textGrob("Deuteranopia", gp = gpar(fontsize = 14))),
  arrangeGrob(gg3, bottom = textGrob("Tritanopia", gp = gpar(fontsize = 14))),
  arrangeGrob(gg4, bottom = textGrob("Protanopia", gp = gpar(fontsize = 14))),
  nrow = 2, ncol = 2,
  top = textGrob("Soil depth uncertainty for different visions", gp = gpar(fontsize = 16, fontface = "bold")))

ggsave("./export/Visualisation_soil_depth_uncertainty.png", grid, width = 20, height = 10)
ggsave("./export/Visualisation_soil_depth_uncertainty.pdf", grid, width = 20, height = 10)

dev.off()

# 06.4 Export final maps =======================================================

# Replace missing values by zero as it occurs only in mountainous area
uncertainty_raster[is.na(uncertainty_raster)] <- 0
predicted_raster[is.na(predicted_raster)] <- 0


# For GeoTiff format
soil_depth_map <- stack(predicted_raster, uncertainty_raster)
soil_depth_crop <- crop(soil_depth_map, survey)
soil_depth_mask <- mask(soil_depth_crop, survey, inverse=FALSE, updatevalue=NA, updateNA=TRUE)
crs(soil_depth_mask) <- "EPSG:32638"
x <- rast(soil_depth_mask)
x_repro <- project(x, "EPSG:4326")
names(x_repro) <- c("Prediction", "Uncertainty")
terra::writeRaster(x_repro,"./export/Soil_depth_prediction_map.tif", overwrite=TRUE)

# For netCDF format
CDF_df <- lapply(1:nlayers(x_repro), function(i) {
  rast(x_repro[[i]])
})

names(CDF_df) <- c("Prediction", "Uncertainty")

soil_to_netcdf <- function(soil_list, output_file, overwrite = FALSE) {
  # Check if file exists and handle overwrite
  if (file.exists(output_file) && !overwrite) {
    stop("File already exists and overwrite = FALSE")
  }
  
  # If file exists and overwrite is TRUE, remove the existing file
  if (file.exists(output_file) && overwrite) {
    file.remove(output_file)
  }
  
  # Get dimensions and CRS from first raster
  r <- soil_list[[1]]
  nx <- ncol(r)
  ny <- nrow(r)
  crs_string <- crs(r)
  
  # Create longitude and latitude vectors
  ext <- ext(r)
  lon <- seq(from = ext[1], to = ext[2], length.out = nx)
  lat <- seq(from = ext[3], to = ext[4], length.out = ny)  # Changed back to ascending order
  
  # Define dimensions
  londim <- ncdim_def("longitude", "degrees_east", lon)
  latdim <- ncdim_def("latitude", "degrees_north", lat)
  
  # Define units for each soil property
  units_list <- list(
    Prediction = "cm",
    Uncertainty  = "cm"
  )
   
  # Create list of variables with appropriate units
  var_list <- list()
  for (var_name in names(soil_list)) {
    var_list[[var_name]] <- ncvar_def(
      name = var_name,
      units = units_list[[var_name]],
      dim = list(londim, latdim),
      missval = NA
    )
  }
  
  # Create netCDF file
  ncout <- nc_create(output_file, var_list)
  
  # Write data
  for (var_name in names(soil_list)) {
    # Convert SpatRaster to matrix and handle orientation
    values_matrix <- t(as.matrix(soil_list[[var_name]], wide=TRUE))
    values_matrix <- values_matrix[,ncol(values_matrix):1]
    ncvar_put(ncout, var_list[[var_name]], vals = values_matrix)
  }
  
  # Add global attributes
  ncatt_put(ncout, 0, "title", "Soil depth until R horizon ")
  ncatt_put(ncout,0,"institution","Tuebingen University, CRC1070 ResourceCultures")
  ncatt_put(ncout, 0, "description", "Soil depth until R horizon for the top 100 cm in the Northern Kurdsitan region of Irak")
  ncatt_put(ncout,0,"author", "Mathias Bellat PhD. candidate at Tuebingen University (mathias.bellat@uni-tuebingen.de)")
  ncatt_put(ncout, 0, "creation_date", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  
  # Add CRS information
  if (!is.na(crs_string)) {
    ncatt_put(ncout, 0, "crs", crs_string)
    ncatt_put(ncout, 0, "spatial_ref", crs_string)
    
    # Add standard CF grid mapping attributes
    ncatt_put(ncout, 0, "Conventions", "CF-1.6")
    ncatt_put(ncout, "longitude", "standard_name", "longitude")
    ncatt_put(ncout, "longitude", "axis", "X")
    ncatt_put(ncout, "latitude", "standard_name", "latitude")
    ncatt_put(ncout, "latitude", "axis", "Y")
  }
  
  # Add variable descriptions
  var_descriptions <- list(
    Prediction = "Soil depth until R horizon up to 100 cm",
    Uncertainty = "Uncertainty of the Quantile Random Forest of soil depth prediction based on 0.05 - 0.95 interval"
  )
  
  # Add variable-specific attributes
  for (var_name in names(soil_list)) {
    ncatt_put(ncout, var_list[[var_name]], "long_name", var_descriptions[[var_name]])
  }
  
  # Close the file
  nc_close(ncout)
}

soil_to_netcdf(CDF_df, "./export/Soil_depth_prediction_map.nc", overwrite = TRUE)



# 06.5 Final visualisations ====================================================
#Change parameters for uncertainy and predictions
raster_resize <- aggregate(soil_depth_map, fact=5, fun=mean)
rasterdf <- raster::as.data.frame(raster_resize, xy = TRUE)
rasterdf <- rasterdf[complete.cases(rasterdf),]
  
bounds <- st_bbox(survey)
xlim_new <- c(bounds["xmin"] - 3000, bounds["xmax"] + 3000)
ylim_new <- c(bounds["ymin"] - 3000, bounds["ymax"] + 3000)
  
palette_wes <- rev(wes_palette("Zissou1", type = "continuous"))
palette_blue <- colorRampPalette(c("lightblue", "blue", "darkblue"))(7)

gg1 <- ggplot() +
 geom_raster(data = rasterdf,
                  aes(x = x, y = y,fill = Prediction )) +
      ggtitle("Soil depth prediction map") +
      scale_fill_gradientn(colors = palette_wes,
                           name = "Soil depth (cm)") +
      annotation_scale(
        location = "br",       
        width_hint = 0.3,     
        height = unit(0.3, "cm"),
        line_col = "black",      
        text_col = "black",    
        bar_cols = c("white", "red") 
      ) +
      annotate("text", label = paste("Projection: WGS84 UTM 38N"), 
               x = Inf, y = -Inf, hjust = 1.5, vjust = -3, size = 3, color = "black") +
      annotation_north_arrow(location = "tr", which_north = "true", height = unit(1, "cm"), width = unit(0.75, "cm")) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
        axis.title = element_blank(),       
        legend.position = "right")+
      coord_equal(ratio = 1) 

ggsave("./export/Visualisation_soil_depth_uncertainty_prediction_maps.pdf", gg1, width = 20, height = 10)
ggsave("./export/Visualisation_soil_depth_uncertainty_prediction_maps.png", gg1, width = 20, height = 10)

# 06.6 Simplified version for publication =======================================
  
gg1 <- ggplot() +
      geom_raster(data = rasterdf, aes(x = x, y = y, fill = Prediction)) +
      ggtitle("Soil depth prediction map") +  
      scale_fill_gradientn(colors = palette_wes, name = "Soil depth in cm") + 
      theme_void() +  
      theme(
        legend.position = "right",
        plot.title = element_text(size = 8),
        legend.title = element_text(size = 10),  
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.7, "cm")) +
      coord_equal(ratio = 1)  
    
gg2 <- ggplot() +
  geom_raster(data = rasterdf, aes(x = x, y = y, fill = Uncertainty)) +
  ggtitle("Soil depth uncertainty prediction map") +  
  scale_fill_gradientn(colors = palette_blue, name = "Uncertainty") + 
  theme_void() +  
  theme(
    legend.position = "right",
    plot.title = element_text(size = 8),
    legend.title = element_text(size = 10),  
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.7, "cm")) +
  coord_equal(ratio = 1)  

combined_plot <- plot_grid(gg1, gg2, ncol = 2) 
 
ggsave("./export/Publication_soil_depth_map.pdf", combined_plot, width = 20, height = 10)
ggsave("./export/Publication_soil_depth_map.png", combined_plot, width = 20, height = 10)

dev.off()
