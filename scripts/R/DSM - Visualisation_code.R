####################################################################
# This script is the digital soil mapping for the                  #
# Northern Kurdistan region (Iraq) soils part 4/4                  #                       
#                                                                  #                                                   
# Author: Mathias Bellat  and Pegah Khosravani                     #
# Affiliation : Tubingen University                                #
# Creation date : 09/10/2024                                       #
# E-mail: mathias.bellat@uni-tubingen.de                           #
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

# Load packages
install.packages("pacman")
library(pacman) #Easier way of loading packages
pacman::p_load(raster, terra, sf, sp, viridis, ggplot2, remotes, RColorBrewer, wesanderson, grid, 
gridExtra, colorspace, mapview, biscale, ncdf4,SpaDES, ggspatial,soiltexture, cowplot, hrbrthemes) # Specify required packages and download it if needed

devtools::install_github("ducciorocchini/cblindplot")
library(cblindplot)

# 0.3 Show session infos =======================================================

sessionInfo()

# 08 Prepare visualisations  ###################################################
# Here we decided to split every run by soil depth to have a better vision
# on the running process.
# 08.1 Prepare data ============================================================

depth <- "0_10"
layers <- list.files(paste0("./export/predictions_DSM/", depth,"/"), pattern = "*tif" , full.names = TRUE)
layers <- layers[1:10] # Remove the normalised stack
raster_stack <- stack(layers)
raster_stack <- raster_stack[[c(8,1,7,4,3,5,9,10,2,6)]] # Re_order the raster position
survey <- st_read("./data/Survey_Area.gpkg", layer = "Survey_Area")

# 08.2 Normalise the texture band on 100% ======================================

tiles.sand <- splitRaster(raster_stack[[7]], nx = 5, ny = 5)
tiles.silt <- splitRaster(raster_stack[[8]], nx = 5, ny = 5) 
tiles.clay <- splitRaster(raster_stack[[9]], nx = 5, ny = 5) 

results <- list()

for (i in 1:length(tiles.sand)) {
x <- stack(tiles.sand[[i]],tiles.silt[[i]], tiles.clay[[i]])
x_df <- raster::as.data.frame(x, xy = TRUE)
texture_df <- x_df[,c(3:5)]
colnames(texture_df) <- c("SAND","SILT", "CLAY")
texture_df <- TT.normalise.sum(texture_df, css.names =  c("SAND","SILT", "CLAY"))
x_df <- cbind(x_df[,c(1:2)], texture_df)
x_normalise <- rasterFromXYZ(x_df)
results[[i]] <- x_normalise
print(paste0(i, " / ", length(tiles.silt)))
}

raster_final <- do.call(merge, results)
crs(raster_final) <- "EPSG:32638"

raster_stack <- stack(raster_stack[[1:6]], raster_stack[[10]])
raster_stack <- stack(raster_stack, raster_final)
raster_stack <- raster_stack[[c(1:6,8:10,7)]]
names(raster_stack) <- c("pH", "CaCO3", "Nt", "Ct", "Corg", "EC", "Sand", "Silt", "Clay", "MWD")

# 08.3 Export final maps =======================================================

# For GeoTiff format
crs(raster_stack) <- "EPSG:32638"
x_croped <- crop(raster_stack, survey)
x_masked <- mask(x_croped, survey)
x_repro <- projectRaster(x_masked, crs = "EPSG:4326")
x <- rast(x_repro)
writeRaster(x, paste0("./export/final_maps/", depth, "_prediction_map.tif"), overwrite=TRUE) #By default already GeoTiff


# For netCDF format
CDF_df <- lapply(1:nlayers(x_repro), function(i) {
  rast(x_repro[[i]])  # Convertir chaque couche en SpatRaster
})

names(CDF_df) <- c("pH", "CaCO3", "Nt", "Ct", "Corg", "EC", "Sand", "Silt", "Clay", "MWD")


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
    pH = "pH units",
    CaCO3 = "%",
    Nt = "%",
    Ct = "%",
    Corg = "%",
    EC = "μS/m",
    Sand = "%",
    Silt = "%",
    Clay = "%",
    MWD = "mm"
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
  ncatt_put(ncout, 0, "title", "Soil properties for 70 - 100 cm depth ")
  ncatt_put(ncout,0,"institution","Tuebingen University, CRC1070 ResourceCultures")
  ncatt_put(ncout, 0, "description", "Soil physicochemical properties in the Northern Kurdsitan region of Irak at 70 - 100 cm depth increment")
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
    pH = "Soil pH (Kcl)",
    CaCO3 = "Calcium carbonate content",
    Nt = "Total nitrogen content",
    Ct = "Total carbon content",
    Corg = "Organic carbon content",
    EC = "Electrical conductivity",
    Sand = "Sand content",
    Silt = "Silt content",
    Clay = "Clay content",
    MWD = "Mean weight diameter"
  )
  
  # Add variable-specific attributes
  for (var_name in names(soil_list)) {
    ncatt_put(ncout, var_list[[var_name]], "long_name", var_descriptions[[var_name]])
  }
  
  # Close the file
  nc_close(ncout)
}

soil_to_netcdf(CDF_df, paste0("./export/final_maps/", depth,"_prediction_map.nc"), overwrite = TRUE)

# 08.4 Visualise for colorblind ===============================================

raster_resize <- aggregate(raster_stack, fact=5, fun=mean)

wes_palette_hcl <- function(palette_name, n = 7) {
  wes_colors <- wes_palette(palette_name, n = n, type = "continuous")
  hex_colors <- as.character(wes_colors)  
  return(hex_colors)
}

palette_wes <- wes_palette_hcl("Zissou1", n = 7)
color_blind_graph <- list()

for (i in 1:nlayers(raster_resize)) {
 
  gg1 <- cblind.plot(raster_resize[[i]], cvd = palette_wes)
  gg2 <- cblind.plot(raster_resize[[i]], cvd = "deuteranopia")
  gg3 <- cblind.plot(raster_resize[[i]], cvd = "tritanopia")
  gg4 <- cblind.plot(raster_resize[[i]], cvd = "protanopia")
  
  plots_with_labels <- arrangeGrob(
    arrangeGrob(gg1, bottom = textGrob("Original", gp = gpar(fontsize = 12))),
    arrangeGrob(gg2, bottom = textGrob("Deuteranopia", gp = gpar(fontsize = 12))),
    arrangeGrob(gg3, bottom = textGrob("Tritanopia", gp = gpar(fontsize = 12))),
    arrangeGrob(gg4, bottom = textGrob("Protanopia", gp = gpar(fontsize = 12))),
    nrow = 2, ncol = 2
  )
  
color_blind_graph[[i]] <-  grid.arrange(plots_with_labels,
               top = textGrob(
                 paste0("Predictions maps of ",names(raster_stack[[i]]) ," at ", depth , " depth"), 
                 gp = gpar(fontsize = 16, fontface = "bold")))
  
  
  png(paste0("./export/visualisations/", depth,"/Prediction_map_",names(raster_resize[[i]]), "_at_",depth,"_depth.png"),    # File name
      width = 1500, height = 1300)
  
  grid.arrange(plots_with_labels,
               top = textGrob(
                 paste0("Predictions maps of ",names(raster_stack[[i]]) ," at ", depth , " depth"), 
                 gp = gpar(fontsize = 16, fontface = "bold")))
  
  dev.off()
  
  pdf(paste0("./export/visualisations/", depth,"/Prediction_map_",names(raster_resize[[i]]), "_at_",depth,"_depth.pdf"),    # File name
      width = 15, height = 13, 
      bg = "white",          
      colormodel = "cmyk") 
  
  grid.arrange(plots_with_labels,
               top = textGrob(
                 paste0("Predictions maps of ",names(raster_stack[[i]]) ," at ", depth , " depth"), 
                 gp = gpar(fontsize = 16, fontface = "bold")))
  
  dev.off()
}

# 08.5 Combined plots of each variables ========================================
  
raster_resize_croped <- crop(raster_resize, survey)
raster_resize_masked <- mask(raster_resize_croped, survey)
rasterdf <- raster::as.data.frame(raster_resize_masked, xy = TRUE)
rasterdf <- rasterdf[complete.cases(rasterdf),]

legend <- c("pH [KCl]", "CaCO3 [%]", "Nt [%]" , "Ct [%]", "Corg [%]", "EC [µS/cm]", "Sand [%]", "Silt [%]", "Clay [%]", "MWD [mm]")
  
bounds <- st_bbox(survey)
xlim_new <- c(bounds["xmin"] - 3000, bounds["xmax"] + 3000)
ylim_new <- c(bounds["ymin"] - 3000, bounds["ymax"] + 3000)

generate_raster_plot <-    function(rasterdf, value, depth, legend) {
       palette_wes <- wes_palette("Zissou1", type = "continuous")
       t <- value +2
    ggplot() +
      geom_raster(data = rasterdf,
                  aes(x = x, y = y,fill = rasterdf[[t]] )) +
      ggtitle(paste0("Prediction map of ", names(rasterdf[t]) , " for ",depth , " cm soil depth")) +
      scale_fill_gradientn(colors = palette_wes,
                           name = legend) +
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
        legend.position = "right",           
        plot.margin = margin(0, 0, 0, 0)
      )+
      coord_equal(ratio = 1) 
}
 
stack_graph <- list()
  
for (i in 1:nlayers(raster_resize_croped)) {
    stack_graph[[i]] <- generate_raster_plot(rasterdf, i, depth, legend[i])
}

  
png(paste0("./export/visualisations/", depth,"/Prediction_maps_at_",depth,"_depth.png"),    # File name
      width = 2000, height = 1700)
  
grid.arrange(grobs = stack_graph, ncol = 3, 
               top = textGrob(paste0("Prediction maps at ", depth , " cm depth"), gp = gpar(fontsize = 12, fontface = "bold")))
  
dev.off()
 
pdf(paste0("./export/visualisations/", depth,"/Prediction_maps_at_",depth,"_depth.pdf"),  
      width = 35, height = 30, 
      bg = "white",          
      colormodel = "cmyk")
  
grid.arrange(grobs = stack_graph, ncol = 3, 
               top = textGrob(paste0("Prediction maps at ", depth , " cm depth"), gp = gpar(fontsize = 12, fontface = "bold")))
  
dev.off()
  
# 08.6 Simplified version for publication =======================================
  
generate_clean_raster_plot <- function(rasterdf, i, legend, depth) {
    t <- i + 2  
    p <- ggplot() +
      geom_raster(data = rasterdf, aes(x = x, y = y, fill = rasterdf[[t]])) +
      ggtitle(paste0("Prediction map of ", colnames(rasterdf)[t])) +  
      scale_fill_gradientn(colors = palette_wes, name = legend) + 
      theme_void() +  
      theme(
        legend.position = "right",
        plot.title = element_text(size = 8),
        legend.title = element_text(size = 6),  
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"),
        plot.margin = margin(0, 0, 0, 0) 
      ) +
      coord_equal(ratio = 1)  
    
    return(p)
  }
  
light_graph <- list()
  
for (i in 1:nlayers(raster_resize_croped)) {
    light_graph[[i]] <- generate_clean_raster_plot(rasterdf, i, legend[i], depth)
}
  
png(paste0("./export/visualisations/", depth,"/Prediction_maps_publication_at_",depth,"_depth.png"),    # File name
    width = 1000, height = 800)

grid.arrange(grobs = light_graph, ncol = 3, 
             top = textGrob(paste0("Prediction maps at ", depth , " cm depth"), gp = gpar(fontsize = 10, fontface = "bold")))

dev.off()

pdf(paste0("./export/visualisations/", depth,"/Prediction_maps_publication_at_",depth,"_depth.pdf"),  
    width = 13, height = 10, 
    bg = "white",          
    colormodel = "cmyk")

grid.arrange(grobs = light_graph, ncol = 3, 
             top = textGrob(paste0("Prediction maps at ", depth , " cm depth"), gp = gpar(fontsize = 10, fontface = "bold")))

dev.off()


# 09 Evaluation with SoilGrid  #################################################
# 09.1 Prepare data ============================================================
rm(list = ls(all.names = TRUE))
# For 0-5 cm increment
files <- list.files("./data/Soil_grid/Values/0_5/", pattern = "*tif", full.names = TRUE)
SoilGrid.zero <- stack(files)
SoilGrid.zero <- projectRaster(SoilGrid.zero, crs = "EPSG:32638")

# For 5-15 cm increment
files <- list.files("./data/Soil_grid/Values/5_15/", pattern = "*tif", full.names = TRUE)
SoilGrid.five <- stack(files)
SoilGrid.five <- projectRaster(SoilGrid.five, crs = "EPSG:32638")

# For 15-30 cm increment
files <- list.files("./data/Soil_grid/Values/15_30/", pattern = "*tif", full.names = TRUE)
SoilGrid.fifteen <- stack(files)
SoilGrid.fifteen <- projectRaster(SoilGrid.fifteen, crs = "EPSG:32638")

# For 30-60 cm increment
files <- list.files("./data/Soil_grid/Values/30_60/", pattern = "*tif", full.names = TRUE)
SoilGrid.thirty <- stack(files)
SoilGrid.thirty <- projectRaster(SoilGrid.thirty, crs = "EPSG:32638")

# For 60-100 cm increment
files <- list.files("./data/Soil_grid/Values/60_100/", pattern = "*tif", full.names = TRUE)
SoilGrid.sixty <- stack(files)
SoilGrid.sixty <- projectRaster(SoilGrid.sixty, crs = "EPSG:32638")

# 09.2 Top soil preparation ====================================================

Prediction.zero <- stack("./export/final_maps/0_10_prediction_map.tif")
Prediction.zero<- projectRaster(Prediction.zero, crs = "EPSG:32638")
Prediction.ten <- stack("./export/final_maps/10_30_prediction_map.tif")
Prediction.ten <- projectRaster(Prediction.ten, crs = "EPSG:32638")

# Resize and resample
SoilGrid.zero_crop <- crop(SoilGrid.zero, Prediction.zero$pH)
SoilGrid.five_crop <- crop(SoilGrid.five, Prediction.zero$pH)
SoilGrid.fifteen_crop <- crop(SoilGrid.fifteen, Prediction.zero$pH)

Prediction.zero_resample <- resample(Prediction.zero, SoilGrid.zero_crop, method = "bilinear")
Prediction.ten_resample <- resample(Prediction.ten, SoilGrid.zero_crop, method = "bilinear")

# Convert into DF
Prediction.zero_df <- raster::as.data.frame(Prediction.zero_resample, xy = TRUE)
Prediction.zero_df <- Prediction.zero_df[complete.cases(Prediction.zero_df),]

Prediction.ten_df <- raster::as.data.frame(Prediction.ten_resample, xy = TRUE)
Prediction.ten_df <- Prediction.ten_df[complete.cases(Prediction.ten_df),]

SoilGrid.zero_df <- raster::as.data.frame(SoilGrid.zero_crop, xy = TRUE)
SoilGrid.zero_df <- SoilGrid.zero_df[complete.cases(SoilGrid.zero_df),]

SoilGrid.five_df <- raster::as.data.frame(SoilGrid.five_crop, xy = TRUE)
SoilGrid.five_df <- SoilGrid.five_df[complete.cases(SoilGrid.five_df),]

SoilGrid.fifteen_df <- raster::as.data.frame(SoilGrid.fifteen_crop, xy = TRUE)
SoilGrid.fifteen_df<- SoilGrid.fifteen_df[complete.cases(SoilGrid.fifteen_df),]

SoilGrid_top_soil <- SoilGrid.zero_df
for (i in 3: length(SoilGrid_top_soil)) {
  SoilGrid_top_soil[i] <- ((SoilGrid.zero_df[i]*5) + (SoilGrid.five_df[i]*10) + (SoilGrid.fifteen_df[i]*15))/30  
}

# Convert to % values and reduce the pH by 10 to fit our values
colnames(SoilGrid_top_soil) <- c("x", "y", "SoilGrid.Clay", "SoilGrid.Corg", "SoilGrid.Nt", "SoilGrid.pH", "SoilGrid.Sand", "SoilGrid.Silt")
SoilGrid_top_soil[,c(3,6:8)] <- SoilGrid_top_soil[,c(3,6:8)]/10

# pH standardisation Aitken and Moody (1991) R2 =  0.78 
SoilGrid_top_soil[6] <- (1.28 * SoilGrid_top_soil[6])  - 0.613

SoilGrid_top_soil[4] <- SoilGrid_top_soil[4]/1000
SoilGrid_top_soil[5] <- SoilGrid_top_soil[5]/10000

# Replace zero value and values under or over the SD from the SoilGrid
for (i in 3:8) {
SoilGrid_top_soil[[i]][SoilGrid_top_soil[[i]] == 0] <- median(SoilGrid_top_soil[[i]])
  
}
sum(SoilGrid_top_soil == 0)
summary(SoilGrid_top_soil[3:8])

replace_sd <- function(x) {
  mean_val <- mean(x, na.rm = TRUE) 
  sd_val <- sd(x, na.rm = TRUE)     
  
  x <- ifelse(x > (mean_val + sd_val), (mean_val + sd_val),  
              ifelse(x < (mean_val - sd_val), (mean_val - sd_val), x))  
  return(x)
}

SoilGrid_top_soil[, 3:8] <- as.data.frame(lapply(SoilGrid_top_soil[, 3:8], replace_sd))

Prediction_top_soil <- Prediction.zero_df
for (i in 3:length(Prediction.zero_df)) {
  Prediction_top_soil[i] <- ((Prediction.zero_df[i]*10) + (Prediction.ten_df[i]*20))/30  
}

# pH standardisation Aitken and Moody (1991) R2 = 0.8
Prediction_top_soil[3] <- (1.175*Prediction_top_soil[3]) - 0.262

# Texture standardisation Minasny and McBratney (2001) 0.063 to 0.05

Texture <- Prediction_top_soil[,c(9:11)]
colnames(Texture) <- c("SAND","SILT", "CLAY")
Texture <- TT.normalise.sum(Texture, css.names =  c("SAND","SILT", "CLAY"))

Texture <- TT.text.transf(
  tri.data = Texture, dat.css.ps.lim = c(0, 0.002, 0.063, 2),  # German system
  base.css.ps.lim = c(0, 0.002, 0.05, 2) # USDA system
)

Prediction_top_soil[,c(9:11)] <- Texture
# 09.3 Sub soil preparation ====================================================

# We decide to match 30 - 70 cm depth increment of our prediction with 30 - 60 cm SoilGrid model

Prediction.thirty <- stack("./export/final_maps/30_50_prediction_map.tif")
Prediction.thirty <- projectRaster(Prediction.thirty, crs = "EPSG:32638")
Prediction.fifty <- stack("./export/final_maps/50_70_prediction_map.tif")
Prediction.fifty <- projectRaster(Prediction.fifty, crs = "EPSG:32638")

# Resize and resample
SoilGrid.thirty_crop <- crop(SoilGrid.thirty, Prediction.thirty$pH)

Prediction.thirty_resample <- resample(Prediction.thirty, SoilGrid.thirty_crop, method = "bilinear")
Prediction.fifty_resample <- resample(Prediction.fifty, SoilGrid.thirty_crop, method = "bilinear")

# Convert into DF
Prediction.thirty_df <- raster::as.data.frame(Prediction.thirty_resample, xy = TRUE)
Prediction.thirty_df <- Prediction.thirty_df[complete.cases(Prediction.thirty_df),]

Prediction.fifty_df <- raster::as.data.frame(Prediction.fifty_resample, xy = TRUE)
Prediction.fifty_df <- Prediction.fifty_df[complete.cases(Prediction.fifty_df),]

SoilGrid.thirty_df <- raster::as.data.frame(SoilGrid.thirty_crop, xy = TRUE)
SoilGrid.thirty_df <- SoilGrid.thirty_df[complete.cases(SoilGrid.thirty_df),]

SoilGrid_sub_soil <- SoilGrid.thirty_df

# Convert to % values and reduce the pH by 10 to fit our values
colnames(SoilGrid_sub_soil) <- c("x", "y", "SoilGrid.Clay", "SoilGrid.Corg", "SoilGrid.Nt", "SoilGrid.pH", "SoilGrid.Sand", "SoilGrid.Silt")
SoilGrid_sub_soil[,c(3,6:8)] <- SoilGrid_sub_soil[,c(3,6:8)]/10

# pH standardisation Aitken and Moody (1991) R2 =  0.78 
SoilGrid_sub_soil[6] <- (1.28 * SoilGrid_sub_soil[6])  - 0.613

SoilGrid_sub_soil[4] <- SoilGrid_sub_soil[4]/1000
SoilGrid_sub_soil[5] <- SoilGrid_sub_soil[5]/10000

for (i in 3:8) {
  SoilGrid_sub_soil[[i]][SoilGrid_sub_soil[[i]] == 0] <- median(SoilGrid_sub_soil[[i]])
}
sum(SoilGrid_sub_soil == 0)
summary(SoilGrid_sub_soil[3:8])

replace_sd <- function(x) {
  mean_val <- mean(x, na.rm = TRUE) 
  sd_val <- sd(x, na.rm = TRUE)     
  
  x <- ifelse(x > (mean_val + sd_val), (mean_val + sd_val),  
              ifelse(x < (mean_val - sd_val), (mean_val - sd_val), x))  
  return(x)
}

SoilGrid_sub_soil[, 3:8] <- as.data.frame(lapply(SoilGrid_sub_soil[, 3:8], replace_sd))
summary(SoilGrid_sub_soil[3:8])

Prediction_sub_soil <- Prediction.thirty_df
for (i in 3:length(Prediction.thirty_df)) {
  Prediction_sub_soil[i] <- ((Prediction.thirty_df[i]*20) + (Prediction.fifty_df[i]*20))/40  
}

# pH standardisation Aitken and Moody (1991) R2 = 0.8
Prediction_sub_soil[3] <- (1.175*Prediction_sub_soil[3]) - 0.262

# Texture standardisation Minasny and McBratney (2001) 0.063 to 0.05
Texture <- Prediction_sub_soil[,c(9:11)]
colnames(Texture) <- c("SAND","SILT", "CLAY")
Texture <- TT.normalise.sum(Texture, css.names =  c("SAND","SILT", "CLAY"))

Texture <- TT.text.transf(
  tri.data = Texture, dat.css.ps.lim = c(0, 0.002, 0.063, 2),  # German system
  base.css.ps.lim = c(0, 0.002, 0.05, 2) # USDA system
)

Prediction_sub_soil[,c(9:11)] <- Texture

# 09.4 Lower soil preparation ==================================================
Prediction.seventy <- stack("./export/final_maps/70_100_prediction_map.tif")
Prediction.seventy <- projectRaster(Prediction.seventy, crs = "EPSG:32638")

# Resize and resample
SoilGrid.sixty_crop <- crop(SoilGrid.sixty, Prediction.seventy$pH)

Prediction.seventy_resample <- resample(Prediction.seventy, SoilGrid.sixty_crop, method = "bilinear")

# Convert into DF
Prediction.seventy_df <- raster::as.data.frame(Prediction.seventy_resample, xy = TRUE)
Prediction.seventy_df <- Prediction.seventy_df[complete.cases(Prediction.seventy_df),]

SoilGrid.sixty_df <- raster::as.data.frame(SoilGrid.sixty_crop, xy = TRUE)
SoilGrid.sixty_df <- SoilGrid.sixty_df[complete.cases(SoilGrid.sixty_df),]

SoilGrid_lower_soil <- SoilGrid.sixty_df

# Convert to % values and reduce the pH by 10 to fit our values
colnames(SoilGrid_lower_soil) <- c("x", "y", "SoilGrid.Clay", "SoilGrid.Corg", "SoilGrid.Nt", "SoilGrid.pH", "SoilGrid.Sand", "SoilGrid.Silt")
SoilGrid_lower_soil[,c(3,6:8)] <- SoilGrid_lower_soil[,c(3,6:8)]/10

# pH standardisation Aitken and Moody (1991) R2 =  0.78 
SoilGrid_lower_soil[6] <- (1.28 * SoilGrid_lower_soil[6])  - 0.613

SoilGrid_lower_soil[4] <- SoilGrid_lower_soil[4]/1000
SoilGrid_lower_soil[5] <- SoilGrid_lower_soil[5]/10000

for (i in 3:8) {
  SoilGrid_lower_soil[[i]][SoilGrid_lower_soil[[i]] == 0] <- median(SoilGrid_lower_soil[[i]])
}
sum(SoilGrid_lower_soil == 0)
summary(SoilGrid_lower_soil[3:8])

replace_sd <- function(x) {
  mean_val <- mean(x, na.rm = TRUE) 
  sd_val <- sd(x, na.rm = TRUE)     
  
  x <- ifelse(x > (mean_val + sd_val), (mean_val + sd_val),  
              ifelse(x < (mean_val - sd_val), (mean_val - sd_val), x))  
  return(x)
}

SoilGrid_lower_soil[, 3:8] <- as.data.frame(lapply(SoilGrid_lower_soil[, 3:8], replace_sd))
summary(SoilGrid_lower_soil[3:8])

Prediction_lower_soil <- Prediction.seventy_df

# pH standardisation Aitken and Moody (1991) R2 = 0.8
Prediction_lower_soil[3] <- (1.175*Prediction_lower_soil[3]) - 0.262

# Texture standardisation Minasny and McBratney (2001) 0.063 to 0.05
Texture <- Prediction_lower_soil[,c(9:11)]
colnames(Texture) <- c("SAND","SILT", "CLAY")
Texture <- TT.normalise.sum(Texture, css.names =  c("SAND","SILT", "CLAY"))

Texture <- TT.text.transf(
  tri.data = Texture, dat.css.ps.lim = c(0, 0.002, 0.063, 2),  # German system
  base.css.ps.lim = c(0, 0.002, 0.05, 2) # USDA system
)

Prediction_lower_soil[,c(9:11)] <- Texture

# 09.5 Explore the relations  ==================================================
# Here we decided to split every run by soil depth to have a better vision
# on the running process.
#===============================================================================
depth <- "lower_soil"
increment <- "lower soil"
compared_map <- merge(SoilGrid_lower_soil, Prediction_lower_soil[,c(1:3,5,7,9:11)], by =c("x", "y"))
colnames(compared_map)
compared_map <- compared_map[, c("x", "y", "SoilGrid.pH", "SoilGrid.Nt", "SoilGrid.Corg","SoilGrid.Sand",    
                                 "SoilGrid.Silt", "SoilGrid.Clay", "pH" ,"Nt", "Corg", "Sand", "Silt" ,"Clay")]
colnames(compared_map)
summary(compared_map[3:8])
summary(compared_map[9:14])

legend <- c("pH", "Nt [%]",  "Corg [%]", "Sand [%]",  "Silt [%]",  "Clay [%]")

two_raster_plot <- function(df, value1, value2, variable, increment) {
gg1 <- ggplot() +
  geom_raster(data = df,  aes(x = x, y = y, fill = value1), show.legend = TRUE) +
  scale_fill_gradientn(colors = hcl.colors(100, "Blues"), name = variable) +
  ggtitle(paste0("SoilGrid 250m map of ", variable, " at the ", increment )) +  
   theme_void() +  
   theme(
     legend.position = "right",
     plot.title = element_text(size = 8),
     legend.title = element_text(size = 6),  
     legend.text = element_text(size = 6),
     legend.key.size = unit(0.3, "cm"),
     plot.margin = margin(0, 0, 0, 0) 
   ) +
   coord_equal(ratio = 1) 


gg2 <- ggplot() +
  geom_raster(data = df,  aes( x = x, y = y, fill = value2), show.legend = TRUE) +
  scale_fill_gradientn(colors = hcl.colors(100, "Blues"), name = variable) +
  ggtitle(paste0("Predicted map of ", variable, " at the ", increment)) +  
  theme_void() +  
  theme(
    legend.position = "right",
    plot.title = element_text(size = 8),
    legend.title = element_text(size = 6),  
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.3, "cm"),
    plot.margin = margin(0, 0, 0, 0) 
  ) +
  coord_equal(ratio = 1) 

two_map <- plot_grid(gg1, gg2, ncol = 2)
return(two_map) 
}

comparaison_maps <- list()

for (i in 3:8) {
  z <- i - 2
  t <- i + 6
  variable <- legend[z]
  value1 <- compared_map[[i]]
  value2 <- compared_map[[t]]
  comparaison_maps[[paste0(variable)]] <-  two_raster_plot(compared_map, value1, value2, variable, increment)
 
ggsave(paste0("./export/visualisations/",depth,"/Two_maps_",names(compared_map[t]),"_",depth,".pdf"),comparaison_maps[[paste0(variable)]], width = 20, height = 10)
ggsave(paste0("./export/visualisations/",depth,"/Two_maps_",names(compared_map[t]),"_",depth,".png"),comparaison_maps[[paste0(variable)]], width = 20, height = 10)
  
}

# It is not possible to automatise selection of column with 'bi_class'
data.pH <- bi_class(compared_map, x = pH, y = SoilGrid.pH , style = "quantile", dim = 3)
data.Nt <- bi_class(compared_map, x = Nt, y = SoilGrid.Nt , style = "quantile", dim = 3)
data.Corg <- bi_class(compared_map, x = Corg, y = SoilGrid.Corg , style = "quantile", dim = 3)
data.Sand <- bi_class(compared_map, x = Sand, y = SoilGrid.Sand , style = "quantile", dim = 3)
data.Silt <- bi_class(compared_map, x = Silt, y = SoilGrid.Silt , style = "quantile", dim = 3)
data.Clay <- bi_class(compared_map, x = Clay, y = SoilGrid.Clay , style = "quantile", dim = 3)

data <- cbind(compared_map,data.pH[15],data.Nt[15],data.Corg[15],data.Sand[15],data.Silt[15],data.Clay[15])
names(data)[15:ncol(data)] <- c("bi_class_pH", "bi_class_Nt", "bi_class_Corg", "bi_class_Sand", "bi_class_Silt", "bi_class_Clay")

comparaison_stats <- list()

for (i in 3:8) {
  t <- i + 6
  z <- i + 12
  value <- names(compared_map[t])
  
  
map <- ggplot() +
  geom_raster(data = data, mapping = aes( x = x, y = y, fill = data[[z]]), show.legend = FALSE) +
  bi_scale_fill(pal = "GrPink", dim = 3) +
  labs(title = paste0("Prediction maps vs. SoilGrid model for ", increment),
    subtitle = paste0("Bivariate map comparison for ", value)) +
  bi_theme() +
    coord_equal(ratio = 1) +  
  theme(
    plot.title = element_text(size = 12),
    plot.subtitle = element_text(size = 10),  
    axis.title = element_blank(),  
    axis.text = element_blank(),   
    plot.margin = margin(0, 0, 0, 0)
  )

legend <- bi_legend(pal = "GrPink",
                    dim = 3,
                    xlab = "Higher values from SoilGrid",
                    ylab = "Higher values from pred. map",
                    size = 7)

finalPlot <- ggdraw() +
  draw_plot(map, 0, 0, 1, 1) +
  draw_plot(legend, 0.01, 0.78, 0.2, 0.2)


# Calculate density of predicted values
dens_predicted <- density(data[[t]])
dens_predicted$y <- dens_predicted$y / sum(dens_predicted$y)

dens_SoilGrid <- density(data[[i]])
dens_SoilGrid$y <- dens_SoilGrid$y / sum(dens_SoilGrid$y)

# Check density
plot(dens_predicted, main = "Density Plot with predicted Normalization", xlab = paste0(value), ylab = "Density")
plot(dens_SoilGrid, main = "Density Plot with SoilGrid Normalization", xlab = paste0(value), ylab = "Density")


gg <- ggplot(data) +
  geom_density(aes(x = data[[t]], fill = "Predicted values"), alpha = 0.5) +  
  geom_density(aes(x = data[[i]], fill = "SoilGrid values"), alpha = 0.5) +  
  scale_fill_manual(values = c("Predicted values" = "#404080", "SoilGrid values" = "#69b3a2")) +
  labs(title = paste0("Density Plot of ", value,  " variable between predicted map and SoilGrid model"),
       x = paste0(value),
       y = "Density") +
  theme_minimal() + 
  theme(legend.title = element_blank(), 
        legend.position = "top") 

combined_plot <- plot_grid(gg, finalPlot, ncol = 2)
comparaison_stats[[paste0(value)]] <- combined_plot 

comparaison_stats[[paste0(value)]]

ggsave(paste0("./export/visualisations/",depth,"/SoilGrid_vs_prediction_for_",value,"_",depth,".pdf"),comparaison_stats[[paste0(value)]], width = 20, height = 10)
ggsave(paste0("./export/visualisations/",depth,"/SoilGrid_vs_prediction_for_",value,"_",depth,".png"),comparaison_stats[[paste0(value)]], width = 20, height = 10)
}

save(compared_map, file= paste0("./export/save/",depth,"_SoilGrid.RData"))

# 09.6 Plot variation of values for both maps  =================================
load("./export/save/top_soil_SoilGrid.RData")
top.soil <- compared_map
load("./export/save/sub_soil_SoilGrid.RData")
sub.soil <- compared_map
load("./export/save/lower_soil_SoilGrid.RData")
lower.soil <- compared_map
final.map <- top.soil[,c(1:2)]


compared.list <- list(top.soil = top.soil,
                      sub.soil = sub.soil,
                      lower.soil = lower.soil)

# Repeat for each variables names. The command 'bi_class' does not work with columns number.
for (i in 1:length(compared.list)) {
  compared.list[[i]] <- bi_class(compared.list[[i]], x = Silt, y = SoilGrid.Silt , style = "quantile", dim = 3)
  split_result <- stringr::str_split_fixed(compared.list[[i]][[length(compared.list[[i]])]], "-", n = 2)
  compared.list[[i]][[paste0(length(compared.list[[i]]), "_1")]] <- as.numeric(split_result[,1])
  compared.list[[i]][[paste0(length(compared.list[[i]]) + 1, "_2")]] <- as.numeric(split_result[,2])
}

df <- top.soil[,c(1:2)]
df$predicted.Silt <- (compared.list[[1]][[length(compared.list[[1]])- 1]] + compared.list[[2]][[length(compared.list[[2]])- 1]] + compared.list[[3]][[length(compared.list[[3]])- 1]])/3
df$SoilGrid.Silt <- (compared.list[[1]][[length(compared.list[[1]])]] + compared.list[[2]][[length(compared.list[[2]])]] + compared.list[[3]][[length(compared.list[[3]])]])/3
df[,c(3:4)] <- round(df[, c(3:4)], digit =0)

# Combine columns
final.map$Silt <- paste(df[[3]], df[[4]], sep = "-")


light_graph <- list()
for (i in 3:length(final.map)) {
  value <- names(final.map[i])
  map <- ggplot() +
    geom_raster(data = final.map, mapping = aes(x = x, y = y, fill = final.map[[i]]), show.legend = FALSE) +
    bi_scale_fill(pal = "GrPink", dim = 3) +
    labs(title = paste0("Bivariate map comparison for ", value)) +
    bi_theme() +
    coord_equal(ratio = 1) +  
    theme(
      plot.title = element_text(size = 12),  
      axis.title = element_blank(),  
      axis.text = element_blank()
    )
  
  if (i == 3) {
    legend <- bi_legend(pal = "GrPink",
                        dim = 3,
                        xlab = "Higher values from SoilGrid",
                        ylab = "Higher values from pred. map",
                        size = 7)
    
    finalPlot <- ggdraw() +
      draw_plot(map, 0, 0, 1, 1) +
      draw_plot(legend, 0.01, 0.78, 0.2, 0.2)
    
    t <- i - 2 
    light_graph[[t]] <- finalPlot
  }
  if (i > 3){

t <- i - 2 
light_graph[[t]] <- map
  }  
}

png("./export/visualisations/combinned/SoilGrid_vs_prediction_values.png",    # File name
    width = 1800, height = 1200)

grid.arrange(grobs = light_graph, ncol = 3, 
             top = textGrob("Bivariate maps of predicted vs. SoilGrid values", gp = gpar(fontsize = 10, fontface = "bold")))

dev.off()

pdf("./export/visualisations/combinned/SoilGrid_vs_prediction_values.pdf",  
    width = 18, height = 12, 
    bg = "white",          
    colormodel = "cmyk")

grid.arrange(grobs = light_graph, ncol = 3, 
             top = textGrob("Bivariate maps of predicted vs. SoilGrid values", gp = gpar(fontsize = 10, fontface = "bold")))

dev.off()

rm(list = ls(all.names = TRUE))
# 10 Uncertainty  ##############################################################
# 10.1 Prepare data ============================================================

# To repeat for every soil depth
depth <- "0_10"
layers <- list.files(paste0("./export/uncertainty_DSM/", depth,"/"), pattern = "*tif" , full.names = TRUE)
raster_stack <- stack(layers)
raster_stack <- raster_stack[[c(8,1,7,4,3,5,9,10,2,6)]]
survey <- st_read("./data/Survey_Area.gpkg", layer = "Survey_Area")
names(raster_stack) <- c("pH", "CaCO3", "Nt", "Ct", "Corg", "EC", "Sand", "Silt", "Clay", "MWD")
crs(raster_stack) <- "EPSG:32638"
x <- rast(raster_stack)
x_croped <- crop(x, survey)
x_masked <- mask(x_croped, survey)
terra::writeRaster(x_masked, paste0("./export/uncertainty_DSM/", depth, "_uncertainty_map.tif"), overwrite=TRUE)

rm(list = ls(all.names = TRUE))

# 10.2 Prepare SG ==============================================================

# For 0-5 cm increment
files <- list.files("./data/Soil_grid/Uncertainty/0_5/", pattern = "*tif", full.names = TRUE)
SoilGrid.zero <- stack(files)
SoilGrid.zero <- projectRaster(SoilGrid.zero, crs = "EPSG:32638")

# For 5-15 cm increment
files <- list.files("./data/Soil_grid/Uncertainty/5_15/", pattern = "*tif", full.names = TRUE)
SoilGrid.five <- stack(files)
SoilGrid.five <- projectRaster(SoilGrid.five, crs = "EPSG:32638")

# For 15-30 cm increment
files <- list.files("./data/Soil_grid/Uncertainty/15_30/", pattern = "*tif", full.names = TRUE)
SoilGrid.fifteen <- stack(files)
SoilGrid.fifteen <- projectRaster(SoilGrid.fifteen, crs = "EPSG:32638")

# For 30-60 cm increment
files <- list.files("./data/Soil_grid/Uncertainty/30_60/", pattern = "*tif", full.names = TRUE)
SoilGrid.thirty <- stack(files)
SoilGrid.thirty <- projectRaster(SoilGrid.thirty, crs = "EPSG:32638")

# For 60-100 cm increment
files <- list.files("./data/Soil_grid/Uncertainty/60_100/", pattern = "*tif", full.names = TRUE)
SoilGrid.sixty <- stack(files)
SoilGrid.sixty <- projectRaster(SoilGrid.sixty, crs = "EPSG:32638")

# 10.3 Top soil preparation ====================================================

Prediction.zero <- stack("./export/uncertainty_DSM/0_10_uncertainty_map.tif")
crs(Prediction.zero) <- "EPSG:32638"
Prediction.ten <- stack("./export/uncertainty_DSM/10_30_uncertainty_map.tif")
crs(Prediction.ten) <- "EPSG:32638"

# Resize and resample
SoilGrid.zero_crop <- crop(SoilGrid.zero, Prediction.zero$pH)
SoilGrid.five_crop <- crop(SoilGrid.five, Prediction.zero$pH)
SoilGrid.fifteen_crop <- crop(SoilGrid.fifteen, Prediction.zero$pH)

Prediction.zero_resample <- resample(Prediction.zero, SoilGrid.zero_crop, method = "bilinear")
Prediction.ten_resample <- resample(Prediction.ten, SoilGrid.zero_crop, method = "bilinear")

# Convert into DF
Prediction.zero_df <- raster::as.data.frame(Prediction.zero_resample, xy = TRUE)
Prediction.zero_df <- Prediction.zero_df[complete.cases(Prediction.zero_df),]

Prediction.ten_df <- raster::as.data.frame(Prediction.ten_resample, xy = TRUE)
Prediction.ten_df <- Prediction.ten_df[complete.cases(Prediction.ten_df),]

SoilGrid.zero_df <- raster::as.data.frame(SoilGrid.zero_crop, xy = TRUE)
SoilGrid.zero_df <- SoilGrid.zero_df[complete.cases(SoilGrid.zero_df),]

SoilGrid.five_df <- raster::as.data.frame(SoilGrid.five_crop, xy = TRUE)
SoilGrid.five_df <- SoilGrid.five_df[complete.cases(SoilGrid.five_df),]

SoilGrid.fifteen_df <- raster::as.data.frame(SoilGrid.fifteen_crop, xy = TRUE)
SoilGrid.fifteen_df<- SoilGrid.fifteen_df[complete.cases(SoilGrid.fifteen_df),]

SoilGrid_top_soil <- SoilGrid.zero_df
for (i in 3: length(SoilGrid_top_soil)) {
  SoilGrid_top_soil[i] <- ((SoilGrid.zero_df[i]*5) + (SoilGrid.five_df[i]*10) + (SoilGrid.fifteen_df[i]*15))/30  
}

# Convert to % values and reduce the pH by 10 to fit our values
colnames(SoilGrid_top_soil) <- c("x", "y", "SoilGrid.Clay", "SoilGrid.Corg", "SoilGrid.Nt", "SoilGrid.pH", "SoilGrid.Sand", "SoilGrid.Silt")
SoilGrid_top_soil[,c(3,6:8)] <- SoilGrid_top_soil[,c(3,6:8)]/10
SoilGrid_top_soil[4] <- SoilGrid_top_soil[4]/1000
SoilGrid_top_soil[5] <- SoilGrid_top_soil[5]/10000

# Replace zero value and values under or over the SD from the SoilGrid
for (i in 3:8) {
  SoilGrid_top_soil[[i]][SoilGrid_top_soil[[i]] == 0] <- median(SoilGrid_top_soil[[i]])
  
}
sum(SoilGrid_top_soil == 0)
summary(SoilGrid_top_soil[3:8])

replace_sd <- function(x) {
  mean_val <- mean(x, na.rm = TRUE) 
  sd_val <- sd(x, na.rm = TRUE)     
  
  x <- ifelse(x > (mean_val + sd_val), (mean_val + sd_val),  
              ifelse(x < (mean_val - sd_val), (mean_val - sd_val), x))  
  return(x)
}

SoilGrid_top_soil[, 3:8] <- as.data.frame(lapply(SoilGrid_top_soil[, 3:8], replace_sd))
summary(SoilGrid_top_soil[3:8])

Prediction_top_soil <- Prediction.zero_df
for (i in 3:length(Prediction.zero_df)) {
  Prediction_top_soil[i] <- ((Prediction.zero_df[i]*10) + (Prediction.ten_df[i]*20))/30  
}

# 10.4 Sub soil preparation ====================================================

# We decide to match 30 - 70 cm depth increment of our prediction with 30 - 60 cm SoilGrid model

Prediction.thirty <- stack("./export/uncertainty_DSM/30_50_uncertainty_map.tif")
crs(Prediction.thirty) <- "EPSG:32638"
Prediction.fifty <- stack("./export/uncertainty_DSM/50_70_uncertainty_map.tif")
crs(Prediction.fifty) <- "EPSG:32638"

# Resize and resample
SoilGrid.thirty_crop <- crop(SoilGrid.thirty, Prediction.thirty$pH)

Prediction.thirty_resample <- resample(Prediction.thirty, SoilGrid.thirty_crop, method = "bilinear")
Prediction.fifty_resample <- resample(Prediction.fifty, SoilGrid.thirty_crop, method = "bilinear")

# Convert into DF
Prediction.thirty_df <- raster::as.data.frame(Prediction.thirty_resample, xy = TRUE)
Prediction.thirty_df <- Prediction.thirty_df[complete.cases(Prediction.thirty_df),]

Prediction.fifty_df <- raster::as.data.frame(Prediction.fifty_resample, xy = TRUE)
Prediction.fifty_df <- Prediction.fifty_df[complete.cases(Prediction.fifty_df),]

SoilGrid.thirty_df <- raster::as.data.frame(SoilGrid.thirty_crop, xy = TRUE)
SoilGrid.thirty_df <- SoilGrid.thirty_df[complete.cases(SoilGrid.thirty_df),]

SoilGrid_sub_soil <- SoilGrid.thirty_df

# Convert to % values and reduce the pH by 10 to fit our values
colnames(SoilGrid_sub_soil) <- c("x", "y", "SoilGrid.Clay", "SoilGrid.Corg", "SoilGrid.Nt", "SoilGrid.pH", "SoilGrid.Sand", "SoilGrid.Silt")
SoilGrid_sub_soil[,c(3,6:8)] <- SoilGrid_sub_soil[,c(3,6:8)]/10
SoilGrid_sub_soil[4] <- SoilGrid_sub_soil[4]/1000
SoilGrid_sub_soil[5] <- SoilGrid_sub_soil[5]/10000

for (i in 3:8) {
  SoilGrid_sub_soil[[i]][SoilGrid_sub_soil[[i]] == 0] <- median(SoilGrid_sub_soil[[i]])
}
sum(SoilGrid_sub_soil == 0)
summary(SoilGrid_sub_soil[3:8])

replace_sd <- function(x) {
  mean_val <- mean(x, na.rm = TRUE) 
  sd_val <- sd(x, na.rm = TRUE)     
  
  x <- ifelse(x > (mean_val + sd_val), (mean_val + sd_val),  
              ifelse(x < (mean_val - sd_val), (mean_val - sd_val), x))  
  return(x)
}

SoilGrid_sub_soil[, 3:8] <- as.data.frame(lapply(SoilGrid_sub_soil[, 3:8], replace_sd))
summary(SoilGrid_sub_soil[3:8])

Prediction_sub_soil <- Prediction.thirty_df
for (i in 3:length(Prediction.thirty_df)) {
  Prediction_sub_soil[i] <- ((Prediction.thirty_df[i]*20) + (Prediction.fifty_df[i]*20))/40  
}

# 10.5 Lower soil preparation ==================================================
Prediction.seventy <- stack("./export/uncertainty_DSM/70_100_uncertainty_map.tif")
crs(Prediction.seventy) <- "EPSG:32638"

# Resize and resample
SoilGrid.sixty_crop <- crop(SoilGrid.sixty, Prediction.seventy$pH)

Prediction.seventy_resample <- resample(Prediction.seventy, SoilGrid.sixty_crop, method = "bilinear")

# Convert into DF
Prediction.seventy_df <- raster::as.data.frame(Prediction.seventy_resample, xy = TRUE)
Prediction.seventy_df <- Prediction.seventy_df[complete.cases(Prediction.seventy_df),]

SoilGrid.sixty_df <- raster::as.data.frame(SoilGrid.sixty_crop, xy = TRUE)
SoilGrid.sixty_df <- SoilGrid.sixty_df[complete.cases(SoilGrid.sixty_df),]

SoilGrid_lower_soil <- SoilGrid.sixty_df

# Convert to % values and reduce the pH by 10 to fit our values
colnames(SoilGrid_lower_soil) <- c("x", "y", "SoilGrid.Clay", "SoilGrid.Corg", "SoilGrid.Nt", "SoilGrid.pH", "SoilGrid.Sand", "SoilGrid.Silt")
SoilGrid_lower_soil[,c(3,6:8)] <- SoilGrid_lower_soil[,c(3,6:8)]/10
SoilGrid_lower_soil[4] <- SoilGrid_lower_soil[4]/1000
SoilGrid_lower_soil[5] <- SoilGrid_lower_soil[5]/10000

for (i in 3:8) {
  SoilGrid_lower_soil[[i]][SoilGrid_lower_soil[[i]] == 0] <- median(SoilGrid_lower_soil[[i]])
}
sum(SoilGrid_lower_soil == 0)
summary(SoilGrid_lower_soil[3:8])

replace_sd <- function(x) {
  mean_val <- mean(x, na.rm = TRUE) 
  sd_val <- sd(x, na.rm = TRUE)     
  
  x <- ifelse(x > (mean_val + sd_val), (mean_val + sd_val),  
              ifelse(x < (mean_val - sd_val), (mean_val - sd_val), x))  
  return(x)
}

SoilGrid_lower_soil[, 3:8] <- as.data.frame(lapply(SoilGrid_lower_soil[, 3:8], replace_sd))
summary(SoilGrid_lower_soil[3:8])

Prediction_lower_soil <- Prediction.seventy_df

# 10.6 Ensemble model ==========================================================

depth <- "lower_soil"
load(paste0("./export/save/",depth,"_SoilGrid.RData"))
prediction_map <- rasterFromXYZ(compared_map)
crs(prediction_map) <- "EPSG:32638"

uncertainty_df <- merge(SoilGrid_lower_soil, Prediction_lower_soil[,c(1:3,5,7,9:11)], by =c("x", "y"))
colnames(uncertainty_map)
uncertainty_df <- uncertainty_df[, c("x", "y", "SoilGrid.pH", "SoilGrid.Nt", "SoilGrid.Corg","SoilGrid.Sand",    
                                 "SoilGrid.Silt", "SoilGrid.Clay", "pH" ,"Nt", "Corg", "Sand", "Silt" ,"Clay")]
uncertainty_map <- rasterFromXYZ(uncertainty_df)
crs(uncertainty_map) <- "EPSG:32638"
uncertainty_map <- resample(uncertainty_map, prediction_map, method= "bilinear")

# Calculate MAE from all residuals
SG_ERROR_ABS<-(abs(uncertainty_map$SoilGrid.pH)+ abs(uncertainty_map$SoilGrid.Nt)+ abs(uncertainty_map$SoilGrid.Corg)+ 
                       abs(uncertainty_map$SoilGrid.Sand)+ abs(uncertainty_map$SoilGrid.Silt)+ abs(uncertainty_map$SoilGrid.Clay))/6
Prediction_ERROR_ABS<-(abs(uncertainty_map$pH)+ abs(uncertainty_map$Nt)+ abs(uncertainty_map$Corg)+ 
                         abs(uncertainty_map$Sand)+ abs(uncertainty_map$Silt)+ abs(uncertainty_map$Clay))/6

residuals <-stack(Prediction_ERROR_ABS, SG_ERROR_ABS)
names(residuals)<-c("ERROR_Prediction","ERROR_SG")

ensemble <- function(predvalues, serrors, basemap = 1){
  serrors <- round(serrors, 2)
  result <- predvalues[[basemap]]
  names(result) <- 'result'
  model <- raster(predvalues[[1]])
  values(model) <- basemap
  model[is.na(result)] <- NA
  minerror <- min(stack(serrors))
  names(model) <- "model"
  names(minerror) <- "error"
  result[serrors[[2]] == minerror] <- predvalues[[2]][serrors[[2]] == minerror]
  model[serrors[[2]] == minerror] <- 2
  minerror <- mask(minerror, result)
  model <- mask(model, result)
  return(stack(result, minerror, model))
}

predictions <- list()
for (i in 1:(nlayers(prediction_map)/2)) {
  x <- stack(prediction_map[[i+6]],prediction_map[[i]])

  # Run the comparison model
  start <- Sys.time()
  predictions[[i]] <- ensemble(predvalues= x , serrors=abs(residuals))
  print(Sys.time() - start)
  
}

# 10.7 Export the evaluation of SG =============================================
r_stack <- stack(predictions[[1]][[3]], predictions[[2]][[3]],predictions[[3]][[3]],predictions[[4]][[3]],
                      predictions[[5]][[3]],predictions[[6]][[3]])

r_stack[is.na(r_stack)] <- 1 
model_raster <- calc(r_stack, fun = function(x) {
  if (any(x == 2)) {
    return(2)  
  } else {
    return(1)  
  }
})

save(model_raster, file= paste0("./export/save/",depth,"_model_selection.RData"))

rm(list = ls(all.names = TRUE))

# 10.8 Plot the models maps  ===================================================
load("./export/save/top_soil_model_selection.RData")
top.soil <- model_raster
load("./export/save/sub_soil_model_selection.RData")
sub.soil <- model_raster
load("./export/save/lower_soil_model_selection.RData")
lower.soil <- model_raster
survey <- st_read("./data/Survey_Area.gpkg", layer = "Survey_Area")

model.map <- stack(top.soil, sub.soil, lower.soil)
names(model.map) <- c("top_soil", "sub_soil", "lower_soil")
crs(model.map) <- "EPSG:32638"
x <- rast(model.map)
x_croped <- crop(x, survey)
x_masked <- mask(x_croped, survey)
x_repro <- project(x_masked, "EPSG:4326")
terra::writeRaster(x_repro, paste0("./export/visualisations/Model_accuracy.tif"), overwrite=TRUE)

model_df <- raster::as.data.frame(x_repro, xy = TRUE)
model_df <- model_df[complete.cases(model_df),]

raster_plot <- function(rasterdf, depth) {
  p <- ggplot() +
    geom_raster(data = rasterdf, aes(x = x, y = y, fill = Model)) +
    ggtitle(paste0("Models for the ", depth)) +  
    scale_fill_manual(
      values = c("Prediction model" = "grey", "SoilGrid model" = "red"),  
      name = "Best model based on residuals"  
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 12),  
      axis.title = element_blank(),  
      axis.text = element_blank()
    )  +
    coord_equal(ratio = 1) 
  
  return(p)
}

graph <- list()
depth <- c("Top soil", "Sub soil", "Lower soil")
for (i in 3:length(model_df)) {
  t <- i -2
  model_df$Model <- ifelse(model_df[[i]] <= 1, "Prediction model", "SoilGrid model")
  graph[[t]] <-raster_plot(model_df, depth[t])
}


png("./export/visualisations/combinned/Model_accuracy.png",
    width = 1800, height = 1200)

grid.arrange(grobs = graph, ncol = 3, 
             top = textGrob("Comparison of models residuals", gp = gpar(fontsize = 10, fontface = "bold")))

dev.off()

pdf("./export/visualisations/combinned/Model_accuracy.pdf",  
    width = 18, height = 12, 
    bg = "white",          
    colormodel = "cmyk")

grid.arrange(grobs = graph, ncol = 3, 
             top = textGrob("Comparison of models residuals", gp = gpar(fontsize = 10, fontface = "bold")))

dev.off()

rm(list = ls(all.names = TRUE))
