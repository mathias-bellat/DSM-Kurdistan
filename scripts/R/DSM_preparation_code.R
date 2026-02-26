####################################################################
# This script is the digital soil mapping for the                  #
# Northern Kurdistan region (Iraq) soils part 1/4                  #                       
#                                                                  #                                                   
# Author: Mathias Bellat  and Pegah Khosravani                     #
# Affiliation : Tubingen University                                #
# Creation date : 09/10/2024                                       #
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

install.packages("pacman")        
#Install and load the "pacman" package (allow easier download of packages)
library(pacman)
pacman::p_load(dplyr, tidyr,ggplot2, mapview, sf, cli, terra,  corrplot, doParallel, viridis, Boruta,  caret,
               quantregForest, readr, rpart, reshape2, usdm, soiltexture, compositions, patchwork)

# 0.3 Show session infos =======================================================

sessionInfo()

# 01 Import data sets ##########################################################

# 01.1 Import soils infos ======================================================

# From the data accessible at the https://doi.org/10.1594/PANGAEA.973700
soil_infos <- read.csv("./data/MIR_spectra_prediction.csv", sep=";")
soil_infos$Depth..bot <- as.factor(soil_infos$Depth..bot)
soil_infos <- soil_infos[,-c(5,7,8)]
soil_list <- split(soil_infos, soil_infos$Depth..bot)

depths <- c("0_10", "10_30", "30_50", "50_70", "70_100")
names(soil_list) <- depths 

soil_infos <- soil_list
for (i in 1:length(soil_infos)) {
  soil_infos[[i]] <- soil_infos[[i]][,-c(2,5)]
  colnames(soil_infos[[i]]) <- c("Site_name","Latitude","Longitude","pH","CaCO3","Nt","Ct","Corg","EC","Sand","Silt","Clay","MWD")
  write.csv(soil_infos[[i]], paste0("./data/Infos_",names(soil_infos[i]),"_soil.csv"))
}

# 01.2 Plot the soil information ===============================================

# Short overview of the values
head(soil_infos[[1]])
head(soil_infos[[2]])
head(soil_infos[[3]])
head(soil_infos[[4]])
head(soil_infos[[5]])

# Histogramm ploting with normal and sqrt values
for (i in 1:length(soil_infos)) {
  
  windows(width = 12, height = 9)
  par(mfrow = c(3, 4))
  for (j in 4:length(soil_infos[[1]])) {
    hist(sqrt(soil_infos[[i]][, j]), main = paste0("Distribution of ", names(soil_infos[[i]][j]), " for ", names(soil_infos[i]) ," soil"), 
         xlab = paste0("Transformed Square Root of ", names(soil_infos[[i]][j])), col = "skyblue", border = "white")
  }
  savePlot(paste0("./export/preprocess/Histogram_sqrt_", names(soil_infos[i]), "_soil.png"), type = "png")
  par(mfrow = c(1, 1))
  dev.off()
  
  windows(width = 12, height = 9)
  par(mfrow = c(3, 4))
  for (j in 4:length(soil_infos[[1]])) {
    hist(soil_infos[[i]][, j], main = paste0("Distribution of ", names(soil_infos[[i]][j]), " for ", names(soil_infos[i]) ," soil"), 
         xlab = paste0("Distribution of ", names(soil_infos[[i]][j])), col = "skyblue", border = "white")
  }
  savePlot(paste0("./export/preprocess/Histogram_", names(soil_infos[i]), "_soil.png"), type = "png")
  
  par(mfrow = c(1, 1))
  dev.off()
  
}  
dev.off()  

# 01.3 Set coordinates =========================================================

# Create a spatial dataframe and convert to WGS84 UTM 38 N coordinates

soil_infos_sp <- soil_infos

for (i in 1:length(soil_infos)) {
  soil_infos_sp[[i]] <- st_as_sf(soil_infos_sp[[i]], coords = c("Longitude", "Latitude"), crs = 4326)
  soil_infos_sp[[i]] <-st_transform(soil_infos_sp[[i]], crs = 32638)
}

mapview(soil_infos_sp[[1]]) + mapview(soil_infos_sp[[2]], col.regions = "red") + mapview(soil_infos_sp[[3]], col.regions = "green") +
  mapview(soil_infos_sp[[4]], col.regions = "pink") + mapview(soil_infos_sp[[5]], col.regions = "darkgrey")

# 01.4 Import covariates raster ================================================

Landsat <- list.files("./data/Landsat/", full.names = TRUE)
Landsat <- Landsat[!grepl("raw", Landsat, ignore.case = TRUE)]
Landsat <- rast(Landsat)

Sentinel <- list.files("./data/Sentinel/", full.names = TRUE)
Sentinel <- Sentinel[!grepl("raw", Sentinel, ignore.case = TRUE)]
Sentinel <- rast(Sentinel)

Terrain <- list.files("./SAGA/", pattern = "*sg-grd-z" , full.names = TRUE)
Terrain <- rast(Terrain[-c(2,4,7,8,9,11,15,25)])
names(Terrain)[names(Terrain) == "DEM_raw [no sinks] [Smoothed]"] <- "DEM"
names(Terrain)[names(Terrain) == "Landforms"] <- "Surface Landform"
names(Terrain)[names(Terrain) == "LS Factor"] <- "Total Catchment Area"
names(Terrain)[names(Terrain) == "Vertical Distance to Channel Network"] <- "Channel Network Distance"
Terrain <- Terrain[[order(names(Terrain))]]

Sentinel <- list.files("./data/Sentinel/", full.names = TRUE)
Sentinel <- Sentinel[!grepl("raw", Sentinel, ignore.case = TRUE)]
Sentinel <- rast(Sentinel)

Others <- list.files("./data/Others/", full.names = TRUE)
Others <- Others[!grepl("raw", Others, ignore.case = TRUE)]
Others <- Others[!grepl("DEM", Others, ignore.case = TRUE)]
Others <- rast(Others)

Others_names <- list.files("./data/Others/")
Others_names <- Others_names[!grepl("raw", Others_names, ignore.case = TRUE)]
Others_names <- Others_names[!grepl("DEM", Others_names, ignore.case = TRUE)]
names(Others) <- gsub("\\.tif$", "", Others_names)

Modis <- list.files("./data/MODIS/", full.names = TRUE)
Modis <- Modis[!grepl("raw", Modis, ignore.case = TRUE)]
Modis <- rast(Modis)
names(Modis) <- gsub("\\_raw$", "", names(Modis))

# RS Landsat 8
Landsat$Landsat8_NDWI <- (Landsat$Landsat8_green_2021_MedianComposite - Landsat$Landsat8_NIR_2021_MedianComposite)/(Landsat$Landsat8_green_2021_MedianComposite + Landsat$Landsat8_NIR_2021_MedianComposite)
Landsat$Landsat8_NDVI <- (Landsat$Landsat8_NIR_2021_MedianComposite - Landsat$Landsat8_red_2021_MedianComposite)/(Landsat$Landsat8_NIR_2021_MedianComposite + Landsat$Landsat8_red_2021_MedianComposite)
Landsat$EVI <- 2.5 * ((Landsat$Landsat8_NIR_2021_MedianComposite - Landsat$Landsat8_red_2021_MedianComposite)/((Landsat$Landsat8_NIR_2021_MedianComposite + 6 * Landsat$Landsat8_red_2021_MedianComposite) - (7.5*Landsat$Landsat8_blue_2021_MedianComposite) + 1))
Landsat$SAVI <- 1.5 * ((Landsat$Landsat8_NIR_2021_MedianComposite - Landsat$Landsat8_red_2021_MedianComposite) / (Landsat$Landsat8_NIR_2021_MedianComposite + Landsat$Landsat8_red_2021_MedianComposite + 0.5)) #Enhanced Vegetation Index
Landsat$TVI <- sqrt(Landsat$Landsat8_NDVI + 0.5)
Landsat$NDMI <- (Landsat$Landsat8_NIR_2021_MedianComposite - Landsat$Landsat8_SWIR1_2021_MedianComposite)/(Landsat$Landsat8_NIR_2021_MedianComposite + Landsat$Landsat8_SWIR1_2021_MedianComposite) # normilized difference moisture index
Landsat$COSRI <- Landsat$Landsat8_NDVI * ((Landsat$Landsat8_blue_2021_MedianComposite + Landsat$Landsat8_green_2021_MedianComposite)/(Landsat$Landsat8_red_2021_MedianComposite + Landsat$Landsat8_NIR_2021_MedianComposite))  # Combined Specteral Response Index
Landsat$LSWI <- (Landsat$Landsat8_NIR_2021_MedianComposite - Landsat$Landsat8_SWIR1_2021_MedianComposite) / (Landsat$Landsat8_NIR_2021_MedianComposite + Landsat$Landsat8_SWIR1_2021_MedianComposite)
Landsat$BrightnessIndex <- sqrt((Landsat$Landsat8_red_2021_MedianComposite^2) + (Landsat$Landsat8_NIR_2021_MedianComposite^2))
Landsat$ClayIndex <- Landsat$Landsat8_SWIR1_2021_MedianComposite / Landsat$Landsat8_SWIR2_2021_MedianComposite
Landsat$SalinityIndex <- (Landsat$Landsat8_SWIR1_2021_MedianComposite - Landsat$Landsat8_SWIR2_2021_MedianComposite) / (Landsat$Landsat8_SWIR1_2021_MedianComposite - Landsat$Landsat8_NIR_2021_MedianComposite)
Landsat$CarbonateIndex <- Landsat$Landsat8_red_2021_MedianComposite / Landsat$Landsat8_green_2021_MedianComposite
Landsat$GypsumIndex <- (Landsat$Landsat8_SWIR1_2021_MedianComposite - Landsat$Landsat8_SWIR2_2021_MedianComposite) / (Landsat$Landsat8_SWIR1_2021_MedianComposite + Landsat$Landsat8_SWIR2_2021_MedianComposite)

# RS Sentinel 2
Sentinel$Sentinel2_NDWI <- (Sentinel$Sentinel2_green_2021_MedianComposite - Sentinel$Sentinel2_NIR_2021_MedianComposite) / (Sentinel$Sentinel2_green_2021_MedianComposite + Sentinel$Sentinel2_NIR_2021_MedianComposite)
Sentinel$Sentinel2_NDVI <- (Sentinel$Sentinel2_NIR_2021_MedianComposite - Sentinel$Sentinel2_red_2021_MedianComposite) / (Sentinel$Sentinel2_NIR_2021_MedianComposite + Sentinel$Sentinel2_red_2021_MedianComposite)
Sentinel$EVI <- 2.5 * ((Sentinel$Sentinel2_NIR_2021_MedianComposite - Sentinel$Sentinel2_red_2021_MedianComposite) / ((Sentinel$Sentinel2_NIR_2021_MedianComposite + 6 * Sentinel$Sentinel2_red_2021_MedianComposite) - (7.5 * Sentinel$Sentinel2_blue_2021_MedianComposite) + 1))   
Sentinel$SAVI <- 1.5 * ((Sentinel$Sentinel2_NIR_2021_MedianComposite - Sentinel$Sentinel2_red_2021_MedianComposite) / (Sentinel$Sentinel2_NIR_2021_MedianComposite + Sentinel$Sentinel2_red_2021_MedianComposite + 0.5))
Sentinel$TVI <- sqrt(Sentinel$Sentinel2_NDVI + 0.5) 
Sentinel$NDMI <- (Sentinel$Sentinel2_NIR_2021_MedianComposite - Sentinel$Sentinel2_SWIR1_2021_MedianComposite) / (Sentinel$Sentinel2_NIR_2021_MedianComposite + Sentinel$Sentinel2_SWIR1_2021_MedianComposite) # normilized difference moisture index
Sentinel$COSRI <- Sentinel$Sentinel2_NDVI * ((Sentinel$Sentinel2_blue_2021_MedianComposite + Sentinel$Sentinel2_green_2021_MedianComposite)/(Sentinel$Sentinel2_red_2021_MedianComposite + Sentinel$Sentinel2_NIR_2021_MedianComposite))  # Combined Specteral Response Index
Sentinel$LSWI <- (Sentinel$Sentinel2_NIR_2021_MedianComposite - Sentinel$Sentinel2_SWIR1_2021_MedianComposite) / (Sentinel$Sentinel2_NIR_2021_MedianComposite + Sentinel$Sentinel2_SWIR1_2021_MedianComposite)
Sentinel$BrightnessIndex <- sqrt((Sentinel$Sentinel2_red_2021_MedianComposite^2) + (Sentinel$Sentinel2_NIR_2021_MedianComposite^2)) 
Sentinel$ClayIndex <- Sentinel$Sentinel2_SWIR1_2021_MedianComposite / Sentinel$Sentinel2_SWIR2_2021_MedianComposite
Sentinel$SalinityIndex <- (Sentinel$Sentinel2_SWIR1_2021_MedianComposite - Sentinel$Sentinel2_SWIR2_2021_MedianComposite) / (Sentinel$Sentinel2_SWIR1_2021_MedianComposite - Sentinel$Sentinel2_NIR_2021_MedianComposite)
Sentinel$CarbonateIndex <- Sentinel$Sentinel2_red_2021_MedianComposite / Sentinel$Sentinel2_green_2021_MedianComposite
Sentinel$GypsumIndex <- (Sentinel$Sentinel2_SWIR1_2021_MedianComposite - Sentinel$Sentinel2_SWIR2_2021_MedianComposite) / (Sentinel$Sentinel2_SWIR1_2021_MedianComposite + Sentinel$Sentinel2_SWIR2_2021_MedianComposite)

# RS MODIS
Modis$SAVI <- 1.5 * ((Modis$MODIS_NIR - Modis$MODIS_Red) / (Modis$MODIS_NIR + Modis$MODIS_Red + 0.5))
Modis$TVI <- sqrt(Modis$MODIS_NDVI + 0.5) 
Modis$BrightnessIndex <- sqrt((Modis$MODIS_Red^2) + (Modis$MODIS_NIR^2)) 

df_names <- data.frame()
for (i in 1:length(names(Terrain))) {
  c <- paste0("TE.",i)
  df_names[i,1] <- c
  df_names[i,2] <- names(Terrain)[i]
}

t <- nrow(df_names)
for (i in 1:length(names(Landsat))) {
  c <- paste0("LA.",i)
  df_names[i+t,1] <- c
  df_names[i+t,2] <- names(Landsat)[i]
}

t <- nrow(df_names)
for (i in 1:length(names(Sentinel))) {
  c <- paste0("SE.",i)
  df_names[i+t,1] <- c
  df_names[i+t,2] <- names(Sentinel)[i]
}

t <- nrow(df_names)
for (i in 1:length(names(Modis))) {
  c <- paste0("MO.",i)
  df_names[i+t,1] <- c
  df_names[i+t,2] <- names(Modis)[i]
}

t <- nrow(df_names)
for (i in 1:length(names(Others))) {
  c <- paste0("OT.",i)
  df_names[i+t,1] <- c
  df_names[i+t,2] <- names(Others)[i]
}

write.table(df_names,"./data/Covariates_names_DSM.txt")
x <- c(Terrain, Landsat, Sentinel, Modis, Others)
names(x) <- df_names[,1]


writeRaster(x, "./data/Stack_layers_DSM.tif", overwrite = TRUE)

# 01.5 Plot the covariates maps ================================================
covariates <- rast("./data/Stack_layers_DSM.tif")
reduce <- aggregate(covariates, fact=10, fun=modal)
writeRaster(reduce, "./data/Stack_layers_DSM_reduce.tif", overwrite = TRUE)
plot(reduce)

# 01.6 Extract the values ======================================================

# Extract the values of each band for the sampling location

df_cov <- soil_infos_sp

for (i in 1:length(df_cov)) { 
  df_cov[[i]] <- extract(covariates, df_cov[[i]], method='simple')
  df_cov[[i]] <- as.data.frame(df_cov[[i]])
  write.csv(df_cov[[i]], paste0("./data/df_",names(df_cov[i]),"_cov_DSM.csv"))
}

# 01.7 Export and save data ====================================================

save(df_cov, soil_infos_sp, file = "./export/save/Pre_process.RData")
rm(list = ls())

# 02 Check the data ############################################################

# 02.1 Import the data and merge ===============================================
make_subdir <- function(parent_dir, subdir_name) {
  path <- file.path(parent_dir, subdir_name)
  
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("✅  Folder created", path)
  } else {
    message("ℹ️ Existing folder : ", path)
  }
  
  return(path)
}

make_subdir("./export", "preprocess")


load(file = "./export/save/Pre_process.RData")

SoilCov <- df_cov
for (i in 1:length(SoilCov)) {
  ID <- 1:nrow(SoilCov[[i]])
  SoilCov[[i]] <- cbind(df_cov[[i]], ID, st_drop_geometry(soil_infos_sp[[i]]))
  cat("There is ", sum(is.na(SoilCov[[i]])== TRUE), "Na values in ", names(SoilCov[i])," soil list \n")
}


# 02.2 Plot and export the correlation matrix ==================================

for (i in 1:length(df_cov)) {
  pdf(paste0("./export/preprocess/Correlation_",names(df_cov[i]), ".pdf"),    # File name
      width = 40, height = 40,  # Width and height in inches
      bg = "white",          # Background color
      colormodel = "cmyk")   # Color model 
  
  
  # Correlation of the data (remove discrete data)
  corrplot(cor(df_cov[[i]][,-c(1,79:81)]),  method = "color", col = viridis(200), 
           type = "upper", 
           addCoef.col = "black", # Add coefficient of correlation
           tl.col = "black", tl.srt = 45, # Text label color and rotation
           number.cex = 0.7, # Size of the text labels
           cl.cex = 0.7, # Size of the color legend text
           cl.lim = c(-1, 1)) # Color legend limits
  
  dev.off()
  
}

# 02.3 Select with VIF correlation =============================================

vif <- df_cov
vif_plot <- df_cov

for (i in 1:length(df_cov)) {
  vif[[i]] <-vifcor(df_cov[[i]], th=0.8)
  vif_df <- as.data.frame(vif[[i]]@results)
  write.table(vif_df, paste0("./export/VIF/vif_results_",names(df_cov[i]) ,"_soil.txt"))
  
  vif_plot[[i]] <- ggplot(vif_df, aes(x = reorder(Variables, VIF), y = VIF)) +
    geom_bar(stat = "identity", fill = "lightblue") +
    coord_flip() +
    theme_minimal() +
    labs(title = paste0("VIF Values for ", names(df_cov[i]) ," soil"), x = "Variables", y = "VIF") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(paste0("./export/VIF/VIF_", names(df_cov[i]),"_soil.png"), vif_plot[[i]], width = 12, height = 8)
  ggsave(paste0("./export/VIF/VIF_", names(df_cov[i]),"_soil.pdf"), vif_plot[[i]], width = 12, height = 8)
}

# 02.4 Statistics ==============================================================
increments <- c("0_10", "10_30", "30_50", "50_70", "70_100")

make_subdir("./export", "Boruta") 

for (depth in increments) {
  make_subdir("./export/Boruta", depth)
  make_subdir("./export/RFE", depth)
  make_subdir("./export/preprocess", depth)
# Basic statistics
  print(names(SoilCov[[depth]]))
  print(nrow(SoilCov[[depth]]))
# If necessary
  SoilCov[[depth]] <- na.omit(SoilCov[[depth]])
  print(head(SoilCov[[depth]]))
  sapply(SoilCov[[depth]], class)
  print(summary(SoilCov[[depth]]))
}


# 03 Check covariates influences  ###############################################

Cov <- list()
nzv_vars <- list()

cl <- makeCluster(6)
registerDoParallel(cl)

# 03.1 Transform the soiltexture and remove nzv ================================

for (depth in increments) {
  depth_name <- gsub("_", " - ", depth)

  SoilCovMLCon <- SoilCov[[depth]][,-c(1,87:88)] # Remove ID and site name

  NumCovLayer = 85 # define number of covariate layer after hot coding
  StartTargetCov = NumCovLayer + 1 # start column after all covariates
  NumDataCol= ncol(SoilCovMLCon) # number of column in all data set

  # Remove and save NZV
  nzv <- nearZeroVar(SoilCovMLCon[,1:NumCovLayer], saveMetrics = TRUE)
  nzv_vars[[depth]] <- rownames(nzv)[nzv$nzv == TRUE]
  print(rownames(nzv)[nzv$nzv == TRUE])
  
  # Realise the PSF transformation
  texture_df <- SoilCovMLCon[,c((NumDataCol-3):(NumDataCol-1))]
  colnames(texture_df) <- c("SAND","SILT", "CLAY")
  texture_df <- TT.normalise.sum(texture_df, css.names =  c("SAND","SILT", "CLAY"))
  colnames(texture_df) <- colnames(SoilCovMLCon[,c((NumDataCol-3):(NumDataCol-1))])
  alr_df <- as.data.frame(alr(texture_df)) # Additive-log ratio with Sand/Clay and Silt/Clay (last column is taken)
  colnames(alr_df) <- c("alr.Sand", "alr.Silt")
  SoilCovMLCon <- SoilCovMLCon[,-c((NumDataCol-3):(NumDataCol-1))] 
  SoilCovMLCon <- cbind(SoilCovMLCon, alr_df, texture_df)
  NumDataCol <- (NumDataCol -1)
  Cov[[depth]] <- SoilCovMLCon
}

# 03.3 Boruta selection =========================================================

FormulaMLCon <- list()
Preprocess <- list()

for (i in 1:(ncol(Cov[[1]]) - NumCovLayer - 3)) { 
  FormulaMLCon[[i]] = as.formula(paste(names(Cov[[1]])[NumCovLayer+i]," ~ ",paste(names(Cov[[1]])[1:NumCovLayer],collapse="+")))
}

# Define traincontrol
TrainControl <- trainControl(method="repeatedcv", 10, 3, allowParallel = TRUE, savePredictions=TRUE)
TrainControlRF <- rfeControl(functions = rfFuncs, method = "repeatedcv", 10, 3, allowParallel = TRUE)
seed=1070

for (depth in increments) {
  depth_name <- gsub("_", " - ", depth)
  Boruta = list() 
  BorutaLabels = list() 
  Boruta_covariates = list()
  SoilCovMLCon <- Cov[[depth]]
  
   cli_progress_bar(
    format = "Boruta {.val {depth}} {.val {var}} {cli::pb_bar} {cli::pb_percent} [{cli::pb_current}/{cli::pb_total}] | \ ETA: {cli::pb_eta} - Time elapsed: {cli::pb_elapsed_clock}",
    total = length(FormulaMLCon), 
    clear = FALSE)
  
  # Individual plots
  for (i in 1:length(FormulaMLCon)) {
    var <- names(SoilCovMLCon)[NumCovLayer+i]
    set.seed(seed)
    Boruta[[i]] <- Boruta(FormulaMLCon[[i]], data = SoilCovMLCon, rfeControl = TrainControl)
    BorutaBank <- TentativeRoughFix(Boruta[[i]])
    
    pdf(paste0("./export/boruta/", depth,"/Boruta_",names(SoilCovMLCon)[NumCovLayer+i], "_for_",depth,"_soil.pdf"),    # File name
        width = 8, height = 8,  # Width and height in inches
        bg = "white",          # Background color
        colormodel = "cmyk")   # Color model 
    
    plot(BorutaBank, xlab = "", xaxt = "n",
         main=paste0("Feature Importance - Boruta ",names(SoilCovMLCon)[NumCovLayer+i]," for ", depth_name ," cm increment"))
    lz <- lapply(1:ncol(BorutaBank$ImpHistory),
                 function(j)BorutaBank$ImpHistory[is.finite(BorutaBank$ImpHistory[,j]),j])
    names(lz) <- c(names(SoilCovMLCon)[1:NumCovLayer],c("sh_Max","sh_Mean","sh_Min"))
    Labels <- sort(sapply(lz,median))
    axis(side = 1,las=2,labels = names(Labels),at = 1:ncol(BorutaBank$ImpHistory),
         cex.axis = 1)
    
    dev.off()  # Close the device
    
    BorutaLabels[[i]] <- sapply(lz,median)
    confirmed_features <- getSelectedAttributes(BorutaBank, withTentative = FALSE)
    Boruta_covariates[[i]] <- cbind(SoilCovMLCon[, confirmed_features], SoilCovMLCon[, NumCovLayer + i])
    colnames(Boruta_covariates[[i]]) <- c(confirmed_features,names(SoilCovMLCon)[NumCovLayer+i])
    write.csv(data.frame(Boruta_covariates[[i]]), paste0("./export/boruta/", depth,"/Boruta_results_",names(SoilCovMLCon)[NumCovLayer+i], "_for_",depth,"_soil.csv"))
    cli_progress_update()
  }
  cli_progress_done()
  
  # Combinned plot
  BorutaResultCon = data.frame();BorutaResultCon = data.frame(BorutaLabels[[1]])
  for (i in 2:length(FormulaMLCon)) { 
    BorutaResultCon[i] = data.frame(BorutaLabels[[i]])}
  BorutaResultCon    = BorutaResultCon[c(1:NumCovLayer),]
  names(BorutaResultCon) <- names(SoilCovMLCon)[c(StartTargetCov:NumDataCol)]
  
  BorutaResultConT = data.frame();BorutaResultConT = data.frame(Boruta[[1]]$finalDecision == "Confirmed")
  for (i in 2:length(FormulaMLCon)) { 
    BorutaResultConT[i] = data.frame(Boruta[[i]]$finalDecision == "Confirmed")}
  names(BorutaResultConT) <- names(SoilCovMLCon)[c(StartTargetCov:NumDataCol)]
  
  BorutaCovPlot = gather(BorutaResultCon,key,value);BorutaCovPlotT = gather(BorutaResultConT,key,value)
  BorutaCovPlot$cov = rep(row.names(BorutaResultCon), (NumDataCol-NumCovLayer))
  names(BorutaCovPlot) = c("Y","Z","X");BorutaCovPlot$Z.1 = BorutaCovPlotT$value
  BorutaCovPlot$Y = factor(BorutaCovPlot$Y);BorutaCovPlot$X = factor(BorutaCovPlot$X)
  BorutaCovPlot$Z.1 <- as.logical(BorutaCovPlot$Z.1)
  BorutaCovPlot$Y = factor(BorutaCovPlot$Y, levels = rev(unique(BorutaCovPlot$Y)))
  
  # Split into two groups
  df1 <- BorutaCovPlot[BorutaCovPlot$X %in% unique(BorutaCovPlot$X)[1:42], ]
  df2 <- BorutaCovPlot[BorutaCovPlot$X %in% unique(BorutaCovPlot$X)[(43):length(unique(BorutaCovPlot$X))], ]
  
  # First part
  p1 <- ggplot(df1, aes(x = X, y = Y)) + 
    geom_tile(aes(fill = Z, colour = Z.1), size = 1) + 
    labs(x = "Variables importance", y = "Soil properties", fill = "Importance") +
    theme_classic() +
    scale_fill_viridis_c(option = "E") +  
    scale_color_manual(values = c('#00000000', 'red')) +
    theme(axis.text.x = element_text(colour = "black", size = 10, angle = 90, hjust = 1),
          axis.text.y = element_text(colour = "black", size = 10),
          legend.position = "right") +  
    geom_text(aes(label = round(Z, 1)), cex = 3) +
    coord_flip() + 
    coord_equal()
  
  # Second part
  p2 <- ggplot(df2, aes(x = X, y = Y)) + 
    geom_tile(aes(fill = Z, colour = Z.1), size = 1, show.legend = FALSE) + 
    labs(x = "Variables importance", y = "Soil properties", fill = "Importance") +
    theme_classic() +
    scale_fill_viridis_c(option = "E") + 
    scale_color_manual(values = c('#00000000', 'red')) +
    theme(axis.text.x = element_text(colour = "black", size = 10, angle = 90, hjust = 1),
          axis.text.y = element_text(colour = "black", size = 10)) +
    geom_text(aes(label = round(Z, 1)), cex = 3) +
    coord_flip() + 
    coord_equal()
  
  # Combinne poth plot
  FigCovImpoBr <- (p1 / p2)  + 
    plot_annotation(title = paste0("Boruta selection of features for the ", depth_name, " cm depth interval"),
                    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
  
  ggsave(paste0("./export/boruta/", depth,"/Boruta_final_combinned_plot_", depth,"_soil.png"), FigCovImpoBr, width = 15, height = 8.5)
  ggsave(paste0("./export/boruta/", depth,"/Boruta_final_combinned_plot_", depth,"_soil.pdf"), FigCovImpoBr, width = 15, height = 8.5)
  plot(FigCovImpoBr) 
  
# 03.4 RFE covariate influence =================================================

  cli_progress_bar(
    format = "RFE {.val {depth}} {.val {var}} {cli::pb_bar} {cli::pb_percent} [{cli::pb_current}/{cli::pb_total}] | \ ETA: {cli::pb_eta} - Time elapsed: {cli::pb_elapsed_clock}",
    total = length(FormulaMLCon), 
    clear = FALSE)
    
  ResultRFECon <- list()
  subsets= c(seq(1,NumCovLayer,5))
  
  for (i in 1:length(FormulaMLCon)) {
    var <- names(SoilCovMLCon)[NumCovLayer+i]
    set.seed(seed)
    ResultRFECon[[i]] <- rfe(FormulaMLCon[[i]], data=SoilCovMLCon, 
                             sizes = subsets, rfeControl = TrainControlRF)
    cli_progress_update()
  }
  cli_progress_done() 
  
  RFE_covariates <- list()
  for (i in 1:length(ResultRFECon)) { 
    RFEpredictors <- predictors(ResultRFECon[[i]])
    RFE_covariates[[i]] <- cbind(SoilCovMLCon[, RFEpredictors], SoilCovMLCon[, NumCovLayer + i])
    colnames(RFE_covariates[[i]]) <- c(RFEpredictors, colnames(SoilCovMLCon[NumCovLayer+i]))
    write.table(data.frame(RFEpredictors), paste0("./export/RFE/", depth,"/RFE_results_",names(SoilCovMLCon)[NumCovLayer+i], "_for_",depth,"_soil.txt"))
    
  }
  
  PlotResultRFE =list()
  for (i in 1:length(ResultRFECon)) { 
    trellis.par.set(caretTheme())
    
    PlotResultRFE[[i]] <- plot(ResultRFECon[[i]],                            
                              type = c("g", "o"),
                              main=paste0("RFE of ",names(SoilCovMLCon)[NumCovLayer+i]," for ",depth_name, " cm increment"),
                              xlab="Optimal variables number") 
    
    pdf(paste0("./export/RFE/", depth,"/RFE_",names(SoilCovMLCon)[NumCovLayer+i], "_for_",depth,"_soil.pdf"),    # File name
        width = 12, height = 12,  # Width and height in inches
        bg = "white",          # Background color
        colormodel = "cmyk")   # Color model 
    
    plot(PlotResultRFE[[i]]) 
    dev.off()
  }
  
# 03.5 Export results ==========================================================
  
  Preprocess[[depth]] <- list(
    NZV = nzv_vars,
    Cov = Cov,
    Selected_cov_boruta = Boruta_covariates,
    Boruta_full_fig = FigCovImpoBr,
    Boruta = Boruta,
    Selected_cov_RFE = RFE_covariates,
    RFE_fig = PlotResultRFE,
    RFE = ResultRFECon
  )
}

save(Preprocess, file = paste0("./export/save/Preprocess.RData"))

# 03.6 Covariates selection table ==============================================

cov_names <- read.table("./data/Covariates_names_DSM.txt")

conv <- data.frame(source = cov_names[,1],
  target  = paste(cov_names[,1], cov_names[,2], sep = "_"),
  stringsAsFactors = FALSE)

conversion <- setNames(conv$target, conv$source)
cov_list <- list()

# For Boruta
for (i in 1:length(Preprocess[[1]]$Selected_cov_boruta)) {
  cov <- list()
  for (depth in increments) {
    cov[[depth]] <- colnames(Preprocess[[depth]]$Selected_cov_boruta[[i]][1:length(Preprocess[[depth]]$Selected_cov_boruta[[i]])-1])
  }
  cov_combinned <- do.call(c, cov)
  cov_table <- conversion[cov_combinned]
  cov_table <- table(cov_table)
  
  df_top5 <- as.data.frame(
    head(sort(cov_table, decreasing = TRUE), 5)
  )
  
  colnames(df_top5) <- c("value", "occurrence")
  cov_list[[i]] <- data.frame(variable = colnames(Preprocess[[depth]]$Selected_cov_boruta[[i]][length(Preprocess[[depth]]$Selected_cov_boruta[[i]])]),
                       num.cov = c(paste0(length(Preprocess[["0_10"]]$Selected_cov_boruta[[i]])-1, "; ", length(Preprocess[["10_30"]]$Selected_cov_boruta[[i]])-1, "; ",
                                          length(Preprocess[["30_50"]]$Selected_cov_boruta[[i]])-1, "; ", length(Preprocess[["50_70"]]$Selected_cov_boruta[[i]])-1, "; ",
                                          length(Preprocess[["70_100"]]$Selected_cov_boruta[[i]])-1)),
                       top.cov = c(paste0(df_top5[1,1], "; ", df_top5[2,1], "; ", df_top5[3,1], "; ", df_top5[4,1], "; ", df_top5[5,1], "; "))
                       )
}

cov_combinned <- do.call(rbind, cov_list)
write.table(cov_combinned, "./export/boruta/Factors_selection.txt")

# For RFE
for (i in 1:length(Preprocess[[1]]$Selected_cov_RFE)) {
  cov <- list()
  for (depth in increments) {
    cov[[depth]] <- colnames(Preprocess[[depth]]$Selected_cov_RFE[[i]][1:length(Preprocess[[depth]]$Selected_cov_RFE[[i]])-1])
  }
  cov_combinned <- do.call(c, cov)
  cov_table <- conversion[cov_combinned]
  cov_table <- table(cov_table)
  
  df_top5 <- as.data.frame(
    head(sort(cov_table, decreasing = TRUE), 5)
  )
  
  colnames(df_top5) <- c("value", "occurrence")
  cov_list[[i]] <- data.frame(variable = colnames(Preprocess[[depth]]$Selected_cov_RFE[[i]][length(Preprocess[[depth]]$Selected_cov_RFE[[i]])]),
                              num.cov = c(paste0(length(Preprocess[["0_10"]]$Selected_cov_RFE[[i]])-1, "; ", length(Preprocess[["10_30"]]$Selected_cov_RFE[[i]])-1, "; ",
                                                 length(Preprocess[["30_50"]]$Selected_cov_RFE[[i]])-1, "; ", length(Preprocess[["50_70"]]$Selected_cov_RFE[[i]])-1, "; ",
                                                 length(Preprocess[["70_100"]]$Selected_cov_RFE[[i]])-1)),
                              top.cov = c(paste0(df_top5[1,1], "; ", df_top5[2,1], "; ", df_top5[3,1], "; ", df_top5[4,1], "; ", df_top5[5,1], "; "))
  )
}

cov_combinned <- do.call(rbind, cov_list)
write.table(cov_combinned, "./export/RFE/Factors_selection.txt")
