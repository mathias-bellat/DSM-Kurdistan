####################################################################
# This script is the digital soil mapping for the                  #
# Northern Kurdistan region (Iraq) soils part 1/4                  #                       
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

install.packages("pacman")        
#Install and load the "pacman" package (allow easier download of packages)
library(pacman)
pacman::p_load(dplyr, tidyr,ggplot2, mapview, sf, sp, terra, raster,  corrplot, viridis, Boruta,  caret,
               quantregForest, readr, rpart, Cubist, reshape2, usdm)

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

Landsat <- raster::stack(list.files("./data/Landsat/", full.names = TRUE))
names(Landsat)

Sentinel <- raster::stack(list.files("./data/Sentinel/", full.names = TRUE))
names(Sentinel)

Terrain <- raster::stack(list.files("./data/Terrain/", full.names = TRUE))
names(Terrain)

Others <- raster::stack(list.files("./data/Others/", full.names = TRUE))
names(Others)

Modis <- raster::stack(list.files("./data/MODIS/", full.names = TRUE))
names(Modis)

# RS Landsat 8
Landsat$EVI <- 2.5 * ((Landsat$Landsat8_NIR_2020_Median - Landsat$Landsat8_red_2020_Median)/(Landsat$Landsat8_NIR_2020_Median + (6*Landsat$Landsat8_red_2020_Median) - (7.5*Landsat$Landsat8_blue_2020_Median + 1)))
Landsat$SAVI <- ((Landsat$Landsat8_NIR_2020_Median - Landsat$Landsat8_red_2020_Median)/(Landsat$Landsat8_NIR_2020_Median + Landsat$Landsat8_red_2020_Median + 0.5)) * (1.5) #Enhanced Vegetation Index
Landsat$NDMI <- (Landsat$Landsat8_NIR_2020_Median - Landsat$Landsat8_SIR1_2020_Median)/(Landsat$Landsat8_NIR_2020_Median + Landsat$Landsat8_SIR1_2020_Median) # normilized difference moisture index
Landsat$COSRI <- ((Landsat$Landsat8_blue_2020_Median - Landsat$Landsat8_green_2020_Median)/(Landsat$Landsat8_red_2020_Median + Landsat$Landsat8_NIR_2020_Median)) * (Landsat$Landsat8_NDVI_2020_Median)  # Combined Specteral Response Index
Landsat$BrightnessIndex <- ((Landsat$Landsat8_NIR_2020_Median)^2 - (Landsat$Landsat8_red_2020_Median)^2) 
Landsat$ClayIndex <- (Landsat$Landsat8_SIR1_2020_Median / Landsat$Landsat8_SIR2_2020_Median)
Landsat$SalinityIndex <- (Landsat$Landsat8_red_2020_Median - Landsat$Landsat8_NIR_2020_Median)/(Landsat$Landsat8_green_2020_Median + Landsat$Landsat8_NIR_2020_Median)
Landsat$CarbonateIndex <- (Landsat$Landsat8_red_2020_Median / Landsat$Landsat8_green_2020_Median)
Landsat$GypsumIndex <- (Landsat$Landsat8_SIR1_2020_Median - Landsat$Landsat8_NIR_2020_Median)/(Landsat$Landsat8_SIR1_2020_Median + Landsat$Landsat8_NIR_2020_Median)

# RS Sentinel 2
Sentinel$EVI <- ((Sentinel$Sentinel2_NIR_2021_MedianComposite - Sentinel$Sentinel2_red_2021_MedianComposite)/((Sentinel$Sentinel2_NIR_2021_MedianComposite + 6 * Sentinel$Sentinel2_red_2021_MedianComposite) - (7.5 * Sentinel$Sentinel2_blue_2021_MedianComposite + 1))) * 2.5   
Sentinel$TVI <- ((((Sentinel$Sentinel2_NIR_2021_MedianComposite - Sentinel$Sentinel2_red_2021_MedianComposite)/(Sentinel$Sentinel2_NIR_2021_MedianComposite + Sentinel$Sentinel2_red_2021_MedianComposite)) + 0.5) ^ 0.5) *100 
Sentinel$SAVI <- ((Sentinel$Sentinel2_NIR_2021_MedianComposite - Sentinel$Sentinel2_red_2021_MedianComposite) * 0.5) /(Sentinel$Sentinel2_NIR_2021_MedianComposite + Sentinel$Sentinel2_red_2021_MedianComposite + 0.5)
Sentinel$LSWI <- (Sentinel$Sentinel2_NIR_2021_MedianComposite - (Sentinel$Sentinel2_SWIR1_2021_MedianComposite+Sentinel$Sentinel2_SWIR2_2021_MedianComposite))/(Sentinel$Sentinel2_NIR_2021_MedianComposite - (Sentinel$Sentinel2_SWIR1_2021_MedianComposite+Sentinel$Sentinel2_SWIR2_2021_MedianComposite))
Sentinel$BrightnessIndex <- (((Sentinel$Sentinel2_NIR_2021_MedianComposite * Sentinel$Sentinel2_red_2021_MedianComposite) + (Sentinel$Sentinel2_green_2021_MedianComposite * Sentinel$Sentinel2_green_2021_MedianComposite))^0.5) / 2 
Sentinel$ClayIndex <- (Sentinel$Sentinel2_redEdge1_2021_MedianComposite / Sentinel$Sentinel2_redEdge3_2021_MedianComposite)
Sentinel$SalinityIndex <- (Sentinel$Sentinel2_red_2021_MedianComposite - Sentinel$Sentinel2_NIR_2021_MedianComposite)/(Sentinel$Sentinel2_blue_2021_MedianComposite + Sentinel$Sentinel2_NIR_2021_MedianComposite)
Sentinel$CarbonateIndex <- (Sentinel$Sentinel2_red_2021_MedianComposite / Sentinel$Sentinel2_blue_2021_MedianComposite)
Sentinel$GypsumIndex <- (Sentinel$Sentinel2_SWIR2_2021_MedianComposite - Sentinel$Sentinel2_redEdge1_2021_MedianComposite)/(Sentinel$Sentinel2_SWIR2_2021_MedianComposite + Sentinel$Sentinel2_redEdge1_2021_MedianComposite)

# RS MODIS
Modis$SAVI <- ((Modis$MODIS_NIR_band - Modis$MODIS_Red_band)/(Modis$MODIS_NIR_band + Modis$MODIS_Red_band + 0.5)) * (1.5)
Modis$BrightnessIndex <- ((Modis$MODIS_Red_band)^2 - (Modis$MODIS_NIR_band)^2)

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
x <- raster::stack(Terrain, Landsat, Sentinel, Modis, Others)
names(x) <- df_names[,1]
x <- rast(x)

terra::writeRaster(x, "./data/Stack_layers_DSM.tif", overwrite = TRUE)
covariates <- stack("./data/Stack_layers_DSM.tif")

# 01.6 Plot the covariates maps ================================================

reduce <- aggregate(covariates, fact=10, fun=mean)

plot(reduce)

# 01.6 Extract the values ======================================================

# Extract the values of each band for the sampling location

df_cov <- soil_infos_sp

for (i in 1:length(df_cov)) { 
  df_cov[[i]] <- raster::extract(covariates, df_cov[[i]], method='simple')
  df_cov[[i]] <- as.data.frame(df_cov[[i]])
  write.csv(df_cov[[i]], paste0("./data/df_",names(df_cov[i]),"_cov_DSM.csv"))
}

# 01.7 Export and save data ====================================================

save(df_cov, soil_infos_sp, file = "./export/save/Pre_process.RData")
rm(list = ls())

# 02 Check the data ############################################################

# 02.1 Import the data and merge ===============================================

load(file = "./export/save/Pre_process.RData")

SoilCov <- df_cov
for (i in 1:length(SoilCov)) {
  ID <- 1:nrow(SoilCov[[i]])
  SoilCov[[i]] <- cbind(df_cov[[i]], ID, st_drop_geometry(soil_infos_sp[[i]]))
  cat("There is ", sum(is.na(SoilCov[[i]])== TRUE), "Na values in ", names(SoilCov[i])," soil list")
}

# 02.2 Plot and export the correlation matrix ==================================

for (i in 1:length(df_cov)) {
  pdf(paste0("./export/preprocess/Correlation_",names(df_cov[i]), ".pdf"),    # File name
      width = 40, height = 40,  # Width and height in inches
      bg = "white",          # Background color
      colormodel = "cmyk")   # Color model 
  
  
  # Correlation of the data
  corrplot(cor(df_cov[[i]]),  method = "color", col = viridis(200), 
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
# Here we decided to split every run by soil depth to have a better vision
# on the running process.
#===============================================================================

x <- "0_10"
depth <- names(SoilCov[x])

# Basic statistics
names(SoilCov[[x]])
nrow(SoilCov[[x]])
# If necessary
SoilCov[[x]] <- na.omit(SoilCov[[x]])
head(SoilCov[[x]])
sapply(SoilCov[[x]], class)
summary(SoilCov[[x]])


# 03 Check covariates influences  ###############################################

SoilCovMLCon <- SoilCov[[x]][,-c(81,82)] # Remove ID and site name

NumCovLayer = 80 # define number of covariate layer after hot coding
StartTargetCov = NumCovLayer + 1 # start column after all covariates
NumDataCol= ncol(SoilCovMLCon) # number of column in all data set

preproc <- preProcess(SoilCovMLCon[,1:NumCovLayer], method=c("range"))
SoilCovMLConTrans <- predict(preproc, SoilCovMLCon[,1:NumCovLayer])
SoilCovMLConTrans <- cbind(SoilCovMLConTrans, SoilCovMLCon[,c(StartTargetCov:NumDataCol)])

# 03.1 Develop models ========================================================
FormulaMLCon  = list()
for (i in 1:(NumDataCol - NumCovLayer)) { 
  FormulaMLCon[[i]] = as.formula(paste(names(SoilCovMLConTrans)[NumCovLayer+i]," ~ ",paste(names(SoilCovMLConTrans)[1:NumCovLayer],collapse="+")))
}

# Define traincontrol

TrainControl <- trainControl(method="repeatedcv", 10, 3, allowParallel = TRUE, savePredictions=TRUE)
seed=1070

# Train different ML algorithms

#rpart (CART)
FitRpartCon  = list()
start_time <- proc.time()
for (i in 1:length(FormulaMLCon)) {
  set.seed(seed)
  FitRpartCon[[i]] <- train(FormulaMLCon[[i]], data=SoilCovMLConTrans, 
                            method="rpart", metric="RMSE", trControl=TrainControl)
  print(names(SoilCovMLConTrans)[i+NumCovLayer])
} 
end_time <- proc.time()
print(end_time - start_time)
print("CART done")


#Knn
FitKnnCon  = list()
start_time <- proc.time()
for (i in 1:length(FormulaMLCon)) {
  set.seed(seed)
  FitKnnCon[[i]] <- train(FormulaMLCon[[i]], data=SoilCovMLConTrans, 
                          method="knn", metric="RMSE", trControl=TrainControl)
  print(names(SoilCovMLConTrans)[i+NumCovLayer])
} 
end_time <- proc.time()
print(end_time - start_time)
print("Knn done")


# SVM
FitSvrCon  = list()
start_time <- proc.time()
for (i in 1:length(FormulaMLCon)) {
  set.seed(seed)
  FitSvrCon [[i]] <- train(FormulaMLCon[[i]], data=SoilCovMLConTrans, 
                           method="svmRadial", metric="RMSE", trControl=TrainControl)
  print(names(SoilCovMLConTrans)[i+NumCovLayer])
}
end_time <- proc.time()
print(end_time - start_time)
print("SVM done")


# Cubist
FitCubCon  = list()
start_time <- proc.time()
for (i in 1:length(FormulaMLCon)) {
  set.seed(seed)
  FitCubCon [[i]] <- train(FormulaMLCon[[i]], data=SoilCovMLConTrans, 
                           method="cubist", metric="RMSE", trControl=TrainControl)
  print(names(SoilCovMLConTrans)[i+NumCovLayer])
} 
end_time <- proc.time()
print(end_time - start_time)
print("Cubist done")


# QRF
FitQRaFCon  = list()
start_time <- proc.time()
for (i in 1:length(FormulaMLCon)) {
  set.seed(seed)
  FitQRaFCon [[i]] <- train(FormulaMLCon[[i]], data=SoilCovMLConTrans, 
                            method="qrf", metric="RMSE", trControl=TrainControl)
  print(names(SoilCovMLConTrans)[i+NumCovLayer])
}
end_time <- proc.time()
print(end_time - start_time)
print("QRF done")

# 03.2 Combine models statistics ===============================================

# Look at the primary results of ML
ModelConList = list()
for (i in 1:length(FormulaMLCon)) {  
  ModelConList[[i]] <- list(CART=FitRpartCon[[i]], Knn=FitKnnCon[[i]],SVM=FitSvrCon[[i]], Cubist=FitCubCon[[i]], QRF=FitQRaFCon[[i]])
}

ResultsModelCon = list()
for (i in 1:length(ModelConList)) {
  ResultsModelCon[[i]] <- resamples(ModelConList[[i]])
}

SummaryModelCon = list()
for (i in 1:length(ResultsModelCon)) {
  SummaryModelCon[[i]] <- summary(ResultsModelCon[[i]])
}


# Scale the models 
ScalesMolel <- list(x=list(relation="free"), y=list(relation="free"))
BwplotModelCon = list()
for (i in 1:length(ResultsModelCon)){ 
  BwplotModelCon[[i]] <- bwplot(ResultsModelCon[[i]], scales=ScalesMolel, main = paste0("Comparative models of ",names(SoilCovMLCon)[NumCovLayer+i], " for ", depth, " soil"))
  
  png(paste0("./export/preprocess/", depth,"/Boxplot_first_run_model_",names(SoilCovMLConTrans)[NumCovLayer+i], "_for_",depth,"_soil.png"),    # File name
      width = 800, height = 800)
  plot(BwplotModelCon[[i]]) 
  dev.off()
}

# Calculate Error indices
Error1Con = list()
for (i in 1:length(FormulaMLCon)) { 
  Error1Con[[i]] <- NaN*seq(length(FormulaMLCon))
  for(j in 1:(3 * length(ModelConList[[i]]))) { 
    Error1Con[[i]][j] <- mean(SummaryModelCon[[i]]$values[[j]])
  }}

ErrorIndex2Con <- data.frame(NaN)
for (i in 1:length(Error1Con)) { 
  ErrorIndexCon <- data.frame(matrix(Error1Con[[i]], nrow = length(ModelConList[[1]]), ncol = 3, byrow=T))
  colnames(ErrorIndexCon) <- c(paste("MAE",names(SoilCovMLCon)[NumCovLayer+i]),
                               paste("RMSE",names(SoilCovMLCon)[NumCovLayer+i]),
                               paste("R2",names(SoilCovMLCon)[NumCovLayer+i]))
  ErrorIndex2Con <- cbind(ErrorIndex2Con,ErrorIndexCon)
  rownames(ErrorIndex2Con) <- names(ModelConList[[1]])
}
write.csv(data.frame(ErrorIndex2Con), paste0("./export/preprocess/", depth,"/First_run_models_results_for_",depth,"_soil.csv"))
  

# 03.3 Look at models covariates influences ====================================

ModelsPlots = list()
for (i in 1:length(ModelConList)) {
  AllVarImportance <- data.frame()
  
  # Cart does not have a variables influence
  for (j in 2:5) {
    var_importance <- varImp(ModelConList[[i]][[j]], scale = TRUE)
    importance_df <- as.data.frame(var_importance$importance)
    importance_df$Variable <- rownames(importance_df)
    importance_df$Model <- names(ModelConList[[i]][j])
    AllVarImportance <- rbind(AllVarImportance, importance_df)
  }
  
  AvgVarImportance <- AllVarImportance %>%
    group_by(Variable) %>%
    summarise(AvgImportance = mean(Overall, na.rm = TRUE)) %>%
    arrange(desc(AvgImportance))
  
  # Select top 20 variables
  Top20Var <- AvgVarImportance %>%
    top_n(20, wt = AvgImportance)
  
  
  AllVarImportanceTop20 <- AllVarImportance %>%
    filter(Variable %in% Top20Var$Variable)
  AllVarImportanceLong <- melt(AllVarImportanceTop20, id.vars = c("Variable", "Model"), 
                               variable.name = "Metric", value.name = "Importance")
  
  
  ModelsPlots[[i]] <- ggplot(AllVarImportanceLong, aes(x = reorder(Variable, Importance), y = Importance, fill = Model)) +
    geom_bar(stat = "identity", position = "dodge") +  
    coord_flip() +  
    labs(title = paste0("Top 20 covariates influence accros all models of ", names(SoilCovMLConTrans)[NumCovLayer+i], " for ", depth, " soil"), 
         x = "Covariates", 
         y = "Importance") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_brewer(palette = "Set3")  
  
  ggsave(paste0("./export/preprocess/", depth,"/First_run_model_top_20_covariates_influence_of_",names(SoilCovMLConTrans)[NumCovLayer+i], "_for_",depth,"_soil.png"), ModelsPlots[[i]], width = 30, height = 10)
  ggsave(paste0("./export/preprocess/", depth,"/First_run_model_top_20_covariates_influence_of_",names(SoilCovMLConTrans)[NumCovLayer+i], "_for_",depth,"_soil.pdf"), ModelsPlots[[i]], width = 30, height = 10)
  plot(ModelsPlots[[i]]) 
}

# 03.5 Boruta selction =========================================================

Boruta = list() 
BorutaLabels = list() 
Boruta_covariates = list()

# Individual plots
for (i in 1:length(FormulaMLCon)) {
  set.seed(seed)
  Boruta[[i]] <- Boruta(FormulaMLCon[[i]], data = SoilCovMLConTrans)
  BorutaBank <- TentativeRoughFix(Boruta[[i]])
  
  pdf(paste0("./export/boruta/", depth,"/Boruta_",names(SoilCovMLConTrans)[NumCovLayer+i], "_for_",depth,"_soil.pdf"),    # File name
      width = 8, height = 8,  # Width and height in inches
      bg = "white",          # Background color
      colormodel = "cmyk")   # Color model 
  
  plot(BorutaBank, xlab = "", xaxt = "n",
       main=paste0("Feature Importance - Boruta ",names(SoilCovMLConTrans)[NumCovLayer+i]," for ", depth ," cm depth"))
  lz <- lapply(1:ncol(BorutaBank$ImpHistory),
               function(j)BorutaBank$ImpHistory[is.finite(BorutaBank$ImpHistory[,j]),j])
  names(lz) <- c(names(SoilCovMLConTrans)[1:NumCovLayer],c("sh_Max","sh_Mean","sh_Min"))
  Labels <- sort(sapply(lz,median))
  axis(side = 1,las=2,labels = names(Labels),at = 1:ncol(BorutaBank$ImpHistory),
       cex.axis = 1)
  
  dev.off()  # Close the device
  
  BorutaLabels[[i]] <- sapply(lz,median)
  confirmed_features <- getSelectedAttributes(BorutaBank, withTentative = FALSE)
  Boruta_covariates[[i]] <- cbind(SoilCovMLConTrans[, confirmed_features], SoilCovMLConTrans[, NumCovLayer + i])
  colnames(Boruta_covariates[[i]]) <- c(confirmed_features,names(SoilCovMLConTrans)[NumCovLayer+i])
  write.csv(data.frame(Boruta_covariates[[i]]), paste0("./export/boruta/", depth,"/Boruta_results_",names(SoilCovMLConTrans)[NumCovLayer+i], "_for_",depth,"_soil.csv"))
  print(names(SoilCovMLConTrans)[NumCovLayer+i])
}  

# Combinned plot
BorutaResultCon = data.frame();BorutaResultCon = data.frame(BorutaLabels[[1]])
for (i in 2:length(FormulaMLCon)) { 
  BorutaResultCon[i] = data.frame(BorutaLabels[[i]])}
BorutaResultCon    = BorutaResultCon[c(1:NumCovLayer),]
names(BorutaResultCon) <- names(SoilCovMLConTrans)[c(StartTargetCov:NumDataCol)]

BorutaResultConT = data.frame();BorutaResultConT = data.frame(Boruta[[1]]$finalDecision == "Confirmed")
for (i in 2:length(FormulaMLCon)) { 
  BorutaResultConT[i] = data.frame(Boruta[[i]]$finalDecision == "Confirmed")}
names(BorutaResultConT) <- names(SoilCovMLConTrans)[c(StartTargetCov:NumDataCol)]

BorutaCovPlot = gather(BorutaResultCon,key,value);BorutaCovPlotT = gather(BorutaResultConT,key,value)
BorutaCovPlot$cov = rep(row.names(BorutaResultCon), (NumDataCol-NumCovLayer))
names(BorutaCovPlot) = c("Y","Z","X");BorutaCovPlot$Z.1 = BorutaCovPlotT$value
BorutaCovPlot$Y = factor(BorutaCovPlot$Y);BorutaCovPlot$X = factor(BorutaCovPlot$X)
BorutaCovPlot$Z.1 <- as.logical(BorutaCovPlot$Z.1)
BorutaCovPlot$Y = factor(BorutaCovPlot$Y, levels = rev(unique(BorutaCovPlot$Y)))

FigCovImpoBr = ggplot(BorutaCovPlot, aes(x = X, y = Y)) + 
  geom_tile(aes(fill = Z, colour = Z.1),  size = 1,show.legend=F) + 
  labs(title = paste0("Boruta combinned plot for ", depth ," depth") , x = "Covariates", y = "Soil properties")+
  theme_classic() +
  scale_fill_gradient(limits = c(min(BorutaCovPlot$Z), max(BorutaCovPlot$Z)),
                      low="#ffffd9", high="#081d58") +
  theme(axis.text.x = element_text(colour = "black", size=10, angle = 90, hjust = 1), 
        axis.text.y = element_text(colour = "black", size=10)) +
  geom_text(aes(label = round(Z, 1)),cex=3) +
  scale_color_manual(values = c('#00000000', 'red')) +
  coord_flip() + 
  coord_equal() 


ggsave(paste0("./export/boruta/", depth,"/Boruta_final_combinned_plot_", depth,"_soil.png"), FigCovImpoBr, width = 30, height = 10)
ggsave(paste0("./export/boruta/", depth,"/Boruta_final_combinned_plot_", depth,"_soil.pdf"), FigCovImpoBr, width = 30, height = 10)
plot(FigCovImpoBr) 


# 03.6 RFE covariate influence =================================================

TrainControlRFE <- rfeControl(functions = rfFuncs,method = "repeatedcv",
                              repeats = 3,verbose = FALSE)
ResultRFECon =list()
subsets= c(seq(1,NumCovLayer,5))
for (i in 1:length(FormulaMLCon)) {
  set.seed(seed)
  ResultRFECon[[i]] <- rfe(FormulaMLCon[[i]], data=SoilCovMLConTrans, 
                           sizes = subsets,rfeControl = TrainControlRFE)
  print(names(SoilCovMLConTrans)[i+NumCovLayer])
}

RFE_covariates <- list()
for (i in 1:length(ResultRFECon)) { 
  RFEpredictors=predictors(ResultRFECon[[i]])
  RFE_covariates[[i]] <- cbind(SoilCovMLConTrans[, RFEpredictors], SoilCovMLConTrans[, NumCovLayer + i])
  colnames(RFE_covariates[[i]]) < c(RFEpredictors,names(SoilCovMLConTrans)[NumCovLayer+i])
  write.table(data.frame(RFEpredictors), paste0("./export/RFE/", depth,"/RFE_results_",names(SoilCovMLConTrans)[NumCovLayer+i], "_for_",depth,"_soil.txt"))
  
}

PlotResultRFE =list()
for (i in 1:length(ResultRFECon)) { 
  trellis.par.set(caretTheme())
  
  PlotResultRFE[[i]] = plot(ResultRFECon[[i]],
                            type = c("g", "o"),
                            main=paste0("RFE of ",names(SoilCovMLConTrans)[NumCovLayer+i]," for ",depth, " soil"),
                            xlab="Optimal variables number")
  
  pdf(paste0("./export/RFE/", depth,"/RFE_",names(SoilCovMLConTrans)[NumCovLayer+i], "_for_",depth,"_soil.pdf"),    # File name
      width = 12, height = 12,  # Width and height in inches
      bg = "white",          # Background color
      colormodel = "cmyk")   # Color model 
  
  plot(PlotResultRFE[[i]]) 
  dev.off()
  
}

# 03.7 Export results ==========================================================
# Change the name of list regarding each soil depth: first, second, third, 
# fourth and fifth.
#===============================================================================


First_depth_preprocess <- list(
  Cov = SoilCovMLConTrans,
  Cov_original = SoilCovMLCon,
  Models = ModelConList,
  Models_plots = ModelsPlots,
  Selected_cov_boruta = Boruta_covariates,
  Selected_cov_RFE = RFE_covariates,
  Boruta_full_fig = FigCovImpoBr,
  Boruta = Boruta,
  RFE = ResultRFECon,
  RFE_fig = PlotResultRFE
)

save(First_depth_preprocess, file = paste0("./export/save/Selected_cov_", depth,".RData"))
rm(list = ls())
