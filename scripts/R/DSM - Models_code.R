####################################################################
# This script is the digital soil mapping for the                  #
# Northern Kurdistan region (Iraq) soils part 2/4                  #                       
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

pacman::p_load(ggplot2, caret, patchwork, quantregForest, Cubist, caretEnsemble, readr, grid, gridExtra, dplyr)

# 0.3 Show session infos =======================================================

sessionInfo()

# 04 Tuning all the models  ####################################################
# Here we decided to split every run by soil depth to have a better vision
# on the running process.
# 04.1 Prepare data ============================================================
depth <- "0_10"
load(file = paste0("./export/save/Selected_cov_",depth,".RData"))
SoilCovML <- First_depth_preprocess$Selected_cov_boruta
seed=1070

# Formula for the selected covariates with Boruta
FormulaML <- list()
for (i in 1:length(SoilCovML)) { 
  SoilCovMLBoruta <- SoilCovML[[i]]
  StartTargetCov = ncol(SoilCovMLBoruta) 
  NumCovLayer = StartTargetCov - 1 
FormulaML[[i]] <- as.formula(paste(names(SoilCovMLBoruta)[StartTargetCov]," ~ ",paste(names(SoilCovMLBoruta)[1:NumCovLayer],collapse="+")))
}

# 04.2 Split the test and train data ===========================================
# Split the data with a 80/20 ration
Train_data <- list()
Test_data <- list()

for (i in 1:length(FormulaML)) { 
  set.seed(seed)  
split <- createDataPartition(SoilCovML[[i]][,length(SoilCovML[[i]])], p = 0.8, list = FALSE, times = 1)
Train_data[[i]] <- SoilCovML[[i]][ split,]
Test_data[[i]]  <- SoilCovML[[i]][-split,]
}

# 04.3 Define train control ====================================================
set.seed(seed)
TrainControl <- trainControl(method="repeatedcv", 10, 3, allowParallel = TRUE, savePredictions=TRUE)
SVMrgrid <- expand.grid(sigma = seq(0.01, 0.5, by = 0.05), C = c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512))
Cubistgrid <- expand.grid(committees = c(1, 5, 10, 15, 20), neighbors = c(0, 1, 2, 3, 5, 7, 9))

# 04.4 Run the models ==========================================================

#rpart (CART)
FitRpartCon  = list()
start_time <- proc.time()
for (i in 1:length(FormulaML)) {
  set.seed(seed)
  FitRpartCon[[i]] <- train(FormulaML[[i]], data=Train_data[[i]], 
                            method="rpart", tuneLength = 20, metric="RMSE", trControl=TrainControl)
  print(names(Train_data[[i]])[length(Train_data[[i]])])
} 
end_time <- proc.time()
print(end_time - start_time)
print("CART done")


#Knn
FitKnnCon  = list()
start_time <- proc.time()
for (i in 1:length(FormulaML)) {
  set.seed(seed)
  FitKnnCon[[i]] <- train(FormulaML[[i]], data=Train_data[[i]], 
                          method="knn", tuneLength = 30, metric="RMSE", trControl=TrainControl)
  print(names(Train_data[[i]])[length(Train_data[[i]])])
} 
end_time <- proc.time()
print(end_time - start_time)
print("Knn done")

# SVM
FitSvrCon  = list()
start_time <- proc.time()
for (i in 1:length(FormulaML)) {
  set.seed(seed)
  FitSvrCon [[i]] <- train(FormulaML[[i]], data=Train_data[[i]], 
                           method="svmRadial", tuneGrid = SVMrgrid, metric="RMSE", trControl=TrainControl)
  print(names(Train_data[[i]])[length(Train_data[[i]])])
}
end_time <- proc.time()
print(end_time - start_time)
print("SVM done")


# Cubist
FitCubCon  = list()
start_time <- proc.time()
for (i in 1:length(FormulaML)) {
  set.seed(seed)
  FitCubCon [[i]] <- train(FormulaML[[i]], data=Train_data[[i]], 
                           method="cubist", tuneGrid = Cubistgrid, metric="RMSE", trControl=TrainControl)
  print(names(Train_data[[i]])[length(Train_data[[i]])])
} 
end_time <- proc.time()
print(end_time - start_time)
print("Cubist done")


# QRF
FitQRaFCon  = list()
start_time <- proc.time()
for (i in 1:length(FormulaML)) {
  set.seed(seed)
  FitQRaFCon [[i]] <- train(FormulaML[[i]], data=Train_data[[i]], 
                            method="qrf", tuneGrid = expand.grid(mtry = seq(1, length(Train_data[[i]]-1), by = 1)), 
                            metric="RMSE", trControl=TrainControl)
  print(names(Train_data[[i]])[length(Train_data[[i]])])
}
end_time <- proc.time()
print(end_time - start_time)
print("QRF done")


# 04.5 Combine models statistics ===============================================

# Look at the primary results of ML
ModelList = list()
for (i in 1:length(FormulaML)) {  
  ModelList[[i]] <- list(CART=FitRpartCon[[i]], Knn=FitKnnCon[[i]],SVM=FitSvrCon[[i]], Cubist=FitCubCon[[i]], QRF=FitQRaFCon[[i]])
}

# Look at the primary results of ML
ResultsModelCon = list()
for (i in 1:length(ModelConList)) {
  ResultsModelCon[[i]] <- resamples(ModelConList[[i]])
  write.csv(as.data.frame(modelCor(ResultsModelCon[[i]])), paste0("./export/models/", depth, "/Models_correlations_of_",names(Train_data[[i]])[length(Train_data[[i]])], "_for_", depth, "_soil.csv"))
  
}

# 04.6 Check the correlation between the models ================================

for (i in 1:length(ModelConList)) {
png(paste0("./export/models/", depth,"/Correlation between the models of ", names(Train_data[[i]])[length(Train_data[[i]])], " for ", depth, " soil.png"),  
    width = 1600, height = 1600)
splom(ResultsModelCon[[i]])
dev.off()

pdf(paste0("./export/models/", depth,"/Correlation between the models of ", names(Train_data[[i]])[length(Train_data[[i]])], " for ", depth, " soil.pdf"),    
width = 15, height = 15, 
bg = "white",          
colormodel = "cmyk") 
splom(ResultsModelCon[[i]])
dev.off()
}

# 05 Create an ensemble learning ###############################################

# 05.1 Run the model with an RF method =========================================
EnsembleModel <- list()
set.seed(seed)
for (i in 1:length(ModelConList)) {
EnsembleModel[[i]] <- caretStack(ModelSecondList[[i]], 
                                 method="rf", 
                                 metric="RMSE", 
                                 trControl=TrainControl, 
                                 tuneGrid = expand.grid(mtry = seq(1, length(Train_data[[i]]-1), by = 1)))
}

# 05.2 Plot models best tunning ================================================
for (i in 1:length(ModelConLisn)){ 
gg1 <- plot(ModelConList[[i]]$CART, main = "CART tuning parameters")
gg2 <- plot(ModelConList[[i]]$Knn, main = "Knn tuning parameters")
gg3 <- plot(ModelConList[[i]]$SVM, main = "SVMr tuning parameters")
gg4 <- plot(ModelConList[[i]]$Cubist, main = "Cubist tuning parameters")
gg5 <- plot(ModelConList[[i]]$QRF, main = "QRF tuning parameters")
gg6 <- plot(EnsembleModel[[i]]$ens_model, main = "Ensemble tuning parameters")

png(paste0("./export/models/", depth,"/Models_tuning_parameters_",names(Train_data[[i]])[length(Train_data[[i]])], "_for_",depth,"_soil.png"),    # File name
   width = 1900, height = 1200)
grid.arrange(arrangeGrob(gg1, gg2, gg3, gg4, gg5, gg6, nrow = 2, ncol = 3),
             top = textGrob(paste0("Models tuning parameters for ",names(Train_data[[i]])[length(Train_data[[i]])] ," at ", depth , " depth"), gp = gpar(fontsize = 16, fontface = "bold")))
dev.off()

pdf(paste0("./export/models/", depth,"/Models_tuning_parameters_",names(Train_data[[i]])[length(Train_data[[i]])], "_for_",depth,"_soil.pdf"),    # File name
    width = 19, height = 12, 
    bg = "white",          
    colormodel = "cmyk") 
grid.arrange(arrangeGrob(gg1, gg2, gg3, gg4, gg5, gg6, nrow = 2, ncol = 3),
             top = textGrob(paste0("Models tuning parameters for ",names(Train_data[[i]])[length(Train_data[[i]])] ," at ", depth , " depth"), gp = gpar(fontsize = 16, fontface = "bold")))
dev.off()
}

# 05.3 Combine all models ======================================================

# Look at the primary results of ML
ModelSecondList = list()
for (i in 1:length(FormulaML)) {  
  ModelSecondList[[i]] <- list(CART=FitRpartCon[[i]], Knn=FitKnnCon[[i]],SVM=FitSvrCon[[i]], Cubist=FitCubCon[[i]], QRF=FitQRaFCon[[i]], Ensemble = EnsembleModel[[i]]$ens_model)
}

ResultsSecondList = list()
for (i in 1:length(ModelSecondList)) {
  ResultsSecondList[[i]] <- resamples(ModelSecondList[[i]])
}

SummarySecondList = list()
for (i in 1:length(ModelSecondList)) {
  SummarySecondList[[i]] <- summary(ResultsSecondList[[i]])
}

# Scale and plot the models performances
ScalesModel <- list(x=list(relation="free"), y=list(relation="free"))
BwplotModelCon = list()
for (i in 1:length(ResultsSecondList)){ 
  BwplotModelCon[[i]] <- bwplot(ResultsSecondList[[i]], scales=ScalesModel, main = paste0("Boxplot of the different models for ",names(Train_data[[i]])[length(Train_data[[i]])], " at ", depth , " cm interval"))
  
  png(paste0("./export/models/", depth,"/Boxplot_final_model_",names(Train_data[[i]])[length(Train_data[[i]])], "_for_",depth,"_soil.png"),    
      width = 1600, height = 1300)
  plot(BwplotModelCon[[i]]) 
  dev.off()
  
  pdf(paste0("./export/models/", depth,"/Boxplot_final_models_",names(Train_data[[i]])[length(Train_data[[i]])], "_for_",depth,"_soil.pdf"),    
      width = 10, height = 8, 
      bg = "white",          
      colormodel = "cmyk") 
  plot(BwplotModelCon[[i]]) 
  dev.off()
}

# 05.4 Calculate error of all models ===========================================

# Calculate Error indices
Error1Con <- list()
for (i in 1:length(FormulaML)) { 
  Error1Con[[i]] <- NaN*seq(length(FormulaML))
  for(j in 1:(3 * length(ModelSecondList[[i]]))) { 
    Error1Con[[i]][j] <- mean(SummarySecondList[[i]]$values[[j]])
  }
}

ErrorIndex2Con <- data.frame(NaN)
for (i in 1:length(Error1Con)) { 
  ErrorIndexCon <- data.frame(matrix(Error1Con[[i]], nrow = length(ModelSecondList[[1]]), ncol = 3, byrow=T))
  colnames(ErrorIndexCon) <- c(paste("MAE",names(Train_data[[i]])[length(Train_data[[i]])]),
                               paste("RMSE",names(Train_data[[i]])[length(Train_data[[i]])]),
                               paste("R2",names(Train_data[[i]])[length(Train_data[[i]])]))
  ErrorIndex2Con <- cbind(ErrorIndex2Con,ErrorIndexCon)
  rownames(ErrorIndex2Con) <- c(names(ModelSecondList[[1]]))
}
write.csv(data.frame(ErrorIndex2Con), paste0("./export/models/", depth,"/Models_results_for_",depth,"_soil.csv"))


# 05.5 Plot variable importance ================================================

ModelsPlots = list()
for (i in 1:length(ModelSecondList)) {
  AllVarImportance <- data.frame()
  
  # Cart does not have a variables influence and Ensemble does not use the same variables
  for (j in 2:5) {
    var_importance <- varImp(ModelSecondList[[i]][[j]], scale = TRUE)
    importance_df <- as.data.frame(var_importance$importance)
    importance_df$Variable <- rownames(importance_df)
    importance_df$Model <- names(ModelSecondList[[i]][j])
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
    labs(title = paste0("Top 20 covariates influence accros all models of ", names(Train_data[[i]])[length(Train_data[[i]])], " for ", depth, " soil"), 
         x = "Covariates", 
         y = "Importance") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_brewer(palette = "Set3")  
  
  ggsave(paste0("./export/models/", depth,"/Final_models_top_20_covariates_influence_of_",names(Train_data[[i]])[length(Train_data[[i]])], "_for_",depth,"_soil.png"), ModelsPlots[[i]], width = 30, height = 10)
  ggsave(paste0("./export/models/", depth,"/Final_models_top_20_covariates_influence_of_",names(Train_data[[i]])[length(Train_data[[i]])], "_for_",depth,"_soil.pdf"), ModelsPlots[[i]], width = 30, height = 10)

}

# 05.6 Export results ==========================================================
# Change the name of list regarding each soil depth: first, second, third, 
# fourth and fifth.
#===============================================================================

First_depth_models <- list(
  Cov = SoilCovML,
  Models = ModelSecondList,
  Models_plots = ModelsPlots,
  Boxplot_models = BwplotModelCon, 
  Train_data = Train_data,
  Test_data = Test_data,
  Results = ErrorIndex2Con
)

save(First_depth_models, file = paste0("./export/save/Models_",depth,"_DSM.RData"))
rm(list= ls())
