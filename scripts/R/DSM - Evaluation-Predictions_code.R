####################################################################
# This script is the digital soil mapping for the                  #
# Northern Kurdistan region (Iraq) soils part 3/4                  #                       
#                                                                  #                                                   
# Author: Mathias Bellat  and Pegah Khosravani                     #
# Affiliation : Tubingen University                                #
# Creation date : 09/10/2024                                       #
# E-mail: mathias.bellat@uni-tuebingen.de                           #
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

pacman::p_load(ggplot2, DescTools, caret, grid, gridExtra, raster, readr, dplyr,terra, quantregForest, Cubist, caretEnsemble)


# 0.3 Show session infos =======================================================

sessionInfo()

# 06 Evaluation process of the models ##########################################
# Here we decided to split every run by soil depth to have a better vision
# on the running process.
# 06.1 Preparation =============================================================
depth <- "0_10"
load(paste0("./export/save/Models_",depth,"_DSM.RData"))
Evaluation <- First_depth_models


# 06.2 Run the loop for predictions ============================================

model_preds <- list()
Q1.Q3 <- list()
for (i in 1:length(Evaluation$Models)) {
  model_preds[[i]] <- list()
  X_test <- Evaluation$Test_data[[i]][,1:length(Evaluation$Test_data[[i]])-1]
  for (j in 1:5) {
    model_preds[[i]][[j]]<- as.vector(predict(Evaluation$Models[[i]][[j]], X_test))
  }
  Q1.Q3[[i]] <- predict(Evaluation$Models[[i]]$QRF$finalModel, X_test, what = c(0.05, 0.5, 0.95))
  Ensemble.pred <- as.data.frame(do.call(cbind,model_preds[[i]]))
  colnames(Ensemble.pred) <- c("CART", "Knn", "SVM", "Cubist", "QRF")
  model_preds[[i]][[6]]<- as.vector(predict(Evaluation$Models[[i]][[6]], Ensemble.pred))
}

Metrics <- list()
for (i in 1:length(model_preds)) {
empty_matrix <- matrix(NA, nrow = 6, ncol = 5) 
Final_stats <- data.frame(empty_matrix)
colnames(Final_stats) <- c("RMSE", "R2", "MAE", "CCC", "PICP")
row.names(Final_stats) <- c("CART", "Knn", "SVM", "Cubist", "QRF", "Ensemble")
  for (j in 1:length(model_preds[[i]])) {
    RMSE <- postResample(model_preds[[i]][[j]], Evaluation$Test_data[[i]][[length(Evaluation$Test_data[[i]])]])
    ccc <- CCC(model_preds[[i]][[j]], Evaluation$Test_data[[i]][[length(Evaluation$Test_data[[i]])]])
    PICP <- NA
    Final_stats[j,] <- as.data.frame(t(c(RMSE, ccc$rho.c$est, PICP)))
  }
PCIP <- (Evaluation$Test_data[[i]][[length(Evaluation$Test_data[[i]])]] >= Q1.Q3[[i]][,1]) & (Evaluation$Test_data[[i]][[length(Evaluation$Test_data[[i]])]] <= Q1.Q3[[i]][,3])  
Final_stats[[5,5]] <- mean(PCIP)*100
Metrics[[i]] <- Final_stats
}



MetricsDF <- data.frame(NaN)
for (i in 1:length(Metrics)) { 
  ErrorIndexCon <- data.frame(Metrics[[i]])
  colnames(ErrorIndexCon) <- c(paste("RMSE",names(Evaluation$Cov[[i]][length(Evaluation$Cov[[i]])])),
                               paste("R2",names(Evaluation$Cov[[i]][length(Evaluation$Cov[[i]])])),
                               paste("MAE",names(Evaluation$Cov[[i]][length(Evaluation$Cov[[i]])])),
                               paste("CCC",names(Evaluation$Cov[[i]][length(Evaluation$Cov[[i]])])),
                               paste("PICP",names(Evaluation$Cov[[i]][length(Evaluation$Cov[[i]])])))
  MetricsDF <- cbind(MetricsDF ,ErrorIndexCon)
  rownames(MetricsDF) <- rownames(Metrics[[1]])
}
write.csv(data.frame(MetricsDF), paste0("./export/evaluation/", depth,"/Models_metrics_for_",depth,"_soil.csv"))


# 06.3 Visualisation of the predictions ========================================
print(Metrics[[10]])

# Selection on the best model based on the metrics. Mainly lowest RMSE but also an appreciation of the R2 and CCC if a gap was visible
Selected_model <- data.frame(t(c("QRF","CART", "CART", "QRF", "CART", "QRF", "Knn", "QRF", "QRF",  "QRF")))
b <- c("QRF","Cubist", "QRF", "QRF", "QRF", "Cubist", "Ensemble", "SVM", "Knn",  "CART") 
c <- c("QRF","Ensemble", "Cubist", "Ensemble", "Knn", "Cubist", "Ensemble", "SVM", "QRF", "Knn")
d <- c("QRF","SVM", "QRF", "SVM", "QRF", "Ensemble", "Knn", "SVM", "QRF",  "Cubist")
e <- c("QRF","CART", "Cubist", "Ensemble", "Cubist", "SVM", "QRF", "QRF", "QRF",  "Ensemble")
Selected_model <- rbind(Selected_model, b, c, d, e)
colnames(Selected_model) <- c("pH", "CaCO3", "Nt", "Ct", "Corg", "EC", "Sand", "Silt", "Clay",  "MWD")
row.names(Selected_model) <- c("0_10", "10_30", "30_50", "50_70", "70_100")
print(Selected_model)

write.table(Selected_model, "./export/evaluation/Final_models_selected.txt", row.names = TRUE, col.names = TRUE, sep = ";")
Selected_model <- read.delim("./export/evaluation/Final_models_selected.txt", sep = ";")

plot_graph <- function(Obs, Pred, legend) {
  model <- lm(Pred  ~ Obs)
  min_value <- min(c(Obs, Pred))
  max_value <- max(c(Obs, Pred))
  
  
  ggplot(data.frame(Obs, Pred), aes(Obs, Pred)) + geom_point(color = "blue") +
    coord_fixed() +
    scale_x_continuous(limits = c(min_value, max_value)) +
    scale_y_continuous(limits = c(min_value, max_value)) +
    geom_abline(aes(slope = 1, intercept = 0, color = "myline1")) +
    geom_abline(aes(slope= model$coefficients[2],intercept = model$coefficients[1], color = "myline2")) +
    scale_colour_manual(name='Lines',labels = c("1:1 Line", "Linear Regression"),
                        values=c(myline1="black", myline2="red")) +
    ggtitle( paste(legend, " Linear regression model at ", depth, " for ", model_name, " model")) +
    xlab(paste("Observed values from", legend)) +
    ylab(paste("Predicted values for", legend)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    annotate("text", x = -Inf, y = Inf, label = round(Metrics[[i]][model_name,2], digits = 4), hjust = -0.7, vjust = 1.5) +
    annotate("text", x = -Inf, y = Inf, label = "R²", hjust = -0.5, vjust = 1.5) +
    annotate("text", x = -Inf, y = Inf, label = round(Metrics[[i]][model_name,1], digits = 4), hjust = -1.2, vjust = 3.5) +
    annotate("text", x = -Inf, y = Inf, label = "RMSE", hjust = -0.2, vjust = 3.5) +
    annotate("text", x = -Inf, y = Inf, label = round(Metrics[[i]][model_name,2], digits = 4), hjust = -1, vjust = 5.5) +
    annotate("text", x = -Inf, y = Inf, label = "MAE", hjust = -0.2, vjust = 5.5) +
    annotate("text", x = -Inf, y = Inf, label = round(Metrics[[i]][model_name,4], digits = 4), hjust = -1, vjust = 7.5) +
    annotate("text", x = -Inf, y = Inf, label = "CCC", hjust = -0.2, vjust = 7.5) 
}

regression_plot <- list()
for (i in 1:length(Evaluation$Models)) {
  names(model_preds[[i]]) <- c("CART", "Knn", "SVM", "Cubist", "QRF", "Ensemble")
  model_name <- Selected_model[depth,i]
  Obs <- model_preds[[i]][[model_name]]
  Pred <- Evaluation$Test_data[[i]][[length(Evaluation$Test_data[[i]])]]
  legend <- names(Evaluation$Cov[[i]][length(Evaluation$Cov[[i]])])
  
  p <- plot_graph(Obs, Pred, legend)
  regression_plot[[i]] <- p
}


png(paste0("./export/evaluation/", depth,"/Linear_regression_of_best_models_for_",depth,"_soil.png"),
    width = 1900, height = 1900)
grid.arrange(do.call(arrangeGrob, c(regression_plot, nrow = 4, ncol = 3)),
             top = textGrob(paste0("Linear regression for best models at ", depth , " depth"), gp = gpar(fontsize = 16, fontface = "bold")))
dev.off()

pdf(paste0("./export/evaluation/", depth,"/Linear_regression_of_best_models_for_",depth,"_soil.pdf"),
    width = 19, height = 19, 
    bg = "white",          
    colormodel = "cmyk") 
grid.arrange(do.call(arrangeGrob, c(regression_plot, nrow = 4, ncol = 3)),
             top = textGrob(paste0("Linear regression for best models at ", depth , " depth"), gp = gpar(fontsize = 16, fontface = "bold")))
dev.off()

# 06.4 Export results ==========================================================
# Change the name of list regarding each soil depth: first, second, third, 
# fourth and fifth.
#===============================================================================

First_depth_evaluation = list(
  Metrics = Metrics,
  Plots = regression_plot
)

save(First_depth_evaluation, file = paste0("./export/save/Evaluation_",depth,"_DSM.RData"))

# 07 Prediction map for DSM ####################################################

# 07.1 Normalise the values of the predictors ==================================

rm(list = ls(all.names = TRUE))
depth_list <- c("0_10", "10_30", "30_50", "50_70", "70_100")

for (t in 1:length(depth_list)) {
  depth <- depth_list[t]
  start_time <- proc.time()
  raster_stack <- stack("./data/Stack_layers_DSM.tif")
  cov <- read.delim(paste0("./data/df_", depth,"_cov_DSM.csv"), sep = ",")
  cov <- cov[2:81]
  cov[] <- lapply(cov , as.numeric)
  

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
  terra::writeRaster(stack_scaled, paste0("./export/predictions_DSM/", depth, "/Stack_raster_normalised_", depth, "_DSM.tif"), overwrite = TRUE)
  end_time <- proc.time()
  print(paste0(depth[t], ((end_time[3] - start_time[3])/60)*t, " /", ((end_time[3] - start_time[3])/60)*5))
}

# 07.2 Prepare predictions  ====================================================
# Here we decided to split every run by soil depth to have a better vision
# on the running process.
rm(list = ls(all.names = TRUE))
depth <- "0_10"
load(paste0("./export/save/Models_",depth,"_DSM.RData"))
raster_stack_normalised <- stack(paste0("./export/predictions_DSM/", depth,"/Stack_raster_normalised_",depth,"_DSM.tif"))
selected_model <- read.delim("./export/evaluation/Final_models_selected.txt", sep =";")
Prediction <- First_depth_models

# 07.3 Prediction loop  ========================================================

for (i in i:length(Prediction$Cov)){
  variable <- names(Prediction$Cov[[i]][length(Prediction$Cov[[i]])])
  model <- selected_model[depth,variable]
  if(model != "Ensemble"){
    # Create an empty raster for storing predicted values
    predicted_raster <- raster_stack_normalised[[1]]  # Use the first layer as a template
    predicted_raster <- writeStart(predicted_raster, paste0("./export/predictions_DSM/",depth,"/DSM_",variable,"_from_",model,"_for_",depth,"_soil.tif"), overwrite = TRUE)
    
    # Subset the stack, keeping only the desired layers
    block_info <- blockSize(raster_stack_normalised)
    names <- names(Prediction$Cov[[i]][1:length(Prediction$Cov[[i]])-1])
    raster_subset <- subset(raster_stack_normalised, names)
    
    # Loop through each block of the raster stack
    for (g in 1:block_info$n) {
      start_time <- proc.time()
      block <- getValuesBlock(raster_subset , row = block_info$row[g], nrows = block_info$nrows[g])
      block <- as.data.frame(block)
      
      # Process the block
      predicted_block <- predict(Prediction$Models[[i]][[model]], newdata = block)
      
      # Write the predicted block to the output raster
      predicted_raster <- writeValues(predicted_raster, predicted_block, block_info$row[g])
      end_time <- proc.time()
      cat(round((g/block_info$n)*100, 1),"% ", ((end_time[3] - start_time[3])*g)/60, " /", ((end_time[3] - start_time[3])*block_info$n/60), " min \n")
    }
    predicted_raster <- writeStop(predicted_raster)
    print(variable)
  }
  else{
    # Create a loop to produce every_models
    dir.create(paste0("./export/predictions_DSM/",depth,"/Ensemble_",variable,"/"))
    for (j in 1:5) {
      model <- names(Prediction$Models[[i]][j])
      # Create an empty raster for storing predicted values
      predicted_raster <- raster_stack_normalised[[1]]  # Use the first layer as a template
      predicted_raster <- writeStart(predicted_raster, paste0("./export/predictions_DSM/",depth,"/Ensemble_",variable,"/",model,"_DSM_of_",variable,"_for_",depth,"_soil.tif"), overwrite = TRUE)
      
      # Loop through each block of the raster stack
      block_info <- blockSize(raster_stack_normalised)
      names <- names(Prediction$Cov[[i]][1:length(Prediction$Cov[[i]])-1])
      raster_subset <- subset(raster_stack_normalised, names)
      
      # Loop through each block of the raster stack
      for (g in 1:block_info$n) {
        start_time <- proc.time()
        block <- getValuesBlock(raster_subset , row = block_info$row[g], nrows = block_info$nrows[g])
        block <- as.data.frame(block)
        # Process the block (apply your model's predictions here)
        
        predicted_block <- predict(Prediction$Models[[i]][[model]], newdata = block)
        
        # Write the predicted block to the output raster
        predicted_raster <- writeValues(predicted_raster, predicted_block, block_info$row[g])
        end_time <- proc.time()
        cat(round((g/block_info$n)*100, 1),"% ", ((end_time[3] - start_time[3])*g)/60, " /", ((end_time[3] - start_time[3])*block_info$n/60), " min \n")
        
      }
      predicted_raster <- writeStop(predicted_raster)
      cat(variable, model, "\n")
    }
    model <- "Ensemble"
    #Select the files
    all_files <- list.files(paste0("./export/predictions_DSM/",depth,"/Ensemble_",variable), full.names = TRUE)
    raster_pred_ensemble <- stack(all_files)
    
    predicted_raster <- raster_pred_ensemble[[1]]  # Use the first layer as a template
    predicted_raster <- writeStart(predicted_raster, paste0("./export/predictions_DSM/",depth,"/DSM_",variable,"_from_",model,"_for_",depth,"_soil.tif"), overwrite = TRUE)
    block_info <- blockSize(raster_pred_ensemble)
    
    # Loop through each block of the raster stack
    for (g in 1:block_info$n) {
      start_time <- proc.time()
      # Read block of raster data
      block <- getValuesBlock(raster_pred_ensemble, row = block_info$row[g], nrows = block_info$nrows[g])
      block <- as.data.frame(block)
      
      # Reorganise the order of the columns
      colnames(block) <- c("CART", "Cubist", "Knn", "QRF", "SVM")
      block <- block %>% select("CART", "Knn", "SVM", "Cubist", "QRF")
      predicted_block <- predict(Prediction$Models[[i]][[model]], newdata = block)
      
      # Write the predicted block to the output raster
      predicted_raster <- writeValues(predicted_raster, predicted_block, block_info$row[g])
      end_time <- proc.time()
      cat(round((g/block_info$n)*100, 1),"% ", ((end_time[3] - start_time[3])*g)/60, " /", ((end_time[3] - start_time[3])*block_info$n/60), " min \n")
    }
    predicted_raster <- writeStop(predicted_raster)
    cat(variable, model, "\n")
  }
}

# 07.4 Prediction of uncertainy  ===============================================

for (i in i:length(Prediction$Cov)){
  variable <- names(Prediction$Cov[[i]][length(Prediction$Cov[[i]])])
  model <- selected_model[depth,variable]

uncertainty_raster <- raster_stack_normalised[[1]]
uncertainty_raster <- writeStart(uncertainty_raster, paste0("./export/uncertainy_DSM/",depth,"/Uncertainty_QRF_",variable,"_for_",depth,"_soil.tif"), overwrite = TRUE)

# Subset the stack, keeping only the desired layers
block_info <- blockSize(raster_stack_normalised)
names <- names(Prediction$Cov[[i]][1:length(Prediction$Cov[[i]])-1])
raster_subset <- subset(raster_stack_normalised, names)

# Loop through each block of the raster stack
for (g in 1:block_info$n) {
  start_time <- proc.time()
  block <- getValuesBlock(raster_subset , row = block_info$row[g], nrows = block_info$nrows[g])
  block <- as.data.frame(block)
  
  # Process the block
  predicted_block <- predict(Prediction$Models[[i]]$QRF$finalModel, newdata = block,  what = c(0.05, 0.5, 0.95))
  
  # Write the predicted block to the output raster
  values <- (predicted_block[,3] - predicted_block[,1]) /predicted_block[,2]
  uncertainty_raster <- writeValues(uncertainty_raster, values, block_info$row[g])
  end_time <- proc.time()
  cat(round((g/block_info$n)*100, 1),"% ", ((end_time[3] - start_time[3])*g)/60, " /", ((end_time[3] - start_time[3])*block_info$n/60), " min \n")
}
uncertainty_raster<- writeStop(uncertainty_raster)
print(variable)
}
