####################################################################
# This script is the digital soil mapping for the                  #
# Northern Kurdistan region (Iraq) soils part 3/4                  #                       
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

pacman::p_load(ggplot2, DescTools, caret, grid, gridExtra, raster, readr, dplyr,  
               terra, quantregForest, compositions, cli, kableExtra, tidyr)


# 0.3 Show session infos =======================================================

sessionInfo()

# 06 Evaluation process of the models ##########################################
# Here we decided to split every run by soil depth to have a better vision
# on the running process.
# 06.1 Preparation =============================================================
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
make_subdir("./export", "evaluation") 

calc_cv <- function(x) {
  (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) * 100
}

increments <- c("0_10", "10_30", "30_50", "50_70", "70_100")
load(file = "./export/save/Models_DSM_Boruta.RData")
Models_Boruta <- Models
load(file = "./export/save/Models_DSM_RFE.RData")
Models_RFE <- Models

# 06.2 Compare Boruta and RFE ==================================================

Metrics_diff <- data.frame()
Boruta_metrics <- Models_Boruta[[1]]$Metrics[,c(1:4)]

RFE_metrics <- Models_RFE[[1]]$Metrics[,c(1:4)]

Metrics_compared <- cbind(Boruta_metrics[,1], RFE_metrics[,1])
for (i in 2:4) {
  Metrics_compared <- cbind(Metrics_compared, Boruta_metrics[,i], RFE_metrics[,i])
}
  
Metrics_compared <- as.data.frame(Metrics_compared)
  
metrics_names <- c("ME", "ME", "RMSE", "RMSE", "R2", "R2", "PICP", "PICP")
for (i in seq(1,7, by = 2)) {
  colnames(Metrics_compared)[i] <- paste0("Boruta_",metrics_names[i])
}
  
for (i in seq(2,8, by = 2)) {
  colnames(Metrics_compared)[i] <- paste0("RFE_",metrics_names[i])
}
  
Metrics_compared[,c(7:8)] <- Metrics_compared[,c(7:8)]*100
print(Metrics_compared)
write.csv(Metrics_compared, "./export/evaluation/Boruta_vs_RFE_0_10_metrics.csv")  

# 06.3 Code for the table format (LaTEX) =======================================

cond1 <- abs(Metrics_compared[[1]]) < abs(Metrics_compared[[2]])
cond3 <- abs(Metrics_compared[[3]]) < abs(Metrics_compared[[4]])
cond5 <- abs(Metrics_compared[[5]] - 1) < abs(Metrics_compared[[6]] - 1)
cond7 <- abs(Metrics_compared[[7]] - 90) < abs(Metrics_compared[[8]] - 90)


Metrics_compared[[1]] <- ifelse(cond1,
                 cell_spec(Metrics_compared[[1]], "html", background = "lightgreen"),
                 as.character(Metrics_compared[[1]]))

Metrics_compared[[3]] <- ifelse(cond3,
                 cell_spec(Metrics_compared[[3]], "html", background = "lightgreen"),
                 as.character(Metrics_compared[[3]]))

Metrics_compared[[5]] <- ifelse(cond5,
                 cell_spec(Metrics_compared[[5]], "html", background = "lightgreen"),
                 as.character(Metrics_compared[[5]]))

Metrics_compared[[7]] <- ifelse(cond7,
                 cell_spec(Metrics_compared[[7]], "html", background = "lightgreen"),
                 as.character(Metrics_compared[[7]]))


kable(Metrics_compared, "html", escape = FALSE) %>%
  kable_styling(full_width = FALSE) %>%
  scroll_box(width = "100%")

# 06.4 Full metrics for table format (LaTEX) ===================================

Metrics <- list()
for (depth in increments) {
  Metrics[[depth]] <- Models_Boruta[[depth]]$Metrics[,c(1:4)]
  Metrics[[depth]][,4] <- Metrics[[depth]][,4]*100
  colnames(Metrics[[depth]]) <- c("ME", "RMSE", "R2", "PCIP")
  
}

Metrics <- lapply(Metrics, function(df){
  df %>%
    mutate(Variable = rownames(.))
})

Metrics <- lapply(names(Metrics), function(name){
  
  Metrics[[name]] %>%
    pivot_longer(
      cols = -Variable,
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    mutate(Profondeur = name)
  
})

final_metrics <- bind_rows(Metrics)



final_metrics <- final_metrics %>%
  pivot_wider(
    names_from = Profondeur,
    values_from = Value
  ) %>%
  arrange(Variable, Metric)

final_metrics[,c(3:7)] <- round(final_metrics[,c(3:7)], digit = 2)

# Change the "html" to "latex"
kable(final_metrics, "latex", escape = FALSE) %>%
  kable_styling(full_width = FALSE) %>%
  scroll_box(width = "100%")

# 06.5 Best tune ===============================================================

Tune_table <- list()
for (depth in increments) {
  depth_name <- gsub("_", "-", depth)
  depth_name <- paste0(depth_name, " cm")
  
  variables <- names(Models_Boruta[[depth]]$Models)
  for (variable in variables) {
    df <- cbind(variable, depth_name)
    df <- cbind(df, Models_Boruta[[depth]]$Models[[variable]]$bestTune)
    Tune_table[[depth]][[variable]] <- rbind(Tune_table[[depth]][[variable]], df)
  }
}

Tune_df <- lapply(1:9, function(i) {
  do.call(rbind, lapply(Tune_table, `[[`, i))
})
  
Tune_df <- do.call(rbind, Tune_df)

row.names(Tune_df) <- NULL

# Change the "html" to "latex"
kable(Tune_df, "latex", escape = FALSE) %>%
  kable_styling(full_width = FALSE) %>%
  scroll_box(width = "100%")

# 07 Prediction map for DSM ####################################################


# 07.1 Prediction loop  ========================================================

cl <- makeCluster(6)
registerDoParallel(cl)

for (depth in increments) {
  make_subdir("./export/prediction_DSM", depth)
  model_depth <- Models[[depth]]$Models
  
  cli_progress_bar(
    format =  " {.val {depth}}, {.val {variable}} {cli::pb_bar} {cli::pb_percent} [{cli::pb_current}/{cli::pb_total}] | ETA: {cli::pb_eta} - Elapsed: {cli::pb_elapsed_clock}",
    total = length(Models[[1]]$Models),
    clear = TRUE
  )
  
  for (i in 1:length(model_depth)) {
    
    variable <- colnames(Models[[depth]]$Cov[[i]][ncol(Models[[depth]]$Cov[[i]])])
    cov_names <- colnames(model_depth[[i]]$trainingData[,-1])
    rast <- covariates[[cov_names]]
    cov_df <- as.data.frame(rast, xy = TRUE)
    cov_df[is.infinite(as.matrix(cov_df))] <- NA
    cov_df <- na.omit(cov_df)
    
    n_blocks <- 10
    cov_df$block <- cut(1:nrow(cov_df), breaks = n_blocks, labels = FALSE)
    cov_df_blocks <- split(cov_df, cov_df$block)
    
    pred_model <- model_depth[[i]]$finalModel
    predicted_blocks <- list()
    
    # Create a loop for blocks
    for (j in 1:10) {
      block <- cov_df_blocks[[j]][,cov_names]
      prediction <- predict(pred_model, newdata = block,  what = c(0.05, 0.5, 0.95))
      
      # Write the predicted block to the output raster
      values <- prediction[,2]
      uncertainty <- (prediction[,3] - prediction[,1]) /prediction[,2]
      predicted_df <- cov_df_blocks[[j]][,c(1:2)]
      predicted_df$prediction <- values 
      predicted_df$uncertainty <- uncertainty
      predicted_blocks[[j]] <- predicted_df
    }
    
    predicted_rast <- do.call(rbind,predicted_blocks)
    predicted_rast <- rast(predicted_rast, type ="xyz")
    crs(predicted_rast) <- crs(rast)
    writeRaster(predicted_rast, paste0("./export/prediction_DSM/",depth,"/Prediction_map_",variable,"_for_",depth,"_soil.tif"), overwrite = TRUE)
    cli_progress_update()
  }
  cli_progress_done()
}