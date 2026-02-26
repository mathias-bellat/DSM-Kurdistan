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

pacman::p_load(ggplot2, caret, patchwork, quantregForest, readr, grid, gridExtra, 
               dplyr, reshape2, cli, doParallel, compositions)

# 0.3 Show session infos =======================================================

sessionInfo()

# 04 Model tuning and run  #####################################################
# Here we decided to split every run by soil depth to have a better vision
# on the running process.
# 04.1 Load the data ===========================================================

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
make_subdir("./export", "models") 

increments <- c("0_10", "10_30", "30_50", "50_70", "70_100")
load(file = "./export/save/Preprocess.RData")

# 04.2 Prepare the data ========================================================

seed <- 1070
Models <- list()

cl <- makeCluster(6)
registerDoParallel(cl)
source("./script/QRF_models.R")

for (depth in increments) {
  SoilCovML <- Preprocess[[depth]]$Selected_cov_boruta
  depth_name <- gsub("_", " - ", depth)
  make_subdir("./export/models", depth) 
  
  FormulaML <- list()
  for (i in 1:length(SoilCovML)) { 
    SoilCovMLBoruta <- SoilCovML[[i]]
    StartTargetCov = ncol(SoilCovMLBoruta) 
    NumCovLayer = StartTargetCov - 1 
    FormulaML[[i]] <- as.formula(paste(names(SoilCovMLBoruta)[StartTargetCov]," ~ ",paste(names(SoilCovMLBoruta)[1:NumCovLayer],collapse="+")))
  }
  
  # 04.3 Run the models ========================================================
  
  cli_progress_bar(
    format = "Models {.val {depth}} {.val {variable}} {cli::pb_bar} {cli::pb_percent} [{cli::pb_current}/{cli::pb_total}] | \ ETA: {cli::pb_eta} - Time elapsed: {cli::pb_elapsed_clock}",
    total = length(SoilCovML), 
    clear = FALSE
  )
  
  # Train control
  set.seed(seed)
  TrainControl <- trainControl(method = "repeatedcv", 10, 3, savePredictions = TRUE, verboseIter = FALSE, allowParallel = FALSE)
  FitQRaFCon  <- list()
  plot_list <- list()
  VarPlots <- list()
  alr_save <- data.frame()
  variables_df <- data.frame()
  metrics_final <- data.frame()
  
  for (i in 1:length(SoilCovML)) {
    
    SoilCovMLBoruta <- SoilCovML[[i]]
    variable <- colnames(SoilCovMLBoruta[ncol(SoilCovMLBoruta)])
    
    qrf_model <- train(FormulaML[[i]], SoilCovMLBoruta,
                       method = qrf_caret, 
                       trControl = TrainControl,
                       metric = "RMSE", 
                       tuneGrid = expand.grid(mtry = seq(1,ncol(SoilCovMLBoruta-1), by = 3), nodesize = seq(1,21, by = 5)),
                       ntree = 500)

    
    FitQRaFCon[[variable]] <- qrf_model

    # 04.4 Compute metrics =====================================================
        
    best_params <- qrf_model$bestTune
    
    pred_best <- qrf_model$pred %>%
      dplyr::filter(
        mtry == best_params$mtry,
        nodesize == best_params$nodesize
      )
    
        if (variable == "alr.Sand") {
          alr_save <- pred_best
          
        } else if (variable == "alr.Silt") {
           x <- pred_best
           
           alr_tex <- list()
           metrics_alr <- data.frame(NULL)
           for (j in 1:4){
             y <- cbind(alr_save[,j], x[,j])
             alr_pred <- as.data.frame(alrInv(y))
             colnames(alr_pred) <- c("Sand", "Silt", "Clay")
             alr_tex[[j]] <- 100*(alr_pred)
              }
            
           texture <- list()
           names(alr_tex) <- colnames(x)[1:4]
           alr_df <- do.call(cbind, alr_tex)
           for (j in 1:3){ 
             texture[[j]] <- alr_df[,c(j, (j+3), (j+6), (j+9))]
             colnames(texture[[j]]) <- colnames(x)[1:4]
             texture[[j]]$Resample <- x$Resample
           }
           
           names(texture) <- c("Sand", "Silt", "Clay")
           for (j in 1:3){ 
             
             metrics_by_best <- texture[[j]] %>%
               group_by(Resample) %>%
               summarise(
                 ME   = mean(pred - obs, na.rm = TRUE),
                 RMSE = sqrt(mean((pred - obs)^2, na.rm = TRUE)),
                 R2   = cor(pred, obs, use = "complete.obs")^2,
                 PICP = mean(obs >= pred.quantile..0.05 & obs <= pred.quantile..0.95, na.rm = TRUE))
             
             metrics_df <- as.data.frame(metrics_by_best)
             
             metrics_summary <- metrics_df %>%
               summarise(
                 ME_mean   = mean(ME),
                 RMSE_mean = mean(RMSE),
                 R2_mean   = mean(R2),
                 PICP_mean = mean(PICP),
                 ME_sd     = sd(ME),
                 RMSE_sd   = sd(RMSE),
                 R2_sd     = sd(R2),
                 PICP_sd   = sd(PICP)
               )
             
             print(metrics_summary)
             metrics_summary$variable <- names(texture[j])
             metrics_alr <- rbind(metrics_alr, metrics_summary)
            }

          } else { NA
            
    }
    

    metrics_by_best <- pred_best %>%
      group_by(Resample) %>%
      summarise(
        ME   = mean(pred - obs, na.rm = TRUE),
        RMSE = sqrt(mean((pred - obs)^2, na.rm = TRUE)),
        R2   = cor(pred, obs, use = "complete.obs")^2,
        PICP = mean(obs >= pred.quantile..0.05 & obs <= pred.quantile..0.95, na.rm = TRUE))
    
    metrics_df <- as.data.frame(metrics_by_best)
    
    metrics_summary <- metrics_df %>%
      summarise(
        ME_mean   = mean(ME),
        RMSE_mean = mean(RMSE),
        R2_mean   = mean(R2),
        PICP_mean = mean(PICP),
        ME_sd     = sd(ME),
        RMSE_sd   = sd(RMSE),
        R2_sd     = sd(R2),
        PICP_sd   = sd(PICP)
      )
    
    print(metrics_summary)
    metrics_summary$variable <- variable
    metrics_final <- rbind(metrics_final, metrics_summary)
    
    
    plot_list[[variable]] <- plot(qrf_model, main = paste0("QRF tuning parameters for ", variable, " at ", depth_name, "cm"))
    
    # 04.5 Plot variable importance ================================================
    
    var_importance <- varImp(FitQRaFCon[[i]], scale = TRUE)
    importance_df <- as.data.frame(var_importance$importance)
    importance_df$Variable <- rownames(importance_df)
    
    AvgVarImportance <- importance_df %>%
      group_by(Variable) %>%
      summarise(importance_df = mean(Overall, na.rm = TRUE)) %>%
      arrange(desc(importance_df))
    
    # Select top 20 variables
    Top20Var <- AvgVarImportance %>%
      top_n(20, wt = importance_df)
    
    
    AllVarImportanceTop20 <- AvgVarImportance %>%
      filter(Variable %in% Top20Var$Variable)
    
    colnames(AllVarImportanceTop20) <- c("Covariable", "Importance")
    
    VarPlots[[i]] <- ggplot(AllVarImportanceTop20, aes(x = reorder(Covariable, Importance), y = Importance)) +
      geom_bar(stat = "identity", position = "dodge", fill = "lightblue") +  
      coord_flip() +  
      labs(title = paste0("Top 20 covariates influence accros all models of ", variable, " for ", depth_name , " cm increment"), 
           x = "Covariates", 
           y = "Importance") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    AllVarImportanceTop20$Property <- variable
    variables_df <- rbind(variables_df , AllVarImportanceTop20)
    
    cli_progress_update()
  }
  cli_progress_done()
  
  metrics_final <- rbind(metrics_final, metrics_alr)
  row.names(metrics_final) <- metrics_final$variable
  write.csv(data.frame(metrics_final[,-9]), paste0("./export/models/", depth,"/Models_metrics_for_",depth,"_soil.csv"))
  write.csv(variables_df, paste0("./export/models/", depth,"/Var_importances_for_",depth,"_soil.csv"))
  
  # 04.5 Plot models best tunning ==============================================
  
  png(paste0("./export/models/", depth,"/Models_tuning_parameters_for_",depth,"_soil.png"),    # File name
      width = 1900, height = 1200)
  grid.arrange(arrangeGrob(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]],
                           plot_list[[6]], plot_list[[7]], plot_list[[8]], plot_list[[9]], nrow = 3, ncol = 3),
               top = textGrob(paste0("Models tuning parameters for ", depth_name , " cm increment"), gp = gpar(fontsize = 14, fontface = "bold")))
  dev.off()
  
  pdf(paste0("./export/models/", depth,"/Models_tuning_parameters_for_",depth,"_soil.pdf"),    # File name
      width = 19, height = 12, 
      bg = "white") 
  
  grid.arrange(arrangeGrob(plot_list[[1]], plot_list[[2]], plot_list[[3]], plot_list[[4]], plot_list[[5]],
                           plot_list[[6]], plot_list[[7]], plot_list[[8]], plot_list[[9]], nrow = 3, ncol = 3),
               top = textGrob(paste0("Models tuning parameters for ", depth_name , " cm increment"), gp = gpar(fontsize = 14, fontface = "bold")))
  dev.off()
  
  png(paste0("./export/models/", depth,"/Variables_importance_for_",depth,"_soil.png"),    # File name
      width = 1900, height = 1200)
  grid.arrange(arrangeGrob(VarPlots[[1]], VarPlots[[2]], VarPlots[[3]], VarPlots[[4]], VarPlots[[5]],
                           VarPlots[[6]], VarPlots[[7]], VarPlots[[8]], VarPlots[[9]], nrow = 3, ncol = 3),
               top = textGrob(paste0("Variables importances for ", depth_name , " cm increment"), gp = gpar(fontsize = 14, fontface = "bold")))
  dev.off()
  
  pdf(paste0("./export/models/", depth,"/Variables_importance_for_",depth,"_soil.pdf"),    # File name
      width = 19, height = 12, 
      bg = "white") 
  
  grid.arrange(arrangeGrob(VarPlots[[1]], VarPlots[[2]], VarPlots[[3]], VarPlots[[4]], VarPlots[[5]],
                           VarPlots[[6]], VarPlots[[7]], VarPlots[[8]], VarPlots[[9]], nrow = 3, ncol = 3),
               top = textGrob(paste0("Variables importances for ", depth_name , " cm increment"), gp = gpar(fontsize = 14, fontface = "bold")))
  dev.off()

  # 04.7 Export results ==========================================================
  
  Models[[depth]] <- list(
    Cov = SoilCovML,
    Models = FitQRaFCon,
    Metrics = metrics_final,
    Var_scores = var_importance,
    Models_plots = plot_list,
    Var_plots = VarPlots)
  
}
stopCluster(cl)  

save(Models, file = "./export/save/Models_DSM.RData")
rm(list= ls())


# 05. Variation of the models ##################################################
# 05.1 RFE =====================================================================

# Just replace Boruta by RFE for each element it is associated with

# 05.1 Split for the SoilGrid ==================================================
# Add this lines


split <- createDataPartition(SoilCovML[[1]][,length(SoilCovML[[1]])], p = 0.8, list = FALSE, times = 1)
Train_data <- lapply(SoilCovML, function(df) df[split, ])
Test_data  <- lapply(SoilCovML, function(df) df[-split, ])

Models[[depth]] <- list(
  Cov = SoilCovML,
  Models = FitQRaFCon,
  Metrics = metrics_final,
  Var_scores = var_importance,
  Models_plots = plot_list,
  Var_plots = VarPlots,
  Train_data = Train_data,
  Test_data = Test_data)
