####################################################################
# This script is a preparation of the FTIR data for the CRC1070    #
#                                                                  #                       
#                                                                  #                                                   
# Author: Mathias Bellat                                           #
# Affiliation : Tubingen University                                #
# Creation date : 19/09/2024                                       #
# E-mail: mathias.bellat@uni-tuebingen.de                          #
####################################################################

# 00 Preparation ###############################################################
# 0.1 Prepare environment ======================================================

# Folder check
getwd()

# Select folder
setwd(getwd())

# Clean up workspace
rm(list = ls(all.names = TRUE))

# 0.2 Install packages =========================================================

# Load packages
install.packages("pacman")        #Install and load the "pacman" package (allow easier download of packages)
library(pacman)
pacman::p_load(prospectr, remotes, caret , dplyr, readr)    ### Install required packages
remotes::install_github("philipp-baumann/simplerspec") #Install package from Baumann for spectral analysis
library(simplerspec)

# 0.3 Show session infos =======================================================

sessionInfo()

# 01 Prepare and import data ---------------------------------------------------

# 01.1 Import the raw spectra files ############################################
lfMIR <- read_opus_univ(fnames = dir("./data/spectra/", full.names = TRUE), extract = c("spc"))

# 01.2 Preparing the MIR Data ##################################################

MIRspec_tbl <- lfMIR %>%
  gather_spc() %>%    # Gather list of spectra data into tibble data frame
  resample_spc(wn_lower = 375, wn_upper = 4500, wn_interval = 4) %>%     # Resample spectra to new wavenumber interval
  average_spc(by = "sample_id") # Average replicate scans per sample_id

MIRspec_tbl_rs <- MIRspec_tbl[seq(1, nrow(MIRspec_tbl)), c("sample_id", "metadata", "wavenumbers_rs", "spc_mean")]
MIRspec_tbl_rs <- MIRspec_tbl_rs[,-c(2,3)]  

# Create a data frame and change the columns names
MIRspec_wn <- data.frame(matrix(unlist(MIRspec_tbl_rs$spc_mean), nrow = nrow(MIRspec_tbl_rs), byrow = TRUE), stringsAsFactors = FALSE)
rownames(MIRspec_wn) <- substring(MIRspec_tbl$sample_id, 1)[seq(1, nrow(MIRspec_tbl))]
wn <- list(names(MIRspec_tbl_rs[[2]][[1]]))
wn <- unlist(wn)
names(MIRspec_wn) <- wn
MIRspec_wn <- MIRspec_wn[, -c(1)] # Remove NA values from the 4499  cm-1

# 01.3 Export the MIR raw data #################################################
write.table(MIRspec_wn, "./export/spectra/Spectra_raw_Wavenumber.txt", dec = ".", sep = ";", row.names = TRUE, col.names = TRUE, append = FALSE)

# 02 Prepare the spectra regarding the state of the art ------------------------

# 2.1 First cleaning of the spectra ############################################

MIRspec_tbl <- lfMIR %>%
  gather_spc() %>%    # Gather list of spectra data into tibble data frame
  resample_spc(wn_lower = 499, wn_upper = 4500, wn_interval = 4) %>%     # Resample spectra to new wavenumber interval
  average_spc(by = "sample_id") # Average replicate scans per sample_id

MIRspec_tbl_rs <- MIRspec_tbl[seq(1, nrow(MIRspec_tbl)), c("sample_id", "metadata", "wavenumbers_rs", "spc_mean")]
MIRspec_tbl_rs <- MIRspec_tbl_rs[,-c(2,3)]  

# Create a data frame and change the columns names
MIRspec_wn <- data.frame(matrix(unlist(MIRspec_tbl_rs$spc_mean), nrow = nrow(MIRspec_tbl_rs), byrow = TRUE), stringsAsFactors = FALSE)
rownames(MIRspec_wn) <- substring(MIRspec_tbl$sample_id, 1)[seq(1, nrow(MIRspec_tbl))]
wn <- list(names(MIRspec_tbl_rs[[2]][[1]]))
wn <- unlist(wn)
names(MIRspec_wn) <- wn
MIRspec_wn <- MIRspec_wn[, -c(1)] # Remove NA values from the 4499  cm-1

# 2.2 Remove interference ######################################################

MIRspec_wn <- MIRspec_wn[, -c(500:511)] # Remove 2,451 - 2500 cm-1 range because of low signal interferances

# 2.3 Remove outlayers values ##################################################

max(MIRspec_wn)  # Remove value higher than + 2 (Ng et al., 2018; Curran et al., 1996)
MIRspec_wn_remove_up <- MIRspec_wn[apply(MIRspec_wn, 1, function(row) any(row > 2)), ] 

# Check lower value and remove it
min(MIRspec_wn)  # Remove value lower than - 2 (Ng et al., 2018; Curran et al., 1996)
MIRspec_wn <- MIRspec_wn[!rowSums(MIRspec_wn < -2),]
MIRspec_wn_remove_low <- MIRspec_wn[apply(MIRspec_wn, 1, function(row) any(row < - 2)), ]

MIRspec_wn_remove <- cbind(MIRspec_wn_remove_low, MIRspec_wn_remove_up)

write.table(MIRspec_wn_remove, "./export/Interference_spectra.txt", dec = ".", sep = ";", row.names = TRUE, col.names = TRUE, append = FALSE)

# 03 Convert into different spectra variation ----------------------------------

# 3.1 Convert in Wavelength ####################################################
MIRspec_wn <- anti_join(MIRspec_wn, MIRspec_wn_remove, by = names(MIRspec_wn))
MIRspec_wl <- MIRspec_wn
wn <- (1 / as.numeric(names(MIRspec_wn))) * 1e7
wn <- round(wn, digits=0)
names(MIRspec_wl) <- wn 
MIRspec_wl <- na.omit(MIRspec_wl)

# 3.2 Prepare the functions ####################################################

# Remove near zero variable
remove_nzv <- function(x){
  y <- nearZeroVar(x, saveMetrics = TRUE)
  ifelse(sum(y$nzv == TRUE) == 0, x <- x, x <- x[-nearZeroVar(x)])
  return(x)
}       

# Remove variable with high correlation
remove_hcd <- function(x){
  y <- findCorrelation(x, cutoff = .98)
  ifelse(sum(y) == 0, x <- x, x <- x[,-y])
  return(x)
}

# 3.3 Convert in other spectra #################################################

Spectra <- MIRspec_wn
#Different transformation of the Spectra according to Ludwig et al. 2023
IRspectraList_Raw <- list("SG 1.5" = as.data.frame(prospectr::savitzkyGolay(Spectra, m = 1, p = 1, w = 5)),
                          "SG 1.11" = as.data.frame(prospectr::savitzkyGolay(Spectra, m = 1, p = 1, w = 11)),
                          "SG 1.17" = as.data.frame(prospectr::savitzkyGolay(Spectra, m = 1, p = 1, w = 17)),
                          "SG 1.23" = as.data.frame(prospectr::savitzkyGolay(Spectra, m = 1, p = 1, w = 23)),
                          "SG 2.5" = as.data.frame(prospectr::savitzkyGolay(Spectra, m = 1, p = 2, w = 5)),
                          "SG 2.11" = as.data.frame(prospectr::savitzkyGolay(Spectra, m = 1, p = 2, w = 11)),
                          "SG 2.17" = as.data.frame(prospectr::savitzkyGolay(Spectra, m = 1, p = 2, w = 17)),
                          "SG 2.23" = as.data.frame(prospectr::savitzkyGolay(Spectra, m = 1, p = 2, w = 23)),
                          "moving averages 5" = as.data.frame(prospectr::movav(Spectra, w = 5)),
                          "moving averages 11" = as.data.frame(prospectr::movav(Spectra, w = 11)),
                          "moving averages 17" = as.data.frame(prospectr::movav(Spectra, w = 17)),
                          "moving averages 23" = as.data.frame(prospectr::movav(Spectra, w = 23)),
                          "SNV-SG" = as.data.frame(prospectr::standardNormalVariate(prospectr::savitzkyGolay(Spectra, m = 1, p = 2, w = 11))), #Best option for machine learning treatment (See Ng et al. 2018)
                          "continuum removal" = as.data.frame(prospectr::continuumRemoval(Spectra, as.numeric(colnames(Spectra)), type = "R")))

IRspectraList <- IRspectraList_Raw

# Remove the near zero and highly correlated values
for (i in 1:length(IRspectraList)) {
  IRspectraList[[i]] <- remove_hcd(IRspectraList[[i]])
  IRspectraList[[i]] <- remove_nzv(IRspectraList[[i]])
}

IRspectraList <- c(IRspectraList, list("raw" = as.data.frame(Spectra)))

# 3.4 Export all the spectra ###################################################

#Export csv of files
for (i in names(IRspectraList)) {
  write.table(IRspectraList[i], file = paste0("./export/spectra/Full_spectra_",names(IRspectraList[i]),".txt"), dec = ".", sep = ";", row.names = TRUE, col.names = TRUE, append = FALSE, fileEncoding = "UTF-8")
}

write.table(MIRspec_wn, "./export/spectra/Spectra_Wavenumber.txt", dec = ".", sep = ";", row.names = TRUE, col.names = TRUE, append = FALSE)
write.table(MIRspec_wl, "./export/spectra/Spectra_Wavelenght.txt", dec = ".", sep = ";", row.names = TRUE, col.names = TRUE, append = FALSE)

# 04 Plot the Spectrum of the spectra ------------------------------------------

# 4.1 Plot the absorbance ######################################################

IRSpectra <- as.data.frame(row.names(MIRspec_wl))
IRSpectra$wn <- MIRspec_wn
IRSpectra$wnA<- log(1/MIRspec_wn) # Convert to wave number Absorbance
IRSpectra$wl <- MIRspec_wl
IRSpectra$wlA <- log(1/MIRspec_wl) # Convert to wavelength Absorbance
rownames(IRSpectra) <- IRSpectra$`row.names(MIRspec_wl)`

png("Soil spectra in wavenumbers absorbance.png", width = 297, height = 210, units = "mm", res = 300)
pdf("Soil spectra in wavenumbers absorbance.pdf", width = 10*2, height = 6*2)

matplot(x = colnames(IRSpectra$wn), y = t(IRSpectra$wnA),
        xlab = expression(paste("Wavenumber ", "(cm"^"-1", ")")),
        ylab = "Absorbance",
        type = "l",   #"l" = ligne
        lty = 1,
        col = 1:nrow(IRSpectra$wn))
dev.off()


png("Soil spectra in wavelenght absorbance.png", width = 297, height = 210, units = "mm", res = 300)
pdf("Soil spectra in wavelenght absorbance.pdf", width = 10*2, height = 6*2)

matplot(x = colnames(IRSpectra$wl), y = t(IRSpectra$wlA),
        xlab = "Wavelenght /nm ",
        ylab = "Absorbance",
        type = "l",   #"l" = ligne
        lty = 1,
        col = 1:nrow(IRSpectra$wl))
dev.off()

# 4.2 Plot the reflectance #####################################################

png("Soil spectra in wavenumbers reflectance.png", width = 297, height = 210, units = "mm", res = 300)
pdf("Soil spectra in wavenumbers reflectance.pdf", width = 10*2, height = 6*2)

matplot(x = colnames(IRSpectra$wn), y = t(IRSpectra$wn),
        xlab = expression(paste("Wavenumber ", "(cm"^"-1", ")")),
        ylab = "Reflectance",
        type = "l",   #"l" = ligne
        lty = 1,
        col = 1:nrow(IRSpectra$wn))
dev.off()


png("Soil spectra in wavelenght reflectance.png", width = 297, height = 210, units = "mm", res = 300)
pdf("Soil spectra in wavelenght reflectance.pdf", width = 10*2, height = 6*2)

matplot(x = colnames(IRSpectra$wl), y = t(IRSpectra$wl),
        xlab = "Wavelenght /nm ",
        ylab = "Reflectance",
        type = "l",   #"l" = ligne
        lty = 1,
        col = 1:nrow(IRSpectra$wl))
dev.off()

# 5 Kennard Stone sampling -----------------------------------------------------
# 5.1 Select the samples #######################################################

Samples_info <- read_delim("./data/Samples_info.csv", delim = ";")
selection <- Samples_info[Samples_info$Depth_cm == "0_10" & grepl("B07_2022", Samples_info$Site_name), ] # You can select the year and the depth

# Select the samples to analyse
MIRsample <- MIRspec_wl[rownames(MIRspec_wl) %in% selection$Lab_ID,]

# 5.2 Run the sampling #########################################################
sample <- kenStone(MIRsample, k = 30, metric = "euclid") # Make samples with Kennard Stone in euclidian distance

MIRsample <- MIRsample[sample$model,]
MIRsample <- row.names(MIRsample)
write.table(MIRsample, "./export/samples/B07_2022 Samples 0 - 10 cm.txt", dec = ".", sep = ";", row.names = FALSE, col.names = FALSE, append = FALSE)