###############################################
# This script will run a PCA on the LFS methylation data and plot age of diagnosis.
# existing by gene methylation data we have. 
# This is the 4th step in the pipeline
library(dplyr)
library(stringr)
library(impute)

# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/analyze')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')

setwd(data_folder)

if('methyl_lsa.RData' %in% dir()){
  
  load(paste0(data_folder, '/methyl_lsa.RData'))
  
}else{
  
  results_file <- "file.txt"
  constructPath <- function(intermediate_folders, parent=results_folder,
                            file=results_file) {
    paste(parent, intermediate_folders, file, sep="/")
  }
  
  incomplete_file <- constructPath("Incomplete")
  imputed_file <- constructPath("Imputed")
  jvmGBLimit <- 8
  
  source(paste0(project_folder, '/Code/Functions/lsaImputation.R'))
  
  ################################################################
  # Read in methyl and clinical data and join by ids
  ################################################################
  
  # Read in methylation data
  methyl <- read.csv(paste0(methyl_data, '/methyl.csv'), stringsAsFactors = FALSE)
  # methyl_tumor <- data.matrix(read.csv(paste0(methyl_data, '/methyl_tumor.csv'), stringsAsFactors = FALSE))
  
  # put ids in rownames for imputation
  rownames(methyl) <- methyl[,1]
  methyl <- methyl[, -1]
  
  # put ids in rownames for imputation
  # rownames(methyl_tumor) <- methyl_tumor[,1]
  # methyl_tumor <- methyl_tumor[, -1]
  
  # remove duplicate rownames,so it can impute
  row.names(methyl)
  
  methyl <- as.matrix(methyl)
  
  # run lsaImputaion of methylation data
  methyl_impute_raw <- lsaImputation(incomplete_data = methyl, sample_rows = TRUE)
  #methyl_impute_raw_tumor <- lsaImputation(incomplete_data = methyl_tumor, sample_rows = TRUE)

  # join rownames and methyl_impute and then erase rownames
  methyl_impute_raw <- cbind(id = rownames(methyl_impute_raw), methyl_impute_raw)
  rownames(methyl_impute_raw) <- NULL

  # Save data to be used later
  write.csv(methyl_impute_raw, paste0(data_folder, '/methyl_impute_raw.csv'))
  
  save.image(paste0(data_folder, '/methyl_lsa.RData'))
  
}
