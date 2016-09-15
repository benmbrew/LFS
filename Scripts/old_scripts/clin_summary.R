
#########################################
# This script will give a summary of the clinical data.
# load library
library(dplyr)
library(ggplot2)
library(reshape2)
######################################### Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/classification_template')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')

# Load in clinical data
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = FALSE)
# change structure of column 
columnStructure <- function(data, col_names) {
  
  for (i in 1:length(col_names)) {
    
    data_column <- col_names[i]
    
    if (grepl('age', data_column)) {
      data[,data_column] <- as.numeric(as.character(data[,data_column]))
    } else {
      data[,data_column] <- as.factor(as.character(data[,data_column]))
    }
    
  }
  return(data)
}

clin <- columnStructure(clin, col_names = names(clin))

###########################
# Groupy by variables and get counts.
cancer <- clin %>%
  group_by(cancer_diagnosis_diagnoses, p53_germline) %>%
  summarise(counts = n(),
            mean_age = mean(age_diagnosis, na.rm = T))


