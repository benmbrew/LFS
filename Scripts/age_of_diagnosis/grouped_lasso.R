### 
# this script will run grouped lasso
library(clValid)
library(fpc)
library(cluster) 

# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/regression_template')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')


# read in methyl impute raw
methyl_impute <- read.csv(paste0(methyl_data, '/methyl_impute_raw.csv'), stringsAsFactors = F)
kmeans <- read.csv(paste0(data_folder, '/kmeans_labels.csv'), stringsAsFactors = F)
hier <- read.csv(paste0(data_folder, '/hier_labels.csv'), stringsAsFactors = F)


# remove X
methyl_impute$X <- NULL

# put ids in rownames for imputation
rownames(methyl_impute) <- methyl_impute[,1]
methyl_impute <- methyl_impute[, -1]
