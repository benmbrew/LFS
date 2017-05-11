####### Script will subset probes used based on highly correlated probes and near zero variance probes
# This is 5th step

##########
# initialize libraries
##########
library(caret)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
model_data <- paste0(data_folder, '/model_data')
clin_data <- paste0(data_folder, '/clin_data')

##########
# load cases only
##########
# full
quan_cases_full <- readRDS(paste0(model_data, '/quan_cases_full.rda'))

funnorm_cases_full <- readRDS(paste0(model_data, '/funnorm_cases_full.rda'))

# full
quan_cases_full <- readRDS(paste0(model_data, '/quan_cases_full.rda'))

funnorm_cases_full <- readRDS(paste0(model_data, '/funnorm_cases_full.rda'))


##########
# drop highly correlated (find appropriate threshold)
##########
cor_mat <- quan_cases[, 8:ncol(quan_cases)]
cor_mat <- cor(as.matrix(cor_mat))


highlyCorDescr <- findCorrelation(cor_mat, cutoff = .70)

quan_mat <- cor_mat[,-highlyCorDescr]

saveRDS(quan_mat,paste0(model_data, '/quan_cases.rda'))

