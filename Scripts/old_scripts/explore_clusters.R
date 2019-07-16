###########
# this script will explore clusters 

##########
# initiate libraries
##########
library(tidyverse)
library(reshape2)

registerDoParallel(1)
##########
# initialize folders
##########
home_folder <- '~/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
clin_data <- paste0(data_folder, '/clin_data')
model_data <- paste0(data_folder, '/model_data')

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
method = 'raw'
k = 5

##########
# read in cluster labels
##########
kmeans_lab_scaled <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_kmeans_labs_scaled.rda')))
kmeans_lab <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_kmeans_labs.rda')))

##########
# read in clinical data
##########
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = F)

# clean clinical idss
clin$ids <-  gsub('A|B|_|-', '', clin$blood_dna_malkin_lab_)

##########
# get clinica data for each cluster
##########
cluster_labs <- kmeans_lab
get_clin <- function(cluster_labs) {
  
  for(i in 1:length(unique(cluster_labs$count))){
    
    temp <- subset(cluster_labs, count == i)
    probe_names <- temp$probe
    
    #
    
  }
    
}

