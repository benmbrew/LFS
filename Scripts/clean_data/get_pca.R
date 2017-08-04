## this script will read in batch data from get_cases, get_controls, or get_valid and explore potential batches and outliers

##########
# initialize libraries
##########
library(dplyr)
library(sva)
library(impute)

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
# read in batch data
##########
betaCases <- readRDS(paste0(methyl_data, '/betaCasesBatch.rda'))
betaControls <- readRDS(paste0(methyl_data, '/betaControlsBatch.rda'))
betaValid <- readRDS(paste0(methyl_data, '/betaValidBatch.rda'))

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# impute data - knn needs samples in columns
##########
betaCasesFull <- removeNA(betaCases, probe_start = 8)
betaControlsFull <- removeNA(betaControls, probe_start = 8) #450168
betaValidFull <- removeNA(betaValid, probe_start = 8) #450168

##########
# first check the potential sentrix_id batch effect (probes start at 8)
##########

getPCA(pca_data = betaCasesFull, 
       column_name = 'sentrix_id', 
       gene_start = 8, 
       pca1 = 1, 
       pca2 = 2, 
       name = 'cases sentrix_id', 
       use_legend = F) # 3358

getPCA(pca_data = betaControlsFull, 
       column_name = 'sentrix_id', 
       gene_start = 8, 
       pca1 = 1, 
       pca2 = 2, 
       name = 'controls sentrix_id', 
       use_legend = F)

getPCA(pca_data = betaValidFull, 
       column_name = 'sentrix_id', 
       gene_start = 8, 
       pca1 = 1, 
       pca2 = 2, 
       name = 'controls sentrix_id', 
       use_legend = F)


##########
# remove outliers 
##########

betaCasesFull <- removeOutlier(betaCasesFull, 
                               cases = T, 
                               controls = F, 
                               val =F)

betaControlsFull <- removeOutlier(betaControlsFull, 
                                  cases = F, 
                                  controls = T, 
                                  val = F)

betaValidFull <- removeOutlier(betaValidFull, 
                               cases = F, 
                               controls = F, 
                               val = T)

##########
# first check the potential sentrix_id batch effect (probes start at 8)
##########

getPCA(pca_data = betaCasesFull, 
       column_name = 'sentrix_id', 
       gene_start = 8, 
       pca1 = 1, 
       pca2 = 2, 
       name = 'cases sentrix_id', 
       use_legend = F) # 3010

getPCA(pca_data = betaControlsFull, 
       column_name = 'sentrix_id', 
       gene_start = 8, 
       pca1 = 1, 
       pca2 = 2, 
       name = 'controls sentrix_id', 
       use_legend = F)

getPCA(pca_data = betaValidFull, 
       column_name = 'sentrix_id', 
       gene_start = 8, 
       pca1 = 1, 
       pca2 = 2, 
       name = 'valid sentrix_id', 
       use_legend = F)


##########
# scale data
##########
betaCasesFull <- scaleData(betaCasesFull, probe_start = 8)
betaControlsFull <- scaleData(betaControlsFull, probe_start = 8)
betaValidFull <- scaleData(betaValidFull, probe_start = 8)

##########
# first check the potential sentrix_id batch effect (probes start at 8)
##########

getPCA(pca_data = betaCasesFull, 
       column_name = 'sentrix_id', 
       gene_start = 8, 
       pca1 = 1, 
       pca2 = 2, 
       name = 'cases sentrix_id', 
       use_legend = F) # 3010

getPCA(pca_data = betaControlsFull, 
       column_name = 'sentrix_id', 
       gene_start = 8, 
       pca1 = 1, 
       pca2 = 2, 
       name = 'controls sentrix_id', 
       use_legend = F)

getPCA(pca_data = betaValidFull, 
       column_name = 'sentrix_id', 
       gene_start = 8, 
       pca1 = 1, 
       pca2 = 2, 
       name = 'valid sentrix_id', 
       use_legend = F)

##########
# combine data and get pca
##########



# get intersection
shared_feats <- Reduce(intersect, list(colnames(betaCasesFull[, 8:ncol(betaCasesFull)]),
                                       colnames(betaControlsFull[, 8:ncol(betaControlsFull)]),
                                       colnames(betaValidFull[, 8:ncol(betaValidFull)])))

# add indicator for type of data
betaCasesFull$type <- 'cases'
betaControlsFull$type <- 'controls'
betaValidFull$type <- 'valid'

# subset data and combine
betaCasesFull <- betaCasesFull[, c('ids', 'type' , shared_feats)]
betaControlsFull <- betaControlsFull[, c('ids', 'type' , shared_feats)]
betaValidFull <- betaValidFull[, c('ids', 'type' , shared_feats)]


# full_data
full_data <- rbind(betaCasesFull,
                   betaControlsFull,
                   betaValidFull)

# look at pca
getPCA(pca_data = full_data, 
       column_name = 'type', 
       gene_start = 3, 
       name = 'full_data', 
       pca1 = 1, 
       pca2 = 2, 
       use_legend = F)
