###########################################

###########################################
# this script will check the relationship between the three controls in 450 and 850k
# It will also create and save linear transformed controls
# finally, it will create and save original controls with same probes used in original data

##########
# load libraries
##########
library(dplyr)
library(ggplot2)
library(graphics)
library(epiR)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
scripts_folder <- paste0(project_folder, '/Scripts')
data_folder <- paste0(project_folder, '/Data')
model_data <- paste0(data_folder, '/model_data')

# set fixed variable- this script should only use raw
method = 'raw'
thresh = 0.05

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# load data
##########
betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new.rda')))
betaControls <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_new.rda')))
betaValid <- readRDS(paste0(model_data, paste0('/', method, '_', 'valid_new.rda')))

##########
# combine data
##########
intersect_names <- Reduce(intersect, list(colnames(betaCases)[8:ncol(betaCases)], 
                                          colnames(betaControls)[8:ncol(betaControls)], 
                                          colnames(betaValid)[8:ncol(betaValid)]))
# organize each data set accordling

# cases
betaCases <- betaCases[, c('ids',
                           'p53_germline',
                           'age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender', 
                           'sentrix_id',
                           intersect_names)]
# controls
betaControls <- betaControls[, c('ids',
                                 'p53_germline',
                                 'age_diagnosis', 
                                 'age_sample_collection', 
                                 'cancer_diagnosis_diagnoses', 
                                 'gender', 
                                 'sentrix_id',
                                 intersect_names)]
#validation
betaValid <- betaValid[, c('ids',
                           'p53_germline',
                           'age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender', 
                           'sentrix_id',
                           intersect_names)]


# remove valid samples that are already present in betaCases
betaValidMod <- betaValid[!betaValid$ids %in% betaCases$ids,]

betaCasesMod <- getModData(betaCases)

##########
# find ids that are in both cases and controls 
##########
shared_ids_case_con <- paste(intersect(betaCases$ids, betaControls$ids), collapse = '|')

# subset cases and contorls based on these ids
cases_sub_case_con <- betaCases[grepl(shared_ids_case_con, betaCases$ids),]
controls_sub_case_con <- betaControls[grepl(shared_ids_case_con, betaControls$ids),]

# get probes
sample_list_con <- get_diff_probes_keep(cases_sub_case_con, 
                                        controls_sub_case_con, 
                                        thresh = thresh)

# save data
saveRDS(sample_list_con, paste0(model_data, paste0('/sample_list_con', '_', thresh, '.rda')))
# sample_list_valid <- readRDS(paste0(model_data, paste0('/sample_list_valid', '_', thresh, '.rda')))


##########
# do the same for validation set 
##########

# find ids that are in both cases and controls 
shared_ids_case_valid <- paste(intersect(betaCases$ids, betaValid$ids), collapse = '|')

# subset cases and contorls based on these ids
cases_sub_case_valid <- betaCases[grepl(shared_ids_case_valid, betaCases$ids, useBytes = T),]

# remove 3506, 3503, 3504, 3507, and 3634
remove_case_ids <- '3506|3503|3504|3507|3634'
cases_sub_case_valid <- cases_sub_case_valid[!grepl(remove_case_ids, cases_sub_case_valid$ids),]

# subset
valid_sub_case_valid <- betaValid[grepl(shared_ids_case_valid, betaValid$ids),]

# get subset for validation now
sample_list_valid <- get_diff_probes_keep(cases = cases_sub_case_valid, 
                                          other_850 = valid_sub_case_valid, 
                                          thresh = thresh)

# save 
saveRDS(sample_list_valid, paste0(model_data, paste0('/sample_list_valid', '_', thresh, '.rda')))

##########
# combine controls and valid probes
##########
 
# get the intersection of valid and controls probes lists
length(union(sample_list_con$names, sample_list_valid$names))
sample_list <- union(sample_list_con$names, sample_list_valid$names)

# put temp_probes in right format
remove_index <- colnames(cases_sub_case_valid) %in% sample_list

length(which(remove_index))

# subset controls 
cases_sub_case_con_sub <- cases_sub_case_con[, !remove_index]
controls_sub_case_con_sub <- controls_sub_case_con[, !remove_index]

# subset validation 
cases_sub_case_valid_sub <- cases_sub_case_valid[, !remove_index]
valid_sub_case_valid_sub <- valid_sub_case_valid[, !remove_index]

# Plt
par(mfrow = c(1,2))

for(i in 1:nrow(cases_sub_case_con)) {
  
  plotCaseCon(cases_sub_case_con, controls_sub_case_con, row_index = i)
  plotCaseCon(cases_sub_case_con_sub, controls_sub_case_con_sub, row_index = i)
  
  
}

# Plt
par(mfrow = c(1,2))

for(i in 1:nrow(cases_sub_case_valid)) {
  
  plotCaseCon(cases_sub_case_valid, valid_sub_case_valid, row_index = i)
  plotCaseCon(cases_sub_case_valid_sub, valid_sub_case_valid_sub, row_index = i)
  
}

# remove small data sets 
rm(valid_sub_case_valid, valid_sub_case_valid_sub,
   controls_sub_case_con, controls_sub_case_con_sub,
   cases_sub_case_con, cases_sub_case_con_sub,
   cases_sub_case_valid, cases_sub_case_valid_sub,
   betaValid)

# remove lists
rm(intersect_names, sample_list, sample_list_con, sample_list_valid)

##########
# read in transformed data
##########

# read in transformed controls
betaControlsTransform <- readRDS(paste0(model_data, paste0('/controls_transform', '_' , method, '_' , thresh, '.rda')))

# remove ids in betaValid that are in betaCases
colnames(betaControlsTransform)[1] <- 'ids'

# read in transformed contrls
betaValidTransform <- readRDS(paste0(model_data, paste0('/valid_transform','_' , method,'_' , thresh, '.rda')))

# remove ids in betaValid that are in betaCases
colnames(betaValidTransform)[1] <- 'ids'

# get model data
betaValidTransform <- betaValidTransform[!betaValidTransform$ids %in% betaCasesMod$ids,]

# homogenize data names by creating "Mod"
betaControlsMod <- betaControls
rm(betaControls)

# subset data
betaCasesSub <- betaCasesMod[, !remove_index]
betaControlsSub <- betaControlsMod[, !remove_index]
betaValidSub <- betaValidMod[, !remove_index]

# make betaCasesFullSub
betaCasesFullSub <- betaCases[, !remove_index]

##########
# create column in each data set as an indicator
##########

# add column in
betaCases$batch <- 'betaCases'
betaCasesFullSub$batch <- 'betaCasesFullSub'
betaCasesMod$batch <- 'betaCasesMod'
betaCasesSub$batch <- 'betaCasesSub'
betaControlsMod$batch <- 'betaControlsMod'
betaControlsSub$batch <- 'betaControlsSub'
betaControlsTransform$batch <- 'betaControlsTransform'
betaValidMod$batch <- 'betaValidMod'
betaValidSub$batch <- 'betaValidSub'
betaValidTransform$batch <- 'betaValidTransform'

# get features and subset 
long_feat <- colnames(betaCases)[8:(ncol(betaCases) - 1)]
short_feat <- colnames(betaCasesFullSub)[8:(ncol(betaCasesFullSub) - 1)]

##########
# subset long_feat
##########

# subset 
betaCases <- betaCases[, c('ids', 'p53_germline', 'age_diagnosis', 'age_sample_collection',
                           'cancer_diagnosis_diagnoses', 'gender', 'sentrix_id', 'batch',
                           long_feat)]

# subset 
betaCasesMod <- betaCasesMod[, c('ids', 'p53_germline', 'age_diagnosis', 'age_sample_collection',
                                 'cancer_diagnosis_diagnoses', 'gender', 'sentrix_id', 'batch',
                                 long_feat)]

# subset 
betaControlsMod <- betaControlsMod[, c('ids', 'p53_germline', 'age_diagnosis', 'age_sample_collection',
                                       'cancer_diagnosis_diagnoses', 'gender', 'sentrix_id', 'batch',
                                        long_feat)]

# subset 
betaValidMod <- betaValidMod[, c('ids', 'p53_germline', 'age_diagnosis', 'age_sample_collection',
                                 'cancer_diagnosis_diagnoses', 'gender', 'sentrix_id', 'batch',
                                 long_feat)]

##########
# subset short feat
##########

#subset
betaCasesFullSub <- betaCasesFullSub[, c('ids', 'p53_germline', 'age_diagnosis', 'age_sample_collection',
                                         'cancer_diagnosis_diagnoses', 'gender', 'sentrix_id', 'batch',
                                         short_feat)]

#subset
betaCasesSub <- betaCasesSub[, c('ids', 'p53_germline', 'age_diagnosis', 'age_sample_collection',
                                 'cancer_diagnosis_diagnoses', 'gender', 'sentrix_id', 'batch',
                                 short_feat)]

#subset
betaControlsSub <- betaControlsSub[, c('ids', 'p53_germline', 'age_diagnosis', 'age_sample_collection',
                                       'cancer_diagnosis_diagnoses', 'gender', 'sentrix_id', 'batch',
                                        short_feat)]

#subset
betaControlsTransform <- betaControlsTransform[, c('ids', 'p53_germline', 'age_diagnosis', 'age_sample_collection',
                                                   'cancer_diagnosis_diagnoses', 'gender', 'sentrix_id', 'batch',
                                                    short_feat)]

#subset
betaValidSub <- betaValidSub[, c('ids', 'p53_germline', 'age_diagnosis', 'age_sample_collection',
                                 'cancer_diagnosis_diagnoses', 'gender', 'sentrix_id', 'batch',
                                  short_feat)]

#subset
betaValidTransform <- betaValidTransform[, c('ids', 'p53_germline', 'age_diagnosis', 'age_sample_collection',
                                             'cancer_diagnosis_diagnoses', 'gender', 'sentrix_id', 'batch',
                                              short_feat)]

##########
# combined data
##########

# get full (full)
full_data <- rbind(betaCases, 
                   betaControlsMod,
                   betaValidMod)

# get full (mod)
full_data_mod  <- rbind(betaCasesMod, 
                        betaControlsMod,
                        betaValidMod)

# get sub data (full)
sub_data  <- rbind(betaCasesFullSub,
                   betaControlsSub,
                   betaValidSub)

# get sub data (mod)
sub_data_mod  <- rbind(betaCasesSub,
                      betaControlsSub,
                      betaValidSub)

# get transform data
trans_data <- rbind(betaCasesFullSub,
                    betaControlsTransform,
                    betaValidTransform)

# get transform data
trans_data_mod <- rbind(betaCasesSub,
                        betaControlsTransform,
                        betaValidTransform)

rm(list = ls()[grepl('beta*', ls())])

# save.image('/home/benbrew/Desktop/compare_probes_quan.RData')
# load('/home/benbrew/Desktop/compare_probes.RData')

##########
# pca of larger data sets
##########

# DATA WITH ALL FEATURES

# full data - 200 cases, all features
getPCA(pca_data = full_data,
       column_name = 'batch', 
       name = '200 cases, all features', 
       gene_start = 9, 
       pca1 = 1, 
       pca2 = 2, 
       use_legend = T)

# full data mod- 80 cases, all features
getPCA(pca_data = full_data_mod,
       column_name = 'batch', 
       name = '80 cases, all features', 
       gene_start = 9, 
       pca1 = 1, 
       pca2 = 2, 
       use_legend = T)

# DATA WITH SCRUBBED FEATURES

# sub data - 200 cases, sub features
getPCA(pca_data = sub_data,
       column_name = 'batch', 
       name = '200 cases, removed features', 
       gene_start = 9, 
       pca1 = 1, 
       pca2 = 2, 
       use_legend = T)

# sub data - 200 cases, sub features
getPCA(pca_data = sub_data_mod,
       column_name = 'batch', 
       name = '80 cases, removed features', 
       gene_start = 9, 
       pca1 = 1, 
       pca2 = 2, 
       use_legend = T)

# sub data - 200 cases, sub features
getPCA(pca_data = trans_data,
       column_name = 'batch', 
       name = '200 cases, removed features and transformed', 
       gene_start = 9, 
       pca1 = 1, 
       pca2 = 2, 
       use_legend = T)

# sub data - 200 cases, sub features
getPCA(pca_data = trans_data_mod,
       column_name = 'batch', 
       name = '80 cases, removed features and transformed', 
       gene_start = 9, 
       pca1 = 1, 
       pca2 = 2, 
       use_legend = T)



##########
# read in sample lists and see 
##########

##########
# explore controls list
##########

# load in list of probes that were "bad" for at least one of the 13 controls
sample_list_con <- readRDS(paste0(model_data, paste0('/sample_list_con', '_', thresh, '.rda')))

##########
# explore valid list
##########

# load in list of probes that were "bad" for at least one of the 8 validation set
sample_list_valid <- readRDS(paste0(model_data, paste0('/sample_list_valid', '_', thresh, '.rda')))


##########
# 0 and 11
##########

