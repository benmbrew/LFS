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

# set fixed variable 
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



##########
# find ids that are in both cases and controls 
##########
shared_ids <- paste(intersect(betaCases$ids, betaControls$ids), collapse = '|')
shared_cols <- intersect(colnames(betaCases), colnames(betaControls))

# subset cases and contorls based on these ids
cases_sub <- betaCases[grepl(shared_ids, betaCases$ids),]
controls_sub <- betaControls[grepl(shared_ids, betaControls$ids),]

# subet by columns 
cases_sub <- cases_sub[, shared_cols]
controls_sub <- controls_sub[, shared_cols]

# subset betaControls as well
betaControlsSub <- betaControls[, shared_cols]

# plot data
plotCaseCon(cases_sub, controls_sub, row_index = 1)
plotCaseCon(cases_sub, controls_sub, row_index = 2)
plotCaseCon(cases_sub, controls_sub, row_index = 3)
plotCaseCon(cases_sub, controls_sub, row_index = 4)
plotCaseCon(cases_sub, controls_sub, row_index = 5)
plotCaseCon(cases_sub, controls_sub, row_index = 6)
plotCaseCon(cases_sub, controls_sub, row_index = 7)
plotCaseCon(cases_sub, controls_sub, row_index = 8)
plotCaseCon(cases_sub, controls_sub, row_index = 9)
plotCaseCon(cases_sub, controls_sub, row_index = 10)
plotCaseCon(cases_sub, controls_sub, row_index = 11)
plotCaseCon(cases_sub, controls_sub, row_index = 12)
plotCaseCon(cases_sub, controls_sub, row_index = 13)

# get probes
sample_list <- get_diff_probes(cases_sub, controls_sub, thresh = thresh)

# put temp_probes in right format
remove_index <- colnames(cases_sub) %in% sample_list

length(which(remove_index))

# subset cases, contorls, and controlsub
cases_sub <- cases_sub[, !remove_index]
controls_sub <- controls_sub[, !remove_index]
betaControlsSub <- betaControlsSub[, !remove_index]

# replot
plotCaseCon(cases_sub, controls_sub, row_index = 1)
plotCaseCon(cases_sub, controls_sub, row_index = 2)
plotCaseCon(cases_sub, controls_sub, row_index = 3)
plotCaseCon(cases_sub, controls_sub, row_index = 4)
plotCaseCon(cases_sub, controls_sub, row_index = 5)
plotCaseCon(cases_sub, controls_sub, row_index = 6)
plotCaseCon(cases_sub, controls_sub, row_index = 7)
plotCaseCon(cases_sub, controls_sub, row_index = 8)
plotCaseCon(cases_sub, controls_sub, row_index = 9)
plotCaseCon(cases_sub, controls_sub, row_index = 10)
plotCaseCon(cases_sub, controls_sub, row_index = 11)
plotCaseCon(cases_sub, controls_sub, row_index = 12)
plotCaseCon(cases_sub, controls_sub, row_index = 13)

# apply transformation
controls_transformed <- linearTransform(cases_sub, controls_sub, betaControlsSub)

# save data
saveRDS(controls_transformed, paste0(model_data, paste0('/controls_transform', '_' , method, '_' , thresh, '.rda')))
saveRDS(betaControlsSub, paste0(model_data, paste0('/controls_no_transform', '_' , method, '_' , thresh, '.rda')))

##########
# do the same for validation set 
##########

# find ids that are in both cases and controls 
shared_ids <- paste(intersect(betaCases$ids, betaValid$ids), collapse = '|')
shared_cols <- intersect(colnames(betaCases), colnames(betaValid))

# subset cases and contorls based on these ids
cases_sub <- betaCases[grepl(shared_ids, betaCases$ids, useBytes = T),]

# remove 3506, 3503, 3504, 3507, and 3634
remove_case_ids <- '3506|3503|3504|3507|3634'
cases_sub <- cases_sub[!grepl(remove_case_ids, cases_sub$ids),]

valid_sub <- betaValid[grepl(shared_ids, betaValid$ids),]

# subet by columns 
cases_sub <- cases_sub[, shared_cols]
valid_sub <- valid_sub[, shared_cols]

# subset betaControls as well
betaValidSub <- betaValid[, shared_cols]

plotCaseCon(cases_sub, valid_sub, row_index = 1)
plotCaseCon(cases_sub, valid_sub, row_index = 2)
plotCaseCon(cases_sub, valid_sub, row_index = 3)
plotCaseCon(cases_sub, valid_sub, row_index = 4)
plotCaseCon(cases_sub, valid_sub, row_index = 5)
plotCaseCon(cases_sub, valid_sub, row_index = 6)
plotCaseCon(cases_sub, valid_sub, row_index = 7)
plotCaseCon(cases_sub, valid_sub, row_index = 8)


# get subset for validation now
sample_list <- get_diff_probes(cases = cases_sub, 
                               other_850 = valid_sub, 
                               thresh = thresh)


# put temp_probes in right format
remove_index <- colnames(cases_sub) %in% sample_list

length(which(remove_index))

# subset cases, contorls, and validub
cases_sub <- cases_sub[, !remove_index]
valid_sub <- valid_sub[, !remove_index]
betaValidSub <- betaValidSub[, !remove_index]

# replot
plotCaseCon(cases_sub, valid_sub, row_index = 1)
plotCaseCon(cases_sub, valid_sub, row_index = 2)
plotCaseCon(cases_sub, valid_sub, row_index = 3)
plotCaseCon(cases_sub, valid_sub, row_index = 4)
plotCaseCon(cases_sub, valid_sub, row_index = 5)
plotCaseCon(cases_sub, valid_sub, row_index = 6)
plotCaseCon(cases_sub, valid_sub, row_index = 7)

# apply the linear transformation
valid_transformed <- linearTransform(cases_sub, valid_sub, betaValidSub)

# save data
saveRDS(valid_transformed, paste0(model_data, paste0('/valid_transform','_' , method,'_' , thresh, '.rda')))
saveRDS(betaValidSub, paste0(model_data, paste0('/valid_no_transform','_' , method,'_' , thresh, '.rda')))
