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



##########
# find ids that are in both cases and controls 
##########
shared_ids_case_con <- paste(intersect(betaCases$ids, betaControls$ids), collapse = '|')

# subset cases and contorls based on these ids
cases_sub_case_con <- betaCases[grepl(shared_ids_case_con, betaCases$ids),]
controls_sub_case_con <- betaControls[grepl(shared_ids_case_con, betaControls$ids),]

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

valid_sub_case_valid <- betaValid[grepl(shared_ids_case_valid, betaValid$ids),]

# get probes
sample_list_con <- get_diff_probes(cases_sub_case_con, 
                                   controls_sub_case_con, 
                                   thresh = thresh)
# get subset for validation now
sample_list_valid <- get_diff_probes(cases = cases_sub_case_valid, 
                                     other_850 = valid_sub_case_valid, 
                                     thresh = thresh)



###
# get the intersection of valid and controls probes lists
length(union(sample_list_con, sample_list_valid))
length(intersect(sample_list_con, sample_list_valid))

# get union and intersection of two techs
sample_list_union <- union(sample_list_con, sample_list_valid)
sample_list_intersection <- intersect(sample_list_con, sample_list_valid)

# put temp_probes in right format
remove_index_union <- colnames(cases_sub_case_valid) %in% sample_list_union
remove_index_int <- colnames(cases_sub_case_valid) %in% sample_list_intersection

# removals
length(which(remove_index_int))
length(which(remove_index_union))


# subset controls int
cases_sub_case_con_int <- cases_sub_case_con[, !remove_index_int]
controls_sub_case_con_int <- controls_sub_case_con[, !remove_index_int]

# subset controls union
cases_sub_case_con_union <- cases_sub_case_con[, !remove_index_union]
controls_sub_case_con_union <- controls_sub_case_con[, !remove_index_union]

# subset main data
betaControlsSubInt <- betaControls[, !remove_index_int]
betaControlsSubUnion <- betaControls[, !remove_index_union]


# subset validation int
cases_sub_case_valid_int <- cases_sub_case_valid[, !remove_index_int]
valid_sub_case_valid_int <- valid_sub_case_valid[, !remove_index_int]

# subset validation int
cases_sub_case_valid_union <- cases_sub_case_valid[, !remove_index_union]
valid_sub_case_valid_union <- valid_sub_case_valid[, !remove_index_union]

# #subset main data
betaValidSubInt <- betaValid[, !remove_index_int]
betaValidSubUnion <- betaValid[, !remove_index_union]

# Plt
par(mfrow = c(1,3))

for(i in 1:nrow(cases_sub_case_con)) {
  
  plotCaseCon(cases_sub_case_con, controls_sub_case_con, row_index = i)
  plotCaseCon(cases_sub_case_con_int, controls_sub_case_con_int, row_index = i)
  plotCaseCon(cases_sub_case_con_union, controls_sub_case_con_union, row_index = i)
  
}

# Plt
par(mfrow = c(1,1))

for(i in 1:nrow(cases_sub_case_valid)) {
  
  plotCaseCon(cases_sub_case_valid, valid_sub_case_valid, row_index = i)
  plotCaseCon(cases_sub_case_valid_int, valid_sub_case_valid_int, row_index = i)
  plotCaseCon(cases_sub_case_valid_union, valid_sub_case_valid_union, row_index = i)
  
  
  
}

# # apply transformation
# controls_transformed <- linearTransform(cases_sub_case_con, controls_sub_case_con, betaControlsSub)
# valid_transformed <- linearTransform(cases_sub_case_valid, valid_sub_case_valid, betaValidSub)

# save data controls
# saveRDS(controls_transformed, paste0(model_data, paste0('/controls_transform', '_' , method, '_' , thresh, '.rda')))
saveRDS(betaControlsSubInt, paste0(model_data, paste0('/controls_no_transform_int', '_' , method, '_' , thresh, '.rda')))
saveRDS(betaControlsSubUnion, paste0(model_data, paste0('/controls_no_transform_union', '_' , method, '_' , thresh, '.rda')))


# save data valid
# saveRDS(valid_transformed, paste0(model_data, paste0('/valid_transform','_' , method,'_' , thresh, '.rda')))
saveRDS(betaValidSubInt, paste0(model_data, paste0('/valid_no_transform_int','_' , method,'_' , thresh, '.rda')))
saveRDS(betaValidSubUnion, paste0(model_data, paste0('/valid_no_transform_union','_' , method,'_' , thresh, '.rda')))
