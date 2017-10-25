
##########
# initialize folders
##########
home_folder <- '~/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
model_data <- paste0(data_folder, '/model_data')

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
method = 'noob'
k = 10
combat = F

##########
# load data
##########


if (combat) {
  
  betaFull <-  readRDS(paste0(model_data, paste0('/', method, '_', 'mod_data_combat_m.rda')))
  
  
} else {
  betaFull <- readRDS(paste0(model_data, paste0('/', method, '_', 'mod_data_m.rda')))
  
}


##########
# get column names
##########
intersect_names <- colnames(betaFull)[9:ncol(betaFull)]

##########
# read in all features from feat data
##########

lfs_feats_m <- readRDS(paste0(feat_data, paste0('/', method, '_', 'lfs_m.rda')))
no_cancer_feats_m <- readRDS(paste0(feat_data, paste0('/', method, '_', 'no_cancer_m.rda')))


# setwd(feat_data)
# file_list_names = list.files()
# 
# 
# # store all raw rda feature lists in feat_list
# feat_list <- lapply(file_list_names, function(x) readRDS(x))
# 
# # order feat list
# feat_list  <- feat_list[order(sapply(feat_list, length), decreasing=F)]
# 
# # select first 10 
# feat_list <- feat_list[5]
# 
# # model_names <- c('enet', 'rf', 'lasso')
# # seeds <- c(1, 2, 3)

model_names <- c('enet')
seeds <- c(1,2,3)
feat_list <- list(no_cancer_feats_m, lfs_feats_m)
file_list_names <- list('no_cancer_m', 'lfs_m')

########## 
# get a training and test set with different families 
##########
cases_full <- betaFull[!grepl('Unaffected', betaFull$cancer_diagnosis_diagnoses),]
controls_full <- betaFull[grepl('Unaffected', betaFull$cancer_diagnosis_diagnoses),]


# remove duplicates from each data set
cases_full <- cases_full[!duplicated(cases_full$ids),]
controls_full <- controls_full[!duplicated(controls_full$ids),]

# remove overlapping familes 
length(unique(cases_full$family_name))
length(unique(controls_full$family_name))

# how many overlapping familes (21) 
family_names <- cases_full$family_name[cases_full$family_name %in% controls_full$family_name]
family_names <- family_names[!duplicated(family_names)]

temp_cases <- cases_full[, 1:20] %>%
  group_by(family_name) %>%
  summarise(counts_x = n())

temp_controls <- controls_full[, 1:20] %>%
  group_by(family_name) %>%
  summarise(counts_y = n())

temp_full <- inner_join(temp_cases, temp_controls, by = 'family_name')

temp_full <- temp_full[order(temp_full$family_name),]


##########
# first remove all overlapping families from cases
##########
cases_sub <- cases_full[!cases_full$family_name %in% temp_full$family_name,]

cases_full_sub <- rbind(cases_sub, controls_full)

saveRDS(cases_full_sub, paste0(model_data, paste0('/', method, '_', 'cases_sub_m.rda')))

##########
# get sub of cases sub by beta valid
##########






