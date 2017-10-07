
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
# fixed variables
##########
method = 'raw'
k = 5

##########
# load data
##########
# read in full m value data 
betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new_m_fam.rda')))
betaControls <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_new_m_fam.rda'))) #34 449936
betaValid <- readRDS(paste0(model_data, paste0('/', method, '_', 'valid_new_m_fam.rda')))
#35 449783


###########
# make id into ids
###########
colnames(betaCases)[1] <- 'ids'
colnames(betaControls)[1] <- 'ids'
colnames(betaValid)[1] <- 'ids'

##########
# remove inf
##########
betaCases <- removeInf(betaCases, probe_start = 8)
betaControls <- removeInf(betaControls, probe_start = 8)
betaValid<- removeInf(betaValid, probe_start = 8)


# get old controls - Mut and 'Unaffected'
betaControlsOld <- subset(betaCases, p53_germline == 'Mut' & 
                            cancer_diagnosis_diagnoses == 'Unaffected')

# get p53, not 'Unaffected'
betaCases <- getModData(betaCases)

# get rid of cancer samples in controls 
betaControls <- betaControls[grepl('Unaffected', betaControls$cancer_diagnosis_diagnoses),]

#subset valid
betaValid <- betaValid[!betaValid$ids %in% betaCases$ids,]

##########
# get intersecting colnames and prepare data for modeling
##########

intersect_names <- Reduce(intersect, list(colnames(betaCases)[8:ncol(betaCases)], 
                                          colnames(betaControls)[8:ncol(betaControls)], 
                                          colnames(betaValid)[8:ncol(betaValid)]))
# assign dataframe identifier
betaCases$type <- 'cases_450k'
betaControls$type <- 'controls_850k'
betaControlsOld$type <- 'controls_450k'
betaValid$type <- 'valid_850k'

# cases
betaCases <- betaCases[, c('ids',
                           'p53_germline',
                           'age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender',
                           'type',
                           'family_name',
                           intersect_names)]
# controls
betaControls <- betaControls[, c('ids',
                                 'p53_germline',
                                 'age_diagnosis', 
                                 'age_sample_collection', 
                                 'cancer_diagnosis_diagnoses', 
                                 'gender', 
                                 'type',
                                 'family_name',
                                 intersect_names)]

# controls
betaControlsOld <- betaControlsOld[, c('ids',
                                       'p53_germline',
                                       'age_diagnosis', 
                                       'age_sample_collection', 
                                       'cancer_diagnosis_diagnoses', 
                                       'gender', 
                                       'type',
                                       'family_name',
                                       intersect_names)]

#validation
betaValid <- betaValid[, c('ids', 
                           'p53_germline',
                           'age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender', 
                           'type',
                           'family_name',
                           intersect_names)]



# get controls full
betaControlsFull <- rbind(betaControls,
                          betaControlsOld)


# remove duplicates from betaControlsFull
length(which(duplicated(betaControlsFull$ids)))
betaControlsFull <- betaControlsFull[!duplicated(betaControlsFull$ids),]

#########

# get full data
betaFull <- rbind(betaCases,
                  betaControls,
                  betaControlsOld,
                  betaControlsFull,
                  betaValid)

betaFull$family_name.1 <- NULL

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

saveRDS(cases_full_sub, paste0(model_data, paste0('/', method, '_', 'cases_sub.rda')))

##########
# get sub of cases sub by beta valid
##########






