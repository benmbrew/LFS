
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
method = 'noob'
k = 5
combat = T
##########
# load data
##########
if(combat) {
  
  full_data <- readRDS(paste0(model_data, paste0('/', method, '_', 'full_data_combat.rda')))
  
  # get sub dat
  betaCases <- full_data[full_data$batch == 'cases',]
  betaControls <- full_data[full_data$batch == 'controls',]
  betaValid <- full_data[full_data$batch == 'valid',]
  
  
} else {
  full_data <- readRDS(paste0(model_data, paste0('/', method, '_', 'full_data.rda')))
  
  # get sub dat
  betaCases <- full_data[full_data$batch == 'cases',]
  betaControls <- full_data[full_data$batch == 'controls',]
  betaValid <- full_data[full_data$batch == 'valid',]
}



###########
# make id into ids
###########
colnames(betaCases)[1] <- 'ids'
colnames(betaControls)[1] <- 'ids'
colnames(betaValid)[1] <- 'ids'

##########
# remove inf
##########
betaCases <- removeInf(betaCases, probe_start = 9)
betaControls <- removeInf(betaControls, probe_start = 9)
betaValid<- removeInf(betaValid, probe_start = 9)


# get old controls - Mut and 'Unaffected'
betaControlsOld <- subset(betaCases, p53_germline == 'Mut' & 
                            cancer_diagnosis_diagnoses == 'Unaffected')


betaControlsWt <- subset(betaCases, p53_germline == 'WT' & 
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
betaControlsWt$type <- 'controls_wt_450k'


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
betaControlsWt <- betaControlsWt[, c('ids',
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
clin_data <- betaCases[, 1:8]

# merge with betacobat
inner_join(clin_data, betaFullCombat)


##########

# get full data
betaFull <- rbind(betaCases,
                  betaControls,
                  betaControlsOld,
                  betaControlsFull,
                  betaControlsWt,
                  betaValid)


if(!outlier) {
  saveRDS(betaFull, paste0(model_data, paste0('/', method, '_', 'full_mod.rda')))
  
} else {
  saveRDS(betaFull, paste0(model_data, paste0('/', method, '_', 'full_mod_outlier.rda')))
  
}


























# betaCasesScale <- scale_data(betaCases, probe_start = 7)
# betaValidScale <- scale_data(betaValid, probe_start = 7)
# betaControlsScale <- scale_data(betaControls, probe_start = 7)
# betaControlsOldScale <- scale_data(betaControlsOld, probe_start = 7)
# betaControlsFullScale <- scale_data(betaControlsFull, probe_start = 7)

# save data
saveRDS(betaCases, paste0(model_data, paste0('/', method, '_', 'cases_mod.rda')))
saveRDS(betaValid, paste0(model_data, paste0('/', method, '_', 'valid_mod.rda')))
saveRDS(betaControls, paste0(model_data, paste0('/', method, '_', 'controls_mod.rda')))
saveRDS(betaControlsOld, paste0(model_data, paste0('/', method, '_', 'controls_old_mod.rda')))
saveRDS(betaControlsFull, paste0(model_data, paste0('/', method, '_', 'controls_full_mod.rda')))

# # save data
# saveRDS(betaCasesScale, paste0(model_data, paste0('/', method, '_', 'cases_mod_scaled.rda')))
# saveRDS(betaValidScale, paste0(model_data, paste0('/', method, '_', 'valid_mod_scaled.rda')))
# saveRDS(betaControlsScale, paste0(model_data, paste0('/', method, '_', 'controls_mod_scaled.rda')))
# saveRDS(betaControlsOldScale, paste0(model_data, paste0('/', method, '_', 'controls_old_mod_scaled.rda')))
# saveRDS(betaControlsFullScale, paste0(model_data, paste0('/', method, '_', 'controls_full_mod_scaled.rda')))



# # # save data
# betaCases <-readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new_m_scaled.rda')))
# betaValid <- readRDS(paste0(model_data, paste0('/', method, '_', 'valid_new_m_scaled.rda')))
# betaControls <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_new_m_scaled.rda')))
# betaControlsOld <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_old_new_m_scaled.rda')))
# betaControlsFull <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_full_m_scaled.rda')))
# 
# kmeans_lab_scaled <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_kmeans_labs_scaled.rda')))

