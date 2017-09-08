##########
# initialize libraries
##########
library(caret)
library(glmnet)
library(randomForest)
library(kernlab)
library(pROC)
library(Metrics)
library(doParallel)
library(nnet)
library(dplyr)
library(bumphunter)
library(sqldf)
library(e1071)
library(reshape2)

registerDoParallel(1)
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
betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new_m.rda')))
betaControls <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_new_m.rda'))) #34 449936
betaValid <- readRDS(paste0(model_data, paste0('/', method, '_', 'valid_new_m.rda')))
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
betaCases$type <- '450k'
betaControls$type <- '850k'
betaControlsOld$type <- '450k'
betaValid$type <- '850k'


# cases
betaCases <- betaCases[, c('ids',
                           'age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender',
                           'type',
                           intersect_names)]
# controls
betaControls <- betaControls[, c('ids',
                                 'age_diagnosis', 
                                 'age_sample_collection', 
                                 'cancer_diagnosis_diagnoses', 
                                 'gender', 
                                 'type',
                                 intersect_names)]

# controls
betaControlsOld <- betaControlsOld[, c('ids',
                                       'age_diagnosis', 
                                       'age_sample_collection', 
                                       'cancer_diagnosis_diagnoses', 
                                       'gender', 
                                       'type',
                                       intersect_names)]

#validation
betaValid <- betaValid[, c('ids', 
                           'age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender', 
                           'type',
                           intersect_names)]



# get controls full
betaControlsFull <- rbind(betaControls,
                          betaControlsOld)


# remove duplicates from betaControlsFull
length(which(duplicated(betaControlsFull$ids)))
betaControlsFull <- betaControlsFull[!duplicated(betaControlsFull$ids),]


##########
# read in scaled data
##########

# # read in data
betaCasesScaled <-readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new_m_scaled.rda')))
betaValidScaled <- readRDS(paste0(model_data, paste0('/', method, '_', 'valid_new_m_scaled.rda')))
betaControlsScaled <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_new_m_scaled.rda')))
betaControlsOldScaled <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_old_new_m_scaled.rda')))
betaControlsFullScaled <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_full_m_scaled.rda')))


##########
# get the variation of the probe that is orthogonal to age of sample collection
##########
get_resid <- function(mod_data) 
  
{
  
  mod_data <- mod_data[complete.cases(mod_data),]

  probes <- colnames(mod_data)[7:ncol(mod_data)]
  
  resid <- list()
  
  for (i in 7:ncol(mod_data)){
    
    temp <- mod_data[, i]
    temp1 <- mod_data$age_sample_collection
    
    resid[[i]] <- lm(temp ~ temp1)$residuals
    
    print(i)
    
  }
  
  resid <- do.call('cbind', resid)
  resid <- apply(resid, 2, function(x) as.numeric(x))
  resid <- as.data.frame(resid)
  resid <- cbind(mod_data$ids, 
                 mod_data$age_diagnosis, 
                 mod_data$age_sample_collection, 
                 mod_data$cancer_diagnosis_diagnoses,
                 mod_data$gender, 
                 mod_data$type,
                 resid)
  
  # change colnames
  colnames(resid) <- c('ids',
                       'age_diagnosis', 
                       'age_sample_collection',
                       'cancer_diagnosis_diagnoses',
                       'gender',
                       'type',
                       probes)
  
  
  
  return(resid)
  
}


##########
# unscaled
##########
# get cases resid
##########
beta_cases_resid <- get_resid(mod_data = betaCases)

##########
# beta valid resid
##########
beta_valid_resid <- get_resid(mod_data = betaValid)

##########
# get validation and cases resid combined
##########
beta_full <- rbind(betaCases,
                   betaValid)

beta_full_resid <- get_resid(mod_data = beta_full)


##########
# save data unscaled residuals
##########
saveRDS(beta_cases_resid, paste0(model_data, paste0('/', method, '_', 'cases_new_m_resid.rda')))
saveRDS(beta_valid_resid, paste0(model_data, paste0('/', method, '_', 'valid_new_m_resid.rda')))
saveRDS(beta_full_resid, paste0(model_data, paste0('/', method, '_', 'full_new_m_resid.rda')))


#############
# scaled
##########
# get cases resid
##########
beta_cases_scaled_resid <- get_resid(mod_data = betaCasesScaled)

##########
# beta valid resid
##########
beta_valid_scaled_resid <- get_resid(mod_data = betaValidScaled)

##########
# get validation and cases resid combined
##########
beta_full_scaled <- rbind(betaCasesScaled,
                          betaValidScaled)

beta_full_scaled_resid <- get_resid(mod_data = beta_full_scaled)


##########
# save data unscaled residuals
##########
saveRDS(beta_cases_scaled_resid, paste0(model_data, paste0('/', method, '_', 'cases_new_m_scaled_resid.rda')))
saveRDS(beta_valid_scaled_resid, paste0(model_data, paste0('/', method, '_', 'valid_new_m_scaled_resid.rda')))
saveRDS(beta_full_scaled_resid, paste0(model_data, paste0('/', method, '_', 'full_new_m_scaled_resid.rda')))



