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

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
scripts_folder <- paste0(project_folder, '/Scripts')
data_folder <- paste0(project_folder, '/Data')
model_data <- paste0(data_folder, '/model_data')

method = 'raw'

##########
# load data
##########
betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new.rda')))
betaControls <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_new.rda')))
betaValid <- readRDS(paste0(model_data, paste0('/', method, '_', 'valid_new.rda')))


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

##########
# function to plot each id against the other
##########
plotCaseCon <- function (cases, controls, row_index) 
{
  cases_cg <- as.numeric(cases[row_index, 8:ncol(cases)])
  controls_cg <- as.numeric(controls[row_index, 8:ncol(controls)])
  
  smoothScatter(cases_cg, 
                controls_cg, 
                main = row_index,
                xlab = 'cases', 
                ylab = 'controls')
  
}

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

##########
# estimate linear model and transform
##########

linearTransform <- function(cases_12, 
                            controls_12, 
                            controls_full) {
  
  probe_model <- list()
  probe_control_result <- list()
  
  for (i in 8:ncol(controls_full)) {
    
    control <- as.data.frame(controls_12[, i])
    cases <- as.data.frame(cases_12[, i])
    model_data <- data.frame(control = control, cases = cases)
    names(model_data) <- c('control', 'cases')
    probe_model[[i]] <- lm(cases ~ control, data = model_data)
    control_full <- as.numeric(controls_full[, i])
    model_data_new <- data.frame(control = control_full)
    names(model_data_new) <- 'control'
    probe_control_result[[i]] <- predict(probe_model[[i]], newdata = model_data_new, type = 'response')
    
    print(i) 
  }
  
  # transpose results
  temp <- do.call(rbind, probe_control_result)
  transform_controls <- t(temp)
  
  # add cg sites
  colnames(transform_controls) <- colnames(controls_full[8:ncol(controls_full)])
  
  # add clinical variables
  transform_controls <- as.data.frame(cbind(id = controls_full$id, 
                                            p53_germline = controls_full$p53_germline, 
                                            cancer_diagnosis_diagnoses = controls_full$cancer_diagnosis_diagnoses, 
                                            age_diagnosis = controls_full$age_diagnosis, 
                                            age_sample_collection = controls_full$age_sample_collection,
                                            gender = controls_full$gender, 
                                            sentrix_id = controls_full$sentrix_id, 
                                            transform_controls))
  
  # make numeric
  transform_controls[, 8:ncol(transform_controls)] <- apply(transform_controls[, 8:ncol(transform_controls)], 
                                                            2, 
                                                            function(x) as.numeric(x))
  
  
  return(transform_controls)
  
}

controls_transformed <- linearTransform(cases_sub, controls_sub, betaControlsSub)

saveRDS(controls_transformed, paste0(model_data, '/controls_transformed.rda'))

# scale the transformed data and save 
controls_transformed[, 8:ncol(controls_transformed)] <- scale(controls_transformed[, 8:ncol(controls_transformed)])

saveRDS(controls_transformed, paste0(model_data, '/controls_transformed_scaled.rda'))


##########
# do the same for validation setW
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

valid_transformed <- linearTransform(cases_sub, valid_sub, betaValidSub)

saveRDS(valid_transformed, paste0(model_data, '/valid_transformed.rda'))

# scale the transformed data and save 
valid_transformed[, 8:ncol(valid_transformed)] <- scale(valid_transformed[, 8:ncol(valid_transformed)])

saveRDS(valid_transformed, paste0(model_data, '/valid_transformed_scaled.rda'))



