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

method = 'raw'

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
plotCaseCon(cases_sub, controls_sub, row_index = 13)



##########
# find way to remove outliers on correlation plot 
##########

# loop through each probe and if their diff is bigger than 0.2 flag it
probe_list <- list()
sample_list <- list()

for (i in 1:nrow(cases_sub)){
  
  temp.sample_cases <- cases_sub[i, ]
  temp.sample_controls <- controls_sub[i, ]
  
  for (j in 9:ncol(temp.sample_cases)) {
    
    temp.probe_cases_name <- names(temp.sample_cases)[j]
    
    temp.probe_cases <- temp.sample_cases[, j]
    temp.probe_controls <- temp.sample_controls[, j]
    
    diff <- abs(temp.probe_cases - temp.probe_controls)
    
    if (diff > 0.2) {
      
    probe_list[[j]] <- temp.probe_cases_name
    
    }
    
  }
  
  probe_names <- do.call(rbind, probe_list)
  sample_list[[i]] <- probe_names
  
  print(i)
  
}

load('/home/benbrew/Desktop/temp.sample_list.RData')

###########
# combine list 
###########

# get union of all "bad" probes 
temp_probes <- as.data.frame(do.call(rbind, sample_list))

# remove duplicates 
temp_probes <- as.character(temp_probes[!duplicated(temp_probes$V1),])

# put temp_probes in right format
remove_index <- colnames(cases_sub) %in% temp_probes

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



##########
# estimate linear model and transform
##########

linearTransform <- function (cases_12, 
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

saveRDS(controls_transformed, paste0(model_data, paste0('/controls_transform', '_' , method, '.rda')))
saveRDS(betaControlsSub, paste0(model_data, paste0('/controls_no_transform', '_' , method, '.rda')))

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


##########
# find way to remove outliers on correlation plot 
##########

# loop through each probe and if their diff is bigger than 0.2 flag it
probe_list <- list()
sample_list <- list()

for (i in 1:nrow(cases_sub)){
  
  temp.sample_cases <- cases_sub[i, ]
  temp.sample_valid <- valid_sub[i, ]
  
  for (j in 9:ncol(temp.sample_cases)) {
    
    temp.probe_cases_name <- names(temp.sample_cases)[j]
    
    temp.probe_cases <- temp.sample_cases[, j]
    temp.probe_valid <- temp.sample_valid[, j]
    
    diff <- abs(temp.probe_cases - temp.probe_valid)
    
    if (diff > 0.2) {
      
      probe_list[[j]] <- temp.probe_cases_name
      
    }
    
  }
  
  probe_names <- do.call(rbind, probe_list)
  sample_list[[i]] <- probe_names
  
  print(i)
  
}

load('/home/benbrew/Desktop/temp.sample_list_valid.RData')

###########
# combine list 
###########

# get union of all "bad" probes 
temp_probes <- as.data.frame(do.call(rbind, sample_list))

# remove duplicates 
temp_probes <- as.character(temp_probes[!duplicated(temp_probes$V1),])

# put temp_probes in right format
remove_index <- colnames(cases_sub) %in% temp_probes

# subset cases, contorls, and controlsub
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
plotCaseCon(cases_sub, valid_sub, row_index = 8)

##########
# plot x vs z1, z2, z3
##########

valid_transformed <- linearTransform(cases_sub, valid_sub, betaValidSub)

saveRDS(valid_transformed, paste0(model_data, paste0('/valid_transform','_' , method, '.rda')))
saveRDS(betaValidSub, paste0(model_data, paste0('/valid_no_transform','_' , method, '.rda')))
