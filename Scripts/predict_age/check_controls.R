###########################################
# this script will check the relationship between the three controls in 450 and 850k
library(dplyr)
library(ggplot2)
# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
scripts_folder <- paste0(project_folder, '/Scripts')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
idat_data <- paste0(methyl_data, '/raw_files')
clin_data <- paste0(data_folder, '/clin_data')
bumphunter_data <- paste0(data_folder, '/bumphunter_data')
model_data <- paste0(data_folder, '/model_data')

##########################
# load 450k and 850k
##########################
load(paste0(idat_data, '/imputed_idat_betas_final_control.RData'))
beta_raw_control <- beta_raw
beta_quan_control <- beta_quan
beta_swan_control <- beta_swan
beta_funnorm_control <- beta_funnorm
rm(beta_raw, beta_quan, beta_swan, beta_funnorm)
load(paste0(idat_data, '/imputed_idat_betas_final.RData'))

#########################
# Check distribution of controls age of sample vs original
#########################
hist(beta_raw$age_sample_collection)
hist(beta_raw_control$age_sample_collection)

##########################
# find controls that are both in 450 and 850
###########################

# subset control data and original data by 3 samples in both
subsetBeta <- function(control, orig) {
  # control <- control[control$cancer_diagnosis_diagnoses != 'Unaffected',]
  control_ids <- paste(control$id, collapse = '|')
  orig <- orig[grepl(control_ids, orig$id), ]
  orig_ids <- paste(orig$id, collapse = '|')
  control <- control[grepl(orig_ids, control$id), ]
  return(list(control, orig))
}

beta_data <- subsetBeta(beta_raw_control, beta_raw)
beta_orig <- beta_data[[2]]
beta_control <- beta_data[[1]]

# order the ids 
beta_orig <- beta_orig[order(beta_orig$id),]
beta_control <- beta_control[order(beta_control$id),]

# find overalpping probes between 450k and 850k for those 3 individuals (?)
orig_features <- colnames(beta_orig[, 5:ncol(beta_orig)])
con_features <- colnames(beta_control[, 5:ncol(beta_control)])
overlaps <- intersect(orig_features, con_features)

# subset orig and control by overlaps
beta_orig <- beta_orig[, c('id', 'age_sample_collection', overlaps)]
beta_control <- beta_control[, c('id' ,'age_sample_collection', overlaps)]
beta_raw_control <- beta_raw_control[, c('id' ,'age_sample_collection', overlaps)]

rm(beta_quan, beta_swan, beta_funnorm, beta_quan_control, beta_swan_control, beta_funnorm_control)

########################
# 3 plots, 450k against 850k, ordering methylation values from smallest to largest
##########################
# # plot for first patient
# orig_1 <- as.numeric(beta_orig[1, 3:ncol(beta_orig)])
# control_1 <- as.numeric(beta_control[1, 3:ncol(beta_orig)])
# 
# smoothScatter(orig_1, control_1, main = '1st Sample',
#               xlab = 'Original', ylab = 'Control')
# 
# # plot for 2nd  patient
# orig_2 <- as.numeric(beta_orig[2,3:ncol(beta_orig)])
# control_2 <- as.numeric(beta_control[2,3:ncol(beta_orig)])
# 
# smoothScatter(orig_2, control_2, main = '2nd Sample',
#               xlab = 'Original', ylab = 'Control')
# 
# # plot for 3rd patient
# orig_3 <- as.numeric(beta_orig[3,3:ncol(beta_orig)])
# control_3 <- as.numeric(beta_control[3,3:ncol(beta_orig)])
# 
# smoothScatter(orig_3, control_3, main = '3rd Sample',
#               xlab = 'Original', ylab = 'Control')
# 
# # plot for 4th patient
# orig_4 <- as.numeric(beta_orig[4,3:ncol(beta_orig)])
# control_4 <- as.numeric(beta_control[4,3:ncol(beta_orig)])
# 
# smoothScatter(orig_4, control_4, main = '4th Sample',
#               xlab = 'Original', ylab = 'Control')
# 
# # plot for 5th patient
# orig_5 <- as.numeric(beta_orig[5,3:ncol(beta_orig)])
# control_5 <- as.numeric(beta_control[5,3:ncol(beta_orig)])
# 
# smoothScatter(orig_5, control_5, main = '5th Sample',
#               xlab = 'Original', ylab = 'Control')
# 
# # plot for 6th patient
# orig_6 <- as.numeric(beta_orig[6,3:ncol(beta_orig)])
# control_6 <- as.numeric(beta_control[6,3:ncol(beta_orig)])
# 
# smoothScatter(orig_6, control_6, main = '6th Sample',
#               xlab = 'Original', ylab = 'Control')
# 
# # plot for 7th patient
# orig_7 <- as.numeric(beta_orig[7,3:ncol(beta_orig)])
# control_7 <- as.numeric(beta_control[7,3:ncol(beta_orig)])
# 
# smoothScatter(orig_7, control_7, main = '7th Sample',
#               xlab = 'Original', ylab = 'Control')
# 
# # plot for 8th patient
# orig_8 <- as.numeric(beta_orig[8,3:ncol(beta_orig)])
# control_8 <- as.numeric(beta_control[8,3:ncol(beta_orig)])
# 
# smoothScatter(orig_8, control_8, main = '8th Sample',
#               xlab = 'Original', ylab = 'Control')
# 
# # plot for 9th patient
# orig_9 <- as.numeric(beta_orig[9,3:ncol(beta_orig)])
# control_9 <- as.numeric(beta_control[9,3:ncol(beta_orig)])
# 
# smoothScatter(orig_9, control_9, main = '9th Sample',
#               xlab = 'Original', ylab = 'Control')
# 
# # plot for 10th patient
# orig_10 <- as.numeric(beta_orig[10,3:ncol(beta_orig)])
# control_10 <- as.numeric(beta_control[10,3:ncol(beta_orig)])
# 
# smoothScatter(orig_10, control_10, main = '10th Sample',
#               xlab = 'Original', ylab = 'Control')
# 
# # plot for 11th patient
# orig_11 <- as.numeric(beta_orig[11,3:ncol(beta_orig)])
# control_11 <- as.numeric(beta_control[11,3:ncol(beta_orig)])
# 
# smoothScatter(orig_11, control_11, main = '11th Sample',
#               xlab = 'Original', ylab = 'Control')
# 
# # plot for 12th patient
# orig_12 <- as.numeric(beta_orig[12,3:ncol(beta_orig)])
# control_12 <- as.numeric(beta_control[12,3:ncol(beta_orig)])
# 
# smoothScatter(orig_12, control_12, main = '12th Sample',
#               xlab = 'Original', ylab = 'Control')

# ######################
# estimate a linear model for each probe between the 12 individuals in controls and original
#####################
probe_model <- list()
probe_control_result <- list()

for (i in 3:ncol(beta_control)) {
  control <- as.data.frame(beta_control[, i])
  orig <- as.data.frame(beta_orig[, i])
  model_data <- data.frame(control = control, orig = orig)
  names(model_data) <- c('control', 'orig')
  probe_model[[i]] <- lm(orig ~ control, data = model_data)
  control <- as.numeric(beta_raw_control[, i])
  model_data_new <- data.frame(control = control)
  names(model_data_new) <- 'control'
  probe_control_result[[i]] <- predict(probe_model[[i]], newdata = model_data_new, type = 'response')
  
  print(i)
}

temp <- do.call(rbind, probe_control_result)
transform_controls <- t(temp)

######################
# add the colnames of beta_raw_control to transform controls
######################

# add cg sites
colnames(transform_controls) <- colnames(beta_raw_control[3:ncol(beta_raw_control)])

# add clinical variables
transform_controls <- as.data.frame(cbind(id = beta_raw_control$id, age_sample_collection = beta_raw_control$age_sample_collection, 
              transform_controls))

# save
rm(beta_control, beta_orig, beta_raw, beta_raw_control, clin, orig, temp, temp1, model_data_new, beta_data, 
   con_features, orig_features, overlaps, probe_control_result, probe_model)

save.image(paste0(model_data, '/transform.controls.RData'))
