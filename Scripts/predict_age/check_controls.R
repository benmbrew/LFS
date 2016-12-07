###########################################
# this script will check the relationship between the three controls in 450 and 850k
# It will also create and save linear transformed controls
# finally, it will create and save original controls with same probes used in original data
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
hist(beta_raw$age_sample_collection, xlab = 'Age in months', 
     main = 'Age of sample collection', col = 'lightblue', 
     xlim = c(0,1000), ylim = c(0,70))
hist(beta_raw_control$age_sample_collection, xlab = 'Age in months',
     main = 'Age of sample collection (controls)', col = 'lightblue',
     xlim = c(0,1000), ylim = c(0,20))

##########################
# find controls that are both in 450 and 850
###########################

# subset control data and original data by 12 samples in both
subsetBeta <- function(control, orig) {
  # control <- control[control$cancer_diagnosis_diagnoses != 'Unaffected',]
  control_ids <- paste(control$id, collapse = '|')
  orig <- orig[grepl(control_ids, orig$id), ]
  orig_ids <- paste(orig$id, collapse = '|')
  control <- control[grepl(orig_ids, control$id),]
  # order the ids 
  orig <- orig[order(orig$id),]
  control <- control[order(control$id),]
  # find overalpping probes between 450k and 850k for those 3 individuals (?)
  orig_features <- colnames(orig[, 6:ncol(orig)])
  con_features <- colnames(control[, 5:ncol(control)])
  overlaps <- intersect(orig_features, con_features)
  
  # subset orig and control by overlaps
  orig <- orig[, c('id', 'age_sample_collection', overlaps)]
  control <- control[, c('id' ,'age_sample_collection', overlaps)]
  raw_control <- beta_raw_control[, c('id' ,'age_sample_collection', overlaps)]

  return(list(orig, control, raw_control))
}


# raw
raw <- subsetBeta(beta_raw_control, beta_raw)
beta_raw_12 <- raw[[1]]
beta_raw_control_12 <- raw[[2]]
beta_raw_control_overlap <- raw[[3]]

#quan
quan <- subsetBeta(beta_quan_control, beta_quan)
beta_quan_12 <- quan[[1]]
beta_quan_control_12 <- quan[[2]]
beta_quan_control_overlap <- quan[[3]]

#swan
swan <- subsetBeta(beta_swan_control, beta_swan)
beta_swan_12 <- swan[[1]]
beta_swan_control_12 <- swan[[2]]
beta_swan_control_overlap <- swan[[3]]

#funnorm
funnorm <- subsetBeta(beta_funnorm_control, beta_funnorm)
beta_funnorm_12 <- funnorm[[1]]
beta_funnorm_control_12 <- funnorm[[2]]
beta_funnorm_control_overlap <- funnorm[[3]]
########################
# 3 plots, 450k against 850k, ordering methylation values from smallest to largest
##########################
par(mfrow = c(3,4))

plotRelationship <- function(beta_orig, beta_control) {
  # plot for first patient
  orig_1 <- as.numeric(beta_orig[1, 3:ncol(beta_orig)])
  control_1 <- as.numeric(beta_control[1, 3:ncol(beta_orig)])
  
  smoothScatter(orig_1, control_1, main = '1st Sample',
                xlab = 'Original', ylab = 'Control')
  
  # plot for 2nd  patient
  orig_2 <- as.numeric(beta_orig[2,3:ncol(beta_orig)])
  control_2 <- as.numeric(beta_control[2,3:ncol(beta_orig)])
  
  smoothScatter(orig_2, control_2, main = '2nd Sample',
                xlab = 'Original', ylab = 'Control')
  
  # plot for 3rd patient
  orig_3 <- as.numeric(beta_orig[3,3:ncol(beta_orig)])
  control_3 <- as.numeric(beta_control[3,3:ncol(beta_orig)])
  
  smoothScatter(orig_3, control_3, main = '3rd Sample',
                xlab = 'Original', ylab = 'Control')
  
  # plot for 4th patient
  orig_4 <- as.numeric(beta_orig[4,3:ncol(beta_orig)])
  control_4 <- as.numeric(beta_control[4,3:ncol(beta_orig)])
  
  smoothScatter(orig_4, control_4, main = '4th Sample',
                xlab = 'Original', ylab = 'Control')
  
  # plot for 5th patient
  orig_5 <- as.numeric(beta_orig[5,3:ncol(beta_orig)])
  control_5 <- as.numeric(beta_control[5,3:ncol(beta_orig)])
  
  smoothScatter(orig_5, control_5, main = '5th Sample',
                xlab = 'Original', ylab = 'Control')
  
  # plot for 6th patient
  orig_6 <- as.numeric(beta_orig[6,3:ncol(beta_orig)])
  control_6 <- as.numeric(beta_control[6,3:ncol(beta_orig)])
  
  smoothScatter(orig_6, control_6, main = '6th Sample',
                xlab = 'Original', ylab = 'Control')
  
  # plot for 7th patient
  orig_7 <- as.numeric(beta_orig[7,3:ncol(beta_orig)])
  control_7 <- as.numeric(beta_control[7,3:ncol(beta_orig)])
  
  smoothScatter(orig_7, control_7, main = '7th Sample',
                xlab = 'Original', ylab = 'Control')
  
  # plot for 8th patient
  orig_8 <- as.numeric(beta_orig[8,3:ncol(beta_orig)])
  control_8 <- as.numeric(beta_control[8,3:ncol(beta_orig)])
  
  smoothScatter(orig_8, control_8, main = '8th Sample',
                xlab = 'Original', ylab = 'Control')
  
  # plot for 9th patient
  orig_9 <- as.numeric(beta_orig[9,3:ncol(beta_orig)])
  control_9 <- as.numeric(beta_control[9,3:ncol(beta_orig)])
  
  smoothScatter(orig_9, control_9, main = '9th Sample',
                xlab = 'Original', ylab = 'Control')
  
  # plot for 10th patient
  orig_10 <- as.numeric(beta_orig[10,3:ncol(beta_orig)])
  control_10 <- as.numeric(beta_control[10,3:ncol(beta_orig)])
  
  smoothScatter(orig_10, control_10, main = '10th Sample',
                xlab = 'Original', ylab = 'Control')
  
  # plot for 11th patient
  orig_11 <- as.numeric(beta_orig[11,3:ncol(beta_orig)])
  control_11 <- as.numeric(beta_control[11,3:ncol(beta_orig)])
  
  smoothScatter(orig_11, control_11, main = '11th Sample',
                xlab = 'Original', ylab = 'Control')
  
  # plot for 12th patient
  orig_12 <- as.numeric(beta_orig[12,3:ncol(beta_orig)])
  control_12 <- as.numeric(beta_control[12,3:ncol(beta_orig)])
  
  smoothScatter(orig_12, control_12, main = '12th Sample',
                xlab = 'Original', ylab = 'Control')
  
}

plotRelationship(beta_raw_12, beta_raw_control_12)
plotRelationship(beta_quan_12, beta_quan_control_12)
plotRelationship(beta_swan_12, beta_swan_control_12)
plotRelationship(beta_funnorm_12, beta_funnorm_control_12)



# ######################
# estimate a linear model for each probe between the 12 individuals in controls and original
#####################

linearTransform <- function(beta_orig_12, 
                            beta_control_12, 
                            beta_overlap) {
  
  probe_model <- list()
  probe_control_result <- list()
  
  for (i in 3:ncol(beta_overlap)) {
    control <- as.data.frame(beta_control_12[, i])
    orig <- as.data.frame(beta_orig_12[, i])
    model_data <- data.frame(control = control, orig = orig)
    names(model_data) <- c('control', 'orig')
    probe_model[[i]] <- lm(orig ~ control, data = model_data)
    control <- as.numeric(beta_overlap[, i])
    model_data_new <- data.frame(control = control)
    names(model_data_new) <- 'control'
    probe_control_result[[i]] <- predict(probe_model[[i]], newdata = model_data_new, type = 'response')
    
    print(i) 
  }
  
  # transpose results
  temp <- do.call(rbind, probe_control_result)
  transform_controls <- t(temp)
  
  # add cg sites
  colnames(transform_controls) <- colnames(beta_overlap[3:ncol(beta_overlap)])
  
  # add clinical variables
  transform_controls <- as.data.frame(cbind(id = beta_overlap$id, 
                                            age_sample_collection = beta_overlap$age_sample_collection, 
                                            transform_controls))
  
  return(transform_controls)
  
}

beta_raw_transform <- linearTransform(beta_raw_12, beta_raw_control_12, beta_raw_control_overlap)

beta_quan_transform <- linearTransform(beta_quan_12, beta_quan_control_12, beta_quan_control_overlap)

beta_swan_transform <- linearTransform(beta_swan_12, beta_swan_control_12, beta_swan_control_overlap)

beta_funnorm_transform <- linearTransform(beta_funnorm_12, beta_funnorm_control_12, beta_funnorm_control_overlap)

# save.image(paste0(model_data, '/controls.RData'))




