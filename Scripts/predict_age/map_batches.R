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


##########
# load data
##########
betaCases <- readRDS(paste0(model_data, '/raw_cases_new_quan.rda'))
betaControls <- readRDS(paste0(model_data, '/raw_controls_new_quan.rda'))
betaValid <- readRDS(paste0(model_data, '/raw_valid_new_quan.rda'))


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

##########
# function to plot each id against the other
##########
cases <- cases_sub
controls <- controls_sub
row_index <- 1
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




