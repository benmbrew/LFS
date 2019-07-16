#############
# this script will be used to explore batches.
# sessionInfo()

###########
# initialize libraries
##########
library(tidyverse)
library(stringr)

##########
# initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')

##########
# load data image from data created in clean_data.R
##########
load('/home/benbrew/Desktop/batch_raw.RData')
# load('/home/benbrew/Desktop/batch_quan.RData')

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

###########
# function that takes an argument and returns batch pca plots
###########

findBatch <- function(time, 
                      cases_and_controls, 
                      sentrix,
                      cases, 
                      controls,
                      all_cases,
                      all_850,
                      use_legend,
                      plot_name)
{
  
  
  if (sentrix) {
    
    # get common features
    intersected_feats <- Reduce(intersect, 
                                list(colnames(betaCases), 
                                     colnames(betaControls), 
                                     colnames(betaValid)))
    
    if (cases) {
      
      betaCases$batch <- betaCases$sentrix_id
      combined_dat <- betaCases
      
      
    } else if(controls) {
      
      betaControls$batch <- betaControls$sentrix_id
      combined_dat <- betaControls
      

    } else {
      
      betaValid$batch <- betaValid$sentrix_id
      combined_dat <- betaValid
      
    }
    combined_dat <- combined_dat[, c('batch', intersected_feats)]
    
  }
  if (time) {
    
    betaCases$batch <- 'cases'
    betaControls$batch <- 'controls'
    betaValid$batch <- 'valid'
    
    # get common features
    intersected_feats <- Reduce(intersect, list(colnames(betaCases), colnames(betaControls), colnames(betaValid)))
    
    # subset data by common featrues 
    betaCases <- betaCases[, c('batch', intersected_feats)]
    betaControls <- betaControls[, c('batch', intersected_feats)]
    betaValid <- betaValid[,c('batch', intersected_feats)]
    
    # remove extra data
    betaCases$batch.1 <- betaControls$batch.1 <- betaValid$batch.1 <- NULL
    combined_dat <- rbind(betaCases, betaControls, betaValid)
    
   
  }
  
  if (cases_and_controls) {
    
    betaCases$batch <- 'cases'
    betaControls$batch <- 'controls'

    # get common features
    intersected_feats <- Reduce(intersect, list(colnames(betaCases), colnames(betaControls), colnames(betaValid)))
    
    # subset data by common featrues 
    betaCases <- betaCases[, c('batch', intersected_feats)]
    betaControls <- betaControls[, c('batch', intersected_feats)]
    betaValid <- betaValid[,c('batch', intersected_feats)]
    
    # remove extra data
    betaCases$batch.1 <- betaControls$batch.1 <- betaValid$batch.1 <- NULL
    combined_dat <- rbind(betaCases, betaControls, betaValid)
    
  }
  
  if (all_cases) {
    
    betaCases$batch <- 'cases'
    betaValid$batch <- 'valid'
    
    # get common features
    intersected_feats <- Reduce(intersect, list(colnames(betaCases), colnames(betaControls), colnames(betaValid)))
    
    # subset data by common featrues 
    betaCases <- betaCases[, c('batch', intersected_feats)]
    betaValid <- betaValid[,c('batch', intersected_feats)]
    
    # remove extra data
    betaCases$batch.1 <- betaValid$batch.1 <- NULL
    combined_dat <- rbind(betaCases, betaValid)
  }
  
  if (all_850) {
    
    betaControls$batch <- 'controls'
    betaValid$batch <- 'valid'
    
    # get common features
    intersected_feats <- Reduce(intersect, list(colnames(betaControls), colnames(betaValid)))
    
    # subset data by common featrues 
    betaControls <- betaControls[, c('batch', intersected_feats)]
    betaValid <- betaValid[,c('batch', intersected_feats)]
    
    # remove extra data
    betaControls$batch.1 <- betaValid$batch.1 <- NULL
    combined_dat <- rbind(betaControls, betaValid)
    
  }
  
  # cases, gender
  getPCA(pca_data = combined_dat,
         column_name = 'batch',
         name = plot_name,
         gene_start = 9,
         pca1 = 1,
         pca2 = 2,
         use_legend = use_legend)


}

# cases sentrix id
findBatch(time = F, 
          cases_and_controls = F, 
          sentrix = T, 
          cases = T, 
          controls = F, 
          all_cases = F,
          all_850 = F,
          plot_name = 'cases sentrix_id', 
          use_legend = F)


# controls sentrix id
findBatch(time = F, 
          cases_and_controls = F, 
          sentrix = T, 
          cases = F, 
          controls = T, 
          all_cases = F, 
          all_850 = F, 
          plot_name = 'controls sentrix_id', 
          use_legend = F)

# validation sentrix id
findBatch(time = F, 
          cases_and_controls = F, 
          sentrix = T, 
          cases = F, 
          controls = F, 
          all_cases = F, 
          all_850 = F, 
          plot_name = 'validation set sentrix_id', 
          use_legend = F)

# using 
findBatch(time = F, 
          cases_and_controls = F,  
          sentrix = F, 
          cases = NULL, 
          controls = NULL, 
          all_cases = T, 
          all_850 = F, 
          plot_name = 'cases & validation: 850 vs 450', 
          use_legend = T)

findBatch(time = F, 
          cases_and_controls = F, 
          sentrix = F, 
          cases = NULL, 
          controls = NULL, 
          all_cases = F, 
          all_850 = T, 
          plot_name = 'controls and valid: both 850k', 
          use_legend = T)

# using all data
findBatch(time = T, 
          cases_and_controls = F, 
          sentrix = F, 
          cases = NULL, 
          controls = NULL, 
          all_cases = F, 
          all_850 = F, 
          plot_name = 'All 3', 
          use_legend = T)



