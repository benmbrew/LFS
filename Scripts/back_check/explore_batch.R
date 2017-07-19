#############
# this script will be used to explore batches.

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
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# load data image from data created in clean_data.R
##########
load('/home/benbrew/Desktop/batch_raw.RData')
# load('/home/benbrew/Desktop/batch_quan.RData')

###########
# function that takes an argument and returns batch pca plots
###########

findBatch <- function(time, 
                      tech, 
                      sentrix,
                      cases, 
                      controls,
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
  
  if (tech) {
    
    betaCases$batch <- 'k450'
    betaControls$batch <- 'k850'
    betaValid$batch <- 'k850'
    
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
  
  # cases, gender
  getPCA(pca_data = combined_dat,
         column_name = 'batch',
         name = plot_name,
         gene_start = 9,
         pca1 = 1,
         pca2 = 2)
  
  if(use_legend) {
    legend('bottomright',  legend = unique(combined_dat$batch), 
           col=1:length(combined_dat$batch), 
           pch=16,  
           cex = 0.7)
  }
  
 

}

# cases sentrix id
findBatch(time = F, tech = F, sentrix = T, cases = T, controls = F, plot_name = 'cases sentrix_id', use_legend = F)
findBatch(time = F, tech = F, sentrix = T, cases = F, controls = T, plot_name = 'controls sentrix_id', use_legend = F)
findBatch(time = F, tech = F, sentrix = T, cases = F, controls = F, plot_name = 'validation set sentrix_id', use_legend = F)

# using all data
findBatch(time = T, tech = F, sentrix = F, cases = NULL, controls = NULL, plot_name = 'Timing Batches', use_legend = T)
findBatch(time = F, tech = T,  sentrix = F, cases = NULL, controls = NULL, plot_name = '850 vs 450', use_legend = T)


