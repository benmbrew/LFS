# this script will be used to predict cancer status using both cases and controls in the training data and the label 
# is now cancer status instead of age.

# source functions script
source('all_functions.R')

# create fixed objects to model and pipeline inputs and saving  

methyl_type = 'beta'
combat = TRUE

if(methyl_type == 'beta'){
  which_methyl <- 'beta'
} else {
  which_methyl <- 'm'
}

if(combat){
  which_combat <- 'use_combat'
} else {
  which_combat <- 'no_combat'
}

num_seeds <- 'five_seeds'
num_folds <- 'five_folds'
is_tech <- 'use_tech'
is_gen <- 'use_gen'

##########
# load genomic methyl set (from controls) - you need genetic locations by probe from this object
##########
ratio_set <- readRDS('../../Data/g_ranges.rda')

# get g_ranges
g_ranges <- ratio_set

# get probes from rownames
g_ranges$probe <- rownames(ratio_set)

# remove ch and duplicatee
g_ranges <- g_ranges[!duplicated(g_ranges$start),]
g_ranges <- g_ranges[!grepl('ch', g_ranges$probe),]

# read in all data
if(methyl_type == 'beta'){
  # beta controls wt 
  # con_450_wt <- readRDS('../../Data/con_wt_450_beta.rda')
  # con_450_wt$tech <- 'batch_1'
  
  if(combat){
    all_cases <- readRDS('../../Data/all_cases_beta_combat.rda')
    all_con <- readRDS('../../Data/all_con_beta_combat.rda')
    all_con_wt <-  readRDS('../../Data/all_con_beta_wt_combat.rda')
    
    message('loaded beta values, with combat correction')
    
  } else {
    
    all_cases <- readRDS('../../Data/all_cases_beta.rda')
    all_con <- readRDS('../../Data/all_con_beta.rda')
    all_con_wt <-  readRDS('../../Data/all_con_beta_wt.rda')
    
    
    # recode to homogenize with combat data
    all_cases$tech <- ifelse(all_cases$tech == '450k', 'batch_1', 'batch_2')
    all_con$tech <- ifelse(all_con$tech == '450k', 'batch_1', 'batch_2')
    
    
    
    message('loaded beta values, with no combat')
  }
} else {
  
  # read in wild type controls (healthy non LFS patients, to be compared with healthy (controls) LFS patients)
  # con_450_wt <- readRDS('../../Data/con_wt_450_m.rda')
  # con_450_wt$tech <- 'batch_1'
  if(combat){
    all_cases <- readRDS('../../Data/all_cases_m_combat.rda')
    all_con <- readRDS('../../Data/all_con_m_combat.rda')
    all_con_wt <-  readRDS('../../Data/all_con_beta_wt_combat.rda')
    
    message('loaded m values, with combat correction')
    
  } else {
    all_cases <- readRDS('../../Data/all_cases_m.rda')
    all_con <- readRDS('../../Data/all_con_m.rda')
    all_con_wt <-  readRDS('../../Data/all_con_beta_wt.rda')
    
    # recode to homogenize with combat data
    all_cases$tech <- ifelse(all_cases$tech == '450k', 'batch_1', 'batch_2')
    all_con$tech <- ifelse(all_con$tech == '450k', 'batch_1', 'batch_2')
    
    
    message('loaded m values, with no combat')
  }
}

