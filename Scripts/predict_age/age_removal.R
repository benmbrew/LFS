# this script will be used to reproduce the original held out validation set for true accuracy. 
# can also test on controls here

# source functions script
source('all_functions.R')

# create fixed objects to model and pipeline inputs and saving  
methyl_type = 'beta'
combat = TRUE
gender = TRUE
tech = FALSE 
how_many_seeds = 100
how_many_folds = 5


# condition on fixed objects to get saving identifiers

if(methyl_type == 'beta'){
  which_methyl <- 'beta'
} else {
  which_methyl <- 'm'
}

if(combat){
  which_combat <- 'use_combat'
  beta_thresh = 0.05
} else {
  which_combat <- 'no_combat'
  beta_thresh = 0.5
}

if(gender){
  is_gen = 'use_gen'
} else {
  is_gen = 'no_gen'
}

num_seeds <- paste0('seeds_', how_many_seeds)
num_folds <- paste0('folds_', how_many_folds)
k_folds <- how_many_folds


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



##########
# seperate 450 fron 850
##########

#CASES
cases_450 <- all_cases[all_cases$tech == 'batch_1',]
cases_450 <- cases_450[!is.na(cases_450$age_sample_collection),]
cases_450 <- cases_450[!is.na(cases_450$age_diagnosis),]
cases_450 <- cases_450[!is.na(cases_450$gender),]

cases_850 <- all_cases[all_cases$tech == 'batch_2',]
cases_850 <- cases_850[!is.na(cases_850$age_sample_collection),]
cases_850 <- cases_850[!is.na(cases_850$age_diagnosis),]
cases_850 <- cases_850[!is.na(cases_850$gender),]

rm(all_cases)

# gender
cases_450 <- cbind(as.data.frame(class.ind(cases_450$gender)), 
                   cases_450)

cases_850 <- cbind(as.data.frame(class.ind(cases_850$gender)), 
                   cases_850)

# rempove old tech variable 
cases_450$gender <- cases_850$gender <- NULL


# # tech
# cases_450 <- cbind(as.data.frame(class.ind(cases_450$tech)), 
#                    cases_450)
# 
# cases_850 <- cbind(as.data.frame(class.ind(cases_850$tech)), 
#                    cases_850)

# rempove old tech variable 
cases_450$tech <- cases_850$tech <- NULL


#########
# make ages back to months (times by 12)
#########

# cases

# 450
cases_450$age_diagnosis <- 
  round(cases_450$age_diagnosis*12, 2)
cases_450$age_sample_collection <- 
  round(cases_450$age_sample_collection*12, 2)

# 850
cases_850$age_diagnosis <- 
  round(cases_850$age_diagnosis*12, 2)
cases_850$age_sample_collection <- 
  round(cases_850$age_sample_collection*12, 2)


#CONTROLS
con_450 <- all_con[all_con$tech == 'batch_1',]
con_450 <- con_450[!is.na(con_450$age_sample_collection),]
con_450 <- con_450[!is.na(con_450$gender),]

con_850 <- all_con[all_con$tech == 'batch_2',]
con_850 <- con_850[!is.na(con_850$age_sample_collection),]
con_850 <- con_850[!is.na(con_850$gender),]

rm(all_con)
# gender
con_450 <- cbind(as.data.frame(class.ind(con_450$gender)), 
                 con_450)

con_850 <- cbind(as.data.frame(class.ind(con_850$gender)), 
                 con_850)

# rempove old tech variable 
con_450$gender <- con_850$gender <- NULL

# # tech
# con_450 <- cbind(as.data.frame(class.ind(con_450$tech)), 
#                  con_450)
# 
# con_850 <- cbind(as.data.frame(class.ind(con_850$tech)), 
#                  con_850)


# rempove old tech variable 
con_450$tech <- con_850$tech <- NULL

# apply to wt controls
all_con_wt <- all_con_wt[!is.na(all_con_wt$age_sample_collection),]
all_con_wt <- all_con_wt[!is.na(all_con_wt$gender),]

# ge tgender 
all_con_wt <- cbind(as.data.frame(class.ind(all_con_wt$gender)), 
                    all_con_wt)

# rempove old tech variable 
all_con_wt$gender <- NULL

# # tech
# all_con_wt <- cbind(as.data.frame(class.ind(all_con_wt$tech)), 
#                  all_con_wt)
all_con_wt$tech <- NULL

# subset to get controls lfs and wild type
con_mut <- all_con_wt[all_con_wt$p53_germline == 'MUT',]
con_wt <- all_con_wt[all_con_wt$p53_germline == 'WT',]
rm(all_con_wt)

#########
# make ages back to months (times by 12)
#########

# 450
con_450$age_diagnosis <- 
  round(con_450$age_diagnosis*12, 2)
con_450$age_sample_collection <- 
  round(con_450$age_sample_collection*12, 2)

# 850
con_850$age_diagnosis <- 
  round(con_850$age_diagnosis*12, 2)
con_850$age_sample_collection <- 
  round(con_850$age_sample_collection*12, 2)

# 450 wt
con_mut$age_sample_collection <- 
  round(con_mut$age_sample_collection*12, 2)

con_wt$age_sample_collection <- 
  round(con_wt$age_sample_collection*12, 2)

###################
# extract literature probes
###################


###################
# use pca
###################


###################
# lm extract probes
###################
