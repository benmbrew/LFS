# this script will be used to reproduce the original held out validation set for true accuracy. 
# can also test on controls her

# source functions script
source('all_functions.R')

# create fixed objects to model and pipeline inputs and saving  
methyl_type = 'beta'
combat = TRUE
gender = TRUE
age_cutoff = 72
tech = FALSE
how_many_seeds = 3
how_many_folds = 5
use_offset = TRUE
remove_age = TRUE

# condition on fixed objects to get saving identifiers

if(methyl_type == 'beta'){
  which_methyl <- 'beta'
} else {
  which_methyl <- 'm'
}

if(combat){
  which_combat <- 'use_combat'
  beta_thresh = 0.01
} else {
  which_combat <- 'no_combat'
  beta_thresh = 0.5
}

if(gender){
  is_gen = 'use_gen'
} else {
  is_gen = 'no_gen'
}

if(use_offset){
  is_offset <- 'use_offset'
} else {
  is_offset <- 'no_offset'
}

if(remove_age){
  is_removed <- 'age_removed'
} else {
  is_removed <- 'no_removed'
  
}


num_seeds <- paste0('seeds_', how_many_seeds)
num_folds <- paste0('folds_', how_many_folds)
k_folds <- how_many_folds

##########
# read in age probes
##########
age_probes <- readRDS('../../Data/age_probes.rda')

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
names(g_ranges)[1] <- 'chr'

#########

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


if(remove_age){
  all_cases <- get_age_removal(all_cases)
  all_con <- get_age_removal(all_con)
  all_con_wt <- get_age_removal(all_con_wt)
  
}



# gender
all_cases <- cbind(as.data.frame(class.ind(all_cases$gender)), 
                   all_cases)

# tech
all_cases <- cbind(as.data.frame(class.ind(all_cases$tech)), 
                   all_cases)

all_con <- cbind(as.data.frame(class.ind(all_con$gender)), 
                   all_con)

all_con <- cbind(as.data.frame(class.ind(all_con$tech)), 
                 all_con)
# rempove old tech variable 
all_cases$gender <- NULL
all_con$gender <- NULL

# rempove old tech variable 
all_cases$tech <- NULL
all_con$tech <- NULL

#########
# make ages back to months (times by 12)
#########

# cases
all_cases$age_diagnosis <- 
  round(all_cases$age_diagnosis*12, 2)
all_cases$age_sample_collection <- 
  round(all_cases$age_sample_collection*12, 2)

# con
all_con$age_diagnosis <- 
  round(all_con$age_diagnosis*12, 2)
all_con$age_sample_collection <- 
  round(all_con$age_sample_collection*12, 2)


# apply to wt controls
all_con_wt <- all_con_wt[!is.na(all_con_wt$age_sample_collection),]
all_con_wt <- all_con_wt[!is.na(all_con_wt$gender),]

# ge tgender 
all_con_wt <- cbind(as.data.frame(class.ind(all_con_wt$gender)), 
                    all_con_wt)

# rempove old tech variable 
all_con_wt$gender <- NULL

# # tech
all_con_wt <- cbind(as.data.frame(class.ind(all_con_wt$tech)),
                 all_con_wt)
all_con_wt$tech <- NULL

# subset to get controls lfs and wild type
con_mut <- all_con_wt[all_con_wt$p53_germline == 'MUT',]
con_wt <- all_con_wt[all_con_wt$p53_germline == 'WT',]
rm(all_con_wt)

# 450 wt
con_mut$age_sample_collection <- 
  round(con_mut$age_sample_collection*12, 2)

con_wt$age_sample_collection <- 
  round(con_wt$age_sample_collection*12, 2)


# run bumphunter on LFS healthy patients (LFS no cancer) and LFS cancer patients (LFS cancer)
bh_feats <- bump_hunter(dat_1 = con_wt, 
                        dat_2 = con_mut, 
                        bump = 'lfs', 
                        boot_num = 5, 
                        beta_thresh = beta_thresh,
                        methyl_type = methyl_type,
                        g_ranges = g_ranges)

# save.image('~/Desktop/temp_beta_nocombat.RData')
# cases
all_cases <- join_new_features(all_cases, new_features = bh_feats)

# controls
all_con <- join_new_features(all_con, new_features = bh_feats)

# wt
con_mut <- join_new_features(con_mut, new_features = bh_feats)
con_wt <- join_new_features(con_wt, new_features = bh_feats)

# store the selected featufres
lfs_bump_probes <- colnames(all_cases)[grepl('^cg', colnames(all_cases))]


rm(ratio_set, bh_feats)

# separate 450 from 850
cases_450 <- all_cases[all_cases$batch_1 == 1,]
cases_850 <- all_cases[all_cases$batch_1 == 0,]
con_450 <- all_con[all_con$batch_1 == 1,]
con_850 <- all_con[all_con$batch_1 == 0,]
rm(all_cases, all_con, con_mut, con_wt)

# combine 450 and 850
all_450 <- rbind(cases_450, 
                 con_450)

all_850 <- rbind(cases_850,
                 con_850)
rm(cases_450, con_450, cases_850, con_850)
###########################

# save trainig  set 
saveRDS(all_450, paste0('validation_age_predictions_surv/', 'all_450_',which_methyl, '_', which_combat, '_',
                          num_seeds, '_', num_folds, '_', is_gen, '_', is_offset,'.rda'))

# save validation set 
saveRDS(all_850, paste0('validation_age_predictions_surv/', 'all_850_',which_methyl, '_', which_combat, '_',
                             num_seeds, '_', num_folds, '_', is_gen,'_', is_offset, '.rda'))

# save associated lfs bumps
saveRDS(lfs_bump_probes, paste0('validation_age_predictions_surv/', 'lfs_bumps_', which_methyl, '_', 
                                which_combat, '.rda'))
############################

# create range for random sampling for cross validation
seed_range <- c(1:how_many_seeds)

# create list to store model
all_test_results <- list()

# save image
# save.image('~/Desktop/temp_valid.RData')

for(random_seed in 1:length(seed_range)) {
  # prepare data sets for modelling
  
  run_model <- function(full_450,
                        k_folds = k_folds,
                        tech = tech,
                        gender = gender,
                        beta_thresh = beta_thresh,
                        methyl_type = methyl_type,
                        age_cutoff = age_cutoff,
                        g_ranges = g_ranges) {
    
    
    
    # combinedata sets 
    full_data <- full_450

    set.seed(random_seed)
    
    # get vector of random folds
    fold_vec <- sample(1:k_folds, nrow(full_data), replace = T)
    
    # test_data_results
    test_data_results <- list()
    
    # combine 
    for(i in 1:k_folds) {
      
      # get train and test index
      train_index <- !grepl(i, fold_vec)
      test_index <- !train_index
      
      # subset to get training and test data
      train_cases <- full_data[train_index, ]
      test_cases <- full_data[test_index, ]
      
      # get controls 450k mut
      temp_con <- full_data[full_data$cancer_diagnosis_diagnoses == 'Unaffected',]
      temp_con <- temp_con[!duplicated(temp_con$tm_donor),]
      temp_cases <- full_data[full_data$cancer_diagnosis_diagnoses != 'Unaffected',]
      temp_cases <- temp_cases[!duplicated(temp_cases$tm_donor),]
      # use cases training and controls to get bumphunter features
      bh_feats <- bump_hunter(dat_1 = temp_cases, 
                              dat_2 = temp_con , 
                              bump = 'cancer', 
                              boot_num = 5, 
                              beta_thresh = beta_thresh,
                              methyl_type = methyl_type,
                              g_ranges = g_ranges)
      rm(temp_con, temp_cases)
      
      
      # get intersect_names
      intersect_names <- names(train_cases)[grepl('^cg', names(train_cases))]
      
      # get feature list
      colnames(bh_feats)[1] <- 'chr'
      remove_features <- inner_join(bh_feats, g_ranges)$probe
      
      # take remove features out of colnames 
      bh_features <- intersect_names[!intersect_names %in% remove_features]
      
      mod_result <- run_enet_surv(training_dat = train_cases,
                                 test_dat = test_cases,
                                 age_cutoff = age_cutoff,
                                 gender = gender,
                                 tech = tech,
                                 bh_features = bh_features)
      
      
      
      
      mod_result$seed <- random_seed
      mod_result$fold <- i
      test_data_results[[i]] <- mod_result
      print(i)
    }
    
    
    # combine list of case and control result data frames and return all result objects (two dataframes and 4 lists)
    cv_testing_results <- do.call(rbind, test_data_results)
    
    return(cv_testing_results)
  }
  
  
  # run model with 5 k fold cross validation
  all_test_results[[random_seed]] <- run_model(full_450 = all_450,
                                               k_folds = k_folds,
                                               tech = tech,
                                               gender = gender,
                                               methyl_type = methyl_type,
                                               beta_thresh = beta_thresh,
                                               age_cutoff = age_cutoff,
                                               g_ranges = g_ranges)
  
  message('finished working on random seed = ', random_seed)
  
}


final_dat <- do.call(rbind, all_test_results)

# # save data
saveRDS(final_dat, paste0('validation_age_predictions_surv/', which_methyl, '_', which_combat, '_',
                          num_seeds, '_', num_folds, '_', is_gen, '_no_tech_correct.rda'))

