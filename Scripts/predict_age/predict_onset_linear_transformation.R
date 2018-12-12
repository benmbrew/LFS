# this scrip will predict on data that was linear transformed and age regressed out.

# source functions script
source('all_functions.R')

# create fixed objects to model and pipeline inputs and saving  
methyl_type = 'm'
gender = FALSE
tech = FALSE 
how_many_seeds = 5
how_many_folds = 5
use_offset = FALSE
remove_age  = TRUE
beta_thresh = 0.05
num_seeds = 50
k_folds = 5

# condition on fixed objects to get saving identifiers

if(methyl_type == 'beta'){
  which_methyl <- 'beta'
} else {
  which_methyl <- 'm'
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


# load normal cases
all_cases <- readRDS('../../Data/all_cases_beta.rda')

#CASES
cases_450 <- all_cases[all_cases$tech == '450k',]
cases_450 <- cases_450[!is.na(cases_450$age_sample_collection),]
cases_450 <- cases_450[!is.na(cases_450$age_diagnosis),]
cases_450 <- cases_450[!is.na(cases_450$gender),]

rm(all_cases)

# load transformed controls 
if(methyl_type == 'beta'){
  valid_transform <- readRDS(paste0('../../Data/model_data/valid_transform_beta.rda'))
  con_transform <- readRDS(paste0('../../Data/model_data/controls_transform_beta.rda'))
  all_con_wt <-  readRDS('../../Data/all_con_beta_wt.rda')
  
} else {
  valid_transform <- readRDS(paste0('../../Data/model_data/valid_transform_m.rda'))
  con_transform <- readRDS(paste0('../../Data/model_data/controls_transform_m.rda'))
  all_con_wt <-  readRDS('../../Data/all_con_m_wt.rda')
  
}

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

# load cases
cases_450 <- cbind(as.data.frame(class.ind(cases_450$gender)), 
                   cases_450)

# rempove old tech variable 
cases_450$gender <- NULL

# gender
con_transform <- cbind(as.data.frame(class.ind(con_transform$gender)), 
                   con_transform)

# rempove old tech variable 
con_transform$gender <- NULL

# gender
valid_transform <- cbind(as.data.frame(class.ind(valid_transform$gender)), 
                       valid_transform)

# rempove old tech variable 
valid_transform$gender <- NULL


# convert age to months
cases_450$age_diagnosis <- 
  round(cases_450$age_diagnosis*12, 2)
cases_450$age_sample_collection <- 
  round(cases_450$age_sample_collection*12, 2)

# controls
con_transform$age_diagnosis <- 
  round(con_transform$age_diagnosis*12, 2)
con_transform$age_sample_collection <- 
  round(con_transform$age_sample_collection*12, 2)

# valie
valid_transform$age_diagnosis <- 
  round(valid_transform$age_diagnosis*12, 2)
valid_transform$age_sample_collection <- 
  round(valid_transform$age_sample_collection*12, 2)

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

# subset to get controls lfs and wild type
con_mut <- all_con_wt[all_con_wt$p53_germline == 'MUT',]
con_wt <- all_con_wt[all_con_wt$p53_germline == 'WT',]
rm(all_con_wt)

names(con_mut)[3] <- 'ids'
names(con_wt)[3] <- 'ids'
names(con_transform)[3] <- 'ids'
names(valid_transform)[3] <- 'ids'



if (remove_age){
  clin_names <- names(cases_450)[1:12]
  feats <- names(cases_450)[13:ncol(cases_450)]
  feats <- feats[!age_probes %in% feats]
  cases_450 <- cases_450[, c(clin_names, feats)]
  con_transform <- con_transform[, c(clin_names, feats)]
  valid_transform <- valid_transform[, c(clin_names, feats)]
  con_wt <- con_wt[, c(clin_names, feats)]
  con_mut <- con_mut[, c(clin_names, feats)]

  
}
registerDoParallel(1)

# run bumphunter on LFS healthy patients (LFS no cancer) and LFS cancer patients (LFS cancer)
bh_feats <- bump_hunter(dat_1 = con_wt, 
                        dat_2 = con_mut, 
                        bump = 'lfs', 
                        boot_num = 5, 
                        beta_thresh = beta_thresh,
                        methyl_type = methyl_type,
                        g_ranges = g_ranges)

# cases
cases_450 <- join_new_features(cases_450, new_features = bh_feats)
con_transform <- join_new_features(con_transform, new_features = bh_feats)
valid_transform <- join_new_features(valid_transform, new_features = bh_feats)
con_mut <- join_new_features(con_mut, new_features = bh_feats)
con_wt <- join_new_features(con_wt, new_features = bh_feats)

# lfs probes 
lfs_bump_probes <- colnames(cases_450)[grepl('^cg', colnames(cases_450))]

rm(ratio_set, bh_feats)


# add dummy tech variable for data sets with only one, replace family_name
names(cases_450)[9] <- 'tech'
names(con_transform)[9] <- 'tech'
names(valid_transform)[9] <- 'tech'

# fill them with Zero
cases_450$tech <- '450k'
con_transform$tech <- '850k'
valid_transform$tech <- '850k'

# do the same to con_mut and con_wt
names(con_mut)[9] <- 'tech'
names(con_wt)[9] <- 'tech'

# fill new variable with right tech indication
con_mut$tech <- '450k'
con_wt$tech <- '450k'

# save trainig  set 
saveRDS(cases_450, paste0('transform_data/', 'cases_450',which_methyl, '_',
                          num_seeds, '_', k_folds, '_', is_gen, '.rda'))
# save validation set 
saveRDS(con_transform, paste0('transform_data/', 'con_transform_',which_methyl, '_',
                              num_seeds, '_', k_folds, '_', is_gen, '.rda'))

# save validation set 
saveRDS(valid_transform, paste0('transform_data/', 'valid_transform_',which_methyl, '_',
                             num_seeds, '_', k_folds, '_', is_gen, '.rda'))

# save associated lfs bumps
saveRDS(lfs_bump_probes, paste0('transform_data/', 'lfs_bumps_', which_methyl, '_', '.rda'))
############################


# create range for random sampling for cross validation
seed_range <- c(1:how_many_seeds)

# create list to store model
all_test_results <- list()

# save image
# save.image('~/Desktop/temp_valid.RData')

for(random_seed in 1:length(seed_range)) {
  # prepare data sets for modelling
  
  run_model <- function(cases_full,
                        controls_full,
                        valid_full,
                        k_folds = k_folds,
                        tech = tech,
                        gender = gender,
                        beta_thresh = beta_thresh,
                        methyl_type = methyl_type,
                        g_ranges = g_ranges) {
    
    
    
    set.seed(random_seed)
    
    # get vector of random folds
    fold_vec <- sample(1:k_folds, nrow(cases_full), replace = T)
    
    # test_data_results
    test_data_results <- list()
    
    # combine 
    for(i in 1:k_folds) {
      
      # get train and test index
      train_index <- !grepl(i, fold_vec)
      test_index <- !train_index
      
      # subset to get training and test data
      train_cases <- cases_full[train_index, ]
      test_cases <- cases_full[test_index, ]
      
      # get controls 450k mut
      temp_con <- controls_full[controls_full$tech == '850k' & controls_full$p53_germline == 'MUT',]
      temp_con <- temp_con[!duplicated(temp_con$tm_donor),]
      
      # use cases training and controls to get bumphunter features
      bh_feats <- bump_hunter(dat_1 = train_cases, 
                              dat_2 = temp_con , 
                              bump = 'cancer', 
                              boot_num = 5, 
                              beta_thresh = beta_thresh,
                              methyl_type = methyl_type,
                              g_ranges = g_ranges)
      rm(temp_con)
      
      
      # get intersect_names
      intersect_names <- names(train_cases)[13:ncol(train_cases)]
      
      # get feature list
      colnames(bh_feats)[1] <- 'chr'
      remove_features <- inner_join(bh_feats, g_ranges)$probe
      
      # take remove features out of colnames 
      bh_features <- intersect_names[!intersect_names %in% remove_features]
      
      # HERE
      # function to predict with all test, controls, controls old, and valid
      mod_result  <- run_enet_all_test(training_dat = train_cases,
                                       test_dat = test_cases,
                                       controls_dat = con_transform,
                                       valid_dat = valid_transform,
                                       age_cutoff = 72,
                                       gender = TRUE,
                                       tech = FALSE,
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
  all_test_results[[random_seed]] <- run_model(cases_full = cases_450,
                                               controls_full = con_transform,
                                               valid_full = valid_transform,
                                               k_folds = k_folds,
                                               tech = tech,
                                               gender = gender,
                                               beta_thresh = beta_thresh,
                                               g_ranges = g_ranges)
  
  message('finished working on random seed = ', random_seed)
  
}


final_dat <- do.call(rbind, all_test_results)

# # save data
saveRDS(final_dat, paste0('transform_data/', which_methyl, '_',
                          num_seeds, '_', k_folds, '_', is_gen, '.rda'))



