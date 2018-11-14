# this script will be used to reproduce the original held out validation set for true accuracy. 
# can also test on controls here

# source functions script
source('all_functions.R')

# create fixed objects to model and pipeline inputs
methyl_type <- 'beta'
combat <- TRUE

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
  con_450_wt <- readRDS('../../Data/con_wt_450_beta.rda')
  con_450_wt$tech <- 'batch_1'
  
  if(combat){
    all_cases <- readRDS('../../Data/all_cases_beta_combat.rda')
    all_con <- readRDS('../../Data/all_con_beta_combat.rda')
    message('loaded beta values, with combat correction')
    
  } else {
    all_cases <- readRDS('../../Data/all_cases_beta.rda')
    all_con <- readRDS('../../Data/all_con_beta.rda')
    
    # recode to homogenize with combat data
    all_cases$tech <- ifelse(all_cases$tech == '450k', 'batch_1', 'batch_2')
    all_con$tech <- ifelse(all_con$tech == '450k', 'batch_1', 'batch_2')
   
    
    
    message('loaded beta values, with no combat')
  }
} else {
  
  # read in wild type controls (healthy non LFS patients, to be compared with healthy (controls) LFS patients)
  con_450_wt <- readRDS('../../Data/con_wt_450_m.rda')
  con_450_wt$tech <- 'batch_1'
  if(combat){
    all_cases <- readRDS('../../Data/all_cases_m_combat.rda')
    all_con <- readRDS('../../Data/all_con_m_combat.rda')
    message('loaded m values, with combat correction')
    
  } else {
    all_cases <- readRDS('../../Data/all_cases_m.rda')
    all_con <- readRDS('../../Data/all_con_m.rda')
    
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


# tech
cases_450 <- cbind(as.data.frame(class.ind(cases_450$tech)), 
                   cases_450)

cases_850 <- cbind(as.data.frame(class.ind(cases_850$tech)), 
                   cases_850)

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

# tech
con_450 <- cbind(as.data.frame(class.ind(con_450$tech)), 
                 con_450)

con_850 <- cbind(as.data.frame(class.ind(con_850$tech)), 
                 con_850)

# rempove old tech variable 
con_450$tech <- con_850$tech <- NULL

# apply to wt controls
con_450_wt <- con_450_wt[!is.na(con_450_wt$age_sample_collection),]
con_450_wt <- con_450_wt[!is.na(con_450_wt$gender),]

# ge tgender 
con_450_wt <- cbind(as.data.frame(class.ind(con_450_wt$gender)), 
                 con_450_wt)

# rempove old tech variable 
con_450_wt$gender <- NULL

# tech
con_450_wt <- cbind(as.data.frame(class.ind(con_450_wt$tech)), 
                 con_450_wt)
con_450_wt$tech <- NULL


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
con_450_wt$age_diagnosis <- 
  round(con_450_wt$age_diagnosis*12, 2)


# split into wt and mut data
con_450_mut <- con_450_wt[con_450_wt$p53_germline == 'MUT',]
con_450_wt <- con_450_wt[con_450_wt$p53_germline == 'WT',]

# run bumphunter on LFS healthy patients (LFS no cancer) and LFS cancer patients (LFS cancer)
bh_feats <- bump_hunter(dat_1 = con_450_wt, 
                        dat_2 = con_450_mut, 
                        bump = 'lfs', 
                        boot_num = 5, 
                        beta_thresh = 0.05,
                        methyl_type = methyl_type,
                        g_ranges = g_ranges)

# subset all dady by bh_feats
cases_full <- 

# prepare data sets for modelling
cases_full <- cases_450
controls_full <- con_450
controls_wt <- con_450_wt
valid_full <- cases_850
combine_controls <- FALSE
extra_controls <- NULL
k_folds = 5
tech <- FALSE
gender <- TRUE
beta_thresh = 0.05
g_ranges = g_ranges
i = 1

run_model <- function(cases_full,
                      controls_full,
                      controls_wt,
                      valid_full,
                      combine_controls,
                      extra_controls,
                      k_folds = k_folds,
                      tech = tech,
                      gender = gender,
                      beta_thresh = beta_thresh,
                      methyl_type = methyl_type,
                      g_ranges = g_ranges) {
  
  
  # create lists to store model results
  temp_cases <- list()
  temp_controls <- list()
  temp_valid <- list()
  temp_models <- list()
  temp_lambda_min <- list()
  temp_lambda_1se <- list()
  temp_alpha <- list()
  
  
  
  # get vector of random folds
  fold_vec <- sample(1:k_folds, nrow(cases_full), replace = T)
  
  # combine 
  for(i in 1:k_folds) {
    
    # get train and test index
    train_index <- !grepl(i, fold_vec)
    test_index <- !train_index
    
    
    train_cases <- cases_full[train_index, ]
    test_cases <- cases_full[test_index, ]
    
    
    
    # use cases training and controls to get bumphunter features
    bh_feats <- bump_hunter(dat_1 = train_cases, 
                            dat_2 = controls_full , 
                            bump = 'cancer', 
                            boot_num = 5, 
                            thresh = beta_thresh,
                            methyl_type = methyl_type,
                            g_ranges = g_ranges)
    
    
    
    # get intersect_names
    intersect_names <- names(train_cases)[14:ncol(train_cases)]
    
    # get feature list
    colnames(bh_feats)[1] <- 'chr'
    remove_features <- inner_join(bh_feats, g_ranges)$probe
    
    # take remove features out of colnames 
    bh_features <- intersect_names[!intersect_names %in% remove_features]
    
    # HERE
    # function to predict with all test, controls, controls old, and valid
    mod_result  <- run_enet_all_test(cases_dat = cases_full,
                                    controls_dat = controls_full,
                                    valid_dat = valid_full,
                                    age_cutoff = 72,
                                    gender = gender, 
                                    tech = tech,
                                    bh_features = bh_features)
    
    
    
    # store restults for each result object in its list category
    temp_cases[[i]] <- mod_result[[1]]
    temp_controls[[i]] <- mod_result[[2]]
    temp_models[[i]] <- mod_result[[3]]
    temp_lambda_min[[i]] <- mod_result[[4]]
    temp_lambda_1se[[i]] <- mod_result[[5]]
    temp_alpha[[i]] <- mod_result[[6]]
    
    
    print(i)
  }
  
  # combine list of case and control result data frames and return all result objects (two dataframes and 4 lists)
  full_cases <- do.call(rbind, temp_cases)
  full_controls <- do.call(rbind, temp_controls)
  
  return(list(full_cases, full_controls, temp_models, temp_lambda_min, temp_lambda_1se, temp_alpha))
}


# run model with 5 k fold cross validation
temp_res <- run_model(cases_full = all_cases,
                      controls_full =  all_con,
                      k_folds = k_folds,
                      tech = tech,
                      gender = gender,
                      beta_thresh = beta_thresh,
                      methyl_type = methyl_type,
                      g_ranges = g_ranges)

# store object resuts for further anaylsis (conusion matrix, accuracy, control validation, etc)
temp_1 <- temp_res[[1]]
temp_3 <- temp_res[[3]]
temp_4 <- temp_res[[4]]
temp_5 <- temp_res[[5]]
temp_6 <- temp_res[[6]]




