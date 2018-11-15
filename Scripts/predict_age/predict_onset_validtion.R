# this script will be used to reproduce the original held out validation set for true accuracy. 
# can also test on controls here

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
all_con_wt <- all_con_wt[!is.na(all_con_wt$age_sample_collection),]
all_con_wt <- all_con_wt[!is.na(all_con_wt$gender),]

# ge tgender 
all_con_wt <- cbind(as.data.frame(class.ind(all_con_wt$gender)), 
                 all_con_wt)

# rempove old tech variable 
all_con_wt$gender <- NULL

# tech
all_con_wt <- cbind(as.data.frame(class.ind(all_con_wt$tech)), 
                 all_con_wt)
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


# run bumphunter on LFS healthy patients (LFS no cancer) and LFS cancer patients (LFS cancer)
bh_feats <- bump_hunter(dat_1 = con_wt, 
                        dat_2 = con_mut, 
                        bump = 'lfs', 
                        boot_num = 5, 
                        beta_thresh = 0.05,
                        methyl_type = methyl_type,
                        g_ranges = g_ranges)

# cases
cases_450 <- join_new_features(cases_450, new_features = bh_feats)
cases_850 <- join_new_features(cases_850, new_features = bh_feats)

# controls
con_450 <- join_new_features(con_450, new_features = bh_feats)
con_850 <- join_new_features(con_850, new_features = bh_feats)

# wt
con_mut <- join_new_features(con_mut, new_features = bh_feats)
con_wt <- join_new_features(con_wt, new_features = bh_feats)

# store the selected featufres
lfs_bump_probes <- colnames(cases_450)[grepl('^cg', colnames(cases_450))]


rm(ratio_set, bh_feats)
# save image
# load('~/Desktop/temp_valid.RData')


# add dummy tech variable for data sets with only one, replace family_name
names(con_450)[10] <- 'batch_2'
names(con_850)[10] <- 'batch_1'
names(cases_450)[10] <- 'batch_2'
names(cases_850)[10] <- 'batch_1'

# fill them with Zero
con_450$batch_2 <- 0
con_850$batch_1 <- 0
cases_450$batch_2 <- 0
cases_850$batch_1 <- 0

# remove family_name from con_mut and con_wt to even the column names
con_mut$family_name <- con_wt$family_name <- NULL

# rbind all controls
controls_all <- rbind(con_450,
                      con_850,
                      con_mut,
                      con_wt)

rm(con_450,
   con_850,
   con_mut, 
   con_wt)

# create range for random sampling for cross validation
seed_range <- c(1:5)

# create list to store model
all_test_results <- list()
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
      
      
      train_cases <- cases_full[train_index, ]
      test_cases <- cases_full[test_index, ]
      
      # get controls 450k mut
      temp_con <- controls_all[controls_all$batch_1 == 1 & controls_all$p53_germline == 'MUT',]
      
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
                                       controls_dat = controls_full,
                                       valid_dat = valid_full,
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
                                               controls_full = controls_all,
                                               valid_full = cases_850,
                                               k_folds = 5,
                                               tech = FALSE,
                                               gender = TRUE,
                                               beta_thresh = 0.05,
                                               g_ranges = g_ranges)
  
  message('finished working on random seed = ', random_seed)
  
  
}

final_dat <- do.call(rbind, all_test_results)

# # save data
# saveRDS(final_dat, paste0('validation_age_predictions/', which_methyl, '_', which_combat, '_', 
#                           num_seeds, '_', num_folds, '_', is_tech, '_', is_gen, '.rda'))

# read in data
final_dat <- readRDS(paste0('validation_age_predictions/', which_methyl, '_', which_combat, '_', 
                            num_seeds, '_', num_folds, '_', is_tech, '_', is_gen, '.rda'))
# get the average over seeds
final_dat$seed <- as.character(final_dat$seed)
dat <- final_dat %>% group_by(seed) %>% summarise(mean_acc = mean(accuracy, na.rm = TRUE),
                                                  mean_alpha = mean(alpha, na.rm = TRUE),
                                                  mean_lambda = mean(lambda, na.rm = TRUE))

# rank by accuracy 
final_dat <- final_dat[order(final_dat$accuracy, decreasing = TRUE),]

# get bets preds
best_preds <- final_dat[final_dat$accuracy > .89,]
# best alpha ~.5, best lambda ~.25

##########
# use the alpha and lambda means to train on all cases and test on all controls and valid
##########
cases <- cases_450
controls <- controls_all
valid <- cases_850
age_cutoff <- 72
gender = TRUE
tech = FALSE
train_lambda = TRUE
alpha_value <- 0.5
lambda_value <- 0.25
bh_features <- lfs_bump_probes
test_model <- function(cases, 
                       controls, 
                       valid, 
                       gender,
                       tech,
                       age_cutoff,
                       alpha_value, 
                       lambda_value,
                       bh_features) {
  
  # get intersection of bh features and real data
  bh_features <- as.character(unlist(bh_features))
  
  intersected_feats <- intersect(bh_features, colnames(cases))
  
  if(gender) {
    intersected_feats <- c('Female', 'Male', intersected_feats)
  }
  if (tech) {
    intersected_feats <- c('batch_1', 'batch_2', intersected_feats)
  }
  
  
  
  # get outcomes
  cases_y <- as.factor(ifelse(cases$age_diagnosis < age_cutoff, 'positive', 'negative'))
  valid_y <- as.factor(ifelse(valid$age_diagnosis < age_cutoff, 'positive', 'negative')) 
  controls_y <-  as.factor(ifelse(controls$age_sample_collection < age_cutoff, 'positive', 'negative'))
  controls_clin <- controls[, !grepl('^cg', colnames(controls))]
  valid_clin <- valid[, !grepl('^cg', names(valid))]
  
  # get model data
  cases <- cases[, intersected_feats]
  controls <- controls[, intersected_feats]
  valid <- valid[, intersected_feats]
  
  # store fixed values
  best_alpha <- alpha_value
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  
  # set parameters for training model
  type_family <- 'binomial'
  type_measure <- 'auc'
  nfolds = 5
  
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  while (temp.non_zero_coeff < 1) { 
    elastic_net.cv_model = cv.glmnet(x = as.matrix(cases)
                                     , y =  cases_y
                                     , alpha = alpha_value
                                     , type.measure = type_measure
                                     , family = type_family
                                     , standardize=FALSE
                                     , nlambda = 100
                                     , nfolds = nfolds
                                     , parallel = TRUE
    )
    
    
   
    # get outcome variables and clin variables
    temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) 
    trained_lambda_index_set = elastic_net.cv_model$lambda[elastic_net.cv_model$lambda < 0.26 & elastic_net.cv_model$lambda > 0.24]
    trained_lambda_value <- trained_lambda_index_set[1]
    trained_lambda_index <- which(elastic_net.cv_model$lambda == trained_lambda_value) 
    
    
    # # number of non zero coefficients at that lambda    
    temp.non_zero_coeff = elastic_net.cv_model$nzero[temp.min_lambda_index] 
    temp.loop_count = temp.loop_count + 1
    
    # set seed for next loop iteration
    as.numeric(Sys.time())-> t 
    set.seed((t - floor(t)) * 1e8 -> seed) 
    if (temp.loop_count > 10) {
      print("diverged")
      temp.min_lambda_index = 50 # if it loops more than 5 times, then model did not converge
      break
    }
  }# while loop ends 
  # print(temp.non_zero_coeff)  
  
  model  = glmnet(x = as.matrix(cases)
                  , y =  cases_y
                  ,alpha = best_alpha
                  ,standardize=FALSE
                  ,nlambda = 100
                  ,family = type_family)
  
  
  # Predictions on controls data
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions_con <- predict(model, 
                                       data.matrix(controls),
                                       type = 'response')
  
  if(train_lamda){
    # get predictions with corresponding lambda.
    test.predictions_con <- temp_test.predictions_con[, trained_lambda_index]
    
  } else {
    # get predictions with corresponding lambda.
    test.predictions_con <- temp_test.predictions_con[, temp.min_lambda_index]
    
  }
  
  # combine predictions and real labels 
  temp_dat_con <- as.data.frame(cbind(controls_age_pred = test.predictions_con, controls_age_label = controls_y, controls_clin))
  
  
  
  # Predictions on validation data
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions_valid <- predict(model, 
                                         data.matrix(valid),
                                         type = 'response')
  
  if(train_lamda){
    # get predictions with corresponding lambda.
    test.predictions_valid <- temp_test.predictions_valid[, trained_lambda_index]
    
  } else {
    # get predictions with corresponding lambda.
    test.predictions_valid <- temp_test.predictions_valid[, temp.min_lambda_index]
    
  }

  # combine predictions and real labels 
  temp_dat_valid <- as.data.frame(cbind(valid_age_pred = test.predictions_valid, vallid_age_label = valid_y, valid_clin))
  
  
  return(list(temp_dat_con, temp_dat_valid))
  
  
}





