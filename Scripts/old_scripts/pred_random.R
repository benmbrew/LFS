##########
# initialize libraries
##########
library(caret)
library(glmnet)
library(randomForest)
library(kernlab)
library(pROC)
library(Metrics)
library(doParallel)
library(nnet)
library(dplyr)
library(bumphunter)
library(sqldf)
library(e1071)
library(reshape2)

registerDoParallel(1)
##########
# initialize folders
##########
home_folder <- '~/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
model_data <- paste0(data_folder, '/model_data')

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
method = 'raw'
k = 5

##########
# load data
##########
# read in full m value data 
betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new_m.rda')))
betaControls <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_new_m.rda'))) #34 449936
betaValid <- readRDS(paste0(model_data, paste0('/', method, '_', 'valid_new_m.rda')))
#35 449783

##########
# read in cluster labels
##########
kmeans_lab <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_kmeans_labs.rda')))

###########
# make id into ids
###########
colnames(betaCases)[1] <- 'ids'
colnames(betaControls)[1] <- 'ids'
colnames(betaValid)[1] <- 'ids'

##########
# remove inf
##########
betaCases <- removeInf(betaCases, probe_start = 8)
betaControls <- removeInf(betaControls, probe_start = 8)
betaValid<- removeInf(betaValid, probe_start = 8)


# get old controls - Mut and 'Unaffected'
betaControlsOld <- subset(betaCases, p53_germline == 'Mut' & 
                            cancer_diagnosis_diagnoses == 'Unaffected')

# get p53, not 'Unaffected'
betaCases <- getModData(betaCases)

# get rid of cancer samples in controls 
betaControls <- betaControls[grepl('Unaffected', betaControls$cancer_diagnosis_diagnoses),]

#subset valid
betaValid <- betaValid[!betaValid$ids %in% betaCases$ids,]

##########
# get intersecting colnames and prepare data for modeling
##########

intersect_names <- Reduce(intersect, list(colnames(betaCases)[8:ncol(betaCases)], 
                                          colnames(betaControls)[8:ncol(betaControls)], 
                                          colnames(betaValid)[8:ncol(betaValid)]))
# assign dataframe identifier
betaCases$type <- '450k'
betaControls$type <- '850k'
betaControlsOld$type <- '450k'
betaValid$type <- '850k'


# cases
betaCases <- betaCases[, c('ids',
                           'age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender',
                           'type',
                           intersect_names)]
# controls
betaControls <- betaControls[, c('ids',
                                 'age_diagnosis', 
                                 'age_sample_collection', 
                                 'cancer_diagnosis_diagnoses', 
                                 'gender', 
                                 'type',
                                 intersect_names)]

# controls
betaControlsOld <- betaControlsOld[, c('ids',
                                       'age_diagnosis', 
                                       'age_sample_collection', 
                                       'cancer_diagnosis_diagnoses', 
                                       'gender', 
                                       'type',
                                       intersect_names)]

#validation
betaValid <- betaValid[, c('ids', 
                           'age_diagnosis', 
                           'age_sample_collection', 
                           'cancer_diagnosis_diagnoses', 
                           'gender', 
                           'type',
                           intersect_names)]



# get controls full
betaControlsFull <- rbind(betaControls,
                          betaControlsOld)


# remove duplicates from betaControlsFull
length(which(duplicated(betaControlsFull$ids)))
betaControlsFull <- betaControlsFull[!duplicated(betaControlsFull$ids),]
# #########
# # train and test random
# #########
# 
# scale_data <- function(data_frame, probe_start) {
#   
#   temp_clin <- data_frame[, 1:(probe_start - 1)]
#   temp_scale <- scale(data_frame[, 7:ncol(data_frame)])
#   temp_dat <- cbind(temp_clin, temp_scale)
#   
#   return(temp_dat)
#   
# }
# 
# betaCases <- scale_data(betaCases, probe_start = 7)
# betaValid <- scale_data(betaValid, probe_start = 7)
# betaControls <- scale_data(betaControls, probe_start = 7)
# betaControlsOld <- scale_data(betaControlsOld, probe_start = 7)
# betaControlsFull <- scale_data(betaControlsFull, probe_start = 7)

# # save data
# saveRDS(betaCases, paste0(model_data, paste0('/', method, '_', 'cases_new_m_scaled.rda')))
# saveRDS(betaValid, paste0(model_data, paste0('/', method, '_', 'valid_new_m_scaled.rda')))
# saveRDS(betaControls, paste0(model_data, paste0('/', method, '_', 'controls_new_m_scaled.rda')))
# saveRDS(betaControlsOld, paste0(model_data, paste0('/', method, '_', 'controls_old_new_m_scaled.rda')))
# saveRDS(betaControlsFull, paste0(model_data, paste0('/', method, '_', 'controls_full_m_scaled.rda')))


# # save data
betaCases <-readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new_m_scaled.rda')))
betaValid <- readRDS(paste0(model_data, paste0('/', method, '_', 'valid_new_m_scaled.rda')))
betaControls <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_new_m_scaled.rda')))
betaControlsOld <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_old_new_m_scaled.rda')))
betaControlsFull <- readRDS(paste0(model_data, paste0('/', method, '_', 'controls_full_m_scaled.rda')))

kmeans_lab_scaled <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_kmeans_labs_scaled.rda')))


beta_cases = betaCases
beta_controls = betaControls
beta_controls_old = betaControlsOld
beta_controls_full = betaControlsFull
beta_valid = betaValid
model_method = 'rf'
cluster_feats = kmeans_lab
mod_feats = intersect_names
feature_set = feature_length[l]
seed_num = seeds[j]
subset = 'both'
k = k

trainTest <- function(beta_cases,
                      beta_controls,
                      beta_controls_old,
                      beta_controls_full,
                      beta_valid,
                      model_method,
                      use_clusters,
                      cluster_feats,
                      mod_feats,
                      feature_set,
                      seed_num, 
                      subset, 
                      k) 
{
  
  
  set.seed(seed_num)
  # get a column for each dataset indicating the fold
  beta_cases <- getFolds(beta_cases, seed_number = seed_num, k = k)
  
  # remove samples that dont have an age of sample collection
  beta_cases <- beta_cases[complete.cases(beta_cases),]
  beta_valid <- beta_valid[complete.cases(beta_valid),]
  
  beta_controls <- beta_controls[!is.na(beta_controls$age_sample_collection),]
  beta_controls_old <- beta_controls_old[!is.na(beta_controls_old$age_sample_collection),]
  beta_controls_full <- beta_controls_full[!is.na(beta_controls_full$age_sample_collection),]
  
  # list to store results
  model_results <- list()
  y_results <- list()
  y_results_cases <- list()
  y_results_valid  <- list()
  y_results_controls  <- list()
  y_results_controls_old <- list()
  y_results_controls_full <- list()
  
  # now write forloop to 
  for (i in 1:k) {
    
    # get x 
    train_index <- !grepl(i, beta_cases$folds)
    test_index <- !train_index
    
    if(use_clusters) {
     mod_feats <-  (cluster_feats %>% group_by(count) %>% sample_n(feature_set/10))$probe
    } else {
      mod_feats <- sample(mod_feats, feature_set, replace = T)
      
    }
    
    if(model_method == 'enet') {
      # # get residuals
      # cases_resid <- getResidual(data = cases, 
      #                            bh_features = bh_feat_sig)
      
      
      mod_result <- runEnetRand(training_dat = beta_cases[train_index,], 
                                controls_dat = beta_controls,
                                controls_dat_old = beta_controls_old,
                                controls_dat_full = beta_controls_full,
                                valid_dat = beta_valid,
                                test_dat = beta_cases[test_index,], 
                                bh_features = mod_feats,
                                gender = F)
      
    } 
    
    if(model_method == 'rf') {
      # # get residuals
      # cases_resid <- getResidual(data = cases, 
      #                            bh_features = bh_feat_sig)
      
      
      mod_result <- runRfRand(training_dat = beta_cases[train_index,], 
                                controls_dat = beta_controls,
                                controls_dat_old = beta_controls_old,
                                controls_dat_full = beta_controls_full,
                                valid_dat = beta_valid,
                                test_dat = beta_cases[test_index,], 
                                bh_features = mod_feats,
                                gender = F)
      
    } 
    
    if(model_method == 'svm') {
      # # get residuals
      # cases_resid <- getResidual(data = cases, 
      #                            bh_features = bh_feat_sig)
      
      
      mod_result <- runSvmRand(training_dat = beta_cases[train_index,], 
                                controls_dat = beta_controls,
                                controls_dat_old = beta_controls_old,
                                controls_dat_full = beta_controls_full,
                                valid_dat = beta_valid,
                                test_dat = beta_cases[test_index,], 
                                bh_features = mod_feats,
                                gender = F)
      
    } 
    
    if(model_method == 'lasso') {
      # # get residuals
      # cases_resid <- getResidual(data = cases, 
      #                            bh_features = bh_feat_sig)
      
      
      mod_result <- runLassoRand(training_dat = beta_cases[train_index,], 
                                controls_dat = beta_controls,
                                controls_dat_old = beta_controls_old,
                                controls_dat_full = beta_controls_full,
                                valid_dat = beta_valid,
                                test_dat = beta_cases[test_index,], 
                                bh_features = mod_feats,
                                gender = F)
      
    } 
    
    
    if(model_method == 'ridge') {
      # # get residuals
      # cases_resid <- getResidual(data = cases, 
      #                            bh_features = bh_feat_sig)
      
      
      mod_result <- runRidgeRand(training_dat = beta_cases[train_index,], 
                                controls_dat = beta_controls,
                                controls_dat_old = beta_controls_old,
                                controls_dat_full = beta_controls_full,
                                valid_dat = beta_valid,
                                test_dat = beta_cases[test_index,], 
                                bh_features = mod_feats,
                                gender = F)
      
    } 
    
    
    
    
    # predictions, age of onset, age of sample collection
    temp_cases <- as.data.frame(t(do.call(rbind, list(mod_result[[1]], mod_result[[2]], mod_result[[3]]))))
    temp_valid <- as.data.frame(t(do.call(rbind, list(mod_result[[4]], mod_result[[5]], mod_result[[6]]))))
    temp_controls <- as.data.frame(t(do.call(rbind, list(mod_result[[7]], mod_result[[8]]))))
    temp_controls_full <- as.data.frame(t(do.call(rbind, list(mod_result[[9]], mod_result[[10]]))))
    temp_controls_old <- as.data.frame(t(do.call(rbind, list(mod_result[[11]], mod_result[[12]]))))
    
    # rename 
    colnames(temp_cases) <- 
      colnames(temp_valid) <- 
      c('preds', 'onset', 'age')
    colnames(temp_controls) <- 
      colnames(temp_controls_full) <- 
      colnames(temp_controls_old) <-
      c('preds', 'age')
    
    
    temp_cases$num_feats <- 
      temp_valid$num_feats <-
      temp_controls$num_feats <- 
      temp_controls_old$num_feats <- 
      temp_controls_full$num_feats <- length(mod_feats)
    
    
    temp_cases$model_method <- 
      temp_valid$model_method <-
      temp_controls$model_method <- 
      temp_controls_old$model_method <- 
      temp_controls_full$model_method <- model_method
    
    temp_cases$seed_number <- 
      temp_valid$seed_number <-
      temp_controls$seed_number <- 
      temp_controls_old$seed_number <- 
      temp_controls_full$seed_number <- j
    
    y_results_cases[[i]] <- temp_cases
    y_results_valid[[i]] <-  temp_valid 
    y_results_controls[[i]] <- temp_controls
    y_results_controls_old[[i]] <- temp_controls_old
    y_results_controls_full[[i]] <- temp_controls_full
    
    # mod_result_resid <- runEnet(training_dat = cases_resid[train_index,], 
    #                             test_dat = cases_resid[test_index,], 
    #                             bh_features = bh_feat_sig,
    #                             gender = T)
    
    
    
  }
  
  return(list(y_results_cases, 
              y_results_valid, 
              y_results_controls, 
              y_results_controls_full, 
              y_results_controls_old))
  
}

temp.cases <- list()
temp.valid <- list()
temp.controls <- list()
temp.controls_old <- list()
temp.controls_full <- list()

temp.cases_2 <- list()
temp.valid_2 <- list()
temp.controls_2 <- list()
temp.controls_old_2 <- list()
temp.controls_full_2 <- list()

full.cases <- list()
full.valid <- list()
full.controls <- list()
full.controls_old <- list()
full.controls_full <- list()



# feature_length <- c(100, 200, 300, 400, 500,
#                     600, 700, 800, 900, 1000,
#                     2000, 3000, 4000, 5000, 6000, 
#                     7000, 8000, 9000, 10000,
#                     20000, 30000, 40000, 50000, 100000, 
#                     200000, 300000, 400000)
# seeds <- c(1, 2, 3, 4, 5)
model_names <- c('enet', 'rf')
feature_length <- c(20, 50, 100, 500, 1000, 5000, 10000)
seeds <- c(1, 2, 3)

for(m in 1:length(model_names)) {
  
  for(l in 1:length(feature_length)) {
    
    for (j in 1:length(seeds)) {
      
      mod_results <- trainTest(beta_cases = betaCases,
                               beta_controls = betaControls,
                               beta_controls_old = betaControlsOld,
                               beta_controls_full = betaControlsFull,
                               beta_valid = betaValid,
                               model_method = model_names[m],
                               use_clusters = T,
                               cluster_feats = kmeans_lab_scaled,
                               mod_feats = intersect_names,
                               feature_set = feature_length[l],
                               seed_num = seeds[j],
                               subset = 'both',
                               k = k)
      
      
      temp.cases[[j]] <- do.call(rbind, mod_results[[1]])
      temp.valid[[j]] <- do.call(rbind, mod_results[[2]])
      temp.controls[[j]] <- do.call(rbind, mod_results[[3]])
      temp.controls_full[[j]] <- do.call(rbind, mod_results[[4]])
      temp.controls_old[[j]] <- do.call(rbind, mod_results[[5]])
      
      
      print(paste0('done with ', j, ' seed'))
    }
    temp.cases_2[[l]] <- do.call(rbind, temp.cases)
    temp.valid_2[[l]] <- do.call(rbind, temp.valid)
    temp.controls_2[[l]] <- do.call(rbind, temp.controls)
    temp.controls_full_2[[l]] <- do.call(rbind, temp.controls_full)
    temp.controls_old_2[[l]] <- do.call(rbind, temp.controls_old)
    
    print(paste0('done with ', l, ' random feature set'))
    
  }
  
  full.cases[[m]] <- do.call(rbind, temp.cases_2)
  full.valid[[m]] <- do.call(rbind, temp.valid_2)
  full.controls[[m]] <- do.call(rbind, temp.controls_2)
  full.controls_full[[m]] <- do.call(rbind, temp.controls_full_2)
  full.controls_old[[m]] <- do.call(rbind, temp.controls_old_2)
  
  print(paste0('done with ',model_names[m] , ' models'))
  
}

# get full data for each type
final_cases <- do.call(rbind, full.cases)
final_valid <- do.call(rbind, full.valid)
final_controls <- do.call(rbind, full.controls)
final_controls_old <- do.call(rbind, full.controls_old)
final_controls_full <- do.call(rbind, full.controls_full)

saveRDS(final_cases, paste0(model_data, '/predict_rand_cases_scaled_clust.rda'))
saveRDS(final_valid, paste0(model_data, '/predict_rand_valid_scaled_clust.rda'))
saveRDS(final_controls, paste0(model_data, '/predict_rand_controls_scaled_clust.rda'))
saveRDS(final_controls_old, paste0(model_data, '/predict_rand_controls_old_scaled_clust.rda'))
saveRDS(final_controls_full, paste0(model_data, '/predict_rand_controls_full_scaled_clust.rda'))

# 
final_cases <- readRDS(paste0(model_data, '/predict_rand_cases_scaled.rda'))
final_valid <- readRDS(paste0(model_data, '/predict_rand_valid_scaled.rda'))
final_controls <- readRDS(paste0(model_data, '/predict_rand_controls_scaled.rda'))
final_controls_full <- readRDS(paste0(model_data, '/predict_rand_controls_full_scaled.rda'))
final_controls_old <- readRDS(paste0(model_data, '/predict_rand_controls_old_unscaled.rda'))

final_cases_c <- readRDS(paste0(model_data, '/predict_rand_cases_scaled_clust.rda'))
final_valid_c <- readRDS(paste0(model_data, '/predict_rand_valid_scaled_clust.rda'))
final_controls_c <- readRDS(paste0(model_data, '/predict_rand_controls_scaled_clust.rda'))
final_controls_full_c <- readRDS(paste0(model_data, '/predict_rand_controls_full_scaled_clust.rda'))
final_controls_old_c <- readRDS(paste0(model_data, '/predict_rand_controls_old_scaled_clust.rda'))


# get beta_dim
cases_dim <- nrow(betaCases[complete.cases(betaCases),])
valid_dim <- nrow(betaValid[complete.cases(betaValid),])
controls_dim <- length(betaControls$age_sample_collection[!is.na(betaControls$age_sample_collection)])
controls_old_dim <- length(betaControlsOld$age_sample_collection[!is.na(betaControlsOld$age_sample_collection)])
controls_full_dim <- length(betaControlsFull$age_sample_collection[!is.na(betaControlsFull$age_sample_collection)])


# remove unneccesarry obejcs 
rm(list = ls(pattern = "*beta")) 

##########
# function that takes results and puts in plot form
##########

data_type <- 
get_plot_data <-
  
  function(final_dat, data_type) {
    

  # get row number and y length
  row_num <- nrow(final_dat)
  
  if(data_type == 'cases'){
    y_length <- nrow(subset(final_dat, model_method == 'enet' & seed_number == 1 & num_feats == 100))
    
    #order by seed_number, num_feats, onset
    final_dat <- final_dat[order(final_dat$seed_number, final_dat$model_method, final_dat$num_feats, final_dat$onset, final_dat$age),]
    
    # get row index
    final_dat$row_index <- rep.int(seq(1, y_length, 1), row_num/y_length)
    
    # group by feats, onset, age
    results <- 
      final_dat %>%
      group_by(num_feats, model_method,onset, age, row_index) %>%
      summarise(mean_preds = mean(preds),
                counts = n())
    
  } else {
    
    y_length <- nrow(subset(final_dat, model_method == 'enet' & seed_number == 1 & num_feats == 100))/5
    
    #order by seed_number, num_feats, onset
    final_dat <- final_dat[order(final_dat$seed_number, final_dat$model_method, final_dat$num_feats),]
    
    # get row index
    final_dat$row_index <- rep.int(seq(1, y_length, 1), row_num/y_length)
    
    if (data_type == 'controls') {
      # group by feats, onset, age
      results <- 
        final_dat %>%
        group_by(num_feats,model_method, age, row_index, seed_number) %>%
        summarise(mean_preds = mean(preds),
                  counts = n())
      
      # group by feats, onset, age
      results <- 
        results %>%
        group_by(num_feats, model_method,age, row_index) %>%
        summarise(mean_preds = mean(mean_preds),
                  counts = n())
      

      stopifnot(all(results$counts == length(seeds)))
      
    } else {
      
      # group by feats, onset, age
      results <- 
        final_dat %>%
        group_by(num_feats, model_method, onset, age, row_index, seed_number) %>%
        summarise(mean_preds = mean(preds),
                  counts = n())
      
      # group by feats, onset, age
      results <- 
        results %>%
        group_by(num_feats, model_method, onset, age, row_index) %>%
        summarise(mean_preds = mean(mean_preds),
                  counts = n())
      
      stopifnot(all(results$counts == length(seeds)))
      
      
    }
   
  }
  return(results)
 
}




# get result matrices
results_cases <- get_plot_data(final_cases, data_type = 'cases')
results_valid <- get_plot_data(final_valid, data_type = 'valid')
results_controls <- get_plot_data(final_controls, data_type = 'controls')
results_controls_old <- get_plot_data(final_controls_old, data_type = 'controls')
results_controls_full <- get_plot_data(final_controls_full, data_type = 'controls')

# get result matrices
results_cases_c <- get_plot_data(final_cases_c, data_type = 'cases')
results_valid_c <- get_plot_data(final_valid_c, data_type = 'valid')
results_controls_c <- get_plot_data(final_controls_c, data_type = 'controls')
results_controls_old_c <- get_plot_data(final_controls_old_c, data_type = 'controls')
results_controls_full_c <- get_plot_data(final_controls_full_c, data_type = 'controls')


# make sure they are right length
stopifnot(nrow(results_cases_c)/length(feature_length)/length(model_names) == 77)
stopifnot(nrow(results_valid)/length(feature_length)/length(model_names)  == 37)
stopifnot(nrow(results_controls)/length(feature_length)/length(model_names)  == 30)
stopifnot(nrow(results_controls_old)/length(feature_length)/length(model_names)  == 24)
stopifnot(nrow(results_controls_full)/length(feature_length)/length(model_names)  ==44)

##########
# combine all data
##########

# add indicator for data set 
results_cases$type <- 'cases'
results_valid$type <- 'valid'
results_controls$type <- 'controls_new'
results_controls_old$type <- 'controls_old'
results_controls_full$type <- 'controls_full'

# add age of diagnosis to result_controls*
results_controls$onset <- NA 
results_controls_full$onset <- NA
results_controls_old$onset <- NA


# add indicator for data set 
results_cases_c$type <- 'cases'
results_valid_c$type <- 'valid'
results_controls_c$type <- 'controls_new'
results_controls_old_c$type <- 'controls_old'
results_controls_full_c$type <- 'controls_full'

# add age of diagnosis to result_controls*
results_controls_c$onset <- NA 
results_controls_full_c$onset <- NA
results_controls_old_c$onset <- NA

options(scipen=999)

# comibne
results_all <- rbind(results_cases,
                     results_valid,
                     results_controls,
                     results_controls_full,
                     results_controls_old)

# comibne
results_all_c <- rbind(results_cases_c,
                     results_valid_c,
                     results_controls_c,
                     results_controls_full_c,
                     results_controls_old_c)

##########
# function that takes all results and plot (add in a fitted controls model at some point)
##########

plot_rand_point <- function(results_data, 
                            model_type,
                            remove_feat_over,
                            remove_sample_over,
                            data_type,
                            plot_type) {
  
  results_data <- results_data[results_data$model_method == model_type,]
  results_data <- results_data[results_data$num_feats < remove_feat_over,]
  results_data <- results_data[results_data$age < remove_sample_over,]

  if(data_type == 'cases'){
    
    temp_data <- subset(results_data, type == data_type)
    
    num_samples <- nrow(temp_data)/length(unique(temp_data$num_feats))
    
    
  }else if(data_type == 'valid') {
    
    temp_data <- subset(results_data, type == data_type)
    
    num_samples <- nrow(temp_data)/length(unique(temp_data$num_feats))
    
    
  }else if(data_type == 'controls_new') {
    
    temp_data <- subset(results_data, type == data_type)
    
    num_samples <- nrow(temp_data)/length(unique(temp_data$num_feats))
    
    
  }else if(data_type == 'controls_old') {
    
    temp_data <- subset(results_data, type == data_type)
    
    num_samples <- nrow(temp_data)/length(unique(temp_data$num_feats))
    
    
  }else if(data_type == 'controls_full') {
    
    temp_data <- subset(results_data, type == data_type)
    
    num_samples <- nrow(temp_data)/length(unique(temp_data$num_feats))
    
    
  } else {
    temp_data <- results_data
    
    num_samples <- nrow(temp_data)
    
    # remove counts and row_index
    temp_data$counts <-
      temp_data$row_index <- NULL
    
    # group by features and get means
    temp_mean <- 
      temp_data %>%
      group_by(num_feats, type) %>%
      summarise(mean_onset = mean(onset),
                mean_age = mean(age),
                mean_pred = mean(mean_preds))
    
    # melt it!
    temp_melt <- melt(temp_mean, id.vars = c('num_feats', 'type'))
    temp_melt$num_feats <- as.numeric(temp_mean$num_feats)
    
    # plot difference
    p <- ggplot(temp_melt, aes(num_feats, value, 
                          colour = type,
                          shape = variable,
                          group = interaction(variable, type))) +
      geom_point(alpha = 0.7, size = 5)  + geom_line(size = 2, alpha = 0.6) +
      xlab('Number of features') + ylab('Patient age (in months)') +
      scale_colour_manual(name = 'Dataset',
                          breaks = c('cases', 'valid', 'controls_new', 'controls_old', 'controls_full'),
                          labels = c('Cases', 'Validation', 'Controls', 'Controls 450k', 'Combined Controls'),
                          values = c('darkred', 'darkblue', 'orange', 'grey', 'lightblue')) +
      scale_shape_manual(name = 'Age variable',
                          breaks = c('mean_onset', 'mean_age', 'mean_pred'),
                          labels = c('Avg onset', 'Avg age', 'Avg predictions'),
                          values = c(15,17,18)) +
      ggtitle(paste0(data_type,'_', num_samples)) +
      theme_bw() + theme(axis.text.y = element_text(size = 12),
                         axis.text.x = element_text(size = 12)) 
    

    return(p)
  }
  
  if(plot_type == 'diff') {
    
    # get absolute value of difference
    temp_data$onset_diff <- abs(temp_data$onset - temp_data$mean_preds)
    temp_data$age_diff <- abs(temp_data$age - temp_data$mean_preds)
  
    temp_result <- temp_data[, c('num_feats', 'onset_diff', 'age_diff')]
    # group by features and get means
    temp_mean <- 
      temp_result %>%
      group_by(num_feats) %>%
      summarise(mean_onset_diff = mean(onset_diff),
                mean_age_diff = mean(age_diff))
    
    temp_melt <- melt(temp_mean, id.vars = 'num_feats')
    temp_melt$num_feats <- as.numeric(temp_melt$num_feats)
    
    # plot difference
    ggplot(temp_melt, aes(num_feats, value, group = variable, colour = variable)) +
      geom_point(alpha = 0.7, size = 3)  + geom_line(size = 2, alpha = 0.7) +
      xlab('Number of features') + ylab('Age errors (AV) (in months)') + 
      ggtitle(paste0(data_type,'_', num_samples)) +
      scale_colour_manual(name = '',
                          breaks = c('mean_onset_diff', 'mean_age_diff'),
                          labels = c('Mean onset diff', 'mean age diff'),
                          values = c('darkred', 'darkblue')) +
      theme_bw() + theme(axis.text.y = element_text(size = 12),
                         axis.text.x = element_text(size = 12)) 
    
    
  } else {
    ###########
    # plot raw predicted and onset and age against number of features 
    # group by features and get means
    temp_mean <- 
      temp_data %>%
      group_by(num_feats) %>%
      summarise(mean_onset = mean(onset),
                mean_age = mean(age),
                mean_pred = mean(mean_preds))
    
    # melt it!
    temp_melt <- melt(temp_mean, id.vars = 'num_feats')
    temp_melt$num_feats <- as.numeric(temp_mean$num_feats)
    
    # plot difference
    ggplot(temp_melt, aes(num_feats, value, group = variable, colour = variable)) +
      geom_point(alpha = 0.7, size = 3) + geom_line(size = 2, alpha=0.7) +
      xlab('Number of features') + ylab('Patient age (in months)') +
      scale_colour_manual(name = '',
                          breaks = c('mean_onset', 'mean_age', 'mean_pred'),
                          labels = c('Onset avg', 'Age avg', 'Pred avg'),
                          values = c('darkred', 'darkblue', 'black')) +
      ggtitle(paste0(data_type,'_', num_samples)) +
      theme_bw() + theme(axis.text.y = element_text(size = 12),
                         axis.text.x = element_text(size = 12)) 
    
  }
  
}


##########
# function for within each feature length
##########

# histrogram comparison
# also plot age against pred diff 
# (maybe we predict well for younger patients or maybe avg are driven by very old people.)

plot_by_feat <- 
  function(results_data,
           model_type,
           plot_type, 
           get_feat, 
           keep_type){

  # get model type
  results_data <- results_data[results_data$model_method == model_type,]  
    
  # get feature set
  sub_dat <- results_data[results_data$num_feats == get_feat, ]

  # keep data types, but first paste them together
  sub_dat <- sub_dat[grepl(keep_type, sub_dat$type),] 
  
  # get rid of unnedded columns
  sub_dat$counts <-
    sub_dat$num_feats <-
    sub_dat$row_index <- NULL
  
  if(plot_type == 'hist') {
    
    sub_dat$model_method <- NULL
    title_name <- paste0(get_feat, ' features', ' ', keep_type)
    
    sub_melt <- melt(sub_dat, id.vars = 'type')
    title_name <- paste0(get_feat, ' features')
    
    ggplot(sub_melt, aes(value, fill=variable)) + 
      geom_histogram(bins = 25, alpha=0.6, position="identity", colour = 'grey') +
      scale_fill_manual(name = '',
                        breaks = c('onset', 'age', 'mean_preds'),
                        labels = c('age of diagnosis', 'age', 'predictions'),
                        values = c('darkblue', 'darkred', 'grey')) +
      theme_bw() + theme(axis.text.y = element_text(size = 12),
                         axis.text.x = element_text(size = 12)) + ggtitle(title_name)
    
  } else if(plot_type == 'pred_diff'){
    
    # plot age (sample collection) againse abs difference (error)
    # get abs difference variable (error)
    if(grepl('controls', keep_type)){
      sub_dat$pred_diff <- sub_dat$mean_preds - sub_dat$age
      
    } else {
      sub_dat$pred_diff <- sub_dat$mean_preds - sub_dat$onset
      
    }
    title_name <- paste0(get_feat, ' features', ' ', keep_type)
    
    # plot difference
    ggplot(sub_dat, aes(age, pred_diff)) +
      geom_point(alpha = 0.7, size = 3, colour = 'red') + 
      xlab('Age of patient (months)') + ylab('Error') +
      ggtitle(title_name) +
      theme_bw() + theme(axis.text.y = element_text(size = 12),
                         axis.text.x = element_text(size = 12)) 
    
    
  } else {
    # plot age against preds, color for onset and age
    sub_dat$type <- NULL
    sub_melt <- melt(sub_dat, id.vars = 'mean_preds')
    title_name <- paste0(get_feat, ' features', ' ', keep_type)
    
    x_max <- max(sub_dat$age)
    y_max <- max(sub_dat$mean_preds)
    
    max_num = max(x_max, y_max)
    
    ggplot(sub_melt, aes(value, mean_preds, colour=variable)) + 
      geom_point(alpha=0.6, size = 3) + 
      xlab('Age') + ylab('Predictions') + 
      xlim(c(0, max_num)) + ylim(c(0, max_num)) +
      scale_colour_manual(name = '',
                        breaks = c('onset', 'age'),
                        labels = c('age of diagnosis', 'age'),
                        values = c('darkblue', 'darkred')) +
      theme_bw() + theme(axis.text.y = element_text(size = 12),
                         axis.text.x = element_text(size = 12)) + ggtitle(title_name) +
      geom_abline(intercept = 0, slope = 1)
    
  }
  

}
summary(as.factor(results_all$num_feats))

# hist(results_all$age[results_all$type == 'cases'])
plot_rand_point(results_data = results_all,
                model_type = 'enet',
                remove_feat_over = 500000, 
                remove_sample_over = 1000000,
                data_type = 'controls_new', 
                plot_type = 'diff')

# histograms 
plot_by_feat(results_data = results_all_c, 
             model_type = 'enet',
             plot_type = 'hist', 
             get_feat = 100, 
             keep_type = 'controls_new')

#########
# function to combine cluster and normal and plot correlations
#########

# give each data set an indicator 
results_all$dat <- 'normal'
results_all_c$dat <- 'clust'

# cbind them
results_full <- rbind(results_all,
                      results_all_c)

names(results_full)
#########
# function that plots correlation by type against number of features
#########
# results_data <- results_full
# data_type <- 'cases'
# model_type = 'enet'
# remove_feat_over = 100000
# remove_sample_over = 100000
get_corr_plots <- function(results_data, 
                           data_type,
                           model_type,
                           remove_feat_over,
                           remove_sample_over){
  
  results_data <- results_data[results_data$model_method == model_type,]
  results_data <- results_data[results_data$num_feats < remove_feat_over,]
  results_data <- results_data[results_data$age < remove_sample_over,]
  
  if(data_type == 'cases'){
    
    temp_data <- subset(results_data, type == data_type)
    
    num_samples <- nrow(temp_data)/length(unique(temp_data$num_feats))
    
    
  }else if(data_type == 'valid') {
    
    temp_data <- subset(results_data, type == data_type)
    
    num_samples <- nrow(temp_data)/length(unique(temp_data$num_feats))
    
    
  }else if(data_type == 'controls_new') {
    
    temp_data <- subset(results_data, type == data_type)
    
    num_samples <- nrow(temp_data)/length(unique(temp_data$num_feats))
    
    
  }else if(data_type == 'controls_old') {
    
    temp_data <- subset(results_data, type == data_type)
    
    num_samples <- nrow(temp_data)/length(unique(temp_data$num_feats))
    
    
  }else if(data_type == 'controls_full') {
    
    temp_data <- subset(results_data, type == data_type)
    
    num_samples <- nrow(temp_data)/length(unique(temp_data$num_feats))
    
  }
  
  temp_mean <- 
    temp_data %>%
    group_by(num_feats, dat) %>%
    summarise(onset_cor = cor(onset, mean_preds),
              age_cor = cor(age, mean_preds))
  
  temp_melt <- melt(temp_mean, id.vars = c('num_feats', 'dat'))
  temp_melt$num_feats <- as.numeric(temp_melt$num_feats)
  # plot difference
  p <- ggplot(temp_melt, aes(num_feats, value, 
                             colour = dat,
                             shape = variable,
                             group = interaction(variable, dat))) +
    geom_point(alpha = 0.7, size = 5)  + geom_line(size = 2, alpha = 0.6) +
    xlab('Number of features') + ylab('Patient age (in months)') +
    scale_colour_manual(name = 'Dataset',
                        breaks = c('normal', 'clust'),
                        labels = c('Normal', 'Clusters'),
                        values = c('darkred', 'darkblue')) +
    scale_shape_manual(name = 'Age variable',
                       breaks = c('age_cor', 'onset_cor'),
                       labels = c('Corr w/ onset', 'Corr w/ age'),
                       values = c(15,17,18)) +
    ggtitle(paste0(data_type,'_', num_samples/2)) +
    theme_bw() + theme(axis.text.y = element_text(size = 12),
                       axis.text.x = element_text(size = 12)) 
  
  
  return(p)
  
  
}

get_corr_plots(results_full, 
               data_type = 'cases', 
               model_type = 'enet', 
               remove_feat_over = 100000, 
               remove_sample_over = 100000)


