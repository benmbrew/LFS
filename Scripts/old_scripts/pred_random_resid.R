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
scaled = T

##########
# load data
##########
if(scaled) {
  # read in full m value data 
  betaFull <- readRDS(paste0(model_data, paste0('/', method, '_', 'full_new_m_scaled_resid.rda')))
  betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new_m_scaled_resid.rda')))
  betaValid <- readRDS(paste0(model_data, paste0('/', method, '_', 'valid_new_m_scaled_resid.rda')))
  #35 449783
} else {
  # read in full m value data 
  betaFull <- readRDS(paste0(model_data, paste0('/', method, '_', 'full_new_m_resid.rda')))
  betaCases <- readRDS(paste0(model_data, paste0('/', method, '_', 'cases_new_m_resid.rda')))
  betaValid <- readRDS(paste0(model_data, paste0('/', method, '_', 'valid_new_m_resid.rda')))
  #35 449783
}


# create lists to store model results
temp.cases <- list()
temp.cases_2 <- list()
full.cases <- list()

# set parameters to loop through
model_names <- c('enet', 'rf', 'svm')
feature_length <- c(100, 1000, 5000, 10000, 50000, 100000)
seeds <- c(1)

# beta_cases = betaCases
# model_method = 'rf'
# mod_feats = colnames(betaCases)[7:ncol(betaCases)]
# feature_set = feature_length[l]
# seed_num = seeds[j]
# k = k

trainTestResid <- function(beta_cases,
                      model_method,
                      mod_feats,
                      feature_set,
                      seed_num, 
                      k) 
{
  
  
  set.seed(seed_num)
  # get a column for each dataset indicating the fold
  beta_cases <- getFolds(beta_cases, seed_number = seed_num, k = k)
  
  # remove samples that dont have an age of sample collection
  beta_cases <- beta_cases[complete.cases(beta_cases),]
  
  # list to store results
  y_results_cases <- list()
  
  
  # now write forloop to 
  for (i in 1:k) {
    
    # get x 
    train_index <- !grepl(i, beta_cases$folds)
    test_index <- !train_index
    
    
    if(model_method == 'enet') {
      # # get residuals
      # cases_resid <- getResidual(data = cases, 
      #                            bh_features = bh_feat_sig)
      
      mod_feats <- sample(mod_feats, feature_set, replace = T)
      
      mod_result <- runEnetRandResid(training_dat = beta_cases[train_index,], 
                                     test_dat = beta_cases[test_index,], 
                                     bh_features = mod_feats,
                                     gender = F)
      
    } 
    
    if(model_method == 'rf') {
      # # get residuals
      # cases_resid <- getResidual(data = cases, 
      #                            bh_features = bh_feat_sig)
      
      mod_feats <- sample(mod_feats, feature_set, replace = T)
      
      mod_result <- runRfRandResid(training_dat = beta_cases[train_index,], 
                                   test_dat = beta_cases[test_index,], 
                                   bh_features = mod_feats,
                                   gender = F)
      
    } 
    
    if(model_method == 'svm') {
      # # get residuals
      # cases_resid <- getResidual(data = cases, 
      #                            bh_features = bh_feat_sig)
      
      mod_feats <- sample(mod_feats, feature_set, replace = T)
      
      mod_result <- runSvmRandResid(training_dat = beta_cases[train_index,], 
                                    test_dat = beta_cases[test_index,], 
                                    bh_features = mod_feats,
                                    gender = F)
      
    } 
    
    if(model_method == 'lasso') {
      # # get residuals
      # cases_resid <- getResidual(data = cases, 
      #                            bh_features = bh_feat_sig)
      
      mod_feats <- sample(mod_feats, feature_set, replace = T)
      
      mod_result <- runLassoRandResid(training_dat = beta_cases[train_index,], 
                                      test_dat = beta_cases[test_index,], 
                                      bh_features = mod_feats,
                                      gender = F)
      
    } 
    
    
    if(model_method == 'ridge') {
      # # get residuals
      # cases_resid <- getResidual(data = cases, 
      #                            bh_features = bh_feat_sig)
      
      mod_feats <- sample(mod_feats, feature_set, replace = T)
      
      mod_result <- runRidgeRandResid(training_dat = beta_cases[train_index,], 
                                      test_dat = beta_cases[test_index,], 
                                      bh_features = mod_feats,
                                      gender = F)
      
    } 
    
    
    # predictions, age of onset, age of sample collection
    temp_cases <- as.data.frame(t(do.call(rbind, list(mod_result[[1]], mod_result[[2]], mod_result[[3]]))))
  
    
    # rename 
    colnames(temp_cases) <-  c('preds', 'onset', 'age')
    
    
    temp_cases$num_feats <- length(mod_feats)
    
    
    temp_cases$model_method <- model_method
    
    temp_cases$seed_number <- j
    
    y_results_cases[[i]] <- temp_cases
    
    
    
  }
  
  return(y_results_cases)
  
}


for(m in 1:length(model_names)) {
  
  for(l in 1:length(feature_length)) {
    
    for (j in 1:length(seeds)) {
      
      # or do betaFull here
      mod_results <- trainTestResid(beta_cases = betaCases,
                               model_method = model_names[m],
                               mod_feats = colnames(betaCases)[7:ncol(betaCases)],
                               feature_set = feature_length[l],
                               seed_num = seeds[j],
                               k = k)
      
      
      temp.cases[[j]] <- do.call(rbind, mod_results)
      
      
      print(paste0('done with ', j, ' seed'))
    }
    temp.cases_2[[l]] <- do.call(rbind, temp.cases)
    
    print(paste0('done with ', l, ' random feature set'))
    
  }
  
  full.cases[[m]] <- do.call(rbind, temp.cases_2)
  
  print(paste0('done with ',model_names[m] , ' models'))
  
}

# get full data for each type
final_cases <- do.call(rbind, full.cases)


saveRDS(final_cases, paste0(model_data, '/predict_rand_cases_resid_scaled.rda'))

final_cases <- readRDS(paste0(model_data, '/predict_rand_cases_resid_scaled.rda'))



# get beta_dim
cases_dim <- nrow(betaCases[complete.cases(betaCases),])
full_dim <- nrow(betaFull[complete.cases(betaFull),])

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


# make sure they are right length
stopifnot(nrow(results_cases)/length(feature_length)/length(model_names) == cases_dim)

##########
# combine all data
##########

options(scipen=999)

##########
# function that takes all results and plot (add in a fitted controls model at some point)
##########
# results_data <- results_cases
# model_type = 'rf'
# remove_feat_over = 4000000
# remove_sample_over = 5000000
plot_rand_point <- function(results_data, 
                            model_type,
                            remove_feat_over,
                            remove_sample_over,
                            data_type,
                            plot_type) {
  
  results_data <- results_data[results_data$model_method == model_type,]
  results_data <- results_data[results_data$num_feats < remove_feat_over,]
  results_data <- results_data[results_data$age < remove_sample_over,]
  
  num_samples <- nrow(results_data)/6
  

  if(plot_type == 'diff') {
    
    temp_data <- results_data
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
    
    temp_data <- results_data
    
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


# hist(results_all$age[results_all$type == 'cases'])
plot_rand_point(results_data = results_cases,
                model_type = 'rf',
                remove_feat_over = 500000, 
                remove_sample_over = 500000,
                data_type = 'cases', 
                plot_type = 'diff')


# histograms 
plot_by_feat(results_data = results_cases, 
             model_type = 'rf',
             plot_type = 'hist', 
             get_feat = 100000, 
             keep_type = 'cases')




get_corr_plots <- function(results_data, 
                           model_type,
                           remove_feat_over,
                           remove_sample_over){
  
  results_data <- results_data[results_data$model_method == model_type,]
  results_data <- results_data[results_data$num_feats < remove_feat_over,]
  results_data <- results_data[results_data$age < remove_sample_over,]
  
  temp_data <- results_data
  
  num_samples <- nrow(temp_data)/length(unique(temp_data$num_feats))
  
  temp_mean <- 
    temp_data %>%
    group_by(num_feats) %>%
    summarise(onset_cor = cor(onset, mean_preds),
              age_cor = cor(age, mean_preds))
  
  temp_melt <- melt(temp_mean, id.vars = c('num_feats'))
  temp_melt$num_feats <- as.numeric(temp_melt$num_feats)
  # plot difference
  p <- ggplot(temp_melt, aes(num_feats, value, 
                             colour = variable,
                             group = variable)) +
    geom_point(alpha = 0.7, size = 5)  + geom_line(size = 2, alpha = 0.6) +
    xlab('Number of features') + ylab('Patient age (in months)') +
    scale_shape_manual(name = 'Age variable',
                       breaks = c('age_cor', 'onset_cor'),
                       labels = c('Corr w/ onset', 'Corr w/ age'),
                       values = c('blue', 'red')) +
    ggtitle(paste0(num_samples/2)) +
    theme_bw() + theme(axis.text.y = element_text(size = 12),
                       axis.text.x = element_text(size = 12)) 
  
  
  return(p)
  
  
}

get_corr_plots(results_data = results_cases, 
               model_type = 'enet', 
               remove_feat_over = 100000, 
               remove_sample_over = 100000)


