# source functions script
source('helper_functions.R')

# create fixed objects to model and pipeline inputs and saving  
methyl_type = 'beta'
data_set  = 'valid'
gender = FALSE
tech = FALSE 
how_many_seeds = 10
how_many_folds = 5
# use_offset = FALSE
remove_age  = TRUE
beta_thresh = 0.05


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

# if(use_offset){
#   is_offset <- 'use_offset'
# } else {
#   is_offset <- 'no_offset'
# }

num_seeds <- paste0('seeds_', how_many_seeds)
num_folds <- paste0('folds_', how_many_folds)
k_folds <- how_many_folds

# saveRDS(temp, paste0('transform_data/', 'valid_test_transform',which_methyl, '_', is_gen, '_', num_seeds,'.rda'))
dat <- readRDS(paste0('transform_data_test/', data_set,'_test_transform',which_methyl, '_', is_gen, '.rda'))

if(data_set == 'valid'){
  # subset data
  dat <- dat[, c('tm_donor','cancer_diagnosis_diagnoses','valid_age_pred', 'valid_age_label', 'age_sample_collection',
                  'alpha', 'non_zero')]
  
  # remove duplicates
 
  temp <- get_acc_val(dat, thresh = 0.5)
  temp <- temp[order(temp$acc, decreasing = TRUE),]
  temp_dedup <- temp[!duplicated(temp$tm_donor),]
  max_alpha = temp_dedup$alpha[which(temp_dedup$acc == max(temp_dedup$acc))]
  
  # get dataset of predictions and labels for both small and large data
  pred_short <- prediction(temp_dedup$valid_age_pred, temp_dedup$valid_age_label)
  # pred_long <- prediction($preds, final_dat$real)
  
  # get performace objects
  perf_s <- performance(prediction.obj = pred_short, measure = 'tpr', x.measure = 'fpr')
  # perf_l <- performance(prediction.obj = pred_long, measure = 'tpr', x.measure = 'fpr')
  
  # plot mean preds
  plot(perf_s)
  abline(a = 0, b =1)
  
  # plot all preds
  plot(perf_l)
  abline(a = 0, b =1)
  
  
  # plot optimal cutoff (TPR vs TNR)
  temp_roc <- performance(pred_short, measure = 'tpr', x.measure = 'tnr')
  cutoffs <- data.frame(cut=temp_roc@alpha.values[[1]], tpr=temp_roc@x.values[[1]],
                        tnr=temp_roc@y.values[[1]])
  
  # combine plots 
  plot(cutoffs$cut, cutoffs$tpr,type="l",col="red")
  par(new=TRUE)
  plot(cutoffs$cut, cutoffs$tnr,type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="")
  
  # get precision recall curves 
  temp_pr <- performance(pred_short, measure = 'prec', x.measure = 'rec')
  plot(temp_pr)
  
  # get lift curve
  temp_lc <- performance(pred_short, measure="lift")
  plot(temp_lc)
  
  
  ###############
  # get optimal cutoff and confusion matrix info 
  ###############
  
  # optimal cutoff
  cost.perf = performance(pred_short, "cost", cost.fp = 1, cost.fn = 1)
  optimal_cutoff <- pred_short@cutoffs[[1]][which.min(cost.perf@y.values[[1]])]
  
  # get confusion matrix function for plotting 
  ConfusionMatrixInfo(data = dat_sample, 
                      predict = 'mean_preds', 
                      actual = 'real_label', 
                      cutoff = optimal_cutoff,
                      get_plot = TRUE)
  
  
  
  dev.off()
  dat_sample$real_label <- as.character(dat_sample$real_label)
  # get cost function by specifying cost before
  cost_fp <- 1
  cost_fn <- 1
  roc_info <- ROCInfo(data = dat_sample, 
                      predict = 'mean_preds', 
                      actual = 'real_label', 
                      cost.fp = cost_fp,
                      cost.fn = cost_fn)
  
  grid.draw(roc_info$plot)
  
  ###############
  # get acc optimal cutoff plot
  ###############
  dev.off()
  acc.perf = performance(pred_short, measure = "acc")
  plot(acc.perf)
  abline(v = optimal_cutoff)
  
  ###############
  # Other plots
  ###############
  
  # distribution of preds
  ggplot(dat_sample, aes(mean_preds)) + 
    geom_histogram(aes(y=..density..), bins = 15, colour="black", fill="white")+
    geom_density(adjust = 1, fill= 'black', alpha = 0.4) +
    labs(title = 'Distribution of risk scores',
         x = 'Mean risk',
         y = 'Density')
  
  # plot age as a function of prediction
  ggplot(dat_sample, 
         aes(age, mean_preds)) + 
    geom_point() +
    labs(title = 'Age in months vs risk score',
         x = 'Age (months)',
         y = 'Mean risk score') +
    geom_smooth(method = 'loess')
  
  
  t_dat <- dat[dat$alpha == 0.5,]
  # get confusion matrix function for plotting 
  t <-ConfusionMatrixInfo(data = t_dat, 
                          predict = 'valid_age_pred', 
                          actual = 'valid_age_label', 
                          cutoff = .6,
                          get_plot = FALSE)
  
  val_dat <- as.data.frame(cbind(t_dat, t))
  val_dat$valid_age_pred <- val_dat$vallid_age_label <- NULL
  
 
} else {
  # subset data
  dat <- dat[, c('tm_donor','cancer_diagnosis_diagnoses','controls_age_pred', 'controls_age_label', 'age_sample_collection',
                  'alpha', 'non_zero')]
  
  # create a function for controls that gets the young predictions for each alpha level
  
  
  temp <- get_young_labels(dat, thresh = 0.6, age = 72)
  
  # remove duplicates
  temp <- temp[!duplicated(temp$tm_donor, fromLast = TRUE),]
  
}





# beta and gender .45
# beta and no gender .45
# m and gender .49
# m and no gender .49


