
# set parameter for whether we should remove age related cgs
remove_age_cgs_lm = F
remove_age_cgs_lit = T

# set preprocessing method
method <- 'noob'

# set type of data, beta or m
methyl_type <- 'm'

# source all_functions.R to load libraries and my functions
source('all_functions.R')

# set data directory
data_dir <- '../../Data/'

# get data
if(paste0(methyl_type, '_processed', '.RData') %in% dir(data_dir)) {
  load(methyl_type, '_processed', '.RData') 
} else {
  source('get_data_new.R')
}

##########
# run model
##########
# https://rdrr.io/rforge/glmmixedlasso/man/glmmlasso.html

run_model <- function(cases_full,
                      controls_full,
                      k_folds = k_folds,
                      beta_thresh = beta_thresh) {
  
  
  # get vector of random folds
  fold_vec <- sample(1:k_folds, nrow(beta_cases_full), replace = T)
  
  # combine 
  for(i in 1:k_folds) {
    
    # get train and test index
    train_index <- !grepl(i, fold_vec)
    test_index <- !train_index
    
    
    beta_train_cases <- beta_cases_full[train_index, ]
    beta_test_cases <- beta_cases_full[test_index, ]
    
    
    # use cases training and controls to get bumphunter features
    bh_feats <- bump_hunter(dat_1 = beta_train_cases, 
                            dat_2 = beta_controls_full, 
                            bump = 'cancer', 
                            boot_num = 5, 
                            thresh = beta_thresh,
                            g_ranges = g_ranges)
    
    
    # get feature list
    colnames(bh_feats)[1] <- 'chr'
    remove_features <- inner_join(bh_feats, g_ranges)$probe
    
    # take remove features out of colnames 
    bh_features <- intersect_names[!intersect_names %in% remove_features]
    
    # function to predict with all test, controls, controls old, and valid
    mod_result  <- predCancer(training_dat = beta_train_cases,
                              controls_dat = beta_controls_full,
                              test_dat = beta_test_cases,
                              age_cutoff = age_cutoff,
                              gender = gender,
                              tech = tech,
                              base_change = base_change,
                              exon_intron = exon_intron,
                              bh_features = bh_features)
    
    temp_cases[[i]] <- mod_result[[1]]
    temp_controls[[i]] <- mod_result[[2]]
    temp_models[[i]] <- mod_result[[3]]
    temp_lambda[[i]] <- mod_result[[4]]
    temp_alpha[[i]] <- mod_result[[5]]
    
    
    print(i)
  }
  
  full_cases <- do.call(rbind, temp_cases)
  full_controls <- do.call(rbind, temp_controls)
  
  return(list(full_cases, full_controls, temp_models, temp_lambda, temp_alpha))
}


##########
# fixed variables
##########


# run full pipeline
full_results <- run_model(cases_full,
                          controls_full,
                          k_folds = k_folds,
                          beta_thresh = beta_thresh)


# get results from list 
temp_cases <- full_results[[1]]

temp_cases <- temp_cases[ , c('test_pred', 'test_label', 
                              'age_diagnosis' ,  'age_sample_collection')]
temp_cases$test_pred_label <- ifelse(temp_cases$test_pred > .5, 1, 0)

temp_cases$pred_is <- ifelse(temp_cases$test_pred_label == temp_cases$test_label, 
                             'good',
                             'bad')


temp_controls <- full_results[[2]]

temp_controls <- temp_controls[ , c('controls_age_pred', 'controls_age_label', 
                                    'age_sample_collection')]


# get person identfier 
temp_controls$p_id <- rep.int(seq(1, 35, 1), 4)

# group by fold and get mean 
temp_controls_pred <- temp_controls %>%
  group_by(p_id) %>%
  summarise(mean_pred = mean(controls_age_pred, na.rm =T)) %>%
  cbind(temp_controls[1:35,])

temp_controls$controls_pred_label <- ifelse(temp_controls$controls_age_pred > .5, 1, 0)

temp_controls$pred_is <- ifelse(temp_controls$controls_pred_label == temp_controls$controls_age_label, 
                                'good',
                                'bad')

# remove original prediction 
temp_controls_pred$controls_age_pred <- NULL
rm(temp_controls,temp_results)

temp_controls_pred <- temp_controls_pred[order(temp_controls_pred$age_sample_collection, decreasing = F),]

##########
# examin cases with prediction objects (ROC, TRP, etc)
##########
temp_pred_cases <- prediction(temp_cases$test_pred, temp_cases$test_label)
temp_roc <- performance(temp_pred_cases, measure = 'tpr', x.measure = 'tnr')
cutoffs <- data.frame(cut=temp_roc@alpha.values[[1]], tpr=temp_roc@x.values[[1]], 
                      tnr=temp_roc@y.values[[1]])

plot(cutoffs$cut, cutoffs$tpr,type="l",col="red")
par(new=TRUE)
plot(cutoffs$cut, cutoffs$tnr,type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="")

temp_pr <- performance(temp_pred_cases, measure = 'prec', x.measure = 'rec')
temp_lc <- performance(temp_pred_cases, measure="lift")

plot(temp_roc)
plot(temp_pr)
plot(temp_lc)
caret::confusionMatrix(round(temp_cases$test_pred, 0.1), temp_cases$test_label)

