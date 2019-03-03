
# source functions script
source('all_functions.R')


# create fixed objects to model and pipeline inputs
methyl_type <- 'm'
combat <- TRUE
if(combat){
  tech <- FALSE
} else {
  tech <- TRUE
}
gender <- TRUE
k_folds <- 10
beta_thresh <- 0.1
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

# create tech and gender dummy variables 

# tech
all_cases <- cbind(as.data.frame(class.ind(all_cases$tech)), 
                               all_cases)

all_con <- cbind(as.data.frame(class.ind(all_con$tech)), 
                             all_con)

# rempove old tech variable 
all_cases$tech <- all_con$tech <- NULL

# gender
all_cases <- cbind(as.data.frame(class.ind(all_cases$gender)), 
                               all_cases)

all_con <- cbind(as.data.frame(class.ind(all_con$gender)), 
                             all_con)

# rempove old tech variable 
all_cases$gender <- all_con$gender <- NULL

# make ages back to months (times by 12)

# cases
all_cases$age_diagnosis <- 
  round(all_cases$age_diagnosis*12, 2)
all_cases$age_sample_collection <- 
  round(all_cases$age_sample_collection*12, 2)

# controls
all_con$age_diagnosis <- 
  round(all_con$age_diagnosis*12, 2)
all_con$age_sample_collection <- 
  round(all_con$age_sample_collection*12, 2)

# Remove NA in age of dianogsis for cases
all_cases <- all_cases[!is.na(all_cases$age_sample_collection),]
all_con <- all_con[!is.na(all_con$age_sample_collection),]

# prepare data sets for modelling
# cases_full <- all_cases
# controls_full <- all_con
# k_folds = 5
# beta_thresh = 0.05
# g_ranges = g_ranges
# i = 1

run_model <- function(cases_full,
                      controls_full,
                      k_folds = k_folds,
                      tech = tech,
                      gender = gender,
                      beta_thresh = beta_thresh,
                      methyl_type = methyl_type,
                      g_ranges = g_ranges) {
  
  
  # create lists to store model results
  temp_cases <- list()
  temp_controls <- list()
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
    
    # function to predict with all test, controls, controls old, and valid
    mod_result  <- run_enet_450_850(training_dat = train_cases,
                                    controls_dat = controls_full,
                                    test_dat = test_cases,
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

# save temporary best model image
# rm(all_cases, all_con, g_ranges, ratio_set)
# save.image('../../Data/temp_best_model_beta_nocombat_nogender_yestech.RData')

## ANALYSIS
# creat an extra test_pred label called test_pred_lab that assumes 'positive' if over 0.50
# and negative otherwise. Don't ovewrite test_pred
temp_1$test_pred_label <- as.factor(ifelse(temp_1$test_pred > 0.5, 'positive', 'negative'))
 
# get the confusion matrix
caret::confusionMatrix(temp_1$test_pred_label, temp_1$test_label)

########
temp_2 <- temp_res[[2]]

# get label for controls
temp_2$controls_age_pred_label <- as.factor(ifelse(temp_2$controls_age_pred > 0.5, 'positive', 'negative'))

# look at controls
temp_con <- temp_2 %>%
  group_by(ids) %>%
  summarise(mean_controls_age_pred = mean(controls_age_pred),
            sum_controls_age_pred_label_positive = sum(controls_age_pred_label == 'positive'),
            sum_controls_age_pred_label_negative = sum(controls_age_pred_label == 'negative')) %>%
  left_join(temp_2, by = 'ids') %>% 
  dplyr::select(ids, mean_controls_age_pred, sum_controls_age_pred_label_positive,
                sum_controls_age_pred_label_negative, controls_age_label,
                age_sample_collection) %>%
  distinct()

# get mean pred label
temp_con$mean_controls_age_pred_label <- as.factor(ifelse(temp_con$mean_controls_age_pred > 0.5, 'positive', 'negative'))


# saveRDS(list(temp_1, temp_2, temp_3, temp_4, temp_5, temp_6), '../../Data/temp_8611.rda')

temp_auc <- performance(temp_pred_cases, measure = 'auc')
temp_roc <- performance(temp_pred_cases, measure = 'tpr', x.measure = 'fpr')
temp_pr <- performance(temp_pred_cases, measure = 'prec', x.measure = 'rec')
temp_roc <- as.data.frame(cbind(unlist(temp_pr@x.values),unlist(temp_pr@y.values)))
temp_lc <- performance(temp_pred_cases, measure="lift")
hist(temp_1_new$test_pred, 
     main = 'Predicted Risk Score on Cases', 
     xlab = 'Risk scores', 
     col = 'darkgrey')




# Controls 
temp_2 <- full_results[[2]]

temp_2 <- temp_2[ , c('controls_age_pred', 'controls_age_label',
                                    'age_sample_collection')]

temp_2 <- temp_2[order(temp_2$controls_age_pred, decreasing = T),]

# get person identfier
temp_2$p_id <- rep.int(seq(1, 41, 1), 5)

# group by fold and get mean
temp_2_pred <- temp_2 %>%
  group_by(p_id) %>%
  summarise(mean_pred = mean(controls_age_pred, na.rm =T)) %>%
  cbind(temp_2[1:41,])

temp_2$controls_pred_label <- ifelse(temp_2$controls_age_pred > .5, 1, 0)

temp_2$pred_is <- ifelse(temp_2$controls_pred_label == temp_2$controls_age_label,
                                'good',
                                'bad')
# remove original prediction
temp_2_pred$controls_age_pred <- NULL

temp_2_pred <- temp_2_pred[order(temp_2_pred$age_sample_collection, decreasing = F),]
