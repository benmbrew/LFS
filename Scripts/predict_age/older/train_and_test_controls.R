

# validation or combined data 
data_used <- 'new'

# get cg regions
cg_gene_regions = 'Body'

# set preprocessing method
method <- 'noob'

# set type of data, beta or m
methyl_type <- 'm'

# set data directory
data_dir <- '../../Data/'


temp_data <- full_data_first %>% filter(ids %in% these_ids)
temp_data <- full_data_first[full_data_first$ids %in% these_ids,]
# source all_functions.R to load libraries and my functions
source('all_functions.R')

# read in wt data
data_wt <- readRDS(paste0(data_dir,paste0('new','_',methyl_type, '_wild_type', '.rda')))

full_data_first <- full_data_first[full_data_first$tm_donor_ != '3955',]
full_data_last <- full_data_last[full_data_last$tm_donor_ != '3955',]

full_data_first_combat <- full_data_first_combat[full_data_first_combat$tm_donor_ != '3955',]
full_data_last_combat <- full_data_last_combat[full_data_last_combat$tm_donor_ != '3955',]

temp <- get_diff_dups(full_data_first)
full_data_first <- temp[[1]]
full_data_last <- temp[[2]]
temp_combat <- get_diff_dups(full_data_first_combat)
full_data_first_combat <- temp_combat[[1]]
full_data_last_combat <- temp_combat[[2]]
# get controls - HERE
# get age dummy cat
full_data_first <- get_age_cat_dummy(full_data_first)
full_data_first_combat <- get_age_cat_dummy(full_data_first_combat)

full_data_last <- get_age_cat_dummy(full_data_last)
full_data_last_combat <- get_age_cat_dummy(full_data_last_combat)

data_wt <- get_age_cat_dummy(data_wt)

##########
# run model
##########
# https://rdrr.io/rforge/glmmixedlasso/man/glmmlasso.html
data_full <- full_data_first
wt_data <- data_wt


run_model <- function(data_full,
                      enet,
                      wt_data,
                      bump_type,
                      gender,
                      tech,
                      fam_num,
                      fam_ratio,
                      age_dum,
                      remove_age_cgs_lit,
                      remove_age_cgs_lm,
                      k_folds = k_folds,
                      beta_thresh = beta_thresh) {
  
  
  intersect_names <- colnames(data_full)
  clin_data_names <- names(data_full)[!grepl('^cg', colnames(data_full))]
  # read in probes associated with age
  if (remove_age_cgs_lit) {
    age_cgs <- readRDS('../../Data/age_probes.rda')
    intersect_names <- intersect_names[!intersect_names %in% age_cgs]
    data_full <- data_full[, colnames(data_full) %in% intersect_names]
  }
  
  
  # remove age
  if (remove_age_cgs_lm) {
    age_cgs_lm <- readRDS('../../Data/age_cgs_lm.rda')
    intersect_names <- intersect_names[!intersect_names %in% age_cgs_lm]
    data_full <- data_full[, colnames(data_full) %in% intersect_names]
  }
  
  
  # get two data sets from data_full - cases and controls (containing both technologies)
  cases_full <- data_full[data_full$cancer_diagnosis_diagnoses != 'Unaffected',]
  controls_full <- data_full[data_full$cancer_diagnosis_diagnoses == 'Unaffected',]
  
  
  # get vector of random folds
  fold_vec <- sample(1:k_folds, nrow(cases_full), replace = T)
  
  # define lists
  temp_cases <- list()
  temp_controls <- list()
  temp_min <- list()
  temp_1se <- list()
  temp_alpha <- list()
  temp_models <- list()
  temp_importance <- list()
  temp_size <- list()
  # combine 
  for(i in 1:k_folds) {
    
    # get train and test index
    train_index <- !grepl(i, fold_vec)
    test_index <- !train_index
    
    # subset into training and test set
    train_cases <- cases_full[train_index, ]
    test_cases <- cases_full[test_index, ]
    
    if(bump_type == 'both') {
      # use cases training and controls to get bumphunter features
      bh_feats <- bump_hunter(dat_1 = train_cases, 
                              dat_2 = controls_full, 
                              wild_type = wt_data,
                              bump = 'lfs', 
                              boot_num = 5, 
                              thresh = beta_thresh,
                              g_ranges = g_ranges)
      
      # get feature list
      colnames(bh_feats)[1] <- 'chr'
      these_features <- inner_join(bh_feats, g_ranges)$probe
      
      # subset data by these_features
      train_cases <- train_cases[, c(clin_data_names, these_features)]
      test_cases <- test_cases[, c(clin_data_names, these_features)]
      controls_temp <- controls_full[, c(clin_data_names, these_features)]
      
      # use cases training and controls to get bumphunter features
      bh_feats <- bump_hunter(dat_1 = train_cases, 
                              dat_2 = controls_temp,
                              wild_type = NULL,
                              bump = 'cancer', 
                              boot_num = 5, 
                              thresh = 0.05,
                              g_ranges = g_ranges)
      
      
      # get feature list
      colnames(bh_feats)[1] <- 'chr'
      remove_features <- inner_join(bh_feats, g_ranges)$probe
      
      # get features
      temp_features <- colnames(train_cases)[23:ncol(train_cases)]
      # take remove features out of colnames 
      remaining_features <- temp_features[!temp_features %in% remove_features]
      
    } else {
      
      # use cases training and controls to get bumphunter features
      bh_feats <- bump_hunter(dat_1 = train_cases, 
                              dat_2 = controls_full,
                              wild_type = NULL,
                              bump = 'cancer', 
                              boot_num = 5, 
                              thresh = beta_thresh,
                              g_ranges = g_ranges)
      
      
      # get feature list
      colnames(bh_feats)[1] <- 'chr'
      remove_features <- inner_join(bh_feats, g_ranges)$probe
      
      # get features
      temp_features <- colnames(train_cases)[23:ncol(train_cases)]
      # take remove features out of colnames 
      remaining_features <- temp_features[!temp_features %in% remove_features]
      
    }
      if(enet){
        
        colnames(controls_full)[grepl('gdna', colnames(controls_full))]
        # function to predict with all test, controls, controls old, and valid
        mod_result  <- run_enet_450_850(training_dat = train_cases,
                                        controls_dat = controls_full,
                                        test_dat = test_cases,
                                        age_cutoff = 72,
                                        age_dum = age_dum,
                                        gender = gender,
                                        tech = tech,
                                        fam_num = fam_num,
                                        fam_ratio = fam_ratio,
                                        bh_features = remaining_features)
        
        
        temp_cases[[i]] <- mod_result[[1]]
        temp_controls[[i]] <- mod_result[[2]]
        temp_models[[i]] <- mod_result[[3]]
        temp_min[[i]] <- mod_result[[4]]
        temp_1se[[i]] <- mod_result[[5]]
        temp_alpha[[i]] <- mod_result[[6]]
        
      } else {
        # function to predict with all test, controls, controls old, and valid
        mod_result  <- run_rf(training_dat = train_cases,
                              controls_dat = controls_full,
                              test_dat = test_cases,
                              age_cutoff = 72,
                              age_dum = age_dum,
                              gender = gender,
                              tech = tech,
                              fam_num = fam_num,
                              fam_ratio = fam_ratio,
                              bh_features = remaining_features)
        
        temp_cases[[i]] <- mod_result[[1]]
        temp_controls[[i]] <- mod_result[[2]]
        temp_models[[i]] <- mod_result[[3]]
        temp_importance[[i]] <- mod_result[[4]]
        temp_size[[i]] <- mod_result[[5]]
        
      }
      
    print(i)
  }
  
  if(enet){
    
    
    full_cases <- do.call(rbind, temp_cases)
    full_controls <- do.call(rbind, temp_controls)
    return(list(full_cases, full_controls, temp_models, temp_min, temp_1se))
    
    
    
  } else {
    
    
    # HERE (validation)
    full_cases <- do.call(rbind, temp_cases)
    full_controls <- do.call(rbind, temp_controls)
    return(list(full_cases, full_controls, temp_models, temp_importance, temp_size))
    
    
  }
  
}

# 4 before the age of 6
# 5 after the age of 6 that we misclassify
##########
# fixed variables
##########
# cancer, differet ages
# full_data, combat, 

# run full pipeline
full_results <- run_model(full_data_first,
                          enet = T,
                          wt_data = data_wt,
                          bump_type = 'both',
                          gender = T,
                          tech = T,
                          fam_num = F,
                          fam_ratio = F,
                          age_dum = T,
                          remove_age_cgs_lit = F,
                          remove_age_cgs_lm = F,
                          k_folds = 5,
                          beta_thresh = 0.01)


# temp_size <- full_results[[5]]
# temp_size_1 <- temp_size[[1]]
# temp_size_2 <- temp_size[[2]]
# temp_size_3 <- temp_size[[3]]
# temp_size_4 <- temp_size[[4]]
# temp_size_5 <- temp_size[[5]]
# temp_importance <- full_results[[4]]
# temp_1 <- temp_importance[[1]]
# temp_2 <- temp_importance[[2]]
# temp_3 <- temp_importance[[3]]
# temp_4 <- temp_importance[[4]]
# temp_5 <- temp_importance[[5]]
# temp_all <- Reduce(intersect, list(temp_1, temp_2, temp_3, temp_4, temp_5))
# 

# get results from list
temp_cases <- full_results[[1]]

temp_cases <- temp_cases[ , c('test_pred', 'test_label',
                              'age_diagnosis' ,  'age_sample_collection')]
temp_cases$test_pred_label <- ifelse(temp_cases$test_pred > .5, 1, 0)

temp_cases$pred_is <- ifelse(temp_cases$test_pred_label == temp_cases$test_label,
                             'good',
                             'bad')

temp_bad <- temp_cases[temp_cases$pred_is == 'bad',]
temp_good <- temp_cases[temp_cases$pred_is == 'good',]

temp_cases_before <- temp_cases[temp_cases$age_sample_collection < 72,]
temp_cases_before$test_pred_label <- 1
temp_cases_after <- temp_cases[temp_cases$age_sample_collection >= 72,]

temp_cases_before <- temp_cases_before[order(temp_cases_before$age_sample_collection),]
temp_cases_after <- temp_cases_after[order(temp_cases_after$age_sample_collection),]
write.csv(temp_cases_after, '~/Desktop/after_6.csv')

temp_cases_after6 <- read.csv('~/Desktop/temp_after_6.csv')
temp_cases_after6$X <- NULL
temp_cases_new <- rbind(temp_cases_after6, 
                        temp_cases_before)

# temp_cases <- temp_cases[order(temp_cases$test_pred, decreasing = T),]
length(which(temp_cases$pred_is == 'good'))

temp_cases <- temp_cases[order(temp_cases$age_sample_collection, decreasing = T),]

# temp_cases_young <- temp_cases[temp_cases$age_sample_collection < 216,]
length(which(temp_cases_young$pred_is == 'good'))

# temp_cases <- temp_cases[order(temp_cases$test_pred, decreasing = T),]
save.image('/../../paper_cases.RData')

##########
# examin cases with prediction objects (ROC, TRP, etc)
##########

ggplot(temp_cases_new, aes(test_pred)) + 
  geom_histogram(aes(y=..count..),      # Histogram with density instead of count on y-axis
                 binwidth=.08,
                 colour="black", fill="black", alpha = 0.6) +
  geom_density(aes(y=..density.. *(103*.08)), 
               alpha=.2, fill= 'blue') + 
  labs(title = 'Distribution of risk scores',
       x = 'Model risk score',
       y = 'Frequency') + theme_igray(base_size = 16, base_family = 'Calibri')
  
# histogram for onset, and onset by cases and controls
# Histogram overlaid with kernel density curve
ggplot(clin, aes(x=age_sample_collection)) +
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.5,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot




ggplot(clin, aes(x=age_diagnosis)) + 
  geom_histogram(aes(y=..count..),      # Histogram with density instead of count on y-axis
                 binwidth=2,
                 colour="black", fill="black") +
  geom_density(aes(y = ..density.. *(430*2)),
               alpha=.2, fill="blue") +
  labs(title = '', x = 'Age of sample collection', y = 'Counts') +
  theme_minimal(base_size = 20, base_family = 'Ubuntu')



temp_pred_cases <- prediction(temp_cases_new$test_pred, temp_cases_new$test_label)
temp_auc <- performance(temp_pred_cases, measure = 'auc')
temp_roc <- performance(temp_pred_cases, measure = 'tpr', x.measure = 'fpr')
temp_pr <- performance(temp_pred_cases, measure = 'prec', x.measure = 'rec')
temp_roc <- as.data.frame(cbind(unlist(temp_pr@x.values),unlist(temp_pr@y.values)))
temp_lc <- performance(temp_pred_cases, measure="lift")
hist(temp_cases_new$test_pred, 
     main = 'Predicted Risk Score on Cases', 
     xlab = 'Risk scores', 
     col = 'darkgrey')

cutoffs <- data.frame(cut=temp_roc@alpha.values[[1]], fpr=temp_roc@x.values[[1]],
                      tpr=temp_roc@y.values[[1]])


save.image('plots_for_poster.RData')
plot(cutoffs$cut, cutoffs$tpr,type="l",col="red", xlab ='Cutoffs', ylab ='' ,main = 'TPR vs TNR', lwd = 2)
par(new=TRUE)


ggplot(cutoffs, aes(cut)) +
  geom_line(aes(y = tnr), colour = "#CC6666", size = 2) +
  geom_line(aes(y = tpr), colour = "#9999CC", size = 2) +
  labs(x = 'Cutoff', y ='', title = 'TPR vs TNR') + theme_igray(base_size = 16, base_family = 'Calibri')


ggplot(cutoffs, aes(tpr, tpr)) +
  geom_line(colour = "#CC6666", size = 2) +
  labs(x = 'TNR', y ='TPR', title = 'ROC curve') +
  theme_igray(base_size = 16, base_family = 'Calibri')
  
plot(cutoffs$cut, 
     cutoffs$tnr,
     lwd = 2,
     type="l",
     col="blue",xaxt="n",yaxt="n",xlab="",ylab="")


plot(temp_roc, 
     lwd = 2.5,
     col = 'blue',bty= 'n')

plot(temp_pr,lwd = 2.5, col = 'darkgrey',bty= 'n')
plot(temp_lc)
caret::confusionMatrix(round(temp_cases$test_pred, 0.1), temp_cases$test_label)

caret::confusionMatrix(temp_cases_new$test_pred_label, temp_cases_new$test_label)










temp_controls <- full_results[[2]]

temp_controls <- temp_controls[ , c('controls_age_pred', 'controls_age_label',
                                    'age_sample_collection')]

temp_controls <- temp_controls[order(temp_controls$controls_age_pred, decreasing = T),]

# get person identfier
temp_controls$p_id <- rep.int(seq(1, 41, 1), 5)

# group by fold and get mean
temp_controls_pred <- temp_controls %>%
  group_by(p_id) %>%
  summarise(mean_pred = mean(controls_age_pred, na.rm =T)) %>%
  cbind(temp_controls[1:41,])

temp_controls$controls_pred_label <- ifelse(temp_controls$controls_age_pred > .5, 1, 0)

temp_controls$pred_is <- ifelse(temp_controls$controls_pred_label == temp_controls$controls_age_label,
                                'good',
                                'bad')
# remove original prediction
temp_controls_pred$controls_age_pred <- NULL

temp_controls_pred <- temp_controls_pred[order(temp_controls_pred$age_sample_collection, decreasing = F),]

# ##########
# # examin cases with prediction objects (ROC, TRP, etc)
# ##########
# temp_pred_cases <- prediction(temp_cases$test_pred, temp_cases$test_label)
# temp_roc <- performance(temp_pred_cases, measure = 'tpr', x.measure = 'tnr')
# cutoffs <- data.frame(cut=temp_roc@alpha.values[[1]], tpr=temp_roc@x.values[[1]], 
#                       tnr=temp_roc@y.values[[1]])
# 
# plot(cutoffs$cut, cutoffs$tpr,type="l",col="red")
# par(new=TRUE)
# plot(cutoffs$cut, cutoffs$tnr,type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="")
# 
# temp_pr <- performance(temp_pred_cases, measure = 'prec', x.measure = 'rec')
# temp_lc <- performance(temp_pred_cases, measure="lift")
# 
# plot(temp_roc)
# plot(temp_pr)
# plot(temp_lc)
# caret::confusionMatrix(round(temp_cases$test_pred, 0.1), temp_cases$test_label)
# 

