
temp_clin <- read_csv('~/Desktop/temp_clin.csv')
temp_cases <- temp_clin[!grepl('Unaffected', temp_clin$cancer_diagnosis_diagnoses),]
temp_controls <- temp_clin[grepl('Unaffected', temp_clin$cancer_diagnosis_diagnoses),]

# validation or combined data 
data_used <- 'old'

# get cg regions
cg_gene_regions = 'Body'

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
  load(paste0(data_dir, methyl_type, '_processed', '.RData'))
} else {
  source(paste0('get_data_', data_used, '.R'))
}

# combine temp_cases with data_cases_full 
cases_family <- as.data.frame(cbind(tm_donor_ = temp_cases$tm_donor_, fam_cancer = temp_cases$fam_cancer))
cases_family$tm_donor_ <- as.numeric(as.character(cases_family$tm_donor_))
data_cases_full <- inner_join(data_cases_full, cases_family, by = 'tm_donor_')

# combine temp_controls with data_controls_full 
controls_family <- as.data.frame(cbind(tm_donor_ = temp_controls$tm_donor_, fam_cancer = temp_controls$fam_cancer))
controls_family$tm_donor_ <- as.numeric(as.character(controls_family$tm_donor_))
data_controls_full <- inner_join(data_controls_full, controls_family, by = 'tm_donor_')

# get family cancer for each data set
data_cases_full <- cbind(as.data.frame(class.ind(data_cases_full$fam_cancer)), data_cases_full)
data_controls_full <- cbind(as.data.frame(class.ind(data_controls_full$fam_cancer)), data_controls_full)
data_cases_full$fam_cancer <- data_controls_full$fam_cancer <- NULL

##########
# run model
##########
# https://rdrr.io/rforge/glmmixedlasso/man/glmmlasso.html

cases_full <- data_cases_full
controls_full <- data_controls_full
k_folds <- 5
beta_thresh <- 0.1

run_model <- function(cases_full,
                      controls_full,
                      control_for_family,
                      bump_type,
                      k_folds = k_folds,
                      beta_thresh = beta_thresh) {
  
  # get vector of random folds
  fold_vec <- sample(1:k_folds, nrow(cases_full), replace = T)
  
  # combine 
  for(i in 1:k_folds) {
    
    # read in probes associated with age
    if (remove_age_cgs_lit) {
      age_cgs <- readRDS('../../Data/age_probes.rda')
      intersect_names <- intersect_names[!intersect_names %in% age_cgs]
    }
    
    
    # remove age
    if (remove_age_cgs_lm) {
      age_cgs_lm <- readRDS('../../Data/age_cgs_lm.rda')
      intersect_names <- intersect_names[!intersect_names %in% age_cgs_lm]
    }
    
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
                              bump = 'both', 
                              boot_num = 5, 
                              thresh = beta_thresh,
                              g_ranges = g_ranges)
      
      
      # get feature list
      colnames(bh_feats)[1] <- 'chr'
      these_features <- inner_join(bh_feats, g_ranges)$probe
      
      # subset data by these_features
      train_cases <- train_cases[, these_features]
      controls_full <- controls_full[, these_features]
    }
    
    # use cases training and controls to get bumphunter features
    bh_feats <- bump_hunter(dat_1 = train_cases, 
                            dat_2 = controls_full, 
                            bump = 'cancer', 
                            boot_num = 5, 
                            thresh = beta_thresh,
                            g_ranges = g_ranges)
    
    
    # get feature list
    colnames(bh_feats)[1] <- 'chr'
    remove_features <- inner_join(bh_feats, g_ranges)$probe
    
    # get features
    temp_features <- colnames(train_cases)[18:ncol(train_cases)]
    # take remove features out of colnames 
    bh_features <- temp_features[!temp_features %in% remove_features]
    
    if(enet){
      # function to predict with all test, controls, controls old, and valid
      mod_result  <- run_enet_450_850(training_dat = train_cases,
                                      controls_dat = controls_full,
                                      test_dat = test_cases,
                                      age_cutoff = 72,
                                      gender = T,
                                      tech = T,
                                      fam_cancer = T,
                                      bh_features = bh_features)
      
      temp_cases[[i]] <- mod_result[[1]]
      temp_controls[[i]] <- mod_result[[2]]
      temp_models[[i]] <- mod_result[[3]]
      temp_lambda[[i]] <- mod_result[[4]]
      temp_alpha[[i]] <- mod_result[[5]]
      
      
    } else {
      # function to predict with all test, controls, controls old, and valid
      mod_result  <- run_rf(training_dat = train_cases,
                            controls_dat = controls_full,
                            test_dat = test_cases,
                            age_cutoff = 72,
                            gender = T,
                            tech = T,
                            fam_cancer = T,
                            bh_features = bh_features)
    }
   
    print(i)
  }
  
  if(enet){
    full_cases <- do.call(rbind, temp_cases)
    full_controls <- do.call(rbind, temp_controls)
    
    return(list(full_cases, full_controls, temp_models, temp_lambda, temp_alpha))
  } else {
    full_cases <- do.call(rbind, temp_cases)
    full_controls <- do.call(rbind, temp_controls)
    
    return(list(full_cases, full_controls))
  }

}


##########
# fixed variables
##########


# run full pipeline
full_results <- run_model(cases_full,
                          controls_full,
                          control_for_family,
                          bump_type,
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

