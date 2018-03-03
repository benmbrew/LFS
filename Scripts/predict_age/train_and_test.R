

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

# get data
if(paste0(data_used,'_',methyl_type, '_final', '.RData') %in% dir(data_dir)) {
  load(paste0(data_dir,paste0(data_used,'_',methyl_type, '_final', '.RData')))
} else {
  source(paste0('get_combat.R'))
}

# read in wt data
data_wt <- readRDS(paste0(data_dir,paste0(data_used,'_',methyl_type, '_wild_type', '.rda')))

# source all_functions.R to load libraries and my functions
source('all_functions.R')
##########
# run model
##########
# https://rdrr.io/rforge/glmmixedlasso/man/glmmlasso.html

data_full <- full_data[, c(1:5000, ncol(full_data))]
enet = F
wt_data <- data_wt
bump_type <- 'both'
k_folds <- 5
beta_thresh <- 0.1
remove_age_cgs_lm = F
remove_age_cgs_lit = T

run_model <- function(data_full,
                      enet,
                      wt_data,
                      bump_type,
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
  temp_lambda <- list()
  temp_alpha <- list()
  temp_models <- list()
  temp_importance <- list()
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
                              thresh = beta_thresh,
                              g_ranges = g_ranges)
      
      
      # get feature list
      colnames(bh_feats)[1] <- 'chr'
      remove_features <- inner_join(bh_feats, g_ranges)$probe
      
      # get features
      temp_features <- colnames(train_cases)[18:ncol(train_cases)]
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
      temp_features <- colnames(train_cases)[18:ncol(train_cases)]
      # take remove features out of colnames 
      remaining_features <- temp_features[!temp_features %in% remove_features]
      
    }
      if(enet){
        
        # function to predict with all test, controls, controls old, and valid
        mod_result  <- run_enet_450_850(training_dat = train_cases,
                                        controls_dat = controls_full,
                                        test_dat = test_cases,
                                        age_cutoff = 72,
                                        gender = T,
                                        tech = T,
                                        fam_num = T,
                                        fam_ratio = T,
                                        bh_features = remaining_features)
        
        
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
                              fam_num = T,
                              fam_ratio = T,
                              bh_features = remaining_features)
        
        temp_cases[[i]] <- mod_result[[1]]
        temp_controls[[i]] <- mod_result[[2]]
        temp_models[[i]] <- mod_result[[3]]
        temp_importance[[i]] <- mod_result[[4]]
        
      }
      
    print(i)
  }
  
  if(enet){
    
    
    full_cases <- do.call(rbind, temp_cases)
    full_controls <- do.call(rbind, temp_controls)
    return(list(full_cases, full_controls, temp_models, temp_lambda, temp_alpha))
    
    
    
  } else {
    
    
    # HERE (validation)
    full_cases <- do.call(rbind, temp_cases)
    full_controls <- do.call(rbind, temp_controls)
    return(list(full_cases, full_controls, temp_models, temp_importance))
    
    
  }
  
}

##########
# fixed variables
##########


# run full pipeline
full_results <- run_model(full_data[, 1:5000],
                          enet = F,
                          wt_data = data_wt,
                          bump_type = 'both',
                          remove_age_cgs_lit =T,
                          remove_age_cgs_lm = F,
                          k_folds = 5,
                          beta_thresh = 0.1)


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

