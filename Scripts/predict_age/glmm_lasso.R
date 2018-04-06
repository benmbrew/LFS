
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

# get age dummy cat
full_data$temp_age_var <- as.factor(ntile(full_data$age_sample_collection, 5))
data_wt$temp_age_var <- as.factor(ntile(data_wt$age_sample_collection, 5))

##########
# run model
##########
# https://rdrr.io/rforge/glmmixedlasso/man/glmmlasso.html

data_full <- full_data[, c(1:5000, ncol(full_data))]
wt_data <- data_wt
bump_type <- 'both'
k_folds <- 5
beta_thresh <- 0.1
remove_age_cgs_lm = F
remove_age_cgs_lit = T

run_model <- function(data_full,
                      wt_data,
                      bump_type,
                      gender,
                      tech,
                      fam_num,
                      fam_ratio,
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
  temp_clin <- list()
  temp_con_clin <- list()
  temp_cases_outcome <- list()
  temp_controls_outcome <- list()
  temp_models <- list()
  # combine 
  for(i in 1:k_folds) {
    
    controls_temp <- controls_full
    
    # get train and test index
    train_index <- !grepl(i, fold_vec)
    test_index <- !train_index
    
    # subset into training and test set
    train_cases <- cases_full[train_index, ]
    test_cases <- cases_full[test_index, ]
    
    # get clinical data
    test_clin <- test_cases[, 1:22 ]
    con_clin <- controls_full[, 1:22 ]
    
    
    if(bump_type == 'both') {
      # use cases training and controls to get bumphunter features
      bh_feats <- bump_hunter(dat_1 = train_cases, 
                              dat_2 = controls_temp, 
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
      controls_temp <- controls_temp[, c(clin_data_names, these_features)]
      
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
      temp_features <- colnames(train_cases)[23:ncol(train_cases)]
      # take remove features out of colnames 
      remaining_features <- temp_features[!temp_features %in% remove_features]
      
    } else if(bump_type == 'lfs_only') {
      
      # read in lfs data 
      lfs_features <- readRDS('../../Data/lfs_features_1.rda')
      
      remaining_features <- lfs_features
      
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
    
    # get model data train_cases, test_cases, and controls_full
    train_cases$outcome <- ifelse(train_cases$age_diagnosis < 72, 1, 0)
    test_cases$outcome <- ifelse(test_cases$age_diagnosis < 72, 1, 0)
    controls_temp$outcome <- ifelse(controls_temp$age_sample_collection < 72, 1, 0)
    
    # cg_sites
    cg_cols <- colnames(train_cases)[grepl('^cg', colnames(train_cases))]
    
   
    # make get lfs_features from cg_cols
    cg_cols <- cg_cols[cg_cols %in% remaining_features]
    
    # cg_cols <- c('fam_num_cancer', 'a', 'b')
    
    # get cg_cols for controls and test
    train_cases <- train_cases[, c('outcome', 'temp_age_var', cg_cols)]
    test_cases <- test_cases[, c('outcome', 'temp_age_var', cg_cols)]
    controls_temp <- controls_temp[, c('outcome', 'temp_age_var', cg_cols)]
    
    # change to factor 
    train_cases$temp_age_var <- as.factor(train_cases$temp_age_var)
    test_cases$temp_age_var <- as.factor(test_cases$temp_age_var)
    controls_temp$temp_age_var <- as.factor(controls_temp$temp_age_var)
    
  
    ## linear mixed model
    fmla <- as.formula(paste("outcome ~ ", paste(cg_cols, collapse= "+")))
    model <- glmmLasso(fmla, 
                       rnd = list(temp_age_var=~1),
                       lambda=10, 
                       data = train_cases)
      
    
    
    
    temp_models[[i]] <- model
    temp_cases_outcome[[i]] <- predict(model, as.data.frame(test_cases), type = 'class')
    temp_controls_outcome[[i]] <- predict(model, as.data.frame(controls_temp))
    temp_clin[[i]] <- test_clin
    temp_con_clin[[i]] <- con_clin
    
    
    print(i)
  }
  
  
  
  return(list(temp_models, temp_cases_outcome, temp_controls_outcome, temp_clin, temp_con_clin))
  
}

##########
# fixed variables
##########


# run full pipeline
full_results <- run_model(full_data_combat,
                          wt_data = data_wt,
                          bump_type = 'lfs_only',
                          gender = F,
                          tech = F,
                          fam_num = F,
                          fam_ratio = F,
                          remove_age_cgs_lit =F,
                          remove_age_cgs_lm = F,
                          k_folds = 5,
                          beta_thresh = 0.1)

# get restults 
cases_outcome <- unlist(full_results[[2]])
cases_clin <- do.call(rbind, full_results[[4]])
cases_results <- as.data.frame(cbind(cases_outcome, cases_clin))

colnames(cases_results)
cases_results[, c('cases_outcome', 'age_diagnosis')]

cases_results$test_pred_label <- ifelse(cases_results$cases_outcome > .5, 1, 0)
cases_results$test_label <- ifelse(cases_results$age_diagnosis < 72, 1, 0)


cases_results$pred_is <- ifelse(cases_results$test_pred_label== cases_results$test_label,
                             'good',
                             'bad')

cases_results <- cases_results[order(cases_results$test_pred, decreasing = T),]

length(which(cases_results$pred_is == 'good'))



# get restults 
con_outcome <- unlist(full_results[[3]])
con_clin <- do.call(rbind, full_results[[5]])
con_results <- as.data.frame(cbind(con_outcome, con_clin))

colnames(con_results)
con_results <- con_results[, c('con_outcome', 'age_sample_collection')]

# get person identfier
con_results$p_id <- rep.int(seq(1, 41, 1), 5)

# group by fold and get mean
con_results_pred <- con_results %>%
  group_by(p_id) %>%
  summarise(mean_pred = mean(con_outcome, na.rm =T)) %>%
  cbind(con_results[1:41,])


# remove original prediction
con_results_pred$controls_age_pred <- NULL

con_results_pred <- con_results_pred[order(con_results_pred$age_sample_collection, decreasing = F),]


con_results$test_pred_label <- ifelse(con_results$con_outcome > .5, 1, 0)
con_results$test_label <- ifelse(con_results$age_sample_collection < 72, 1, 0)


con_results$pred_is <- ifelse(con_results$test_pred_label== con_results$test_label,
                                'good',
                                'bad')

con_results <- con_results[order(con_results$test_pred, decreasing = T),]

length(which(con_results$pred_is == 'good'))




