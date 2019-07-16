

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
  load(paste0(data_dir,paste0(data_used,'_',methyl_type, '_final_outlier', '.RData')))
} else {
  source(paste0('get_combat.R'))
}

# read in wt data
data_wt <- readRDS(paste0(data_dir,paste0(data_used,'_',methyl_type, '_wild_type', '.rda')))


full_data <- full_data[full_data$tm_donor_ != '3955',]
full_data_combat <- full_data[full_data_combat$tm_donor_ != '3955',]


# source all_functions.R to load libraries and my functions
source('all_functions.R')
##########
# run model
##########
# https://rdrr.io/rforge/glmmixedlasso/man/glmmlasso.html


run_model <- function(data_full,
                      wt_data,
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
  
  # get vector of random folds
  fold_vec <- sample(1:k_folds, nrow(data_full), replace = T)
  
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
    train_data <- data_full[train_index, ]
    test_data <- data_full[test_index, ]
  
    # use cases training and controls to get bumphunter features
    bh_feats <- bump_hunter_surv(dat = train_data, 
                            wild_type = wt_data,
                            boot_num = 5, 
                            thresh = beta_thresh,
                            g_ranges = g_ranges)
    
    # get feature list
    colnames(bh_feats)[1] <- 'chr'
    these_features <- inner_join(bh_feats, g_ranges)$probe
    
    # subset data by these_features
    train_data <- train_data[, c(clin_data_names, these_features)]
    test_data <- test_data[, c(clin_data_names, these_features)]

       # function to predict with all test, controls, controls old, and valid
    mod_result  <- run_enet_surv(training_dat = train_data,
                                 test_dat = test_data,
                                 age_cutoff = 72,
                                 gender = gender,
                                 tech = tech,
                                 fam_num = fam_num,
                                 fam_ratio = fam_ratio,
                                 bh_features = these_features)
    
    
    temp_cases[[i]] <- mod_result[[1]]
    temp_models[[i]] <- mod_result[[2]]

 
  }
  temp_results <- do.call(rbind, temp_cases)
  
  return(list(temp_results, temp_models))
}

##########
# fixed variables
##########


# run full pipeline
full_results <- run_model(data_full = full_data,
                          wt_data = data_wt,
                          gender = F,
                          tech = T,
                          fam_num = F,
                          fam_ratio = F,
                          remove_age_cgs_lit = T,
                          remove_age_cgs_lm = F,
                          k_folds = 5,
                          beta_thresh = 0.5)

temp <- full_results[[1]]

temp <- temp[ ,c('test_pred', 'age_diagnosis', 'age_sample_collection', 'cancer_diagnosis_diagnoses')]
temp <- temp[order(temp$test_pred),]
