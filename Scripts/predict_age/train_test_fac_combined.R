

if (combined) {
  
  # get training data - rg_train
  # get testing data here
  rg_test_con <- combineArrays(rg_test, rg_controls)
  rg_test_val <- combineArrays(rg_test, rg_valid)
  
  # preprocess rg sets 
  beta_test_con <- preprocessMethod(rg_test_con, preprocess = method, only_beta_values = T)
  beta_test_val <- preprocessMethod(rg_test_val, preprocess = method, only_beta_values = T)
  
  # training data controls
  beta_train_full <- process_rg_set_single(beta_data = beta_train[1:max_columns,], 
                                           id_map = id_map, 
                                           clin = clin)
  
  # get model data on training 
  beta_train_mod <- getModData(beta_train_full)
  
  
  # get controls data
  controls_data_list  <- combine_clean_split(combined_data = beta_test_con_full, 
                                             train_data = beta_train_full, 
                                             controls = T)
  
  # get validation data
  valid_data_list  <- combine_clean_split(combined_data = beta_test_val_full, 
                                          train_data = beta_train_full, 
                                          controls = F)
  
  
  if (combined_data_type == 'controls') {
    
    # get single datasets
    beta_train_cases <- controls_data_list[[1]]
    beta_test_cases <- controls_data_list[[2]]
    beta_controls <- controls_data_list[[3]]
    beta_valid <- controls_data_list[[3]]
    
    
  } else {
    
    # get single datasets
    beta_train_cases <- valid_data_list[[1]]
    beta_test_cases <- valid_data_list[[2]]
    beta_controls <- valid_data_list[[3]]
    beta_valid <- valid_data_list[[3]]
    
    beta_valid <- beta_valid[complete.cases(beta_valid),]
    
    
  }
  
  ##########
  # use cases training and controls to get bumphunter features
  ##########
  
  # get bumphunter cases and controls 
  bh_cases_train <- controls_data_list[[1]]
  bh_controls <- controls_data_list[[3]]
  
  # run bumphunter
  bh_feats <- bump_hunter(dat_1 = bh_cases_train, 
                          dat_2 = bh_controls, 
                          bump = 'cancer', 
                          boot_num = 3, 
                          m_beta_thresh = m_beta_thresh,
                          g_ranges = g_ranges)
  
  
  # get feature list
  colnames(bh_feats)[1] <- 'chr'
  remove_features <- inner_join(bh_feats, g_ranges)$probe
  
  ##########
  # bumphunter 
  ##########
  intersect_names <- Reduce(intersect, list(colnames(m_train_cases)[10:ncol(m_train_cases)],
                                            colnames(m_test_cases)[10:ncol(m_test_cases)],
                                            colnames(m_controls)[10:ncol(m_controls)],
                                            colnames(m_valid)[10:ncol(m_valid)]))
  
  # take remove features out of colnames 
  bh_features <- intersect_names[!intersect_names %in% remove_features]
  
  
  # make sure no Nas
  m_train_cases <- m_train_cases[complete.cases(m_train_cases),]
  m_test_cases <- m_test_cases[complete.cases(m_test_cases),]
  
  # controls and valid
  m_controls <- m_controls[!is.na(m_controls$age_sample_collection),]
  m_controls <- m_controls[!is.na(m_controls$gender),]
  
  # get random features of length bh_features
  rand_feats <- sample(bh_features, length(bh_features))
  
  # remove ages for training 
  m_train_cases <- remove_ages(m_train_cases, 
                               max_age = case_max_age)
  
  # remove ages for testing 
  m_test_cases <- remove_ages(m_test_cases, 
                              max_age = case_max_age)
  
  # remove ages for testing 
  m_controls_mod <- remove_ages(m_controls, 
                                max_age = con_max_age)
  
  # remove ages for testing 
  m_valid_mod <- remove_ages(m_valid, 
                             max_age = val_max_age)
  
  # get gender variable for each data set
  m_train_cases <- cbind(as.data.frame(class.ind(m_train_cases$gender)), m_train_cases)
  m_test_cases <- cbind(as.data.frame(class.ind(m_test_cases$gender)), m_test_cases)
  m_controls_mod <- cbind(as.data.frame(class.ind(m_controls_mod$gender)), m_controls_mod)
  m_valid_mod <- cbind(as.data.frame(class.ind(m_valid_mod$gender)), m_valid_mod)
  
  # function to predict with all test, controls, controls old, and valid
  mod_result <- runEnetRandFac(training_dat = m_train_cases, 
                               controls_dat = m_controls_mod,
                               valid_dat = m_valid_mod,
                               test_dat = m_test_cases, 
                               age_cutoff = 48,
                               bh_features = bh_features,
                               rand_feats = rand_feats,
                               gender = T)
  
  # returns test_stats, test_stats_age, test_stats_controls, test_stats_valid, test_stats_age_valid
  temp.results <- get_class_results(mod_result, 
                                    dims_of_dat = length(bh_features), 
                                    mod_name = 'enet', 
                                    feat_name = 'enet',
                                    seed_number = 1)
  
  results_list[[i]] <- temp.results[[1]]
  clin_results[[i]] <- temp.results[[2]]
  