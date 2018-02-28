



run_model <- function(cases_full,
                      controls_full,
                      k_folds = k_folds,
                      beta_thresh = beta_thresh) {

  
  surv_results <- list()
  
  full_data <- rbind(beta_cases_full,
                     beta_controls_full)
  
  full_data <- full_data[full_data$age_sample_collection < 400,]
  
  fold_vec <- sample(1:k_folds, nrow(full_data), replace = T)
  
  
  
  # combine 
  for(i in 1:k_folds) {
    
    # get train and test index
    train_index <- !grepl(i, fold_vec)
    test_index <- !train_index
    
    # split into training and test 
    full_train_cases <- full_data[train_index, ]
    full_test_cases <- full_data[test_index, ]
    
    # function to predict with all test, controls, controls old, and valid
    surv_results[[i]]  <- run_coxreg(training_dat = full_train_cases,
                                     test_dat = full_test_cases,
                                     age_cutoff = age_cutoff,
                                     random_forest = random_forest,
                                     rf_surv_fac = rf_surv_fac,
                                     rf_surv_con = rf_surv_con,
                                     gender = gender,
                                     tech = tech,
                                     base_change = base_change,
                                     exon_intron = exon_intron,
                                     intersect_names = intersect_names)
    
    
    
    
    print(i)
  }
  
  
  
  surv_final <- do.call(rbind, surv_results)
  return(surv_final)
  
  
}

