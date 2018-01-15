

training_dat = full_train_cases
test_dat = full_test_cases
age_cutoff = age_cutoff
gender = gender
tech = tech
base_change = base_change
exon_intron = exon_intron

run_coxreg <- function(training_dat, 
                       controls_dat,
                       test_dat,
                       age_cutoff,
                       gender,
                       tech,
                       base_change,
                       exon_intron) {
  
  
  # get survival time as days to onset and days to sample collection in one column for training dat
  time_to_event <- training_dat$age_diagnosis
  missing_ind <- is.na(time_to_event)
  time_to_collection <- training_dat$age_sample_collection
  
  time_to_event[missing_ind] <- time_to_collection[missing_ind]
  
  # get cancer status 
  cancer_status <- training_dat$cancer_diagnosis_diagnoses != 'Unaffected'
  
  test_clin <- test_dat[, 14:23]

  intersected_feats <- intersect(intersect_names, colnames(training_dat))
  
  # get test y 
  test_y <- test_dat$age_
  
  if(gender) {
    
    intersected_feats <- append('M', intersected_feats)
    intersected_feats <- append('F', intersected_feats)
  }
  
  if (tech) {
    
    intersected_feats <- append('a', intersected_feats)
    intersected_feats <- append('b', intersected_feats)
  }
  
  if (base_change){
    
    
    intersected_feats <- append('none', intersected_feats)
    intersected_feats <- append('A', intersected_feats)
    intersected_feats <- append('C', intersected_feats)
    intersected_feats <- append('G', intersected_feats)
    intersected_feats <- append('T', intersected_feats)
  }
  
  if(exon_intron) {
    
    intersected_feats <- append('exon', intersected_feats)
    intersected_feats <- append('intron', intersected_feats)
    intersected_feats <- append('not_clear', intersected_feats)
    
  }
  
  # get bumphunter features
  training_dat <- training_dat[, intersected_feats]
  test_dat <- test_dat[, intersected_feats]
  
  # get survival object
  surv_outcome <- Surv(time_to_event,cancer_status, type = 'right')
  
  # fit coxph
  
  # start elastic net tuning
  N_CV_REPEATS = 2
  nfolds = 3
  
  ###### ENET
  # create vector and list to store best alpha on training data. alpha is the parameter that choses the 
  # the optimal proportion lambda, the tuning parameter for L1 (ridge) and L2 (lasso)
  elastic_net.cv_error = vector()
  elastic_net.cv_model = list()
  elastic_net.ALPHA <- c(1:9) / 10 # creates possible alpha values for model to choose from
  
  # set parameters for training model
  type_family <- 'cox'
  type_measure <- 'deviance'
  
  # create error matrix for for opitmal alpha that can run in parraellel if you have bigger data 
  # or if you have a high number fo N_CV_REPEATS
  temp.cv_error_matrix <- foreach (temp = 1:N_CV_REPEATS, .combine=rbind, .errorhandling="stop") %do% {      
    for (alpha in 1:length(elastic_net.ALPHA)) # for i in 1:9 - the model will run 9 times
    {      
      elastic_net.cv_model[[alpha]] = cv.glmnet(x = as.matrix(training_dat)
                                                , surv_outcome
                                                , alpha = elastic_net.ALPHA[alpha] # first time with 0.1 and so on
                                                , type.measure = type_measure
                                                , family = type_family
                                                , standardize = FALSE 
                                                , nfolds = nfolds 
                                                , nlambda = 10
                                                , parallel = TRUE
      )
      elastic_net.cv_error[alpha] = min(elastic_net.cv_model[[alpha]]$cvm)
    }
    elastic_net.cv_error # stores 9 errors    
  }
  
  if (N_CV_REPEATS == 1) {
    temp.cv_error_mean = temp.cv_error_matrix
  } else {
    temp.cv_error_mean = apply(temp.cv_error_matrix, 2, mean) # take the mean of the 5 iterations  
    # as your value for alpha
  }
  
  # stop if you did not recover error for any models 
  stopifnot(length(temp.cv_error_mean) == length(elastic_net.ALPHA))
  
  # get index of best alpha (lowest error) - alpha is values 0.1-0.9
  temp.best_alpha_index = which(min(temp.cv_error_mean) == temp.cv_error_mean)[length(which(min(temp.cv_error_mean) == temp.cv_error_mean))] 
  # print(paste("Best ALPHA:", elastic_net.ALPHA[temp.best_alpha_index])) # print the value for alpha
  best_alpha <- elastic_net.ALPHA[temp.best_alpha_index]
  temp.non_zero_coeff = 0
  temp.loop_count = 0
  # loop runs initially because temp.non_zero coefficient <3 and then stops 
  # usually after one iteration because the nzero variable selected by lambda is greater that 3. if it keeps looping
  # it they are never greater than 1, then the model does not converge. 
  while (temp.non_zero_coeff < 1) { 
    elastic_net.cv_model = cv.glmnet(x = as.matrix(training_dat)
                                     , surv_outcome
                                     , alpha = elastic_net.ALPHA[temp.best_alpha_index]
                                     , type.measure = type_measure
                                     , family = type_family
                                     , standardize=FALSE
                                     , nlambda = 100
                                     , nfolds = nfolds
                                     , parallel = TRUE
    )
    
    # get optimal lambda - the tuning parameter for ridge and lasso
    # THIS IS IMPORTANT BECAUSE WHEN YOU TRAIN THE MODEL ON 100 SEPERATE VALUES OF LAMBDA
    # AND WHEN YOU TEST THE MODEL IT WILL RETURN PREDCITION FOR ALL THOSE VALUES (1-100). YOU NEED TO 
    # GRAB THE PREDICTION WITH SAME LAMBDA THAT YOU TRAINED ON. ITS ALL IN THE CODE, BUT JUST WANTED TO 
    # GIVE YOU REASONS
    temp.min_lambda_index = which(elastic_net.cv_model$lambda == elastic_net.cv_model$lambda.min) 
    
    # # number of non zero coefficients at that lambda    
    temp.non_zero_coeff = elastic_net.cv_model$nzero[temp.min_lambda_index] 
    temp.loop_count = temp.loop_count + 1
    
    # set seed for next loop iteration
    as.numeric(Sys.time())-> t 
    set.seed((t - floor(t)) * 1e8 -> seed) 
    if (temp.loop_count > 10) {
      print("diverged")
      temp.min_lambda_index = 50 # if it loops more than 5 times, then model did not converge
      break
    }
  }# while loop ends 
  # print(temp.non_zero_coeff)  
  
  model  = glmnet(x = as.matrix(training_dat)
                  , surv_outcome
                  ,alpha = elastic_net.ALPHA[temp.best_alpha_index]
                  ,standardize=FALSE
                  ,nlambda = 100
                  ,family = type_family)
  
  # Predictions on test data
  
  # This returns 100 prediction with 1-100 lambdas
  temp_test.predictions <- predict(model, 
                                   data.matrix(test_dat),
                                   type = 'response')
  
  # get predictions with corresponding lambda.
  test.predictions <- temp_test.predictions[, temp.min_lambda_index]
  
  # combine predictions and real labels 
  temp_dat <- as.data.frame(cbind(test_pred = test.predictions,  test_clin))
  
  
  ###########################################################################################
  return(temp_dat
  
  
}
