# source all_functions.R to load libraries and my functions
source('all_functions.R')

library(coxme)
##########
# read in clinical data
##########
clin <- read.csv('../../Data/clin_data/clinical_two.csv', stringsAsFactors = F)

# clean clinical ids
clin$ids <-  gsub('A|B|_|-', '', clin$blood_dna_malkin_lab_)

#########
# classification 
#########

# options should include random effects, survival, normal classification, variables 
clin_dat <- clin

recode_clin <- function(clin_dat) {
  # keep only family_name, tm_donor_, p53_germline, gdna.base.change, gdna.exon.intron,gdna.codon, protein.codon.change, protein.codon.num, 
  # splice.delins.snv, codon72.npro,mdm2.nG 
  clin_dat <- clin_dat %>% select(ids, age_diagnosis, age_sample_collection, gender,family_name, cancer_diagnosis_diagnoses, tm_donor_, p53_germline, gdna.base.change, gdna.exon.intron, 
                                  gdna.codon, protein.codon.change, protein.codon.num, splice.delins.snv, 
                                  codon72.npro, mdm2.nG)
  
  # remove data without ids
  clin_dat <- clin_dat %>% 
    filter(!is.na(ids)) %>% 
    filter(p53_germline == 'Mut') %>% 
    filter(!duplicated(tm_donor_)) %>%
    filter(!duplicated(ids)) %>%
    filter(!is.na(cancer_diagnosis_diagnoses))
  
  # remove rows that are both NA
  clin_dat <- clin_dat[!is.na(clin_dat$age_diagnosis) | !is.na(clin_dat$age_sample_collection),]
  
  # gdna.base.change
  clin_dat$gdna_base_change <- ifelse(clin_dat$gdna.base.change == 'C>T',
                                      'C>T',
                                      ifelse(clin_dat$gdna.base.change == 'G>A', 
                                             'G>A',
                                             ifelse(clin_dat$gdna.base.change == 'no bp change',
                                                    'no_bp_change', 'other')))
  
  # recode gdna.exon.intron
  clin_dat$gdna_exon_intron <- ifelse(grepl('Exon', clin_dat$gdna.exon.intron), 'exon',
                                      ifelse(grepl('Intron', clin_dat$gdna.exon.intron), 'intron', 'not_clear'))
  
  # recode protein.codon.change
  clin_dat$protein_codon_change <- ifelse(clin_dat$protein.codon.change == 'Arg>Gln',
                                          'Arg>Gln', 
                                          ifelse(clin_dat$protein.codon.change == 'Arg>His',
                                                 'Arg>His',
                                                 ifelse(clin_dat$protein.codon.change == 'no_codon_change',
                                                        'none',
                                                        ifelse(clin_dat$protein.codon.change =='Thr>Thr',
                                                               'Thr>Thr',
                                                               ifelse(clin_dat$protein.codon.change == 'Arg>Cys',
                                                                      'Arg>Cys', 
                                                                      ifelse(is.na(clin_dat$protein.codon.change), 
                                                                             NA, 'other'))))))
  
  
  # transition a <-> g and c <-> t
  clin_dat$gdna_base_change_tran <- ifelse(grepl('A>G|G>A|C>T|T>C', clin_dat$gdna.base.change),
                                                   'transition', 
                                                   ifelse(grepl('A>C|C>A|G>T|T>G', clin_dat$gdna.base.change),
                                                          'transversion',
                                                          ifelse(grepl('no', clin_dat$gdna.base.change),
                                                                 'no_base_change', 'other')))

  # recode final 3 
  clin_dat$splice_snv <- clin_dat$splice.delins.snv
  clin_dat$codon72_npro <- clin_dat$codon72.npro
  clin_dat$mdm2_ng <- clin_dat$mdm2.nG
  
  # subset data 
  clin_dat <- clin_dat %>%
    select(age_diagnosis, age_sample_collection, gender,family_name, cancer_diagnosis_diagnoses, 
           gdna_base_change, gdna_base_change_tran,gdna_exon_intron, protein_codon_change, splice_snv, codon72_npro, mdm2_ng)
  
  return(clin_dat)
 
}

clin <- recode_clin(clin)

# subset to get only ACC and unaffectd
clin <- clin %>%
  filter(grepl('Unaffected|ACC', cancer_diagnosis_diagnoses))

# https://www.r-bloggers.com/cox-proportional-hazards-model/

##########
# fit coxph
##########

coxph_pred_clin <- function(clin_dat){
  # break into training and test set 
  k_folds <- 5
  fold_vec <- sample(1:k_folds, nrow(clin_dat), replace = T)
  
  # get list to store results
  results_list <- list()
  
  # combine 
  for(i in 1:k_folds) {
    
    # get train and test index
    train_index <- !grepl(i, fold_vec)
    test_index <- !train_index
    
    # get data sets
    clin_train  <- clin_dat[train_index,]
    clin_test <- clin_dat[test_index,]
    
    # get test outcome 
    test_y <- clin_test$cancer_diagnosis_diagnoses
    
    # get survival time as days to onset and days to sample collection in one column for training dat
    time_to_event <- clin_train$age_diagnosis
    missing_ind <- is.na(time_to_event)
    time_to_collection <- clin_train$age_sample_collection
    
    time_to_event[missing_ind] <- time_to_collection[missing_ind]
    
    # get cancer status 
    cancer_status <- ifelse(clin_train$cancer_diagnosis_diagnoses != 'Unaffected', 1, 0)
    
    # estimate clinph model.
    coxfit <-coxph(formula = Surv(time_to_event, cancer_status, type = 'right') ~ gender + gdna_base_change + gdna_exon_intron, data = clin_train)
    
    # predict on test set 
    pred_y <- predict(coxfit, clin_test, type = 'risk')
    
    # combine with real y
    results_dat <- as.data.frame(cbind(pred_y, test_y))
    results_list[[i]] <- results_dat
   
  }
  final_results <- do.call(rbind,results_list)
  return(final_results)
}

##########
# fit a random forest 
##########
clin_dat <- clin
colnames(clin)
random_forest_pred_clin <- function(clin_dat,
                                    age_cutoff,
                                    features){
  
  # get variables of interest 
  clin_dat <- clin_dat %>%
    select(age_diagnosis, age_sample_collection, gender, cancer_diagnosis_diagnoses, gdna_base_change, 
           gdna_exon_intron)
  
  # get controls data 
  clin_dat_controls <-clin_dat %>%
    filter(cancer_diagnosis_diagnoses == 'Unaffected') 
  clin_dat_controls$cancer_diagnosis_diagnoses <- NULL
  # get cases 
  clin_dat_cases <- clin_dat %>%
    filter(cancer_diagnosis_diagnoses != 'Unaffected')
  clin_dat_cases$cancer_diagnosis_diagnoses <- NULL
  
  
  # break into training and test set 
  k_folds <- 5
  fold_vec <- sample(1:k_folds, nrow(clin_dat_cases), replace = T)
  
  # get list to store results
  results_list <- list()
  
  # combine 
  for(i in 1:k_folds) {
    
    # get train and test index
    train_index <- !grepl(i, fold_vec)
    test_index <- !train_index
    
    # get data sets
    clin_train  <- clin_dat_cases[train_index,]
    clin_test <- clin_dat_cases[test_index,]
    
    # get test outcome 
    train_y <- make.names(ifelse(clin_train$age_diagnosis < age_cutoff, 1, 0))
    test_y <-  make.names(ifelse(clin_test$age_diagnosis < age_cutoff, 1, 0))
    
    # remove age from columns 
    clin_train <- clin_train[,!grepl('age', colnames(clin_train))]
    clin_test <- clin_test[,!grepl('age', colnames(clin_test))]
    
    
    # determines how you train the model.
    fitControl <- trainControl( 
      method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
      number = 4, 
      classProbs = TRUE,     
      repeats = 1,
      allowParallel = TRUE)
    
    # mtry: Number of variables randomly sampled as candidates at each split.
    # ntree: Number of trees to grow.
    mtry <- sqrt(ncol(clin_train))
    tunegrid <- expand.grid(.mtry=mtry)
    
    model  <- train(x = clin_train
                    , y =train_y
                    , method = "rf"
                    , metric = "ROC"
                    , trControl = fitControl
                    , tuneGrid = tunegrid
                    , importance = T
                    , verbose = FALSE)
    
    temp <- varImp(model)[[1]]
    importance <- cbind(rownames(temp), temp$Overall)
    
    # predict on test data
    test.predictions <- predict(model,
                                newdata = clin_test)
    
    # get controls
    test.predictions_controls <- predict(model,
                                         controls_dat)
    
    
    
    # get controls old
    test.predictions_controls_old <- predict(model,
                                             controls_dat_old)
    
    
    # get controls
    test.predictions_controls_full <- predict(model,
                                              controls_dat_full)
    
    
    # estimate coxph model.
    coxfit <-coxph(formula = Surv(time_to_event, cancer_status, type = 'right') ~ gender + gdna_base_change + gdna_exon_intron, data = clin_train)
    
    # predict on test set 
    pred_y <- predict(coxfit, clin_test, type = 'risk')
    
    # combine with real y
    results_dat <- as.data.frame(cbind(pred_y, test_y))
    results_list[[i]] <- results_dat
    
  }
  final_results <- do.call(rbind,results_list)
  return(final_results)
}



