
# source all_functions.R to load libraries and my functions
source('all_functions.R')

##########
# get base directory for 4 batch
##########
path_to_cases_tor <- '../../Data/methyl_data/cases_toronto'
path_to_cases_mon <- '../../Data/methyl_data/cases_montreal'
path_to_controls <- '../../Data/methyl_data/controls'
path_to_valid <- '../../Data/methyl_data/validation'


##########
# read in meth array - Data/methyl_data/cases_toronto, cases_montreal, controls, validation
##########

# cases 
rgCasesT <- read.metharray.exp(path_to_cases_tor, recursive = T)
rgCasesM <- read.metharray.exp(path_to_cases_mon, recursive = T)

# combine cases arrays 
rgCases <- combineArrays(rgCasesT, rgCasesM)
rm(rgCasesT, rgCasesM)

# controls
rgControls <- read.metharray.exp(path_to_controls, recursive = T)

rgValid <- read.metharray.exp(path_to_valid, recursive = T)

##########
# load genomic methyl set (from controls) - you need genetic locations by probe from this object
##########
ratio_set <- readRDS('../../Data/model_data/raw_ratio_set.rda')

# get granges object
g_ranges <- as.data.frame(getLocations(ratio_set))

# get probes from rownames
g_ranges$probe <- rownames(g_ranges)

# remove ch and duplicatee
g_ranges <- g_ranges[!duplicated(g_ranges$start),]
g_ranges <- g_ranges[!grepl('ch', g_ranges$probe),]

##########
# read in clinical data
##########
clin <- read.csv('../../Data/clin_data/clinical_two.csv', stringsAsFactors = F)

# clean clinical ids
clin$ids <-  gsub('A|B|_|-', '', clin$blood_dna_malkin_lab_)


##########
# cases 
##########

# cases batch1
id_map_tor <- read.csv(paste0(path_to_cases_tor, '/SampleSheet.csv'), stringsAsFactors = F)

#cases batch2
id_map_mon <- read.csv(paste0(path_to_cases_mon, '/SampleSheet.csv'), stringsAsFactors = F)
id_map_mon$Project <- NULL

# combine id_map and id_map_other
id_map <- rbind(id_map_tor, id_map_mon)
rm(id_map_tor, id_map_mon)

# clean id map
id_map <- cleanIdMap(id_map)


##########
# Controls batch1
##########
id_map_con <- read.csv(paste0(path_to_controls, '/SampleSheet.csv'), stringsAsFactors = F)

# clean idmap
id_map_con <- cleanIdMap(id_map_con)

##########
# valid
##########
id_map_val <- read.csv(paste0(path_to_valid, '/SampleSheet.csv'))

# homogenize valid map data with cases and controls
id_map_val <- id_map_val[, c('Sample.ID', 'Sample.Group', 'Sentrix.Barcode', 'Sample.Section',
                             'Project', 'Pool_ID', 'Sample_Well')]

# sub out '.' for '_'
colnames(id_map_val) <- gsub('.', '_', colnames(id_map_val), fixed = T)

# change 'Sample_ID' to 'Sample_Name' and 'Sentrix_Barcode' to 'Sentrix_ID'
colnames(id_map_val)[1] <- 'Sample_Name'
colnames(id_map_val)[3] <- 'Sentrix_ID'
colnames(id_map_val)[4] <- 'Sentrix_Position'
colnames(id_map_val)[5] <- 'Sample_Plate'

# clean idmap
id_map_val <- cleanIdMap(id_map_val)



##########
# remove outliers (previously determined) from rgset before normalization
##########
rgControls <- remove_outliers(rgSet = rgControls,
                              id_map = id_map_con,
                              method = 'doesnt_matter',
                              type = 'controls')

rgValid <- remove_outliers(rgSet = rgValid,
                           id_map = id_map_val,
                           method = 'doesnt_matter',
                           type = 'valid')


##########
# subset data - remove controls probes on each data set only if raw preprocessing
##########

# cases
rg_cases <- subset_rg_set(rg_set = rgCases,
                          keep_gender = F,
                          keep_controls = T,
                          keep_snps = T,
                          get_island = NULL,
                          get_chr = NULL,
                          get_type = NULL)

# controls
rg_controls <- subset_rg_set(rg_set = rgControls,
                             keep_gender = F,
                             keep_controls = T,
                             keep_snps = T,
                             get_island = NULL,
                             get_chr = NULL,
                             get_type = NULL)

# valid
rg_valid <- subset_rg_set(rg_set = rgValid,
                          keep_gender = F,
                          keep_controls = T,
                          keep_snps = T,
                          get_island = NULL,
                          get_chr = NULL,
                          get_type = NULL)



rg_cases = rg_cases
rg_controls = rg_controls
rg_valid = rg_valid
k_folds = 5
m_beta_thresh = 0.1
method = 'noob'

full_pipeline <- function(rg_cases, 
                          rg_controls, 
                          rg_valid, 
                          k_folds,
                          m_beta_thresh,
                          method,
                          controls) {
  
  # list to store cv results
  results_list <- list()
  clin_results <- list()
  
  # get vector of random folds
  fold_vec <- sample(1:k_folds, dim(rg_cases)[2], replace = T)
  
  # arguments or functions for subsetting data in different ways 
  # for example, remove chromosome, type I and type II probes, Chr17, remove Ch6
  
  # preprocess controls and valid
  beta_controls <- preprocessMethod(rg_controls, preprocess = method)
  beta_valid <- preprocessMethod(rg_valid, preprocess = method)
  # save.image('~/Desktop/temp_rg_data.RData')
  # load('~/Desktop/temp_rg_data.RData')
  
  # get controls old 
  
  for(i in 1:k_folds) {
    
    # get train and test index
    train_index <- !grepl(i, fold_vec)
    test_index <- !train_index
    
    # get training rgsets 
    rg_train <- rg_cases[,train_index]
    rg_test <- rg_cases[,test_index]
    
    # preprocess cases
    beta_train <- preprocessMethod(rg_train, preprocess = method)
    beta_test <- preprocessMethod(rg_test, preprocess = method)
    
    # do cases first (will return list of 2, second element is old controls)
    beta_train_cases <- process_rg_set_single(beta_data = beta_train, 
                                              id_map = id_map, 
                                              clin = clin)
    
    # get test set 
    beta_test_cases <- process_rg_set_single(beta_data = beta_test, 
                                             id_map = id_map, 
                                             clin = clin)
    
    # get controls
    beta_controls_mod <- process_rg_set_single(beta_data = m_controls, 
                                               id_map = id_map_con, 
                                               clin = clin)
    
    # get valid
    beta_valid_mod <- process_rg_set_single(beta_data = m_valid, 
                                            id_map = id_map_val, 
                                            clin = clin)
    
    ##########
    # function to remove infinite values
    ##########
    
    # cases
    m_train_cases <- removeInf(m_train_cases, probe_start = 10)
    m_test_cases <- removeInf(m_test_cases, probe_start = 10)
    
    # controls 
    m_controls_mod <- removeInf(m_controls_mod, probe_start = 10)
    
    # valid 
    m_valid_mod <- removeInf(m_valid_mod, probe_start = 10)
    
    ##########
    # bumphunter 
    ##########
    intersect_names <- Reduce(intersect, list(colnames(m_train_cases)[10:ncol(m_train_cases)],
                                              colnames(m_test_cases)[10:ncol(m_test_cases)],
                                              colnames(m_controls_mod)[10:ncol(m_controls_mod)],
                                              colnames(m_valid_mod)[10:ncol(m_valid_mod)]))
    
    
    # train data 
    m_train_cases <- m_train_cases[, c('ids',
                                       'p53_germline',
                                       'cancer_diagnosis_diagnoses',
                                       'age_diagnosis',
                                       'age_sample_collection',
                                       'gender',
                                       'sentrix_id',
                                       'family_name',
                                       intersect_names)]
    
    # test data 
    m_test_cases <- m_test_cases[, c('ids',
                                     'p53_germline',
                                     'cancer_diagnosis_diagnoses',
                                     'age_diagnosis',
                                     'age_sample_collection',
                                     'gender',
                                     'sentrix_id',
                                     'family_name',
                                     intersect_names)]
    
    # controls data 
    m_controls_mod <- m_controls_mod[, c('ids',
                                         'p53_germline',
                                         'cancer_diagnosis_diagnoses',
                                         'age_diagnosis',
                                         'age_sample_collection',
                                         'gender',
                                         'sentrix_id',
                                         'family_name',
                                         intersect_names)]
    
    # train data 
    m_valid_mod <- m_valid_mod[, c('ids',
                                   'p53_germline',
                                   'cancer_diagnosis_diagnoses',
                                   'age_diagnosis',
                                   'age_sample_collection',
                                   'gender',
                                   'sentrix_id',
                                   'family_name',
                                   intersect_names)]
    
    # get old controls (mut and no cancer from cases)
    m_controls_mod_old <- rbind(subset(m_train_cases, p53_germline == 'Mut' & 
                                         cancer_diagnosis_diagnoses == 'Unaffected'),
                                subset(m_test_cases, p53_germline == 'Mut' & 
                                         cancer_diagnosis_diagnoses == 'Unaffected'))
    
    # get model data in cases for training and test
    m_train_cases <- getModData(m_train_cases)
    m_test_cases <- getModData(m_test_cases)
    
    # full data 
    m_train_cases <- remove_dups(m_train_cases, m_test_cases)[[1]]
    m_test_cases <- remove_dups(m_train_cases, m_test_cases)[[2]]
    
    # get rid of cancer samples in controls 
    m_controls_mod <- m_controls_mod[grepl('Unaffected', m_controls_mod$cancer_diagnosis_diagnoses),]
    
    #subset valid - get ids from train and test
    case_ids <- append(m_train_cases$ids, m_test_cases$ids)
    m_valid_mod <- m_valid_mod[!m_valid_mod$ids %in% case_ids,]
    
    # determine which controls will be used
    if(controls == 'old') {
      m_controls_mod <- m_controls_mod_old
      
    }
    
    if (controls == 'full') {
      m_controls_mod <- rbind(m_controls_mod, m_controls_mod_old)
      
    } 
    
    ##########
    # remove max ages 
    ##########
    
    # remove ages for training 
    m_train_cases <- remove_ages(m_train_cases, 
                                 max_age = case_max_age)
    
    # remove ages for testing 
    m_test_cases <- remove_ages(m_test_cases, 
                                max_age = case_max_age)
    
    # remove ages for testing 
    m_controls_mod <- remove_ages(m_controls_mod, 
                                  max_age = con_max_age)
    
    # remove ages for testing 
    m_valid_mod <- remove_ages(m_valid_mod, 
                               max_age = val_max_age)
    
    ##########
    # use cases training and controls to get bumphunter features
    ##########
    bh_feats <- bump_hunter(dat_1 = m_train_cases, 
                            dat_2 = m_controls_mod, 
                            bump = 'cancer', 
                            boot_num = 3, 
                            m_beta_thresh = m_beta_thresh,
                            g_ranges = g_ranges)
    
    
    # get feature list
    colnames(bh_feats)[1] <- 'chr'
    remove_features <- inner_join(bh_feats, g_ranges)$probe
    
    # take remove features out of colnames 
    bh_features <- intersect_names[!intersect_names %in% remove_features]
    
    # get random features of length bh_features
    rand_feats <- sample(bh_features, length(bh_features))
    
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
    
    print(i)
    
  }
  
  #combine all the folds 
  results_final <- do.call(rbind, results_list)
  results_final_clin <- do.call(rbind, clin_results)
  
  return(list(results_final, results_final_clin))
}


##########
# fixed variables
##########

method = 'noob'
k_folds = 5
max_age = 1000
max_columns = 10000
m_beta_thresh = 0.5
combined = T
combined_data_type = 'controls'

#### TO FIX AND CHECK
# remove sex chromosomes and just control for gender
# check to see that you are removing cancer probes from bumphunter in the pipeline
# age, cotrols  HERE HHHHEREEREREGHETRRE use all of data (valid with cases and old controls), remove age signature with bumphunter
# run full pipeline 
full_results <- full_pipeline(rg_cases = rg_cases, 
                              rg_controls = rg_controls, 
                              rg_valid = rg_valid,
                              k_folds = k_folds,
                              m_beta_thresh = m_beta_thresh,
                              method = method, 
                              controls = 'normal')

# save results 
saveRDS(full_results, paste0('~/Desktop/', method, '_', k_folds, '_', 
                             max_age, '_', max_columns, '_', 
                             m_beta_thresh, '_', combined, '_', 
                             combined_data_type, '48_normal_no_gender_full_results.rda'))
