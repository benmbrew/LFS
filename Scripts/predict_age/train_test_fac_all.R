
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
id_map_cases <- rbind(id_map_tor, id_map_mon)
rm(id_map_tor, id_map_mon)

# clean id map
id_map_cases <- cleanIdMap(id_map_cases)


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

# save.image('~/Desktop/temp_450_850.RData')
# load('~/Desktop/temp_450_850.RData')

# cases
rg_cases <- subset_rg_set(rg_set = rgCases,
                          keep_gender = F,
                          keep_controls = T,
                          keep_snps = T,
                          get_island = "Island",
                          get_chr = NULL,
                          get_type = NULL)

# controls
rg_controls <- subset_rg_set(rg_set = rgControls,
                             keep_gender = F,
                             keep_controls = T,
                             keep_snps = T,
                             get_island = "Island",
                             get_chr = NULL,
                             get_type = NULL)

# valid
rg_valid <- subset_rg_set(rg_set = rgValid,
                          keep_gender = F,
                          keep_controls = T,
                          keep_snps = T,
                          get_island = "Island",
                          get_chr = NULL,
                          get_type = NULL)



rg_cases = rgCases
rg_controls = rgControls
rg_valid = rgValid
age_cutoff = 72
k_folds = 3
beta_thresh = 0.1
method = 'noob'

full_pipeline <- function(rg_cases, 
                          rg_controls, 
                          rg_valid, 
                          age_cutoff,
                          k_folds,
                          beta_thresh,
                          method,
                          controls) {
  
  # list to store cv results
  temp_cases <- list()
  temp_controls <- list()
  temp_valid <- list()
  
  
  # get vector of random folds
  fold_vec <- sample(1:k_folds, dim(rg_cases)[2], replace = T)
  
  # preprocess controls and valid
  beta_cases <-  preprocessMethod(rg_cases, preprocess = method)
  beta_controls <- preprocessMethod(rg_controls, preprocess = method)
  beta_valid <- preprocessMethod(rg_valid, preprocess = method)
  
  # get controls
  beta_cases <- process_rg_set_single(beta_data = beta_cases[1:50000,], 
                                      id_map = id_map_cases, 
                                      clin = clin)
  # get controls
  beta_controls_mod <- process_rg_set_single(beta_data = beta_controls[1:50000,], 
                                             id_map = id_map_con, 
                                             clin = clin)
  
  # get valid
  beta_valid_mod <- process_rg_set_single(beta_data = beta_valid[1:50000,], 
                                          id_map = id_map_val, 
                                          clin = clin)
  
  
  # get intersecting name
  intersect_names <- Reduce(intersect, list(colnames(beta_cases)[10:ncol(beta_cases)],
                                            colnames(beta_controls_mod)[10:ncol(beta_controls_mod)],
                                            colnames(beta_valid_mod)[10:ncol(beta_valid_mod)]))
  
  # get clin names
  clin_names <- colnames(beta_cases)[1:9]
  
  
  # train data 
  beta_cases <- beta_cases[, c(clin_names,
                               intersect_names)]
  
  # controls data 
  beta_controls_mod <- beta_controls_mod[, c(clin_names,
                                             intersect_names)]
  
  # train data 
  beta_valid_mod <- beta_valid_mod[, c(clin_names,
                                       intersect_names)]
  
  # get old controls (mut and no cancer from cases)
  beta_controls_mod_old <- rbind(subset(beta_cases, p53_germline == 'Mut' & 
                                          cancer_diagnosis_diagnoses == 'Unaffected'))
  
  # get model data in cases for training and test
  beta_cases <- getModData(beta_cases)
  
  # get rid of cancer samples in controls 
  beta_controls_mod <- beta_controls_mod[grepl('Unaffected', beta_controls_mod$cancer_diagnosis_diagnoses),]
  
  #subset valid - get ids from train and test
  case_ids <- beta_cases$ids
  beta_valid_mod <- beta_valid_mod[!beta_valid_mod$ids %in% case_ids,]
  
  # remove NAs 
  beta_cases <-beta_cases[complete.cases(beta_cases),]
  
  # combine beta cases and beta valid
  beta_cases_full <- rbind(beta_cases,
                           beta_valid_mod)
  
  # combine beta_controls and beta_controls_old
  beta_controls_full <- rbind(beta_controls_mod, 
                              beta_controls_mod_old)
  beta_controls_full <- beta_controls_full[!duplicated(beta_controls_full$ids),]
  
  # add an indicator for 450 and 850
  beta_cases_full$tech <- ifelse(grepl('^57|97', beta_cases_full$sentrix_id), 'a', 'b')
  beta_controls_full$tech <- ifelse(grepl('^57|97', beta_controls_full$sentrix_id), 'a', 'b')
  
  # get vector of random folds
  fold_vec <- sample(1:k_folds, nrow(beta_cases_full), replace = T)
  
  # combine 
  for(i in 1:k_folds) {
    
    # get train and test index
    train_index <- !grepl(i, fold_vec)
    test_index <- !train_index
    
    # split into training and test 
    beta_train_cases <- beta_cases_full[train_index, ]
    beta_test_cases <- beta_cases_full[test_index, ]
    
    ##########
    # use cases training and controls to get bumphunter features
    ##########
    bh_feats <- bump_hunter(dat_1 = beta_train_cases[, -ncol(beta_train_cases)], 
                            dat_2 = beta_controls_full[, -ncol(beta_controls_full)], 
                            bump = 'cancer', 
                            boot_num = 5, 
                            thresh = beta_thresh,
                            g_ranges = g_ranges)
    
    
    # get feature list
    colnames(bh_feats)[1] <- 'chr'
    remove_features <- inner_join(bh_feats, g_ranges)$probe
    
    # take remove features out of colnames 
    bh_features <- intersect_names[!intersect_names %in% remove_features]
    
    
    # get gender variable for each data set
    beta_train_cases <- cbind(as.data.frame(class.ind(beta_train_cases$gender)), beta_train_cases)
    beta_test_cases <- cbind(as.data.frame(class.ind(beta_test_cases$gender)), beta_test_cases)
    beta_controls_full <- cbind(as.data.frame(class.ind(beta_controls_full$gender)), beta_controls_full)
    
    # get tech variable for each data set
    beta_train_cases <- cbind(as.data.frame(class.ind(beta_train_cases$tech)), beta_train_cases)
    beta_test_cases <- cbind(as.data.frame(class.ind(beta_test_cases$tech)), beta_test_cases)
    beta_controls_full <- cbind(as.data.frame(class.ind(beta_controls_full$tech)), beta_controls_full)
    
    # function to predict with all test, controls, controls old, and valid
    mod_result  <- run_enet_450_850(training_dat = beta_train_cases,
                                    controls_dat = beta_controls_full,
                                    test_dat = beta_test_cases,
                                    age_cutoff = age_cutoff,
                                    bh_features = bh_features,
                                    gender = T,
                                    tech = T)
    
    temp_cases[[i]] <- mod_result[[1]]
    temp_controls[[i]] <- mod_result[[2]]
    
    print(i)
  }
  full_cases <- do.call(rbind, temp_cases)
  full_controls <- do.call(rbind, temp_controls)
  
  return(list(full_cases, full_controls, full_valid))
}


##########
# fixed variables
##########

age_cutoff = 48
method = 'noob'
k_folds = 5
beta_thresh = 0.1
control_type = 'full'

# run full pipeline 
full_results <- full_pipeline(rg_cases = rg_cases, 
                              rg_controls = rg_controls, 
                              rg_valid = rg_valid,
                              age_cutoff = age_cutoff,
                              k_folds = k_folds,
                              beta_thresh = beta_thresh,
                              method = method, 
                              controls = control_type)

# save results 
saveRDS(full_results, paste0('../../Data/results_data/', method, '_', age_cutoff, '_', control_type, '_', beta_thresh, '.rda'))

