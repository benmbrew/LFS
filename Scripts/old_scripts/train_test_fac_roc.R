
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
# remove outliers (previously determined) from rgset before normalization
##########
rgControls <- remove_outliers(rgSet = rgControls,
                              id_map = id_map_con,
                              method = 'doesnt_matter',
                              type = 'controls')

##########
# subset data - remove controls probes on each data set only if raw preprocessing
##########

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




save.image('~/Desktop/temp_beta_data.RData')
load('~/Desktop/temp_beta_data.RData')
source('all_functions.R')



full_pipeline <- function(rg_cases, 
                          rg_controls, 
                          k_folds,
                          beta_thresh,
                          method,
                          controls) {
  
  # list to store cv results
  results_list <- list()
  clin_results <- list()
  
  # get vector of random folds
  fold_vec <- sample(1:k_folds, dim(rg_cases)[2], replace = T)
  
  # preprocess controls and valid
  beta_controls <- preprocessMethod(rg_controls, preprocess = method)

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
    beta_train_cases <- process_rg_set_single(beta_data = beta_train[1:50000,], 
                                              id_map = id_map_cases, 
                                              clin = clin)
    
    # get test set 
    beta_test_cases <- process_rg_set_single(beta_data = beta_test[1:50000,], 
                                             id_map = id_map_cases, 
                                             clin = clin)
    
    # get controls
    beta_controls_mod <- process_rg_set_single(beta_data = beta_controls[1:50000,], 
                                               id_map = id_map_con, 
                                               clin = clin)
  
    if(method == 'raw'){
      # cases
      beta_train_cases <- removeInf(beta_train_cases, probe_start = 10)
      beta_test_cases <- removeInf(beta_test_cases, probe_start = 10)
      
      # controls 
      beta_controls_mod <- removeInf(beta_controls_mod, probe_start = 10)
      
    }
    
    ##########
    # bumphunter 
    ##########
    intersect_names <- Reduce(intersect, list(colnames(beta_train_cases)[10:ncol(beta_train_cases)],
                                              colnames(beta_test_cases)[10:ncol(beta_test_cases)],
                                              colnames(beta_controls_mod)[10:ncol(beta_controls_mod)]))
    
    # get clin names
    clin_names <- colnames(beta_train_cases)[1:9]
    
    
    # train data 
    beta_train_cases <- beta_train_cases[, c(clin_names,
                                             intersect_names)]
    
    # test data 
    beta_test_cases <- beta_test_cases[, c(clin_names,
                                           intersect_names)]
    
    # controls data 
    beta_controls_mod <- beta_controls_mod[, c(clin_names,
                                               intersect_names)]
  
    # get old controls (mut and no cancer from cases)
    beta_controls_mod_old <- rbind(subset(beta_train_cases, p53_germline == 'Mut' & 
                                            cancer_diagnosis_diagnoses == 'Unaffected'),
                                   subset(beta_test_cases, p53_germline == 'Mut' & 
                                            cancer_diagnosis_diagnoses == 'Unaffected'))
    
    # get model data in cases for training and test
    beta_train_cases <- getModData(beta_train_cases)
    beta_test_cases <- getModData(beta_test_cases)
    
    # full data 
    beta_train_cases <- remove_dups(beta_train_cases, beta_test_cases)[[1]]
    beta_test_cases <- remove_dups(beta_train_cases, beta_test_cases)[[2]]
    
    # get rid of cancer samples in controls 
    beta_controls_mod <- beta_controls_mod[grepl('Unaffected', beta_controls_mod$cancer_diagnosis_diagnoses),]
    
    # #subset valid - get ids from train and test
    # case_ids <- append(beta_train_cases$ids, beta_test_cases$ids)
    # beta_valid_mod <- beta_valid_mod[!beta_valid_mod$ids %in% case_ids,]
    # 
    # remove NAs 
    beta_train_cases <-beta_train_cases[complete.cases(beta_train_cases),]
    beta_test_cases <-beta_test_cases[complete.cases(beta_test_cases),]
    
    # determine which controls will be used
    if(controls == 'old') {
      beta_controls_mod <- beta_controls_mod_old
    }
    
    if (controls == 'full') {
      beta_controls_mod <- rbind(beta_controls_mod, beta_controls_mod_old)
    } 
    
    ##########
    # use cases training and controls to get bumphunter features
    ##########
    bh_feats <- bump_hunter(dat_1 = beta_train_cases, 
                            dat_2 = beta_controls_mod, 
                            bump = 'cancer', 
                            boot_num = 3, 
                            thresh = beta_thresh,
                            g_ranges = g_ranges)
    
    
    # get feature list
    colnames(bh_feats)[1] <- 'chr'
    remove_features <- inner_join(bh_feats, g_ranges)$probe
    
    # take remove features out of colnames 
    bh_features <- intersect_names[!intersect_names %in% remove_features]
    
    # get random features of length bh_features
    rand_feats <- sample(bh_features, length(bh_features))
    
    # get gender variable for each data set
    beta_train_cases <- cbind(as.data.frame(class.ind(beta_train_cases$gender)), beta_train_cases)
    beta_test_cases <- cbind(as.data.frame(class.ind(beta_test_cases$gender)), beta_test_cases)
    # beta_controls_mod <- cbind(as.data.frame(class.ind(beta_controls_mod$gender)), beta_controls_mod)
    # beta_valid_mod <- cbind(as.data.frame(class.ind(beta_valid_mod$gender)), beta_valid_mod)
  
    # function to predict with all test, controls, controls old, and valid
    mod_result <- runEnetRoc(training_dat = beta_train_cases, 
                             test_dat = beta_test_cases, 
                             age_cutoff = 72,
                             bh_features = bh_features,
                             gender = T)
    
  
  
  results_list[[i]] <- mod_result

  print(i)
  
  }
  #combine all the folds 
  results_final <- do.call(rbind, results_list)

return(results_final)
}

outcome_noob <- full_pipeline(rg_cases = rg_cases, 
                      rg_controls = rg_controls, 
                      k_folds = 3, 
                      beta_thresh = 0.1, 
                      method = 'noob', 
                      controls = 'normal')


# saveRDS(temp, '~/Desktop/temp_72_perf.rda')
temp <- readRDS('~/Desktop/temp_72_perf.rda')
# 

# get pred object from predictions and real values 
temp_pred <- prediction(temp$test_pred, temp$test_label)
temp_auc <- performance(temp_pred, measure = 'auc')
temp_roc <- performance(temp_pred, measure = 'tpr', x.measure = 'fpr')
temp_pr <- performance(temp_pred, measure = 'prec', x.measure = 'rec')
# temp_ss <- performance(temp_pred,  measure="spec", x.measure="sens")
temp_lc <- performance(temp_pred, measure="lift")


plot(temp_roc, axes=F, col = adjustcolor('black', alpha.f = 0.8), lty =1, lwd = 2.5)


pdf('~/Desktop/roc_pr.pdf')
# plot roc
plot(temp_roc, axes=F, col = adjustcolor('black', alpha.f = 0.8), lty =1, lwd = 2.5)
legend('bottomright' ,  bty = 'n', legend = paste0('auc', ' = ', round(unlist(temp_auc@y.values), 2)))
abline(a = 0, b = 1, lty = 2)

# plot precision recall
plot(temp_pr, axes=F, col = adjustcolor('black', alpha.f = 0.8), lty =1, lwd = 2.5, 
     xlim = c(0, 1), ylim = c(0, 1))

# plot lift curve
plot(temp_lc, axes=F, col = adjustcolor('black', alpha.f = 0.8), lty =1, lwd = 2.5)




# Assuming you already have a vector of probabilities (called probs) computed with your model and the 
# true class labels are in your data frame as df$label (0 and 1) this code should work:
fg <- temp$test_pred[temp$test_label == 1]
bg <- temp$test_pred[temp$test_label == 0]

pdf('~/Desktop/roc_pr.pdf')

# ROC Curve    
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
plot(roc, col = adjustcolor('black', alpha.f = 0.8), lty =1, lwd = 3, legend = F)
abline(0, 1)

# PR Curve
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
plot(pr, col = adjustcolor('black', alpha.f = 0.8), lty =1, lwd = 3, legend = F)






