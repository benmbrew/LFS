

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
if(paste0(data_used,'_',methyl_type, '_final_outlier', '.RData') %in% dir(data_dir)) {
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

run_model <- function(data_full,
                      enet,
                      gender,
                      tech,
                      fam_num,
                      fam_ratio,
                      k_folds = k_folds,
                      beta_thresh = beta_thresh) {

  
  intersect_names <- colnames(data_full)
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
  temp_full <- list()
  temp_alpha <- list()
  temp_models <- list()
  temp_importance <- list()
  
  # combine 
  for(i in 1:k_folds) {
    
    # get train and test index
    train_index <- !grepl(i, fold_vec)
    test_index <- !train_index
    
    train_full <- data_full[train_index, ]
    test_full <- data_full[test_index, ]
    
    
    # use cases training and controls to get bumphunter features
    bh_feats <- bump_hunter_lfs(dat = train_full, 
                                wild_type = data_wt,
                                boot_num = 5, 
                                thresh = beta_thresh,
                                g_ranges = g_ranges)
    
    
    # get feature list
    colnames(bh_feats)[1] <- 'chr'
    remove_features <- inner_join(bh_feats, g_ranges)$probe
    # get features
    temp_features <- colnames(train_full)[18:ncol(train_full)]
    # take remove features out of colnames 
    remaining_features <- temp_features[!temp_features %in% remove_features]
    # function to predict with all test, controls, controls old, and valid
    
    if(enet){
      mod_result  <- pred_cancer_enet(train_dat = train_full,
                                      test_dat = test_full,
                                      age_cutoff = 72,
                                      gender = gender,
                                      tech = tech,
                                      fam_num = fam_num,
                                      fam_ratio = fam_ratio,
                                      bh_features = remaining_features)
      
      
      temp_models[[i]] <- mod_result[[1]]
      temp_full[[i]] <- mod_result[[2]]
      temp_alpha[[i]] <- mod_result[[3]]
      
      
    } else {
      mod_result  <- pred_cancer_rf(train_dat = train_full,
                                    test_dat = test_full,
                                    age_cutoff = 72,
                                    gender = gender,
                                    tech = tech,
                                    fam_num = fam_num,
                                    fam_ratio = fam_ratio,
                                    bh_features = remaining_features)
      
     
      temp_models[[i]] <- mod_result[[1]]
      temp_full[[i]] <- mod_result[[2]]
      temp_importance[[i]] <- mod_result[[3]]
    }
   
    print(i)
  }
  
  if(enet) {
    test_results <- do.call(rbind, temp_full)

    return(list(test_results, temp_models, temp_alpha))
  } else {
    
    
    test_results <- do.call(rbind, temp_full)
    
    return(list(test_results, temp_models, temp_importance))
    
  }

}
  
  



##########
# fixed variables
##########


# run full pipeline
full_results <- run_model(data_full = full_data,
                          enet = F,
                          gender = T,
                          tech = T,
                          fam_num = F,
                          fam_ratio = F,
                          k_folds = 5,
                          beta_thresh = 1)
# get results from list
temp_cases <- full_results[[1]]

temp_cases <- temp_cases[ , c('test_pred.yes', 'test_label',
                              'age_diagnosis' ,  'age_sample_collection')]
temp_cases$test_pred_label <- ifelse(temp_cases$test_pred > .5, 1, 0)

temp_cases$pred_is <- ifelse(temp_cases$test_pred_label == temp_cases$test_label,
                             'good',
                             'bad')

temp_cases <- temp_cases[order(temp_cases$test_pred.yes, decreasing = T),]

temp_cases_young <- temp_cases[temp_cases$age_sample_collection < 216,]

# temp_cases <- temp_cases[order(temp_cases$test_pred, decreasing = T),]


##########
# examin cases with prediction objects (ROC, TRP, etc)
##########

temp_pred_cases <- prediction(temp_cases$test_pred.yes, temp_cases$test_label)
performance(temp_pred_cases, measure = 'acc')
temp_auc <- performance(temp_pred_cases, measure = 'auc')
temp_roc <- performance(temp_pred_cases, measure = 'tpr', x.measure = 'tnr')
temp_pr <- performance(temp_pred_cases, measure = 'prec', x.measure = 'rec')
temp_lc <- performance(temp_pred_cases, measure="lift")
hist(temp_cases$test_pred, main = 'Predicted Risk Score on Cases', xlab = 'Risk scores', col = 'darkgrey')

cutoffs <- data.frame(cut=temp_roc@alpha.values[[1]], tpr=temp_roc@x.values[[1]],
                      tnr=temp_roc@y.values[[1]])

plot(cutoffs$cut, cutoffs$tpr,type="l",col="red", xlab ='Cutoffs', ylab ='' ,main = 'TPR vs TNR by Cutoff')
par(new=TRUE)
plot(cutoffs$cut, cutoffs$tnr,type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="")


plot(temp_roc, lwd = 2.5, col = 'darkgrey',bty= 'n')
plot(temp_pr,lwd = 2.5, col = 'darkgrey',bty= 'n')
plot(temp_lc)
caret::confusionMatrix(round(temp_cases$test_pred.yes, 0.1), temp_cases$test_label)



