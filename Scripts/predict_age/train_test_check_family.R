
# load data from desktop 
load('~/Desktop/m_values.RData')

# source all_functions.R to load libraries and my functions
source('all_functions.R')

remove_age_cgs_lit <- TRUE
remove_age_cgs_lm <- FALSE

##########
# subset data - remove controls probes on each data set only if raw preprocessing
##########

# get intersecting name
intersect_names <- Reduce(intersect, list(colnames(beta_cases)[12:ncol(beta_cases)],
                                          colnames(beta_controls_mod)[12:ncol(beta_controls_mod)],
                                          colnames(beta_valid_mod)[12:ncol(beta_valid_mod)]))

# get the probes that are associated with genes
intersect_names <- intersect(intersect_names, 
                             gene_probes)

# read in probes associated with age
if (remove_age_cgs_lit) {
  age_cgs <- readRDS('../../Data/age_probes.rda')
  intersect_names <- intersect_names[!intersect_names %in% age_cgs]
}


# remove age
if (remove_age_cgs_lm) {
  age_cgs_lm <- readRDS('../../Data/age_cgs_lm.rda')
  intersect_names <- intersect_names[!intersect_names %in% age_cgs_lm]
}

# get clin names
clin_names <- colnames(beta_cases)[1:11]


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

# remove duplicate tm_donor
beta_cases_full <- beta_cases_full[!duplicated(beta_cases_full$tm_donor_),]

# combine beta_controls and beta_controls_old
beta_controls_full <- rbind(beta_controls_mod_old, 
                            beta_controls_mod)
beta_controls_full <- beta_controls_full[!duplicated(beta_controls_full$tm_donor_),]


# add an indicator for 450 and 850
beta_cases_full$tech <- ifelse(grepl('^57|97', beta_cases_full$sentrix_id), 'a', 'b')
beta_controls_full$tech <- ifelse(grepl('^57|97', beta_controls_full$sentrix_id), 'a', 'b')

# get gender variable for each data set
beta_cases_full <- cbind(as.data.frame(class.ind(beta_cases_full$gender)), beta_cases_full)
beta_controls_full <- cbind(as.data.frame(class.ind(beta_controls_full$gender)), beta_controls_full)

# get tech variable for each data set
beta_cases_full <- cbind(as.data.frame(class.ind(beta_cases_full$tech)), beta_cases_full)
beta_controls_full <- cbind(as.data.frame(class.ind(beta_controls_full$tech)), beta_controls_full)

# remove na in both
beta_cases_full <- beta_cases_full[!is.na(beta_cases_full$age_sample_collection),]
beta_controls_full <- beta_controls_full[!is.na(beta_controls_full$age_sample_collection),]

temp_cases <- beta_cases_full[, c('ids', 'gender', 'family_name')]
temp_controls <- beta_controls_full[, c('ids', 'gender', 'family_name')]

# write_csv(temp_cases, '../../Data/cases_family.csv')
temp_cases <- read_csv('../../Data/cases_family.csv')
temp_cases$ids <- as.character(temp_cases$ids)

# remove unneeded objects
rm(beta_controls, beta_valid, id_map_cases, id_map_con, id_map_val)


# create a training and test set with as much overlap of families as possible 
# assign temp_cases train and test to beta_cases_full
# sort both data sets to ensure correct data
beta_cases_full <- beta_cases_full[order(beta_cases_full$ids),]
temp_cases <- temp_cases[order(temp_cases$ids),]

beta_cases_full$train_test <- temp_cases$same_family_test


# create a training and test set with no overlap 
beta_train <- beta_cases_full[beta_cases_full$train_test == 'train',]
beta_test <- beta_cases_full[beta_cases_full$train_test == 'test',]

# remove last two variables 
beta_train$tech <- beta_train$train_test <- 
  beta_test$tech <- beta_test$train_test <- NULL

beta_controls_full$tech <- NULL

# use cases training and controls to get bumphunter features
bh_feats <- bump_hunter(dat_1 = beta_train, 
                        dat_2 = beta_controls_full, 
                        bump = 'cancer', 
                        boot_num = 5, 
                        thresh = beta_thresh,
                        g_ranges = g_ranges)


# get feature list
colnames(bh_feats)[1] <- 'chr'
remove_features <- inner_join(bh_feats, g_ranges)$probe

# take remove features out of colnames 
bh_features <- intersect_names[!intersect_names %in% remove_features]

# function to predict with all test, controls, controls old, and valid
mod_result  <- run_enet_family(training_dat = beta_train,
                               test_dat = beta_test,
                               age_cutoff = 72,
                               gender = TRUE,
                               tech = TRUE,
                               bh_features = bh_features)


##########
# examin cases with prediction objects (ROC, TRP, etc)
##########
library(ROCR)

temp_pred_cases <- prediction(mod_result$test_pred, as.character(mod_result$test_label))
temp_roc <- performance(temp_pred_cases, measure = 'tpr', x.measure = 'tnr')
cutoffs <- data.frame(cut=temp_roc@alpha.values[[1]], tpr=temp_roc@x.values[[1]], 
                      tnr=temp_roc@y.values[[1]])

plot(cutoffs$cut, cutoffs$tpr,type="l",col="red")
par(new=TRUE)
plot(cutoffs$cut, cutoffs$tnr,type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="")


temp_pr <- performance(temp_pred_cases, measure = 'prec', x.measure = 'rec')
temp_lc <- performance(temp_pred_cases, measure="lift")

plot(temp_roc)
plot(temp_pr)
plot(temp_lc)
caret::confusionMatrix(round(temp_cases$test_pred, 0.1), temp_cases$test_label)




