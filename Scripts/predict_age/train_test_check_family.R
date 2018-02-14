
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

write_csv(temp_cases, '../../Data/cases_family.csv')



# create a training and test set with as much overlap of families as possible 



# create a training and test set with no overlap 




beta_train_cases <- beta_cases_full[train_index, ]
beta_test_cases <- beta_cases_full[test_index, ]


# use cases training and controls to get bumphunter features
bh_feats <- bump_hunter(dat_1 = beta_train_cases, 
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
mod_result  <- run_enet_450_850(training_dat = beta_train_cases,
                                controls_dat = beta_controls_full,
                                test_dat = beta_test_cases,
                                age_cutoff = age_cutoff,
                                gender = gender,
                                tech = tech,
                                base_change = base_change,
                                exon_intron = exon_intron,
                                bh_features = bh_features)



