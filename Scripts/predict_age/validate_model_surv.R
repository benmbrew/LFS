# source functions script
source('all_functions.R')

# create fixed objects to model and pipeline inputs and saving  
methyl_type = 'beta'
combat = TRUE
gender = TRUE
age_cutoff = 72
tech = TRUE
remove_age = TRUE
use_offset = TRUE
num_seeds = 'seeds_3'
num_folds = 'folds_5'

# condition on fixed objects to get saving identifiers

if(methyl_type == 'beta'){
  which_methyl <- 'beta'
} else {
  which_methyl <- 'm'
}

if(combat){
  which_combat <- 'use_combat'
  beta_thresh = 0.01
} else {
  which_combat <- 'no_combat'
  beta_thresh = 0.5
}

if(gender){
  is_gen = 'use_gen'
} else {
  is_gen = 'no_gen'
}

if(use_offset){
  is_offset <- 'use_offset'
} else {
  is_offset <- 'no_offset'
}

if(remove_age){
  is_removed <- 'age_removed'
} else {
  is_removed <- 'no_removed'
  
}

##########
# load genomic methyl set (from controls) - you need genetic locations by probe from this object
##########
ratio_set <- readRDS('../../Data/g_ranges.rda')

# get g_ranges
g_ranges <- ratio_set

# get probes from rownames
g_ranges$probe <- rownames(ratio_set)

# remove ch and duplicatee
g_ranges <- g_ranges[!duplicated(g_ranges$start),]
g_ranges <- g_ranges[!grepl('ch', g_ranges$probe),]
names(g_ranges)[1] <- 'chr'

###########################

# save trainig  set 
all_450 <-readRDS(paste0('validation_age_predictions_surv/', 'all_450_',which_methyl, '_', which_combat, '_',
                         num_seeds, '_', num_folds, '_', is_gen, '_', is_offset,'.rda'))

# save validation set 
all_850 <- readRDS(paste0('validation_age_predictions_surv/', 'all_850_',which_methyl, '_', which_combat, '_',
                          num_seeds, '_', num_folds, '_', is_gen,'_', is_offset, '.rda'))

# save associated lfs bumps
lfs_bump_probes <- readRDS(paste0('validation_age_predictions_surv/', 'lfs_bumps_', which_methyl, '_', 
                                  which_combat, '.rda'))
############################

temp_cases <- all_450[all_450$cancer_diagnosis_diagnoses != 'Unaffected',]
temp_con <- all_450[all_450$cancer_diagnosis_diagnoses == 'Unaffected',]

# run bumphunter on full_450 to remove cancer signature
# use cases training and controls to get bumphunter features
bh_feats <- bump_hunter(dat_1 = temp_cases, 
                        dat_2 = temp_con , 
                        bump = 'cancer', 
                        boot_num = 5, 
                        beta_thresh = beta_thresh,
                        methyl_type = methyl_type,
                        g_ranges = g_ranges)
rm(temp_con, temp_cases)


# get intersect_names
intersect_names <- names(all_450)[grepl('^cg', names(all_450))]

# get feature list
colnames(bh_feats)[1] <- 'chr'
remove_features <- inner_join(bh_feats, g_ranges)$probe

# take remove features out of colnames 
bh_features <- intersect_names[!intersect_names %in% remove_features]

# subset all data by bh_features
all_450 <- remove_cancer_feats(all_450, bh_feats = bh_features)
all_850 <- remove_cancer_feats(all_850, bh_feats = bh_features)

# get intersect_names
intersect_names <- names(all_450)[grepl('^cg', names(all_450))]

# get feature list
colnames(bh_feats)[1] <- 'chr'
remove_features <- inner_join(bh_feats, g_ranges)$probe

# take remove features out of colnames 
bh_features <- intersect_names[!intersect_names %in% remove_features]

# put 450 cases with 850 controls 
cases_450 <- all_450[all_450$cancer_status == 'cases',]
con_450 <- all_450[all_450$cancer_status != 'cases',]

cases_850 <- all_850[all_850$cancer_status == 'cases',]
con_850 <- all_850[all_850$cancer_status != 'cases',]

rm(all_450, all_850)

# rejoin data
all_train <- rbind(cases_450,
                   con_850)

all_test <- rbind(cases_850,
                  con_450)

rm(cases_450, con_450, cases_850, con_850)
## creat list to store results for alpha
result_list <- list()

alpha_values <- (1:10/10)

for(i in 1:length(alpha_values)){
  alpha_num <- alpha_values[i]
  
  message('working on alpha = ', alpha_num)
  result_list[[i]]  <- run_enet_surv(training_dat = all_450,
                                test_dat = all_850,
                                age_cutoff = age_cutoff,
                                gender = gender,
                                tech = tech,
                                alpha_value = alpha_num,
                                bh_features = bh_features)
  
  
}

temp <- do.call('rbind', result_list)


##########
# validation set
##########

# read in cases_450
saveRDS(temp, paste0('validation_age_predictions_surv/', 'valid_test', which_methyl, '_', 
                     which_combat, '_', is_gen, '.rda'))

