
source('helper_functions.R')

# create fixed objects to model and pipeline inputs and saving  
methyl_type = 'beta'
combat = TRUE
gender = TRUE
age_cutoff = 72
tech = FALSE 
how_many_seeds = 3
how_many_folds = 5
use_offset = TRUE
remove_age = TRUE

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


num_seeds <- paste0('seeds_', how_many_seeds)
num_folds <- paste0('folds_', how_many_folds)
k_folds <- how_many_folds

# read in data
final_dat <- readRDS(paste0('validation_age_predictions_surv/', which_methyl, '_', which_combat, '_',
                          num_seeds, '_', num_folds, '_', is_gen, '_tech_correct.rda'))

# get avg over seeds 
sub_dat <- final_dat %>% 
  group_by(tm_donor) %>% 
  mutate(mean_risk = mean(test_pred, na.rm = TRUE)) 

sub_dat <- sub_dat[!duplicated(sub_dat$tm_donor),]
sub_dat$cancer_status <- ifelse(sub_dat$cancer_status == 'cases', 'cancer', 'no_cancer')

# keep only necessary columns
names(sub_dat)
sub_dat <- sub_dat[, c('tm_donor','Female', 'Male', 'cancer_diagnosis_diagnoses', 'age_diagnosis', 
                       'age_sample_collection','cancer_status' ,'mean_risk')]

# bar plot colored by cancer status
ggplot(sub_dat, aes(reorder(tm_donor, -mean_risk), mean_risk, fill = cancer_status)) +
  geom_bar(stat = 'identity', position = 'dodge', alpha = 0.6, width = 1) +
  scale_fill_manual(name = '',
                    values = c('red', 'blue')) +
  labs(x = 'Patients (n = 89)',
       y = 'Mean relative risk score') +
  theme_minimal() +
  theme(axis.text.x = element_text(size =0))
