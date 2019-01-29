source('all_functions.R')

# remove WT 
remove_wild_type <- function(m_or_beta_values){
  m_or_beta_values <- m_or_beta_values[m_or_beta_values$p53_germline == 'MUT',]
  return(m_or_beta_values)
}


# start with over and go through every iteration while control_age =
# set fixed variables
# set fixed variables
model_type = 'surv'
size = 'used_bh'
fit_type = 'rf'
outcome_type = '6_years'
control_age = FALSE # need to rerun control age = FALSE and then redo with coxph
age_combat = TRUE
gender = TRUE
method = 'quan'
combat = 'combat_1'
which_methyl = 'beta'
beta_thresh = 0.01
optimal_cutoff = 0.5

# create objects to indicate method and model details when saving
age_cutoff = 72
trained_lambda = FALSE
tech = FALSE


if(trained_lambda){
  is_lambda <- 'lambda_test'
} else {
  is_lambda <- 'lambda_train'
}

if(control_age){
  age_control <- 'control_age'
} else {
  age_control <- 'no_control_age'
  
}

if(age_combat){
  combat_age <- 'combat_age'
} else {
  combat_age <- 'pc_age'
}

if(gender){
  is_gen = 'use_gen'
} else {
  is_gen = 'no_gen'
}

if(age_combat){
  
  if(size =='used_bh'){
    cases_450 <-  readRDS(paste0('../../Data/', method,'/cases_450_small_cv_age_combat_', combat,'.rda'))
    cases_850 <- readRDS(paste0('../../Data/', method,'/cases_850_small_cv_age_combat_', combat,'.rda'))
    con_mut <- readRDS(paste0('../../Data/', method,'/con_mut_small_cv_age_combat_', combat,'.rda'))
    con_850 <- readRDS( paste0('../../Data/', method,'/con_850_small_cv_age_combat_', combat,'.rda'))
    con_wt <- readRDS(paste0('../../Data/', method,'/con_wt_small_cv_age_combat_', combat,'.rda'))
    
  } else {
    cases_450 <-  readRDS(paste0('../../Data/', method,'/cases_450_cv_age_combat_', combat,'.rda'))
    cases_850 <- readRDS(paste0('../../Data/', method,'/cases_850_cv_age_combat_', combat,'.rda'))
    con_mut <- readRDS(paste0('../../Data/', method,'/con_mut_cv_age_combat_', combat,'.rda'))
    con_850 <- readRDS( paste0('../../Data/', method,'/con_850_cv_age_combat_', combat,'.rda'))
    con_wt <- readRDS(paste0('../../Data/', method,'/con_wt_cv_age_combat_', combat,'.rda'))
    
   
  }
  
} else {
  if(size =='used_bh'){
    cases_450 <-  readRDS(paste0('../../Data/', method,'/cases_450_small_cv', combat,'.rda'))
    cases_850 <- readRDS(paste0('../../Data/', method,'/cases_850_small_cv', combat,'.rda'))
    con_mut <- readRDS(paste0('../../Data/', method,'/con_mut_small_cv', combat,'.rda'))
    con_850 <- readRDS( paste0('../../Data/', method,'/con_850_small_cv', combat,'.rda'))
    con_wt <- readRDS(paste0('../../Data/', method,'/con_wt_small_cv', combat,'.rda'))
    
  } else {
    cases_450 <-  readRDS(paste0('../../Data/', method,'/cases_450_cv', combat,'.rda'))
    cases_850 <- readRDS(paste0('../../Data/', method,'/cases_850_cv', combat,'.rda'))
    con_mut <- readRDS(paste0('../../Data/', method,'/con_mut_cv', combat,'.rda'))
    con_850 <- readRDS( paste0('../../Data/', method,'/con_850_cv', combat,'.rda'))
    con_wt <- readRDS(paste0('../../Data/', method,'/con_wt_cv', combat,'.rda'))
    
   
    # 
  }
  
}




# # load lambda
# model_params <- readRDS(paste0('data_cv/results/model_params_', data_type, '_', base_num, '_',
#                              num_seeds, '_', k_folds, '_', is_gen, '_', model_type, '.rda'))
# 
# mean_lambda <- model_params[1]
####
##########
# load genomic methyl set (from controls) - you need genetic locations by probe from this object
##########
g_ranges <- readRDS('../../Data/g_ranges.rda')

# get probes from rownames
g_ranges$probe <- rownames(g_ranges)

# remove ch and duplicatee
g_ranges <- g_ranges[!duplicated(g_ranges$start),]
g_ranges <- g_ranges[!grepl('ch', g_ranges$probe),]

names(g_ranges)[1] <- 'chr'

# remove cancer signature
bh_feats <- bump_hunter(dat_1 = cases_450, 
                        dat_2 = con_mut , 
                        bump = 'cancer', 
                        boot_num = 5, 
                        beta_thresh = beta_thresh,
                        methyl_type = methyl_type,
                        g_ranges = g_ranges)


# get intersect_names
intersect_names <- names(cases_450)[grepl('^cg', names(cases_450))]

# get feature list
colnames(bh_feats)[1] <- 'chr'
remove_features <- inner_join(bh_feats, g_ranges)$probe

# take remove features out of colnames 
bh_features <- intersect_names[!intersect_names %in% remove_features]

# subset all data by bh_features
cases_450 <- remove_cancer_feats(cases_450, bh_feats = bh_features)
con_850 <- remove_cancer_feats(con_850, bh_feats = bh_features)
con_mut <- remove_cancer_feats(con_mut, bh_feats = bh_features)
cases_850 <- remove_cancer_feats(cases_850, bh_feats = bh_features)

# get age factor to control for age.
get_age_fac <- function(temp_dat){
  temp_dat$age_fac <- ifelse(temp_dat$age_sample_collection >= 140, 'age_1', 'age_2')
  return(temp_dat)
}

# combine data
all_450 <- rbind(cases_450,
                 con_mut)
all_850 <- rbind(cases_850,
                 con_850)

all_450 <- get_age_fac(all_450)
all_850 <- get_age_fac(all_850)

rm(cases_450, cases_850, con_mut, con_wt, con_850)
# age
all_450 <- cbind(as.data.frame(class.ind(all_450$age_fac)), 
                   all_450)

# rempove old agevariable 
all_450$age_fac <- NULL

# age
all_850 <- cbind(as.data.frame(class.ind(all_850$age_fac)), 
                 all_850)

# rempove old agevariable 
all_850$age_fac <- NULL


## creat list to store results for alpha

result_list <- list()
model_list <- list()
alpha_values <- (1:10)/10

for(i in 1:length(alpha_values)){
  alpha_num <- alpha_values[i]
  
  message('working on alpha = ', alpha_num)
  temp_results   <- run_enet_surv(training_dat = all_450,
                                     test_dat = all_850,
                                     age_cutoff = age_cutoff,
                                     fit_type = fit_type,
                                     gender = gender,
                                     tech = tech,
                                     control_age = control_age,
                                     alpha_value = alpha_num,
                                     bh_features = bh_features)
  
  result_list[[i]] <- temp_results[[1]]
  model_list[[i]] <- temp_results[[2]]
  
  
}

temp <- do.call('rbind', result_list)


##########
# validation set
##########
# set fixed variables
# set fixed variables
model_type = 'surv'
size = 'used_bh'
fit_type = 'rf'
control_age = TRUE
age_combat = TRUE
gender = TRUE
method = 'funnorm'
combat = 'combat_sen'
if(control_age){
  age_control <- 'control_age'
} else {
  age_control <- 'no_control_age'
  
}

if(age_combat){
  combat_age <- 'combat_age'
} else {
  combat_age <- 'pc_age'
}

if(gender){
  is_gen = 'use_gen'
} else {
  is_gen = 'no_gen'
}


# options(scipen = 000)

# read in cases_450
temp <- readRDS(paste0('surv_data_test/', 'surv_test_', method,'_',size,'_',is_gen, '_',combat,'_', model_type, '_', age_control,'_',combat_age,'_',fit_type,'.rda'))

# ggplot(temp, aes( alpha,test_pred, fill = cancer_status)) +
#   geom_bar(stat = 'identity', position = 'dodge')
# 
best_alpha = 0.5
temp  = temp[temp$alpha == best_alpha,]

temp$age <- ifelse(temp$cancer_status == 'cases', temp$age_diagnosis, temp$age_sample_collection)
ggplot(temp, aes(age,test_pred, fill = cancer_status)) +
  geom_bar(stat = 'identity', position = 'dodge', width =1) + 
  geom_vline(xintercept = 72)

ggplot(temp, aes(age,test_pred, color = cancer_status)) +
  geom_point(size = 4, alpha = 0.4) + 
  geom_vline(xintercept = 72)

library(ggrepel)
ggplot(temp, aes(age,test_pred, color = cancer_status)) +
  geom_point(size = 2, alpha = 0.4) + 
  geom_text_repel(aes(label = cancer_diagnosis_diagnoses), cex = 3)+
  geom_vline(xintercept = 72)

ggplot(temp, aes(age,test_pred, color = cancer_status)) +
  # geom_point(size = 4, alpha = 0.4) + 
  geom_text(aes(label = cancer_diagnosis_diagnoses))+
  geom_vline(xintercept = 72)
