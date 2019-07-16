# this script will regress age out of each probe

# get functions
source('all_functions.R')

# remove WT 
remove_wild_type <- function(m_or_beta_values){
  m_or_beta_values <- m_or_beta_values[m_or_beta_values$p53_germline == 'MUT',]
  return(m_or_beta_values)
}


data_type = 'm'
combat = FALSE

if(combat){
  used_combat <- 'used_combat'
} else {
  used_combat <- 'no_combat'
}

if(data_type == 'beta'){
  
  if(combat){
    all_cases <- readRDS('all_data/all_cases_beta_combat.rda')
    all_con <- readRDS('all_data/all_con_beta_combat.rda')
    all_con_wt <- readRDS('all_data/all_con_beta_wt_combat.rda')
    
    # recode to 450k and 850k
    all_cases$tech <- ifelse(all_cases$tech == 'batch_1', '450k', '850k')
    all_con$tech <- ifelse(all_con$tech == 'batch_1', '450k', '850k')
    
  } else {
    # cases
    all_cases <- readRDS('all_data/all_cases_beta.rda')
    all_con <- readRDS('all_data/all_con_beta.rda')
    all_con_wt <- readRDS('all_data/all_con_beta_wt.rda')
    
  }
  
} else {

  if(combat){
    
    all_cases  <- readRDS('all_data/all_cases_m_combat.rda')
    all_con  <- readRDS('all_data/all_con_m_combat.rda')
    all_con_wt <- readRDS('all_data/all_con_m_wt_combat.rda')
    
    all_cases$tech <- ifelse(all_cases$tech == 'batch_1', '450k', '850k')
    all_con$tech <- ifelse(all_con$tech == 'batch_1', '450k', '850k')
    
  } else {
    # cases
    all_cases <- readRDS('all_data/all_cases_m.rda')
    all_con <- readRDS('all_data/all_con_m.rda')
    all_con_wt <- readRDS('all_data/all_con_m_wt.rda')
    # subset to get controls lfs and wild type
    

  }
}

# get controls 
all_con_wt <- all_con_wt[!duplicated(all_con_wt$tm_donor),]
con_wt <- all_con_wt[all_con_wt$p53_germline == 'WT',]
con_mut <- all_con_wt[all_con_wt$p53_germline == 'MUT',]

rm(all_con_wt)

#  remove duplcates
get_data <- function(cases, controls) {
  
  # split up by tech
  cases_450 <- cases[cases$tech == '450k',]
  cases_850 <- cases[cases$tech == '850k',]
  
  cases_450 <- cases_450[!duplicated(cases_450$tm_donor),]
  cases_850 <- cases_850[!duplicated(cases_850$tm_donor),]
  cases_450 <- remove_wild_type(cases_450)
  cases_850 <- remove_wild_type(cases_850)
  cases_450 <- cases_450[!is.na(cases_450$age_diagnosis),]
  cases_450 <- cases_450[!is.na(cases_450$age_sample_collection),]
  cases_850 <- cases_850[!is.na(cases_850$age_diagnosis),]
  cases_850 <- cases_850[!is.na(cases_850$age_sample_collection),]
  
  # controls
  con_850 <- controls[controls$tech == '850k',]
  con_450 <- controls[controls$tech == '450k',]
  
  
  con_450 <- remove_wild_type(con_450)
  con_850 <- remove_wild_type(con_850)
  con_450 <- con_450[!is.na(con_450$age_sample_collection),]
  con_850 <- con_850[!is.na(con_850$age_sample_collection),]
  
  return(list(cases_450, cases_850, con_450, con_850))
}

data_list <- get_data(all_cases, all_con)

# unlist data 
cases_450 <- data_list[[1]]
cases_850 <- data_list[[2]]
con_450 <- data_list[[3]]
con_850 <- data_list[[4]]

rm(data_list)

get_risidual_data <- function(temp_dat){
  # subset by mut, and complete cases for age diagnosis and age sample collection

  cg_start <- which(grepl('^cg', names(temp_dat)))[1]
  feature_names <- colnames(temp_dat)[cg_start:ncol(temp_dat)]
  
  resid <- list()
  
  for (i in cg_start:ncol(temp_dat)) {
    
    temp_response <- temp_dat[, i]
    temp_var <- temp_dat$age_sample_collection
    
    resid[[i]] <- lm(temp_response ~ temp_var)$residuals
    
    print(i)
    
  }
  
  resid_data <- as.data.frame(do.call('cbind', resid))
  names(resid_data) <- feature_names
  final_dat <- as.data.frame(cbind(id = temp_dat$id, 
                                  p53_germline = temp_dat$p53_germline, 
                                  cancer_diagnosis_diagnoses = temp_dat$cancer_diagnosis_diagnoses, 
                                  age_diagnosis = temp_dat$age_diagnosis, 
                                  age_sample_collection = temp_dat$age_sample_collection,
                                  gender = temp_dat$gender, 
                                  sentrix_id = temp_dat$sentrix_id, 
                                  family_name = temp_dat$family_name,
                                  tm_donor = temp_dat$tm_donor,
                                  cancer_status = temp_dat$cancer_status,
                                  tech = temp_dat$tech,
                                  resid_data))

  return(final_dat)
  
  
}

# cases 450
cases_450_resid <- get_risidual_data(cases_450)
saveRDS(cases_450_resid, paste0('residual_data/cases_450_resid_', 
                                data_type, 
                                '_', 
                                used_combat, '.rda'))
# cases 850
cases_850_resid <- get_risidual_data(cases_850)
saveRDS(cases_850_resid, paste0('residual_data/cases_850_resid_', 
                                data_type, 
                                '_', 
                                used_combat, '.rda'))

# con 450
con_450_resid <- get_risidual_data(con_450)
saveRDS(con_450_resid, paste0('residual_data/con_450_resid_', 
                                data_type, 
                                '_', 
                                used_combat, '.rda'))

# con 450
con_850_resid <- get_risidual_data(con_850)
saveRDS(con_850_resid, paste0('residual_data/con_850_resid_', 
                              data_type, 
                              '_', 
                              used_combat, '.rda'))

# cases 450
con_wt <- get_risidual_data(con_wt)
saveRDS(con_wt, paste0('residual_data/con_wt_resid_', 
                                data_type, 
                                '_', 
                                used_combat, '.rda'))

# cases 450
con_mut <- get_risidual_data(con_mut)
saveRDS(con_mut, paste0('residual_data/con_mut_resid_', 
                       data_type, 
                       '_', 
                       used_combat, '.rda'))


