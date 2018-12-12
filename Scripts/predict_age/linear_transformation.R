# this script will map beta and m values from 850 to 450 and 450 to 850

# get functions
source('all_functions.R')

data_type = 'm'

if(data_type == 'beta'){
  # read in data
  all_cases_beta <- readRDS('../../Data/all_cases_beta_transform.rda')
  all_con_beta <- readRDS('../../Data/all_con_beta_transform.rda')
  
} else {
  all_cases_m <- readRDS('../../Data/all_cases_m_transform.rda')
  all_con_m <- readRDS('../../Data/all_con_m_transform.rda')
}


# creat function for beta and m

transform_beta_m_values <- function(all_cases, all_controls){
  
  # remove overlapping ids to get a transfromed set for model
  all_cases <- all_cases[order(all_cases$tech),]
  
  # split up by tech
  cases_450 <- all_cases[all_cases$tech == '450k',]
  cases_850 <- all_cases[all_cases$tech == '850k',]
  
  cases_450 <- cases_450[!duplicated(cases_450$tm_donor),]
  cases_850 <- cases_850[!duplicated(cases_850$tm_donor),]
  
  cases_450 <- remove_wild_type(cases_450)
  cases_850 <- remove_wild_type(cases_850)
  
  cases_450 <- cases_450[!is.na(cases_450$age_diagnosis),]
  cases_450 <- cases_450[!is.na(cases_450$age_sample_collection),]
  cases_850 <- cases_850[!is.na(cases_850$age_diagnosis),]
  cases_850 <- cases_850[!is.na(cases_850$age_sample_collection),]
  
  # get overlapping ids for cases across tech
  tm_cases_450 <- unique(cases_450$tm_donor)
  tm_cases_850 <- unique(cases_850$tm_donor)
  
  shared_cases_tm <- intersect(tm_cases_450, tm_cases_850)
  shared_cases_tm <- paste(shared_cases_tm, collapse = '|')
  
  shared_cases_450 <- cases_450[grepl(shared_cases_tm, cases_450$tm_donor),]
  shared_cases_850 <- cases_850[grepl(shared_cases_tm, cases_850$tm_donor),]
  
  cases_valid <- cases_850[!grepl(shared_cases_tm, cases_850$tm_donor),]
  
  rm(all_cases, cases_850)
 
  con_850 <- all_controls[all_controls$tech == '850k',]
  
  # # get controls
  # con_850 <- con_850[grepl('Unaffected', con_850$cancer_diagnosis_diagnoses),]
  # 
  # con_450 <- remove_wild_type(con_450)
  con_850 <- remove_wild_type(con_850)
  
  # con_450 <- con_450[!is.na(con_450$age_sample_collection),]
  con_850 <- con_850[!is.na(con_850$age_sample_collection),]
  
  # get overlapping ids for controls across tech
  # tm_con_450 <- unique(con_450$tm_donor)
  tm_con_850 <- unique(con_850$tm_donor)
  
  shared_450 <- intersect(tm_cases_450, tm_con_850)
  shared_con_tm <- paste(shared_450, collapse = '|')
  
  # shared_con_450 <- con_450[grepl(shared_con_tm, con_450$tm_donor),]
  shared_con_850 <- con_850[grepl(shared_con_tm, con_850$tm_donor),]
  shared_cases_con_450 <- cases_450[grepl(shared_con_tm, cases_450$tm_donor),]
  
  # get con valid 
  con_valid <- con_850[!grepl(shared_con_tm, con_850$tm_donor),]
  
  rm(all_con, cases_450, shared_cases_tm, shared_con_tm, con_850)
  rm(tm_cases_450, tm_cases_850, shared_450, tm_con_850)
  
  transform_valid <- linearTransform(shared_cases_450, 
                                     shared_cases_850, 
                                     cases_valid)
  transform_controls <- linearTransform(shared_cases_con_450, 
                                        shared_con_850, con_valid)
  
  return(list(transform_valid, transform_controls))
  
  
}

temp_transform <- transform_beta_m_values(all_cases_m, all_con_m)
transform_valid <- temp_transform[[1]]
transform_controls <- temp_transform[[2]]


if(data_type == 'beta'){
  # save both data sets 
  saveRDS(transform_valid, '../../Data/model_data/valid_transform_beta.rda')
  saveRDS(transform_controls, '../../Data/model_data/controls_transform_beta.rda')
  
} else {
  saveRDS(transform_valid, '../../Data/model_data/valid_transform_m.rda')
  saveRDS(transform_controls, '../../Data/model_data/controls_transform_m.rda')
  
}



