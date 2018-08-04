# get libraries 
library(tidyverse)

# read in m value 
m_values_valli <- readRDS('~/Desktop/m_value_valli.rds')[, 1:300]

# read in 850k
all_cases_m <- readRDS('../../Data/all_cases_m.rda')
all_con_m<- readRDS('../../Data/all_con_m.rda')

all_dat <- rbind(all_cases_m,
                 all_con_m)

# remove duplicated tm _don
all_dat <- all_dat[, 1:300]

# get mut 
all_dat <- all_dat[all_dat$p53_germline == 'MUT',]

# remove na
all_dat <- all_dat[!is.na(all_dat$ids),]
all_dat <- all_dat[!is.na(all_dat$tm_donor),]
all_dat <- all_dat[!is.na(all_dat$p53_germline),]
all_dat <- all_dat[!is.na(all_dat$age_sample_collection),]

# get dups
dup_ids <- all_dat$tm_donor[duplicated(all_dat$tm_donor)]
temp_dups <- all_dat[all_dat$tm_donor %in% dup_ids,]
temp_dups <- temp_dups[order(temp_dups$tm_donor),]

# loop through dups and get first onset or first age 
result_list <- list()
unique_tm_ids <- unique(all_dat$tm_donor)
i = 1
for(i in 1:length(unique_tm_ids)){
  tm_donor_id <- unique_tm_ids[i]
  sub_dups <- all_dat[all_dat$tm_donor == tm_donor_id,]
  if(all(sub_dups$cancer_diagnosis_diagnoses == 'Unaffected')){
    sub_dups <- sub_dups[order(sub_dups$age_sample_collection, decreasing = F),]
    
  } else {
    sub_dups <- sub_dups[order(sub_dups$age_diagnosis, decreasing = F),]

  }
  result_list[[i]] <- sub_dups
}
all_dat_ordered <- do.call('rbind', result_list)

all_dat_ordered <- all_dat_ordered[!duplicated(all_dat_ordered$tm_donor),]

all_ids <- all_dat$ids
v_ids<- m_values_valli$ids

in_v_not_all <- v_ids[!v_ids %in% all_ids]
in_all_not_v <- all_ids[!all_ids %in% v_ids]

in_v_not_all <- as.data.frame(in_v_not_all)
in_v_not_all$dups <- ifelse(in_v_not_all$in_v_not_all == '2447', 'not_there', 'dups')
in_all_not_v <- as.data.frame(in_all_not_v)

remove_these <- as.character(in_v_not_all$in_v_not_all)
remove_these <- remove_these[remove_these != '2407']

new_v <- m_values_valli[!m_values_valli$ids %in% remove_these,]


# check to see what type of data
unique(m_values_valli$p53_germline) # all Mut
sort(unique(m_values_valli$cancer_diagnosis_diagnoses)) # cases and controls
summary(as.factor(m_values_valli$cancer_diagnosis_diagnoses))

# # check discrepancies between m_values_vallie and all_cases_m, all_con_m - from 155 tp 154
m_values_valli <- m_values_valli[!is.na(m_values_valli$age_sample_collection),]

# get cases and con from m_values_valli
cases_v <- m_values_valli[!grepl('Unaffected', m_values_valli$cancer_diagnosis_diagnoses),]
con_v <- m_values_valli[grepl('Unaffected', m_values_valli$cancer_diagnosis_diagnoses),]
rm(m_values_valli)



# ids in new data but not in valli "4176" "3391" "3879" "3880" "2264"
# ids in vallie but not in new data "2760" "1094" "2407" "2485" "2921" "3804" "2413" "2455" "3164" "3319" "3425" "4472" "2815" "2446" "3273"

# read desktop temp data
load('~/Desktop/temp_450_rg.RData')

# # read in from the beginning and dont remove tm_donor to see how this and m_values_valli diffs
# # # read in current data
# all_cases_m <- readRDS('../../Data/all_cases_m.rda')
# all_con_m <- readRDS('../../Data/all_con_m.rda')
# # 
# # # subset by mut
# all_cases_m <- all_cases_m[all_cases_m$p53_germline == 'MUT',]
# all_con_m <- all_con_m[all_con_m$p53_germline == 'MUT',]
# # 
# # # remove NA from age of sample collection
# all_cases_m <- all_cases_m[!is.na(all_cases_m$age_sample_collection),] # from 105 to 102
# all_con_m <- all_con_m[!is.na(all_con_m$age_sample_collection),] # no change

# 

# 
# # combine two data types
# all <- rbind(all_cases_m,
#              all_con_m)
# 
# v <- rbind(cases_v,
#            con_v)
# rm(all_cases_m, all_con_m,
#    cases_v, con_v)
# 
# # remove dups
# all <- all[!duplicated(all$tm_donor),]
# 
# 
