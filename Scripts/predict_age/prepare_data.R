source('all_functions.R')

# remove WT 
remove_wild_type <- function(m_or_beta_values){
  m_or_beta_values <- m_or_beta_values[m_or_beta_values$p53_germline == 'MUT',]
  return(m_or_beta_values)
}


# set fixed variables
method = 'swan'
combat = 'combat_sen'
beta_thresh = 0.05
methyl_type = 'beta'

# save beta data
# cases

if(combat == 'normal'){
  all_cases <- readRDS(paste0('../../Data/', method,'/all_cases_beta.rda'))
  all_con <- readRDS(paste0('../../Data/', method,'/all_con_beta.rda'))
  all_con_wt <- readRDS(paste0('../../Data/', method,'/all_con_beta_wt.rda'))
  
}

if(combat == 'combat_1'){
  all_cases <- readRDS(paste0('../../Data/', method,'/all_cases_beta_combat.rda'))
  all_con <- readRDS(paste0('../../Data/', method,'/all_con_beta_combat.rda'))
  all_con_wt <- readRDS(paste0('../../Data/', method,'/all_con_beta_wt.rda'))
  
  all_cases$tech <- ifelse(all_cases$tech == 'batch_1', '450k', '850k')
  all_con$tech <- ifelse(all_con$tech == 'batch_1', '450k', '850k')
  all_con_wt$tech <- ifelse(all_con_wt$tech == 'batch_1', '450k', '850k')
}

if(combat == 'combat_gen'){
  all_cases <- readRDS(paste0('../../Data/', method,'/all_cases_beta_combat_gen.rda'))
  all_con <- readRDS(paste0('../../Data/', method,'/all_con_beta_combat_gen.rda'))
  all_con_wt <- readRDS(paste0('../../Data/', method,'/all_con_beta_wt_combat_gen.rda'))
  
  all_cases$tech <- ifelse(all_cases$tech == 'batch_1', '450k', '850k')
  all_con$tech <- ifelse(all_con$tech == 'batch_1', '450k', '850k')
  all_con_wt$tech <- ifelse(all_con_wt$tech == 'batch_1', '450k', '850k')
}

if(combat == 'combat_sen'){
  all_cases <- readRDS(paste0('../../Data/', method,'/all_cases_beta_combat_sen.rda'))
  all_con <- readRDS(paste0('../../Data/', method,'/all_con_beta_combat_sen.rda'))
  all_con_wt <- readRDS(paste0('../../Data/', method,'/all_con_beta_wt_combat_sen.rda'))
  
  all_cases$tech <- ifelse(all_cases$tech == 'batch_1', '450k', '850k')
  all_con$tech <- ifelse(all_con$tech == 'batch_1', '450k', '850k')
  all_con_wt$tech <- ifelse(all_con_wt$tech == 'batch_1', '450k', '850k')
}


# get controls 
all_con_wt <- all_con_wt[!duplicated(all_con_wt$tm_donor),]
con_wt <- all_con_wt[all_con_wt$p53_germline == 'WT',]
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
con_mut <- data_list[[3]]
con_850 <- data_list[[4]]

rm(data_list)

##########
# read in age probes
##########
age_probes <- readRDS('../../Data/age_probes.rda')

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

##########
# create variables
##########

# load cases
cases_450 <- cbind(as.data.frame(class.ind(cases_450$gender)), 
                   cases_450)

# rempove old tech variable 
cases_450$gender <- NULL

# gender
con_850 <- cbind(as.data.frame(class.ind(con_850$gender)), 
                 con_850)

# rempove old tech variable 
con_850$gender <- NULL

# gender
cases_850 <- cbind(as.data.frame(class.ind(cases_850$gender)), 
                   cases_850)

# rempove old tech variable 
cases_850$gender <- NULL


# convert age to months
cases_450$age_diagnosis <- 
  round(cases_450$age_diagnosis*12, 2)
cases_450$age_sample_collection <- 
  round(cases_450$age_sample_collection*12, 2)

# controls
con_850$age_diagnosis <- 
  round(con_850$age_diagnosis*12, 2)
con_850$age_sample_collection <- 
  round(con_850$age_sample_collection*12, 2)

# valie
cases_850$age_diagnosis <- 
  round(cases_850$age_diagnosis*12, 2)
cases_850$age_sample_collection <- 
  round(cases_850$age_sample_collection*12, 2)

# apply to wt controls
con_wt <- con_wt[!is.na(con_wt$age_sample_collection),]
con_wt <- con_wt[!is.na(con_wt$gender),]
con_mut <- con_mut[!is.na(con_mut$age_sample_collection),]
con_mut <- con_mut[!is.na(con_mut$gender),]

# ge tgender 
con_wt <- cbind(as.data.frame(class.ind(con_wt$gender)), 
                con_wt)
con_mut <- cbind(as.data.frame(class.ind(con_mut$gender)), 
                 con_mut)

# rempove old tech variable 
con_wt$gender <- NULL
con_mut$gender <- NULL

# subset to get controls lfs and wild type
names(con_wt)[3] <- 'ids'
names(con_mut)[3] <- 'ids'

names(con_850)[3] <- 'ids'
names(cases_850)[3] <- 'ids'

# remove age from literature
clin_names <- names(cases_450)[!grepl('^cg', names(cases_450))]
feats <- names(cases_450)[grepl('^cg', names(cases_450))]
feats <- feats[!feats %in% age_probes]
cases_450 <- cases_450[, c(clin_names, feats)]
con_850 <- con_850[, c(clin_names, feats)]
cases_850 <- cases_850[, c(clin_names, feats)]
con_wt <- con_wt[, c(clin_names, feats)]
con_mut <- con_mut[, c(clin_names, feats)]


# run bumphunter on LFS healthy patients (LFS no cancer) and LFS cancer patients (LFS cancer)
bh_feats <- bump_hunter(dat_1 = con_wt, 
                        dat_2 = con_mut, 
                        bump = 'lfs', 
                        boot_num = 5, 
                        beta_thresh = beta_thresh,
                        methyl_type = methyl_type,
                        g_ranges = g_ranges)

# cases
cases_450_small <- join_new_features(cases_450, new_features = bh_feats)
con_850_small <- join_new_features(con_850, new_features = bh_feats)
cases_850_small <- join_new_features(cases_850, new_features = bh_feats)
con_mut_small <- join_new_features(con_mut, new_features = bh_feats)
con_wt_small <- join_new_features(con_wt, new_features = bh_feats)

# lfs probes 
lfs_bump_probes <- colnames(cases_450)[grepl('^cg', colnames(cases_450))]

rm(bh_feats)

# add dummy tech variable for data sets with only one, replace family_name
names(cases_450)[9] <- 'tech'
names(con_850)[9] <- 'tech'
names(cases_850)[9] <- 'tech'

# fill them with Zero
cases_450$tech <- '450k'
con_850$tech <- '850k'
cases_850$tech <- '850k'

# do the same to con_mut and con_wt
names(con_mut)[9] <- 'tech'
names(con_wt)[9] <- 'tech'

# fill new variable with right tech indication
con_mut$tech <- '450k'
con_wt$tech <- '450k'

saveRDS(cases_450_small, paste0('../../Data/', method,'/cases_450_small_norm_', combat,'.rda'))
saveRDS(cases_850_small, paste0('../../Data/', method,'/cases_850_small_norm_', combat,'.rda'))
saveRDS(con_mut_small, paste0('../../Data/', method,'/con_mut_small_norm_', combat,'.rda'))
saveRDS(con_850_small, paste0('../../Data/', method,'/con_850_small_norm_', combat,'.rda'))
saveRDS(con_wt_small, paste0('../../Data/', method,'/con_wt_small_norm_', combat,'.rda'))

saveRDS(cases_450, paste0('../../Data/', method,'/cases_450_norm_', combat,'.rda'))
saveRDS(cases_850, paste0('../../Data/', method,'/cases_850_norm_', combat,'.rda'))
saveRDS(con_mut, paste0('../../Data/', method,'/con_mut_norm_', combat,'.rda'))
saveRDS(con_850, paste0('../../Data/', method,'/con_850_norm_', combat,'.rda'))
saveRDS(con_wt, paste0('../../Data/', method,'/con_wt_norm_', combat,'.rda'))

