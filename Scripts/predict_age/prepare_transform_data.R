# this scrip will predict on data that was linear transformed and age regressed out.

# source functions script
source('all_functions.R')

# create fixed objects to model and pipeline inputs and saving  
method = 'quan'

# read both data sets 
all_cases <- readRDS(paste0('../../Data/', method,'/all_cases_beta.rda'))
all_con_wt <- readRDS(paste0('../../Data/', method,'/all_con_beta_wt.rda'))
valid_transform <- readRDS(paste0('../../Data/', method,'/valid_transform_beta.rda'))
con_transform <- readRDS(paste0('../../Data/', method,'/controls_transform_beta.rda'))

#CASES
cases_450 <- all_cases[all_cases$tech == '450k',]
cases_450 <- cases_450[!is.na(cases_450$age_sample_collection),]
cases_450 <- cases_450[!is.na(cases_450$age_diagnosis),]
cases_450 <- cases_450[!is.na(cases_450$gender),]

rm(all_cases)

# WT
all_con_wt <- all_con_wt[!duplicated(all_con_wt$tm_donor),]
con_wt <- all_con_wt[all_con_wt$p53_germline == 'WT',]
con_mut <- all_con_wt[all_con_wt$p53_germline == 'MUT',]

rm(all_con_wt)


# condition on fixed objects to get saving identifiers
# if(data_type == 'beta'){
  which_methyl <- 'beta'
  beta_thresh = 0.05
#   
# } else {
#   which_methyl <- 'm'
#   beta_thresh= 0.5
# }

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
con_transform <- cbind(as.data.frame(class.ind(con_transform$gender)), 
                       con_transform)

# rempove old tech variable 
con_transform$gender <- NULL

# gender
valid_transform <- cbind(as.data.frame(class.ind(valid_transform$gender)), 
                         valid_transform)

# rempove old tech variable 
valid_transform$gender <- NULL


# convert age to months
cases_450$age_diagnosis <- 
  round(cases_450$age_diagnosis*12, 2)
cases_450$age_sample_collection <- 
  round(cases_450$age_sample_collection*12, 2)

# controls
con_transform$age_diagnosis <- 
  round(con_transform$age_diagnosis*12, 2)
con_transform$age_sample_collection <- 
  round(con_transform$age_sample_collection*12, 2)

# valie
valid_transform$age_diagnosis <- 
  round(valid_transform$age_diagnosis*12, 2)
valid_transform$age_sample_collection <- 
  round(valid_transform$age_sample_collection*12, 2)

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

names(con_transform)[3] <- 'ids'
names(valid_transform)[3] <- 'ids'

# remove age from literature
clin_names <- names(cases_450)[!grepl('^cg', names(cases_450))]
clin_names <- clin_names[-c(11,12)]
feats <- names(cases_450)[grepl('^cg', names(cases_450))]
feats <- feats[!feats %in% age_probes]
cases_450 <- cases_450[, c(clin_names, feats)]
con_transform <- con_transform[, c(clin_names, feats)]
valid_transform <- valid_transform[, c(clin_names, feats)]
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
con_transform_small <- join_new_features(con_transform, new_features = bh_feats)
valid_transform_small <- join_new_features(valid_transform, new_features = bh_feats)
con_mut_small <- join_new_features(con_mut, new_features = bh_feats)
con_wt_small <- join_new_features(con_wt, new_features = bh_feats)

# lfs probes 
lfs_bump_probes <- colnames(cases_450)[grepl('^cg', colnames(cases_450))]

rm(bh_feats)


# add dummy tech variable for data sets with only one, replace family_name
names(cases_450)[9] <- 'tech'
names(con_transform)[9] <- 'tech'
names(valid_transform)[9] <- 'tech'

# fill them with Zero
cases_450$tech <- '450k'
con_transform$tech <- '850k'
valid_transform$tech <- '850k'

# do the same to con_mut and con_wt
names(con_mut)[9] <- 'tech'
names(con_wt)[9] <- 'tech'

# fill new variable with right tech indication
con_mut$tech <- '450k'
con_wt$tech <- '450k'


# save trainig  set 
saveRDS(cases_450_small, paste0('../../Data/', method,'/cases_450_small.rda'))

# save validation set 
saveRDS(con_transform_small, paste0('../../Data/', method,'/con_transform_small.rda'))

# save validation set 
saveRDS(valid_transform_small, paste0('../../Data/', method,'/valid_transform_small.rda'))

# save wt controls 
saveRDS(con_wt_small, paste0('../../Data/', method,'/con_wt_small.rda'))

# save mut controls 
saveRDS(con_mut_small, paste0('../../Data/', method,'/con_mut_small.rda'))



# save trainig  set 
saveRDS(cases_450, paste0('../../Data/', method,'/cases_450.rda'))

# save validation set 
saveRDS(con_transform, paste0('../../Data/', method,'/con_transform.rda'))

# save validation set 
saveRDS(valid_transform, paste0('../../Data/', method,'/valid_transform.rda'))

# save wt controls 
saveRDS(con_wt, paste0('../../Data/', method,'/con_wt.rda'))

# save mut controls 
saveRDS(con_mut, paste0('../../Data/', method,'/con_mut.rda'))

