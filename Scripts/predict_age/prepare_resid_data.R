# this scrip will predict on data that was linear transformed and age regressed out.

# source functions script
source('all_functions.R')

# create fixed objects to model and pipeline inputs and saving  
data_type = 'm'
combat = FALSE

if(combat){
  used_combat <- 'used_combat'
} else {
  used_combat <- 'no_combat'
}


# read in all data
cases_450 <- readRDS(paste0('residual_data/', 'cases_450_resid_',data_type, '_',used_combat,'.rda'))
con_850 <- readRDS(paste0('residual_data/', 'con_850_resid_',data_type,'_',used_combat, '.rda'))
cases_850 <- readRDS(paste0('residual_data/', 'cases_850_resid_',data_type,'_',used_combat,'.rda'))
con_wt <-  readRDS(paste0('residual_data/', 'con_wt_resid_',data_type,'_',used_combat, '.rda'))
con_mut <-  readRDS(paste0('residual_data/', 'con_mut_resid_',data_type,'_',used_combat, '.rda'))

lfs_bump_probes <- readRDS(paste0('transform_data_cv/', 'lfs_bumps_', data_type,'_.rda'))

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

# contorls wt and mut
con_wt$age_sample_collection <- 
  round(con_wt$age_sample_collection*12, 2)

con_mut$age_sample_collection <- 
  round(con_mut$age_sample_collection*12, 2)


# remove age from literature
clin_names <- names(cases_450)[!grepl('^cg', names(cases_450))]
feats <- names(cases_450)[grepl('^cg', names(cases_450))]
feats <- feats[!feats %in% age_probes]
cases_450 <- cases_450[, c(clin_names, feats)]
con_850 <- con_850[, c(clin_names, feats)]
cases_850 <- cases_850[, c(clin_names, feats)]
con_wt <- con_wt[, c(clin_names, feats)]
con_mut <- con_mut[, c(clin_names, feats)]


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

# save trainig  set 
saveRDS(cases_450, paste0('residual_data_cv/', 'cases_450_',data_type, '_', used_combat,'.rda'))

# save validation set 
saveRDS(con_850, paste0('residual_data_cv/', 'con_850_',data_type,'_', used_combat, '.rda'))

# save validation set 
saveRDS(cases_850, paste0('residual_data_cv/', 'cases_850_',data_type,'_', used_combat,'.rda'))

# save wt controls 
saveRDS(con_wt, paste0('residual_data_cv/', 'con_wt_',data_type, '_', used_combat,'.rda'))

# save mut controls 
saveRDS(con_mut, paste0('residual_data_cv/', 'con_mut_',data_type,'_', used_combat,'.rda'))



