source('all_functions.R')

# remove WT 
remove_wild_type <- function(m_or_beta_values){
  m_or_beta_values <- m_or_beta_values[m_or_beta_values$p53_germline == 'MUT',]
  return(m_or_beta_values)
}


# set fixed variables
method = 'noob'
combat = 'combat_1_new'
remove_leading_pcs = 'first'

# condition on fixed objects to get saving identifiers
which_methyl = 'beta'
beta_thresh = 0.05

cases_450 <- readRDS(paste0('../../Data/', method,'/cases_450_', combat,'.rda'))
cases_850 <- readRDS(paste0('../../Data/', method,'/cases_850_', combat,'.rda'))
con_850 <- readRDS(paste0('../../Data/', method,'/con_850_', combat,'.rda'))
con_mut <- readRDS(paste0('../../Data/', method,'/con_450_', combat,'.rda'))
con_wt <- readRDS(paste0('../../Data/', method,'/con_wt_', combat,'.rda'))

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
                        boot_num = 50, 
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

saveRDS(cases_450_small, paste0('../../Data/', method,'/cases_450_small_cv', combat,'.rda'))
saveRDS(cases_850_small, paste0('../../Data/', method,'/cases_850_small_cv', combat,'.rda'))
saveRDS(con_mut_small, paste0('../../Data/', method,'/con_mut_small_cv', combat,'.rda'))
saveRDS(con_850_small, paste0('../../Data/', method,'/con_850_small_cv', combat,'.rda'))
saveRDS(con_wt_small, paste0('../../Data/', method,'/con_wt_small_cv', combat,'.rda'))

saveRDS(cases_450, paste0('../../Data/', method,'/cases_450_cv', combat,'.rda'))
saveRDS(cases_850, paste0('../../Data/', method,'/cases_850_cv', combat,'.rda'))
saveRDS(con_mut, paste0('../../Data/', method,'/con_mut_cv', combat,'.rda'))
saveRDS(con_850, paste0('../../Data/', method,'/con_850_cv', combat,'.rda'))
saveRDS(con_wt, paste0('../../Data/', method,'/con_wt_cv', combat,'.rda'))

