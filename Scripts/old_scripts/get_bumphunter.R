##########
# initialize libraries
##########
library(caret)
library(glmnet)
library(randomForest)
library(kernlab)
library(pROC)
library(Metrics)
library(ModelMetrics)
library(doParallel)
library(nnet)
library(dplyr)
library(bumphunter)
library(sqldf)
library(e1071)
library(reshape2)

registerDoParallel(1)

##########
# initialize folders
##########

home_folder <- '~/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
feat_data <- paste0(data_folder, '/feat_data')
model_data <- paste0(data_folder, '/model_data')

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
method = 'noob'
k = 5
combat = F

##########
# load data
##########


if(combat) {
  
  betaFull <-  readRDS(paste0(model_data, paste0('/', method, '_', 'mod_data_combat_m.rda')))
  
  
} else {
  
  betaFull <- readRDS(paste0(model_data, paste0('/', method, '_', 'mod_data_m.rda')))
  
}



# # read features
# island_probes <- readRDS(paste0(feat_data, paste0('/', 'raw', '_', 'island_probes.rda')))
# chr_17_int <- readRDS(paste0(feat_data, paste0('/', method, '_', 'chr_17_int.rda')))
# chr_17_island_int <- readRDS(paste0(feat_data, paste0('/', method, '_', 'chr_17_island_int.rda')))

##########
# load genomic methyl set (from controls) - you need genetic locations by probe from this object
##########
ratio_set <- readRDS(paste0(model_data, paste0('/', 'raw', '_', 'ratio_set.rda')))

# get granges object
g_ranges <- as.data.frame(getLocations(ratio_set))

# get probes from rownames
g_ranges$probe <- rownames(g_ranges)

# remove ch and duplicatee
g_ranges <- g_ranges[!duplicated(g_ranges$start),]
g_ranges <- g_ranges[!grepl('ch', g_ranges$probe),]

# get column names
intersect_names <- colnames(betaFull)[9:ncol(betaFull)]

##########
# function that runs bumphunter on differe pop
##########
# dat <- betaFull
# # subset_feats <- island_probes
# probe_start = 9


bump_hunter_master <-
  function (dat, 
            subset_feats, 
            m_thresh, 
            probe_start, 
            which_bumps) {
  
  if(!is.null(subset_feats)) {
    dat_clin <- dat[, 1:(probe_start -1)]
    int_names <- intersect(intersect_names, subset_feats)
    dat <-cbind(dat_clin, dat[, int_names])
  } else {
    subset_feats <- intersect_names
  }
  
  
  if(which_bumps == 'age') {
    # get all dat
    dat <- dat[grepl('Unaffected', dat$cancer_diagnosis_diagnoses),]
    
    # remove duplicates 
    dat <- dat[!duplicated(dat$ids),]
    
    # split age at 18 (216 months)
    dat_under <- dat[dat$age_sample_collection < 216,]
    dat_over <- dat[dat$age_sample_collection >= 216,]
    
    # remove NAs in ids 
    dat_over <- dat_over[!is.na(dat_over$ids),]
    dat_under <- dat_under[!is.na(dat_under$ids),]
    
    # change type 
    dat_over$type <- 'over_18'
    dat_under$type <- 'under_18'
    
    # Run bumphunter on kids and adults
    temp <- bump_hunter(dat_1 = dat_over, 
                        dat_2 = dat_under, 
                        boot_num = 3, 
                        m_beta_thresh = m_thresh)
    
    age_probes <-inner_join(temp, g_ranges, by = 'start' )
    
    # get significant 
    age_sig <- age_probes$probe[age_probes$p.value < 0.10]
    
    # remove these from the subsetfeats 
    temp_age <- subset_feats[!subset_feats %in% age_sig]

    return(temp_age)
    
  }
   
    if(which_bumps == 'cancer') {
      
      # get caes and controls
      cases <- dat[!grepl('Unaffected',dat$cancer_diagnosis_diagnoses),]
      controls <- dat[grepl('Unaffected',dat$cancer_diagnosis_diagnoses),]
      
      controls <- subset(controls, p53_germline == 'Mut')
      
      subset(beta_cases, 
               p53_germline == 'Mut')
      
      controls_sub <- controls[!duplicated(controls$ids),]
      length(which(duplicated(controls$ids)))
      
      cases$type <- 'cancer'
      controls_sub$type <- 'no_cancer'
     
      
      temp_cancer <- bump_hunter(cases, 
                                 controls_sub, 
                                 boot_num = 3, 
                                 m_beta_thresh = m_thresh)
      
      cancer_probes <- inner_join(temp_cancer, g_ranges, by = 'start')
      
      # get significant 
      cancer_sig <- cancer_probes$probe[cancer_probes$p.value < m_thresh]
      
      # remove these from the subsetfeats 
      temp_cancer <- subset_feats[!subset_feats %in% cancer_sig]
      
      return(temp_cancer)
      
    
    }
    
    if(which_bumps == 'lfs') {
      
      # get all controls
      controls <- dat[grepl('Unaffected', dat$cancer_diagnosis_diagnoses),]
      
      # remove duplicates 
      controls <- controls[!duplicated(controls$ids),]
      
      # split age at 18 (216 months)
      # length(which(controls$age_sample_collection <216))
      controls_mut <- controls[grepl('Mut', controls$p53_germline),]
      controls_wt <- controls[grepl('WT', controls$p53_germline),]
      
      # remove NAs in ids 
      controls_mut <- controls_mut[!is.na(controls_mut$ids),]
      controls_wt <- controls_wt[!is.na(controls_wt$ids),]
      
      # change type 
      controls_mut$type <- 'mut'
      controls_wt$type <- 'wt'
      
      # Run bumphunter on kids and adults
      temp <- bump_hunter(dat_1 = controls_wt, 
                          dat_2 = controls_mut, 
                          boot_num = 3, 
                          m_beta_thresh = m_thresh)
      
      lfs_probes <-inner_join(temp, g_ranges, by = 'start' )
      
      # get significant 
      lfs_probes <- lfs_probes$probe[lfs_probes$p.value < m_thresh]
      
      # remove these from the subsetfeats 
      temp_lfs <- subset_feats[!subset_feats %in% lfs_probes]

      return(temp_lfs)
      
    }

}

# bumps = age, cancer, lfs
# subsets = island_probes, loci27k_probes, loci27k_island_probes



##########
# full data
##########
# age
age <- bump_hunter_master(dat = betaFull, 
                          subset_feats = NULL, 
                          m_thresh = 0.2, 
                          probe_start = 9, 
                          which_bumps = 'age')


# cancer
cancer <- bump_hunter_master(dat = betaFull, 
                             subset_feats = NULL, 
                             m_thresh = 0.2, 
                             probe_start = 9, 
                             which_bumps = 'cancer')

keep_these <- intersect_names[!intersect_names %in% cancer]
saveRDS(keep_these, paste0(feat_data, paste0('/', method, '_', 'no_cancer_m.rda')))


# cancer
lfs <- bump_hunter_master(dat = betaFull, 
                          subset_feats = NULL, 
                          m_thresh = 0.5, 
                          probe_start = 9, 
                          which_bumps = 'lfs')

saveRDS(lfs, paste0(feat_data, paste0('/', method, '_', 'lfs_m.rda')))


##########
# chr17 probes
##########
# age, chr_17_int
no_age_chr17 <- bump_hunter_master(dat = betaFull, 
                                    subset_feats = chr_17_int, 
                                    m_thresh = 0.01, 
                                    probe_start = 9, 
                                    which_bumps = 'age')

# remove from insterscted names and save 
keep_these <- chr_17_int[!chr_17_int %in% no_age_chr17]
saveRDS(keep_these, paste0(feat_data, paste0('/', method, '_', 'no_age_chr17_features.rda')))

# age, chr_17_int
no_cancer_chr17 <- bump_hunter_master(dat = betaFull, 
                                       subset_feats = chr_17_int, 
                                       m_thresh = 0.1, 
                                       probe_start = 9, 
                                       which_bumps = 'cancer')

# remove from intersected names and save
keep_these <- chr_17_int[!chr_17_int %in% no_cancer_chr17]
saveRDS(keep_these, paste0(feat_data, paste0('/', method, '_', 'no_cancer_chr17_features.rda')))

# get union, remove from intersected names and save 
union_island <- union(no_age_chr17, no_cancer_chr17)

# remove from union names and save
keep_these <- chr_17_int[!chr_17_int %in% union_island]
saveRDS(keep_these, paste0(feat_data, paste0('/', method, '_', 'no_age_cancer_union_chr17_features.rda')))

# age, chr_17_int
lfs_chr17 <- bump_hunter_master(dat = betaFull, 
                                 subset_feats = chr_17_int, 
                                 m_thresh = 1, 
                                 probe_start = 9, 
                                 which_bumps = 'lfs')

# directly save, no removal
saveRDS(lfs_chr17 , paste0(feat_data, paste0('/', method, '_', 'lfs_chr17_features.rda')))





##########
# chr17island probes is
##########
# age, chr_17_island_int
no_age_chr17_island <- bump_hunter_master(dat = betaFull, 
                                   subset_feats = chr_17_island_int, 
                                   m_thresh = 0.1, 
                                   probe_start = 9, 
                                   which_bumps = 'age')

# remove from insterscted names and save 
keep_these <- chr_17_island_int[!chr_17_island_int %in% no_age_chr17_island]
saveRDS(keep_these, paste0(feat_data, paste0('/', method, '_', 'no_age_chr17_island_features.rda')))

# age, chr_17_island_int
no_cancer_chr17_island <- bump_hunter_master(dat = betaFull, 
                                      subset_feats = chr_17_island_int, 
                                      m_thresh = 0.1, 
                                      probe_start = 9, 
                                      which_bumps = 'cancer')

# remove from intersected names and save
keep_these <- chr_17_island_int[!chr_17_island_int %in% no_cancer_chr17_island]
saveRDS(keep_these, paste0(feat_data, paste0('/', method, '_', 'no_cancer_chr17_island_features.rda')))

# get union, remove from intersected names and save 
union_island <- union(no_age_chr17_island, no_cancer_chr17_island)

# remove from union names and save
keep_these <- chr_17_island_int[!chr_17_island_int %in% union_island]
saveRDS(keep_these, paste0(feat_data, paste0('/', method, '_', 'no_age_cancer_union_chr17_island_features.rda')))

# age, chr_17_island_int
lfs_chr17_island <- bump_hunter_master(dat = betaFull, 
                                subset_feats = chr_17_island_int, 
                                m_thresh = 1, 
                                probe_start = 9, 
                                which_bumps = 'lfs')

# directly save, no removal
saveRDS(lfs_chr17_island , paste0(feat_data, paste0('/', method, '_', 'lfs_chr17_features.rda')))

##########
# 450k island probes
##########

# age, island_probes
no_age_island <- bump_hunter_master(dat = betaFull, 
                                    subset_feats = island_probes, 
                                    m_thresh = 0.1, 
                                    probe_start = 9, 
                                    which_bumps = 'age')

# remove from insterscted names and save 
keep_these <- island_probes[!island_probes %in% no_age_island]
saveRDS(keep_these, paste0(feat_data, paste0('/', method, '_', 'no_age_island_features.rda')))

# age, island_probes
no_cancer_island <- bump_hunter_master(dat = betaFull, 
                                       subset_feats = island_probes, 
                                       m_thresh = 0.1, 
                                       probe_start = 9, 
                                       which_bumps = 'cancer')

# remove from intersected names and save
keep_these <- island_probes[!island_probes %in% no_cancer_island]
saveRDS(keep_these, paste0(feat_data, paste0('/', method, '_', 'no_cancer_island_features.rda')))

# get union, remove from intersected names and save 
union_island <- union(no_age_island, no_cancer_island)

# remove from union names and save
keep_these <- island_probes[!island_probes %in% union_island]
saveRDS(keep_these, paste0(feat_data, paste0('/', method, '_', 'no_age_cancer_union_island_features.rda')))

# age, island_probes
lfs_island <- bump_hunter_master(dat = betaFull, 
                                 subset_feats = island_probes, 
                                 m_thresh = 0.1, 
                                 probe_start = 9, 
                                 which_bumps = 'lfs')

# directly save, no removal
saveRDS(lfs_island , paste0(feat_data, paste0('/', method, '_', 'lfs_island_features.rda')))


##########
# 27k probes
##########

# age, loci_27ks
no_age_loci27k <- bump_hunter_master(dat = betaFull, 
                                     subset_feats = loci27k_probes, 
                                     m_thresh = 0.1, 
                                     probe_start = 9, 
                                     which_bumps = 'age')

# remove from insterscted names and save 
keep_these <- loci27k_probes[!loci27k_probes %in% no_age_loci27k]
saveRDS(keep_these, paste0(feat_data, paste0('/', method, '_', 'no_age_loci27k_features.rda')))

# age, loci_27ks
no_cancer_loci27k <- bump_hunter_master(dat = betaFull, 
                                        subset_feats = loci27k_probes, 
                                        m_thresh = 0.1, 
                                        probe_start = 9, 
                                        which_bumps = 'cancer')

# remove from insterscted names and save 
keep_these <- loci27k_probes[!loci27k_probes %in% no_cancer_loci27k]
saveRDS(keep_these, paste0(feat_data, paste0('/', method, '_', 'no_cancer_loci27k_features.rda')))

# get union, remove from intersected names and save 
union_loci27k <- union(no_age_loci27k, no_cancer_loci27k)

# remove from union names and save
keep_these <- union_loci27k[!loci27k_probes %in% union_loci27k]
saveRDS(keep_these, paste0(feat_data, paste0('/', method, '_', 'no_age_cancer_union_loci27k_features.rda')))

# age, loci_27ks
lfs_loci27k <- bump_hunter_master(dat = betaFull, 
                                  subset_feats = loci27k_probes,  
                                  m_thresh = 0.1, 
                                  probe_start = 9, 
                                  which_bumps = 'lfs')

# directly save, no removal
saveRDS(lfs_loci27k, paste0(feat_data, paste0('/', method, '_', 'lfs_loci27k_features.rda')))



##########
# 27k island probes
##########


# age, loci_27ks
no_age_loci27k_island <- bump_hunter_master(dat = betaFull, 
                                            subset_feats = loci27k_island_probes, 
                                            m_thresh = 0.1, 
                                            probe_start = 9, 
                                            which_bumps = 'age')

# remove from insterscted names and save 
keep_these <- loci27k_probes[!loci27k_island_probes %in% no_age_loci27k_island]
saveRDS(keep_these, paste0(feat_data, paste0('/', method, '_', 'no_age_loci27k_island_features.rda')))

# age, loci_27ks
no_cancer_loci27k_island <- bump_hunter_master(dat = betaFull, 
                                               subset_feats = loci27k_island_probes, 
                                               m_thresh = 0.1, 
                                               probe_start = 9, 
                                               which_bumps = 'cancer')

# remove from insterscted names and save 
keep_these <- loci27k_island_probes[!loci27k_island_probes %in% no_cancer_loci27k_island]
saveRDS(keep_these, paste0(feat_data, paste0('/', method, '_', 'no_cancer_loci27k_island_features.rda')))

# get union, remove from intersected names and save 
union_loci27k_island <- union(no_age_loci27k_island, no_cancer_loci27k_island)

# remove from union names and save
keep_these <- union_loci27k_island[!loci27k_island_probes %in% union_loci27k]
saveRDS(keep_these, paste0(feat_data, paste0('/', method, '_', 'no_age_cancer_union_loci27k_island_features.rda')))

# age, loci_27ks
lfs_loci27k_island <- bump_hunter_master(dat = betaFull, 
                                         subset_feats = loci27k_island_probes,  
                                         m_thresh = 0.1, 
                                         probe_start = 9, 
                                         which_bumps = 'lfs')

# directly save, no removal
saveRDS(lfs_loci27k_island, paste0(feat_data, paste0('/', method, '_', 'lfs_loci27k_island_features.rda')))






