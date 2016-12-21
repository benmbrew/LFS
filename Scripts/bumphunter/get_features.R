####### Script will map regions to probes and save the final bh features for both cancer and lfs
# this is 6th step in pipeline

##########
# initialize libraries
##########
library(minfi)
library(bumphunter)
library(dplyr)
library(FDb.InfiniumMethylation.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(impute)
library(GenomicRanges)
library(biovizBase)
library(GEOquery)
library(IlluminaHumanMethylation450kmanifest)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
model_data <- paste0(data_folder, '/model_data')

#########
# load bump_hunter_lfs
#########
load(paste0(model_data, '/beta_p53_bh.RData'))

#########
# load bump_hunter_cancer data
#########
load(paste0(model_data, '/beta_cancer_bh.RData'))

#########
# change column names of cg_locations for merging
#########
cg_locations$X <- NULL
colnames(cg_locations) <- paste0(colnames(cg_locations), '_', 'rgSet')

###########
# create function that grabs probe site and gene name for results from bumphunter
###########

# all bh features have some data with region length > 0. 
# cg_locations have no retions with length > 0.

results <- list()
results_data <- list()

getProbe <- function(data) {
  
  # first make seqnames in cg_locations and chr in tab character vectors for identification
  cg_locations$seqnames_rgSet <- as.character(cg_locations$seqnames_rgSet)
  data$chr <- as.character(data$chr)
  
  # loop through chromosomes and match then loop through locations and match
  for (i in unique(data$chr)) {
    
    sub_data <- data[data$chr == i,] # some length >0
    sub_rg <- cg_locations[cg_locations$seqnames_rgSet == i, ] # no length >0
    
    for (j in 1:nrow(sub_data)) {
      
      chr_data <- sub_data[j,]
      
      if (any(sub_rg$start_rgSet >= chr_data$start & sub_rg$start_rgSet <= chr_data$end)) {
        
        results[[j]] <- cbind(sub_rg[sub_rg$start_rgSet >= chr_data$start & sub_rg$start_rgSet <= chr_data$end,], chr_data)
        
      }
      
      print(j)
      
    }
    
    results_data[[i]] <- do.call('rbind', results)
    
    print(i)
  }
  
  # combine into data frame 
  dat_cg <- as.data.frame(do.call('rbind', results_data))
  
  # clean up rownames and create column for cg site 
  dat_cg$probe_rgSet <- as.character(dat_cg$probe_rgSet)
  rownames(dat_cg) <- NULL
  
  # remove chromosome info from probe colum
  dat_cg$probe_rgSet <- substr(dat_cg$probe_rgSet, nchar(dat_cg$probe_rgSet) - 9, nchar(dat_cg$probe))
  
  # for those that have an extra 1 at the end and missing a c, paste a c back on and remove 1 
  dat_cg$probe_rgSet <- ifelse(!grepl('c', dat_cg$probe_rgSet), paste0('c', dat_cg$probe_rgSet), dat_cg$probe_rgSet)
  dat_cg$probe_rgSet <- ifelse(nchar(dat_cg$probe_rgSet) == 11, substr(dat_cg$probe_rgSet, 1, 10), dat_cg$probe_rgSet)
  
  # remove duplicates from dat_cg
  dat_cg <- dat_cg[!duplicated(dat_cg$probe_rgSet),]
  
  # keep only necessary columns
  dat_cg <- dat_cg[, c('chr' , 'start_rgSet','end_rgSet', 'probe_rgSet', 'p.value', 'fwer')]
  
  # rename
  colnames(dat_cg) <- c('chr', 'start', 'end', 'probe', 'p.value', 'fwer')
  
  return(dat_cg)
  
}

##########
# apply function to cancer bh
##########

# beta raw
beta_raw_bal_counts_cancer_features <- getProbe(beta_raw_bal_counts_cancer)
beta_raw_bal_cancer_features <- getProbe(beta_raw_bal_cancer)
beta_raw_unbal_cancer_features <- getProbe(beta_raw_unbal_cancer)

# beta swan
beta_swan_bal_counts_cancer_features <- getProbe(beta_swan_bal_counts_cancer)
beta_swan_bal_cancer_features <- getProbe(beta_swan_bal_cancer)
beta_swan_unbal_cancer_features <- getProbe(beta_swan_unbal_cancer)

# beta quan
beta_quan_bal_counts_cancer_features <- getProbe(beta_quan_bal_counts_cancer)
beta_quan_bal_cancer_features <- getProbe(beta_quan_bal_cancer)
beta_quan_unbal_cancer_features <- getProbe(beta_quan_unbal_cancer)

# beta funnorm
beta_funnorm_bal_counts_cancer_features <- getProbe(beta_funnorm_bal_counts_cancer)
beta_funnorm_bal_cancer_features <- getProbe(beta_funnorm_bal_cancer)
beta_funnorm_unbal_cancer_features <- getProbe(beta_funnorm_unbal_cancer)

# beta raw_cancer_intersection

# beta swan_cancer_intersection

# beta quan_cancer_intersection

# beta funnorm_cancer_intersection

# beta raw_cancer_union

# beta swan_cancer_union

# beta quan_cancer_union

# beta funnorm_cancer_union


##########
# apply function to p53 bh
##########

# beta raw
beta_raw_bal_counts_p53_features <- getProbe(beta_raw_bal_counts_p53)
beta_raw_bal_p53_features <- getProbe(beta_raw_bal_p53)
beta_raw_unbal_p53_features <- getProbe(beta_raw_unbal_p53)

# beta swan
beta_swan_bal_counts_p53_features <- getProbe(beta_swan_bal_counts_p53)
beta_swan_bal_p53_features <- getProbe(beta_swan_bal_p53)
beta_swan_unbal_p53_features <- getProbe(beta_swan_unbal_p53)

# beta quan
beta_quan_bal_counts_p53_features <- getProbe(beta_quan_bal_counts_p53)
beta_quan_bal_p53_features <- getProbe(beta_quan_bal_p53)
beta_quan_unbal_p53_features <- getProbe(beta_quan_unbal_p53)

# beta funnorm
beta_funnorm_bal_counts_p53_features <- getProbe(beta_funnorm_bal_counts_p53)
beta_funnorm_bal_p53_features <- getProbe(beta_funnorm_bal_p53)
beta_funnorm_unbal_p53_features <- getProbe(beta_funnorm_unbal_p53)

# beta raw_p53_intersection

# beta swan_p53_intersection

# beta quan_p53_intersection

# beta funnorm_p53_intersection

# beta raw_p53_union

# beta swan_p53_union

# beta quan_p53_union

# beta funnorm_p53_union


##########
# remove unneeded objects
##########

rm(beta_raw_bal_counts_p53, beta_raw_bal_p53, beta_raw_unbal_p53, 
   beta_swan_bal_counts_p53, beta_swan_bal_p53, beta_swan_unbal_p53,
   beta_quan_bal_counts_p53, beta_quan_bal_p53, beta_quan_unbal_p53,
   beta_funnorm_bal_counts_p53, beta_funnorm_bal_p53, beta_funnorm_unbal_p53,
   beta_raw_bal_counts_cancer, beta_raw_bal_cancer, beta_raw_unbal_cancer, 
   beta_swan_bal_counts_cancer, beta_swan_bal_cancer, beta_swan_unbal_cancer,
   beta_quan_bal_counts_cancer, beta_quan_bal_cancer, beta_quan_unbal_cancer,
   beta_funnorm_bal_counts_cancer, beta_funnorm_bal_cancer, beta_funnorm_unbal_cancer,
   cg_locations, cg_locations)

save.image(paste0(model_data, '/bh_features.RData'))

