####### Script will map regions to probes and save the final bh features for both cancer and lfs
# this is 6th step in pipeline

##########
# initialize libraries
##########
library(minfi)
library(dplyr)
library(GenomicRanges)
library(sqldf)

# # IMPORTANT
# The ‘p.value’ is the percent of candidate regions obtained from the permutations that are as extreme as
# the observed region. These p-values should be interpreted with care as the theoretical proporties
# are not well understood. The ‘fwer’ is the proportion of permutations that had at least one region as
# extreme as the observed region. We compute p.values and FWER for the area of the regions 
# (as opposed to length and value as a pair) as well.


##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
model_data <- paste0(data_folder, '/model_data')

##########
# load features and cg_locations
##########
# load all data
load(paste0(model_data, '/modal_feat_pred_20.RData'))


# remove cases and controls for now
rm(list=ls(pattern="cases"))
rm(list=ls(pattern="controls"))

# load gene locations
cg_locations <- read.csv(paste0(model_data, '/cg_locations.csv'))

##########
# change column names of cg_locations for merging
##########
cg_locations$X <- NULL
colnames(cg_locations) <- paste0(colnames(cg_locations), '_', 'rgSet')

##########
# create function that grabs probe site and gene name for results from bumphunter
##########
# all bh features have some data with region length > 0. 
# cg_locations have no retions with length > 0.

getProbe <- function(data) {
  
  results <- list()
  results_data <- list()
  
  # first make seqnames in cg_locations and chr in tab character vectors for identification
  cg_locations$seqnames_rgSet <- as.character(cg_locations$seqnames_rgSet)
  data$chr <- as.character(data$chr)
  
  # use sql data table to merger validators with model_data based on age range
  result = sqldf("select * from cg_locations
                 inner join data
                 on cg_locations.start_rgSet between data.start and data.end")
  
  # keep only necessary columns
  result <- result[, c('chr' , 'start_rgSet','end_rgSet', 'probe_rgSet', 'p.value', 'fwer', 'run')]
  
  # rename cols
  colnames(result) <- c('chr', 'start', 'end', 'probe', 'p.value', 'fwer', 'run')
  
  # get sig results
  result_sig <- result[result$p.value < 0.05,]
  
  # get fwer results
  result_fwer <- result[result$fwer == 0,]
  
  return(list(result, result_sig, result_fwer))
  
}

##########
# full_full, full_sub, sub_full, sub_sub
##########
# full full 
bal_feat <- getProbe(wt_bal)
bal_feat_all <- bal_feat[[1]]
bal_feat_sig <- bal_feat[[2]]
bal_feat_fwer <- bal_feat[[3]]

# sub full 
bal_full_feat <- getProbe(wt_bal_full)
bal_full_feat_all <- bal_full_feat[[1]]
bal_full_feat_sig <- bal_full_feat[[2]]
bal_full_feat_fwer <- bal_full_feat[[3]]


# remove data
rm(wt_bal, wt_bal_full, bal_feat,bal_full_feat, cg_locations)


##########
# get run from each dataset
##########

getRun <- function(data, run_num)
{
  data <- data[data$run == run_num,]
  data <- data[!duplicated(data$probe),]
  data_feat <- as.data.frame(as.character(data$probe))
  return(data_feat)
}


##########
# assign feat to vectors
##########
# bal_feat
pred_bal_20_all <- getRun(bal_feat_all, .20)
pred_bal_20_sig <- getRun(bal_feat_sig, .20)
pred_bal_20_fwer <- getRun(bal_feat_fwer, .20)

# bal_feat
pred_bal_full_20_all <- getRun(bal_full_feat_all, .20)
pred_bal_full_20_sig <- getRun(bal_full_feat_sig, .20)
pred_bal_full_20_fwer <- getRun(bal_full_feat_fwer, .20)


##############
##########
# remove unneccssary objects
##########
rm(list=ls(pattern="feat"))
rm(bumpHunterBalanced, getProbe, getRun, getBalAge)
# 
# ##########
# remove empty object
##########


rm(clin_data, data_folder, home_folder, methyl_data, model_data,project_folder)


# save.image('/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/LFS/Data/model_data/pred_20.RData')
# load('/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/LFS/Data/model_data/pred_20.RData')



