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
load(paste0(model_data, '/modal_feat_surv_20.RData'))


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
full_full_feat <- getProbe(full_full)
ff_feat_all <- full_full_feat[[1]]
ff_feat_sig <- full_full_feat[[2]]
ff_feat_fwer <- full_full_feat[[3]]

# sub full 
sub_full_feat <- getProbe(sub_full)
sf_feat_all <- sub_full_feat[[1]]
sf_feat_sig <- sub_full_feat[[2]]
sf_feat_fwer <- sub_full_feat[[3]]

# sub full 
full_sub_feat <- getProbe(full_sub)
fs_feat_all <- full_sub_feat[[1]]
fs_feat_sig <- full_sub_feat[[2]]
fs_feat_fwer <- full_sub_feat[[3]]

# sub full 
sub_sub_feat <- getProbe(sub_sub)
ss_feat_all <- sub_sub_feat[[1]]
ss_feat_sig <- sub_sub_feat[[2]]
ss_feat_fwer <- sub_sub_feat[[3]]



# remove data
rm(sub_sub_feat, sub_sub, full_sub_feat, full_sub,
   sub_full_feat, sub_full,  full_full_feat, full_full)


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
# ff
surv_ff_20_all <- getRun(ff_feat_all, .20)
surv_ff_20_sig <- getRun(ff_feat_sig, .20)
surv_ff_20_fwer <- getRun(ff_feat_fwer, .20)

# fs
surv_fs_20_all <- getRun(fs_feat_all, .20)
surv_fs_20_sig <- getRun(fs_feat_sig, .20)
surv_fs_20_fwer <- getRun(fs_feat_fwer, .20)

# sf
surv_sf_20_all <- getRun(sf_feat_all, .20)
surv_sf_20_sig <- getRun(sf_feat_sig, .20)
surv_sf_20_fwer <- getRun(sf_feat_fwer, .20)

# ss
surv_ss_20_all <- getRun(ss_feat_all, .20)
surv_ss_20_sig <- getRun(ss_feat_sig, .20)
surv_ss_20_fwer <- getRun(ss_feat_fwer, .20)



##############
##########
# remove unneccssary objects
##########
rm(list=ls(pattern="feat"))
rm(bumpHunterBalanced, getProbe, getRun, getBalAge)

##########
# remove empty object
##########
rm(surv_ff_20_fwer)
rm(surv_fs_20_fwer)
rm(surv_ss_20_fwer)

rm(clin_data, data_folder, home_folder, methyl_data, model_data,project_folder)


# save.image('/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/LFS/Data/model_data/surv_20.RData')
# load('/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/LFS/Data/model_data/surv_20.RData')



