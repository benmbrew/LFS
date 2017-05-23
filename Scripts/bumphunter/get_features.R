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
load(paste0(model_data, '/modal_feat_surv.RData'))

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
# apply quan, even/uneven, gen, sam
##########
###quan
# quan even full 
quan_even_full_feat <- getProbe(quan_even_full)
quan_even_full_feat_all <- quan_even_full_feat[[1]]
quan_even_full_feat_sig <- quan_even_full_feat[[2]]
quan_even_full_feat_fwer <- quan_even_full_feat[[3]]

# quan even sub 
quan_even_sub_bal_feat <- getProbe(quan_even_sub_bal)
quan_even_sub_bal_feat_all <- quan_even_sub_bal_feat[[1]]
quan_even_sub_bal_feat_sig <- quan_even_sub_bal_feat[[2]]
quan_even_sub_bal_feat_fwer <- quan_even_sub_bal_feat[[3]]

# quan uneven full 
quan_uneven_full_feat <- getProbe(quan_uneven_full)
quan_uneven_full_feat_all <- quan_uneven_full_feat[[1]]
quan_uneven_full_feat_sig <- quan_uneven_full_feat[[2]]
quan_uneven_full_feat_fwer <- quan_uneven_full_feat[[3]]

# quan uneven sub 
quan_uneven_sub_feat <- getProbe(quan_uneven_sub)
quan_uneven_sub_feat_all <- quan_uneven_sub_feat[[1]]
quan_uneven_sub_feat_sig <- quan_uneven_sub_feat[[2]]
quan_uneven_sub_feat_fwer <- quan_uneven_sub_feat[[3]]

###funnorm
# funnorm even full 
funnorm_even_full_feat <- getProbe(funnorm_even_full)
funnorm_even_full_feat_all <- funnorm_even_full_feat[[1]]
funnorm_even_full_feat_sig <- funnorm_even_full_feat[[2]]
funnorm_even_full_feat_fwer <- funnorm_even_full_feat[[3]]

# funnorm even sub 
funnorm_even_sub_bal_feat <- getProbe(funnorm_even_sub_bal)
funnorm_even_sub_bal_feat_all <- funnorm_even_sub_bal_feat[[1]]
funnorm_even_sub_bal_feat_sig <- funnorm_even_sub_bal_feat[[2]]
funnorm_even_sub_bal_feat_fwer <- funnorm_even_sub_bal_feat[[3]]

# funnorm uneven full 
funnorm_uneven_full_feat <- getProbe(funnorm_uneven_full)
funnorm_uneven_full_feat_all <- funnorm_uneven_full_feat[[1]]
funnorm_uneven_full_feat_sig <- funnorm_uneven_full_feat[[2]]
funnorm_uneven_full_feat_fwer <- funnorm_uneven_full_feat[[3]]

# funnorm uneven sub 
funnorm_uneven_sub_feat <- getProbe(funnorm_uneven_sub)
funnorm_uneven_sub_feat_all <- funnorm_uneven_sub_feat[[1]]
funnorm_uneven_sub_feat_sig <- funnorm_uneven_sub_feat[[2]]
funnorm_uneven_sub_feat_fwer <- funnorm_uneven_sub_feat[[3]]

###raw
# raw even full 
raw_even_full_feat <- getProbe(raw_even_full)
raw_even_full_feat_all <- raw_even_full_feat[[1]]
raw_even_full_feat_sig <- raw_even_full_feat[[2]]
raw_even_full_feat_fwer <- raw_even_full_feat[[3]]

# raw even sub 
raw_even_sub_bal_feat <- getProbe(raw_even_sub_bal)
raw_even_sub_bal_feat_all <- raw_even_sub_bal_feat[[1]]
raw_even_sub_bal_feat_sig <- raw_even_sub_bal_feat[[2]]
raw_even_sub_bal_feat_fwer <- raw_even_sub_bal_feat[[3]]

# raw uneven full 
raw_uneven_full_feat <- getProbe(raw_uneven_full)
raw_uneven_full_feat_all <- raw_uneven_full_feat[[1]]
raw_uneven_full_feat_sig <- raw_uneven_full_feat[[2]]
raw_uneven_full_feat_fwer <- raw_uneven_full_feat[[3]]

# raw uneven sub 
raw_uneven_sub_feat <- getProbe(raw_uneven_sub)
raw_uneven_sub_feat_all <- raw_uneven_sub_feat[[1]]
raw_uneven_sub_feat_sig <- raw_uneven_sub_feat[[2]]
raw_uneven_sub_feat_fwer <- raw_uneven_sub_feat[[3]]

# save.image('/home/benbrew/Desktop/temp_bh_sig_fwer.RData')
# load('/home/benbrew/Desktop/temp_bh_sig_fwer.RData')


# remove data
rm(raw_uneven_sub_feat, raw_uneven_full_feat, raw_even_sub_bal_feat, raw_even_full_feat,
   funnorm_uneven_sub_feat, funnorm_uneven_full_feat, funnorm_even_sub_bal_feat, funnorm_even_full_feat,
   quan_uneven_sub_feat, quan_uneven_full_feat, quan_even_sub_bal_feat, quan_even_full_feat)

# remove additional data
rm(raw_uneven_full, raw_uneven_sub, raw_even_full, raw_even_sub_bal,
   quan_uneven_full, quan_uneven_sub, quan_even_full, quan_even_sub_bal,
   funnorm_uneven_full, funnorm_uneven_sub, funnorm_even_full, funnorm_even_sub_bal, cg_locations)
# remove empty ones
rm(raw_even_full_feat_fwer, raw_uneven_full_feat_fwer, raw_uneven_sub_feat_fwer)
# you may not need to run getRun because you only chose one diff setting 
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


##### quan
# even full
quan_even_full <- getRun(quan_even_full_feat_all, .10)
quan_even_full_fwer <- getRun(quan_even_full_feat_fwer, .10)
quan_even_full_sig <- getRun(quan_even_full_feat_sig, .10)

# even sub_bal
quan_even_sub_bal <- getRun(quan_even_sub_bal_feat_all, .10)
quan_even_sub_bal_fwer <- getRun(quan_even_sub_bal_feat_fwer, .10)
quan_even_sub_bal_sig <- getRun(quan_even_sub_bal_feat_sig, .10)


# uneven full
quan_uneven_full <- getRun(quan_uneven_full_feat_all, .10)
quan_uneven_full_fwer <- getRun(quan_uneven_full_feat_fwer, .10)
quan_uneven_full_sig <- getRun(quan_uneven_full_feat_sig, .10)

# uneven sub_bal
quan_uneven_sub <- getRun(quan_uneven_sub_feat_all, .10)
quan_uneven_sub_fwer <- getRun(quan_uneven_sub_feat_fwer, .10)
quan_uneven_sub_sig <- getRun(quan_uneven_sub_feat_sig, .10)

##### funnorm
# even full
funnorm_even_full <- getRun(funnorm_even_full_feat_all, .10)
funnorm_even_full_fwer <- getRun(funnorm_even_full_feat_fwer, .10)
funnorm_even_full_sig <- getRun(funnorm_even_full_feat_sig, .10)

# even sub_bal
funnorm_even_sub_bal <- getRun(funnorm_even_sub_bal_feat_all, .10)
funnorm_even_sub_bal_fwer <- getRun(funnorm_even_sub_bal_feat_fwer, .10)
funnorm_even_sub_bal_sig <- getRun(funnorm_even_sub_bal_feat_sig, .10)

# uneven full
funnorm_uneven_full <- getRun(funnorm_uneven_full_feat_all, .10)
funnorm_uneven_full_fwer <- getRun(funnorm_uneven_full_feat_fwer, .10)
funnorm_uneven_full_sig <- getRun(funnorm_uneven_full_feat_sig, .10)

# uneven sub_bal
funnorm_uneven_sub <- getRun(funnorm_uneven_sub_feat_all, .10)
funnorm_uneven_sub_fwer <- getRun(funnorm_uneven_sub_feat_fwer, .10)
funnorm_uneven_sub_sig <- getRun(funnorm_uneven_sub_feat_sig, .10)

##### raw
# even full
raw_even_full <- getRun(raw_even_full_feat_all, .10)
# raw_even_full_fwer <- getRun(raw_even_full_feat_fwer, .10)
raw_even_full_sig <- getRun(raw_even_full_feat_sig, .10)

# even sub_bal
raw_even_sub_bal <- getRun(raw_even_sub_bal_feat_all, .10)
raw_even_sub_bal_fwer <- getRun(raw_even_sub_bal_feat_fwer, .10)
raw_even_sub_bal_sig <- getRun(raw_even_sub_bal_feat_sig, .10)


# uneven full
raw_uneven_full <- getRun(raw_uneven_full_feat_all, .10)
# raw_uneven_full_fwer <- getRun(raw_uneven_full_feat_fwer, .10)
raw_uneven_full_sig <- getRun(raw_uneven_full_feat_sig, .10)

# uneven sub_bal
raw_uneven_sub <- getRun(raw_uneven_sub_feat_all, .10)
raw_uneven_sub_fwer <- getRun(raw_uneven_sub_feat_fwer, .10)
raw_uneven_sub_sig <- getRun(raw_uneven_sub_feat_sig, .10)



##############
##########
# remove unneccssary objects
##########
rm(list=ls(pattern="feat"))
rm(bumpHunterBalanced, getProbe, getRun, getBalAge)

# save.image(paste0(model_data, '/bh_feat_surv.RData'))
# load(paste0(model_data, '/bh_feat_surv.RData'))



