####### Script will map regions to probes and save the final bh features for both cancer and lfs
# this is 8th step in pipeline

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
# load features
load(paste0(model_data, '/modal_feat.RData'))

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
# apply raw/quan, even/uneven, batch/unbatch
##########
##RAW
#raw even
raw_even <- getProbe(raw_even)
raw_even_bh <- raw_even[[1]]
raw_even_bh_sig <- raw_even[[2]]
raw_even_bh_fwer <- raw_even[[3]]

#raw uneven
raw_uneven <- getProbe(raw_uneven)
raw_uneven_bh <- raw_uneven[[1]]
raw_uneven_bh_sig <- raw_uneven[[2]]
raw_uneven_bh_fwer <- raw_uneven[[3]]

#raw even batch
raw_batch_even <- getProbe(raw_batch_even)
raw_batch_even_bh <- raw_batch_even[[1]]
raw_batch_even_bh_sig <- raw_batch_even[[2]] 
raw_batch_even_bh_fwer <- raw_batch_even[[3]]

#raw uneven
raw_batch_uneven <- getProbe(raw_batch_uneven)
raw_batch_uneven_bh <- raw_batch_uneven[[1]]
raw_batch_uneven_bh_sig <- raw_batch_uneven[[2]]
raw_batch_uneven_bh_fwer <- raw_batch_uneven[[3]]

rm(raw_even, raw_uneven, raw_batch_even, raw_batch_uneven)
##quan
#quan even
quan_even <- getProbe(quan_even)
quan_even_bh <- quan_even[[1]]
quan_even_bh_sig <- quan_even[[2]]
quan_even_bh_fwer <- quan_even[[3]]

#quan uneven
quan_uneven <- getProbe(quan_uneven)
quan_uneven_bh <- quan_uneven[[1]]
quan_uneven_bh_sig <- quan_uneven[[2]]
quan_uneven_bh_fwer <- quan_uneven[[3]]

#quan even batch
quan_batch_even <- getProbe(quan_batch_even)
quan_batch_even_bh <- quan_batch_even[[1]]
quan_batch_even_bh_sig <- quan_batch_even[[2]] 
quan_batch_even_bh_fwer <- quan_batch_even[[3]]

#quan uneven
quan_batch_uneven <- getProbe(quan_batch_uneven)
quan_batch_uneven_bh <- quan_batch_uneven[[1]]
quan_batch_uneven_bh_sig <- quan_batch_uneven[[2]]
quan_batch_uneven_bh_fwer <- quan_batch_uneven[[3]]

rm(quan_even, quan_uneven, quan_batch_even, quan_batch_uneven, cg_locations)

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
# raw 10
##########
#NOT BATCH
# even
raw_even_10 <- getRun(raw_even_bh, run_num = 0.10)
raw_even_sig_10 <- getRun(raw_even_bh_sig, run_num = 0.10)
raw_even_fwer_10 <- getRun(raw_even_bh_fwer, run_num = 0.10)

# uneven
raw_uneven_10 <- getRun(raw_uneven_bh, run_num = 0.10)
raw_uneven_sig_10 <- getRun(raw_uneven_bh_sig, run_num = 0.10)
raw_uneven_fwer_10 <- getRun(raw_uneven_bh_fwer, run_num = 0.10)

#BATCH
# even
raw_batch_even_batch_10 <- getRun(raw_batch_even_bh, run_num = 0.10)
raw_batch_even_batch_sig_10 <- getRun(raw_batch_even_bh_sig, run_num = 0.10)
raw_batch_even_batch_fwer_10 <- getRun(raw_batch_even_bh_fwer, run_num = 0.10)

# uneven
raw_batch_uneven_10 <- getRun(raw_batch_uneven_bh, run_num = 0.10)
raw_batch_uneven_sig_10 <- getRun(raw_batch_uneven_bh_sig, run_num = 0.10)
raw_batch_uneven_fwer_10 <- getRun(raw_batch_uneven_bh_fwer, run_num = 0.10)

##########
# raw 20
##########
#NOT BATCH
# even
raw_even_20 <- getRun(raw_even_bh, run_num = 0.20)
raw_even_sig_20 <- getRun(raw_even_bh_sig, run_num = 0.20)
raw_even_fwer_20 <- getRun(raw_even_bh_fwer, run_num = 0.20)

# uneven
raw_uneven_20 <- getRun(raw_uneven_bh, run_num = 0.20)
raw_uneven_sig_20 <- getRun(raw_uneven_bh_sig, run_num = 0.20)
raw_uneven_fwer_20 <- getRun(raw_uneven_bh_fwer, run_num = 0.20)

#BATCH
# even
raw_batch_even_batch_20 <- getRun(raw_batch_even_bh, run_num = 0.20)
raw_batch_even_batch_sig_20 <- getRun(raw_batch_even_bh_sig, run_num = 0.20)
raw_batch_even_batch_fwer_20 <- getRun(raw_batch_even_bh_fwer, run_num = 0.20)

# uneven
raw_batch_uneven_20 <- getRun(raw_batch_uneven_bh, run_num = 0.20)
raw_batch_uneven_sig_20 <- getRun(raw_batch_uneven_bh_sig, run_num = 0.20)
raw_batch_uneven_fwer_20 <- getRun(raw_batch_uneven_bh_fwer, run_num = 0.20)

##########
# raw 30
##########
#NOT BATCH
# even
raw_even_30 <- getRun(raw_even_bh, run_num = 0.30)
raw_even_sig_30 <- getRun(raw_even_bh_sig, run_num = 0.30)
raw_even_fwer_30 <- getRun(raw_even_bh_fwer, run_num = 0.30)

# uneven
raw_uneven_30 <- getRun(raw_uneven_bh, run_num = 0.30)
raw_uneven_sig_30 <- getRun(raw_uneven_bh_sig, run_num = 0.30)
raw_uneven_fwer_30 <- getRun(raw_uneven_bh_fwer, run_num = 0.30)

#BATCH
# even
raw_batch_even_batch_30 <- getRun(raw_batch_even_bh, run_num = 0.30)
raw_batch_even_batch_sig_30 <- getRun(raw_batch_even_bh_sig, run_num = 0.30)
raw_batch_even_batch_fwer_30 <- getRun(raw_batch_even_bh_fwer, run_num = 0.30)

# uneven
raw_batch_uneven_30 <- getRun(raw_batch_uneven_bh, run_num = 0.30)
raw_batch_uneven_sig_30 <- getRun(raw_batch_uneven_bh_sig, run_num = 0.30)
raw_batch_uneven_fwer_30 <- getRun(raw_batch_uneven_bh_fwer, run_num = 0.30)

##########
# raw 40
##########
#NOT BATCH
# even
raw_even_40 <- getRun(raw_even_bh, run_num = 0.40)
raw_even_sig_40 <- getRun(raw_even_bh_sig, run_num = 0.40)
raw_even_fwer_40 <- getRun(raw_even_bh_fwer, run_num = 0.40)

# uneven
raw_uneven_40 <- getRun(raw_uneven_bh, run_num = 0.40)
raw_uneven_sig_40 <- getRun(raw_uneven_bh_sig, run_num = 0.40)
raw_uneven_fwer_40 <- getRun(raw_uneven_bh_fwer, run_num = 0.40)

#BATCH
# even
raw_batch_even_batch_40 <- getRun(raw_batch_even_bh, run_num = 0.40)
raw_batch_even_batch_sig_40 <- getRun(raw_batch_even_bh_sig, run_num = 0.40)
raw_batch_even_batch_fwer_40 <- getRun(raw_batch_even_bh_fwer, run_num = 0.40)

# uneven
raw_batch_uneven_40 <- getRun(raw_batch_uneven_bh, run_num = 0.40)
raw_batch_uneven_sig_40 <- getRun(raw_batch_uneven_bh_sig, run_num = 0.40)
raw_batch_uneven_fwer_40 <- getRun(raw_batch_uneven_bh_fwer, run_num = 0.40)

##########
# raw 50
##########
#NOT BATCH
# even
raw_even_50 <- getRun(raw_even_bh, run_num = 0.50)
raw_even_sig_50 <- getRun(raw_even_bh_sig, run_num = 0.50)
raw_even_fwer_50 <- getRun(raw_even_bh_fwer, run_num = 0.50)

# uneven
raw_uneven_50 <- getRun(raw_uneven_bh, run_num = 0.50)
raw_uneven_sig_50 <- getRun(raw_uneven_bh_sig, run_num = 0.50)
raw_uneven_fwer_50 <- getRun(raw_uneven_bh_fwer, run_num = 0.50)

#BATCH
# even
raw_batch_even_batch_50 <- getRun(raw_batch_even_bh, run_num = 0.50)
raw_batch_even_batch_sig_50 <- getRun(raw_batch_even_bh_sig, run_num = 0.50)
raw_batch_even_batch_fwer_50 <- getRun(raw_batch_even_bh_fwer, run_num = 0.50)

# uneven
raw_batch_uneven_50 <- getRun(raw_batch_uneven_bh, run_num = 0.50)
raw_batch_uneven_sig_50 <- getRun(raw_batch_uneven_bh_sig, run_num = 0.50)
raw_batch_uneven_fwer_50 <- getRun(raw_batch_uneven_bh_fwer, run_num = 0.50)


##########
# quan 10
##########
#NOT BATCH
# even
quan_even_10 <- getRun(quan_even_bh, run_num = 0.10)
quan_even_sig_10 <- getRun(quan_even_bh_sig, run_num = 0.10)
quan_even_fwer_10 <- getRun(quan_even_bh_fwer, run_num = 0.10)

# uneven
quan_uneven_10 <- getRun(quan_uneven_bh, run_num = 0.10)
quan_uneven_sig_10 <- getRun(quan_uneven_bh_sig, run_num = 0.10)
quan_uneven_fwer_10 <- getRun(quan_uneven_bh_fwer, run_num = 0.10)

#BATCH
# even
quan_batch_even_batch_10 <- getRun(quan_batch_even_bh, run_num = 0.10)
quan_batch_even_batch_sig_10 <- getRun(quan_batch_even_bh_sig, run_num = 0.10)
quan_batch_even_batch_fwer_10 <- getRun(quan_batch_even_bh_fwer, run_num = 0.10)

# uneven
quan_batch_uneven_10 <- getRun(quan_batch_uneven_bh, run_num = 0.10)
quan_batch_uneven_sig_10 <- getRun(quan_batch_uneven_bh_sig, run_num = 0.10)
quan_batch_uneven_fwer_10 <- getRun(quan_batch_uneven_bh_fwer, run_num = 0.10)


##########
# quan 20
##########
#NOT BATCH
# even
quan_even_20 <- getRun(quan_even_bh, run_num = 0.20)
quan_even_sig_20 <- getRun(quan_even_bh_sig, run_num = 0.20)
quan_even_fwer_20 <- getRun(quan_even_bh_fwer, run_num = 0.20)

# uneven
quan_uneven_20 <- getRun(quan_uneven_bh, run_num = 0.20)
quan_uneven_sig_20 <- getRun(quan_uneven_bh_sig, run_num = 0.20)
quan_uneven_fwer_20 <- getRun(quan_uneven_bh_fwer, run_num = 0.20)

#BATCH
# even
quan_batch_even_batch_20 <- getRun(quan_batch_even_bh, run_num = 0.20)
quan_batch_even_batch_sig_20 <- getRun(quan_batch_even_bh_sig, run_num = 0.20)
quan_batch_even_batch_fwer_20 <- getRun(quan_batch_even_bh_fwer, run_num = 0.20)

# uneven
quan_batch_uneven_20 <- getRun(quan_batch_uneven_bh, run_num = 0.20)
quan_batch_uneven_sig_20 <- getRun(quan_batch_uneven_bh_sig, run_num = 0.20)
quan_batch_uneven_fwer_20 <- getRun(quan_batch_uneven_bh_fwer, run_num = 0.20)

##########
# quan 30
##########
#NOT BATCH
# even
quan_even_30 <- getRun(quan_even_bh, run_num = 0.30)
quan_even_sig_30 <- getRun(quan_even_bh_sig, run_num = 0.30)
quan_even_fwer_30 <- getRun(quan_even_bh_fwer, run_num = 0.30)

# uneven
quan_uneven_30 <- getRun(quan_uneven_bh, run_num = 0.30)
quan_uneven_sig_30 <- getRun(quan_uneven_bh_sig, run_num = 0.30)
quan_uneven_fwer_30 <- getRun(quan_uneven_bh_fwer, run_num = 0.30)

#BATCH
# even
quan_batch_even_batch_30 <- getRun(quan_batch_even_bh, run_num = 0.30)
quan_batch_even_batch_sig_30 <- getRun(quan_batch_even_bh_sig, run_num = 0.30)
quan_batch_even_batch_fwer_30 <- getRun(quan_batch_even_bh_fwer, run_num = 0.30)

# uneven
quan_batch_uneven_30 <- getRun(quan_batch_uneven_bh, run_num = 0.30)
quan_batch_uneven_sig_30 <- getRun(quan_batch_uneven_bh_sig, run_num = 0.30)
quan_batch_uneven_fwer_30 <- getRun(quan_batch_uneven_bh_fwer, run_num = 0.30)

##########
# quan 40
##########
#NOT BATCH
# even
quan_even_40 <- getRun(quan_even_bh, run_num = 0.40)
quan_even_sig_40 <- getRun(quan_even_bh_sig, run_num = 0.40)
quan_even_fwer_40 <- getRun(quan_even_bh_fwer, run_num = 0.40)

# uneven
quan_uneven_40 <- getRun(quan_uneven_bh, run_num = 0.40)
quan_uneven_sig_40 <- getRun(quan_uneven_bh_sig, run_num = 0.40)
quan_uneven_fwer_40 <- getRun(quan_uneven_bh_fwer, run_num = 0.40)

#BATCH
# even
quan_batch_even_batch_40 <- getRun(quan_batch_even_bh, run_num = 0.40)
quan_batch_even_batch_sig_40 <- getRun(quan_batch_even_bh_sig, run_num = 0.40)
quan_batch_even_batch_fwer_40 <- getRun(quan_batch_even_bh_fwer, run_num = 0.40)

# uneven
quan_batch_uneven_40 <- getRun(quan_batch_uneven_bh, run_num = 0.40)
quan_batch_uneven_sig_40 <- getRun(quan_batch_uneven_bh_sig, run_num = 0.40)
quan_batch_uneven_fwer_40 <- getRun(quan_batch_uneven_bh_fwer, run_num = 0.40)

##########
# quan 50
##########
#NOT BATCH
# even
quan_even_50 <- getRun(quan_even_bh, run_num = 0.50)
quan_even_sig_50 <- getRun(quan_even_bh_sig, run_num = 0.50)
quan_even_fwer_50 <- getRun(quan_even_bh_fwer, run_num = 0.50)

# uneven
quan_uneven_50 <- getRun(quan_uneven_bh, run_num = 0.50)
quan_uneven_sig_50 <- getRun(quan_uneven_bh_sig, run_num = 0.50)
quan_uneven_fwer_50 <- getRun(quan_uneven_bh_fwer, run_num = 0.50)

#BATCH
# even
quan_batch_even_batch_50 <- getRun(quan_batch_even_bh, run_num = 0.50)
quan_batch_even_batch_sig_50 <- getRun(quan_batch_even_bh_sig, run_num = 0.50)
quan_batch_even_batch_fwer_50 <- getRun(quan_batch_even_bh_fwer, run_num = 0.50)

# uneven
quan_batch_uneven_50 <- getRun(quan_batch_uneven_bh, run_num = 0.50)
quan_batch_uneven_sig_50 <- getRun(quan_batch_uneven_bh_sig, run_num = 0.50)
quan_batch_uneven_fwer_50 <- getRun(quan_batch_uneven_bh_fwer, run_num = 0.50)

##########
# remove unneccssary objects
##########
rm(list=ls(pattern="bh"))
rm(bumpHunterBalanced, getProbe, getRun)


##########
# remove any object filled with NAs
##########
rm(quan_batch_even_batch_sig_30, quan_batch_even_batch_sig_40, 
   quan_batch_even_batch_sig_50 , quan_batch_uneven_sig_20, 
   quan_batch_uneven_sig_30,quan_batch_uneven_sig_40,
   quan_batch_uneven_sig_50, quan_even_sig_20, 
   quan_even_sig_30, quan_even_sig_40, quan_even_sig_50,
   quan_uneven_sig_30, quan_uneven_sig_40, quan_uneven_sig_50,
   raw_batch_even_batch_sig_30, raw_batch_even_batch_sig_40, 
   raw_batch_even_batch_sig_50, raw_batch_uneven_sig_20, 
   raw_batch_uneven_sig_30,raw_batch_uneven_sig_40,
   raw_batch_uneven_sig_50,  raw_even_sig_20, 
   raw_even_sig_30, raw_even_sig_40, raw_even_sig_50,
   raw_uneven_sig_30, raw_uneven_sig_40, i)

save.image(paste0(model_data, '/bh_feat.RData'))

