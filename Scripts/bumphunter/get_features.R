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
# apply quan, even/uneven, gen, sam
##########
###quan

## even

#quan even
quan_even <- getProbe(quan_even)
quan_even_bh <- quan_even[[1]]
quan_even_bh_sig <- quan_even[[2]]
quan_even_bh_fwer <- quan_even[[3]]

#quan even gen
quan_even_gen <- getProbe(quan_even_gen)
quan_even_gen_bh <- quan_even_gen[[1]]
quan_even_gen_bh_sig <- quan_even_gen[[2]]
quan_even_gen_bh_fwer <- quan_even_gen[[3]]

#quan even gen sen
quan_even_gen_sen <- getProbe(quan_even_gen_sen)
quan_even_gen_sen_bh <- quan_even_gen_sen[[1]]
quan_even_gen_sen_bh_sig <- quan_even_gen_sen[[2]]
quan_even_gen_sen_bh_fwer <- quan_even_gen_sen[[3]]

#quan even gen sam
quan_even_gen_sam <- getProbe(quan_even_gen_sam)
quan_even_gen_sam_bh <- quan_even_gen_sam[[1]]
quan_even_gen_sam_bh_sig <- quan_even_gen_sam[[2]]
quan_even_gen_sam_bh_fwer <- quan_even_gen_sam[[3]]

# save.image('/home/benbrew/Desktop/feat_Temp_temp.RData')

## uneven

#quan uneven
quan_uneven <- getProbe(quan_uneven)
quan_uneven_bh <- quan_uneven[[1]]
quan_uneven_bh_sig <- quan_uneven[[2]]
quan_uneven_bh_fwer <- quan_uneven[[3]]

#quan uneven gen
quan_uneven_gen <- getProbe(quan_uneven_gen)
quan_uneven_gen_bh <- quan_uneven_gen[[1]]
quan_uneven_gen_bh_sig <- quan_uneven_gen[[2]]
quan_uneven_gen_bh_fwer <- quan_uneven_gen[[3]]

#quan uneven gen sen
quan_uneven_gen_sen <- getProbe(quan_uneven_gen_sen)
quan_uneven_gen_sen_bh <- quan_uneven_gen_sen[[1]]
quan_uneven_gen_sen_bh_sig <- quan_uneven_gen_sen[[2]]
quan_uneven_gen_sen_bh_fwer <- quan_uneven_gen_sen[[3]]

#quan uneven gen sam
quan_uneven_gen_sam <- getProbe(quan_uneven_gen_sam)
quan_uneven_gen_sam_bh <- quan_uneven_gen_sam[[1]]
quan_uneven_gen_sam_bh_sig <- quan_uneven_gen_sam[[2]]
quan_uneven_gen_sam_bh_fwer <- quan_uneven_gen_sam[[3]]

rm(quan_even, quan_even_gen, quan_even_gen_sam, quan_even_gen_sen,
   quan_uneven, quan_uneven_gen, quan_uneven_gen_sam, quan_uneven_gen_sen)

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
# quan 10
##########
###NOT GEN
# even
quan_even_10 <- getRun(quan_even_bh, run_num = 0.10)
quan_even_sig_10 <- getRun(quan_even_bh_sig, run_num = 0.10)
quan_even_fwer_10 <- getRun(quan_even_bh_fwer, run_num = 0.10)

# uneven
quan_uneven_10 <- getRun(quan_uneven_bh, run_num = 0.10)
quan_uneven_sig_10 <- getRun(quan_uneven_bh_sig, run_num = 0.10)
quan_uneven_fwer_10 <- getRun(quan_uneven_bh_fwer, run_num = 0.10)

###GEN

#even
quan_even_gen_10 <- getRun(quan_even_gen_bh, run_num = 0.10)
quan_even_gen_sig_10 <- getRun(quan_even_gen_bh_sig, run_num = 0.10)
quan_even_gen_fwer_10 <- getRun(quan_even_gen_bh_fwer, run_num = 0.10)

# uneven
quan_uneven_gen_10 <- getRun(quan_uneven_gen_bh, run_num = 0.10)
quan_uneven_gen_sig_10 <- getRun(quan_uneven_gen_bh_sig, run_num = 0.10)
quan_uneven_gen_fwer_10 <- getRun(quan_uneven_gen_bh_fwer, run_num = 0.10)

###GEN SAM

#even
quan_even_gen_sam_10 <- getRun(quan_even_gen_sam_bh, run_num = 0.10)
quan_even_gen_sam_sig_10 <- getRun(quan_even_gen_sam_bh_sig, run_num = 0.10)
quan_even_gen_sam_fwer_10 <- getRun(quan_even_gen_sam_bh_fwer, run_num = 0.10)

# uneven
quan_uneven_gen_sam_10 <- getRun(quan_uneven_gen_sam_bh, run_num = 0.10)
quan_uneven_gen_sam_sig_10 <- getRun(quan_uneven_gen_sam_bh_sig, run_num = 0.10)
quan_uneven_gen_sam_fwer_10 <- getRun(quan_uneven_gen_sam_bh_fwer, run_num = 0.10)


###GEN SEN

#even
quan_even_gen_sen_10 <- getRun(quan_even_gen_sen_bh, run_num = 0.10)
quan_even_gen_sen_sig_10 <- getRun(quan_even_gen_sen_bh_sig, run_num = 0.10)
quan_even_gen_sen_fwer_10 <- getRun(quan_even_gen_sen_bh_fwer, run_num = 0.10)

# uneven
quan_uneven_gen_sen_10 <- getRun(quan_uneven_gen_sen_bh, run_num = 0.10)
quan_uneven_gen_sen_sig_10 <- getRun(quan_uneven_gen_sen_bh_sig, run_num = 0.10)
quan_uneven_gen_sen_fwer_10 <- getRun(quan_uneven_gen_sen_bh_fwer, run_num = 0.10)

##############

##########
# quan 20
##########
###NOT GEN
# even
quan_even_20 <- getRun(quan_even_bh, run_num = 0.20)
quan_even_sig_20 <- getRun(quan_even_bh_sig, run_num = 0.20)
quan_even_fwer_20 <- getRun(quan_even_bh_fwer, run_num = 0.20)

# uneven
quan_uneven_20 <- getRun(quan_uneven_bh, run_num = 0.20)
quan_uneven_sig_20 <- getRun(quan_uneven_bh_sig, run_num = 0.20)
quan_uneven_fwer_20 <- getRun(quan_uneven_bh_fwer, run_num = 0.20)

###GEN

#even
quan_even_gen_20 <- getRun(quan_even_gen_bh, run_num = 0.20)
quan_even_gen_sig_20 <- getRun(quan_even_gen_bh_sig, run_num = 0.20)
quan_even_gen_fwer_20 <- getRun(quan_even_gen_bh_fwer, run_num = 0.20)

# uneven
quan_uneven_gen_20 <- getRun(quan_uneven_gen_bh, run_num = 0.20)
quan_uneven_gen_sig_20 <- getRun(quan_uneven_gen_bh_sig, run_num = 0.20)
quan_uneven_gen_fwer_20 <- getRun(quan_uneven_gen_bh_fwer, run_num = 0.20)

###GEN SAM

#even
quan_even_gen_sam_20 <- getRun(quan_even_gen_sam_bh, run_num = 0.20)
quan_even_gen_sam_sig_20 <- getRun(quan_even_gen_sam_bh_sig, run_num = 0.20)
quan_even_gen_sam_fwer_20 <- getRun(quan_even_gen_sam_bh_fwer, run_num = 0.20)

# uneven
quan_uneven_gen_sam_20 <- getRun(quan_uneven_gen_sam_bh, run_num = 0.20)
quan_uneven_gen_sam_sig_20 <- getRun(quan_uneven_gen_sam_bh_sig, run_num = 0.20)
quan_uneven_gen_sam_fwer_20 <- getRun(quan_uneven_gen_sam_bh_fwer, run_num = 0.20)


###GEN SEN

#even
quan_even_gen_sen_20 <- getRun(quan_even_gen_sen_bh, run_num = 0.20)
quan_even_gen_sen_sig_20 <- getRun(quan_even_gen_sen_bh_sig, run_num = 0.20)
quan_even_gen_sen_fwer_20 <- getRun(quan_even_gen_sen_bh_fwer, run_num = 0.20)

# uneven
quan_uneven_gen_sen_20 <- getRun(quan_uneven_gen_sen_bh, run_num = 0.20)
quan_uneven_gen_sen_sig_20 <- getRun(quan_uneven_gen_sen_bh_sig, run_num = 0.20)
quan_uneven_gen_sen_fwer_20 <- getRun(quan_uneven_gen_sen_bh_fwer, run_num = 0.20)

##############




##########
# quan 30
##########
###NOT GEN
# even
quan_even_30 <- getRun(quan_even_bh, run_num = 0.30)
quan_even_sig_30 <- getRun(quan_even_bh_sig, run_num = 0.30)
quan_even_fwer_30 <- getRun(quan_even_bh_fwer, run_num = 0.30)

# uneven
quan_uneven_30 <- getRun(quan_uneven_bh, run_num = 0.30)
quan_uneven_sig_30 <- getRun(quan_uneven_bh_sig, run_num = 0.30)
quan_uneven_fwer_30 <- getRun(quan_uneven_bh_fwer, run_num = 0.30)

###GEN

#even
quan_even_gen_30 <- getRun(quan_even_gen_bh, run_num = 0.30)
quan_even_gen_sig_30 <- getRun(quan_even_gen_bh_sig, run_num = 0.30)
quan_even_gen_fwer_30 <- getRun(quan_even_gen_bh_fwer, run_num = 0.30)

# uneven
quan_uneven_gen_30 <- getRun(quan_uneven_gen_bh, run_num = 0.30)
quan_uneven_gen_sig_30 <- getRun(quan_uneven_gen_bh_sig, run_num = 0.30)
quan_uneven_gen_fwer_30 <- getRun(quan_uneven_gen_bh_fwer, run_num = 0.30)

###GEN SAM

#even
quan_even_gen_sam_30 <- getRun(quan_even_gen_sam_bh, run_num = 0.30)
quan_even_gen_sam_sig_30 <- getRun(quan_even_gen_sam_bh_sig, run_num = 0.30)
quan_even_gen_sam_fwer_30 <- getRun(quan_even_gen_sam_bh_fwer, run_num = 0.30)

# uneven
quan_uneven_gen_sam_30 <- getRun(quan_uneven_gen_sam_bh, run_num = 0.30)
quan_uneven_gen_sam_sig_30 <- getRun(quan_uneven_gen_sam_bh_sig, run_num = 0.30)
quan_uneven_gen_sam_fwer_30 <- getRun(quan_uneven_gen_sam_bh_fwer, run_num = 0.30)


###GEN SEN

#even
quan_even_gen_sen_30 <- getRun(quan_even_gen_sen_bh, run_num = 0.30)
quan_even_gen_sen_sig_30 <- getRun(quan_even_gen_sen_bh_sig, run_num = 0.30)
quan_even_gen_sen_fwer_30 <- getRun(quan_even_gen_sen_bh_fwer, run_num = 0.30)

# uneven
quan_uneven_gen_sen_30 <- getRun(quan_uneven_gen_sen_bh, run_num = 0.30)
quan_uneven_gen_sen_sig_30 <- getRun(quan_uneven_gen_sen_bh_sig, run_num = 0.30)
quan_uneven_gen_sen_fwer_30 <- getRun(quan_uneven_gen_sen_bh_fwer, run_num = 0.30)

##############
##########
# remove unneccssary objects
##########
rm(list=ls(pattern="bh"))
rm(bumpHunterBalanced, getProbe, getRun)


##########
# remove any object filled with NAs
##########
rm(quan_even_gen_sam_sig_30, quan_even_gen_sen_sig_20,
   quan_even_gen_sen_sig_30, quan_even_gen_sig_30,
   quan_even_sig_30, quan_uneven_gen_sam_sig_20,
   quan_uneven_gen_sam_sig_30, quan_uneven_gen_sen_sig_20,
   quan_uneven_gen_sen_sig_30, quan_uneven_gen_sig_20,
   quan_uneven_gen_sig_30, quan_uneven_sig_30, cg_locations)

save.image(paste0(model_data, '/bh_feat.RData'))

