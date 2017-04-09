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
even <- getProbe(even)
even_bh <- even[[1]]
even_bh_sig <- even[[2]]
even_bh_fwer <- even[[3]]

#quan even gen
even_gen <- getProbe(even_gen)
even_gen_bh <- even_gen[[1]]
even_gen_bh_sig <- even_gen[[2]]
even_gen_bh_fwer <- even_gen[[3]]

#quan even gen sen
even_gen_sen_type <- getProbe(even_gen_sen_type)
even_gen_sen_type_bh <- even_gen_sen_type[[1]]
even_gen_sen_type_bh_sig <- even_gen_sen_type[[2]]
even_gen_sen_type_bh_fwer <- even_gen_sen_type[[3]]

#quan even gen sam
even_gen_sam_type <- getProbe(even_gen_sam_type)
even_gen_sam_type_bh <- even_gen_sam_type[[1]]
even_gen_sam_type_bh_sig <- even_gen_sam_type[[2]]
even_gen_sam_type_bh_fwer <- even_gen_sam_type[[3]]

#quan even gen type
even_gen_type <- getProbe(even_type)
even_gen_type_bh <- even_gen_type[[1]]
even_gen_type_bh_sig <- even_gen_type[[2]]
even_gen_type_bh_fwer <- even_gen_type[[3]]

save.image('/home/benbrew/Desktop/march_22.RData')

## uneven

#quan uneven
uneven <- getProbe(uneven)
uneven_bh <- uneven[[1]]
uneven_bh_sig <- uneven[[2]]
uneven_bh_fwer <- uneven[[3]]

#quan uneven gen
uneven_gen <- getProbe(uneven_gen)
uneven_gen_bh <- uneven_gen[[1]]
uneven_gen_bh_sig <- uneven_gen[[2]]
uneven_gen_bh_fwer <- uneven_gen[[3]]

#quan uneven gen sen
uneven_gen_sen_type <- getProbe(uneven_gen_sen_type)
uneven_gen_sen_type_bh <- uneven_gen_sen_type[[1]]
uneven_gen_sen_type_bh_sig <- uneven_gen_sen_type[[2]]
uneven_gen_sen_type_bh_fwer <- uneven_gen_sen_type[[3]]

#quan uneven gen sam
uneven_gen_sam_type <- getProbe(uneven_gen_sam_type)
uneven_gen_sam_type_bh <- uneven_gen_sam_type[[1]]
uneven_gen_sam_type_bh_sig <- uneven_gen_sam_type[[2]]
uneven_gen_sam_type_bh_fwer <- uneven_gen_sam_type[[3]]

#quan uneven gen type
uneven_gen_type <- getProbe(uneven_type)
uneven_gen_type_bh <- uneven_gen_type[[1]]
uneven_gen_type_bh_sig <- uneven_gen_type[[2]]
uneven_gen_type_bh_fwer <- uneven_gen_type[[3]]

rm(even, even_gen, even_gen_sam, even_gen_sen, even_type,
   uneven, uneven_gen, uneven_gen_sam, uneven_gen_sen, uneven_type)

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
even_10 <- getRun(even_bh, run_num = 0.10)
even_sig_10 <- getRun(even_bh_sig, run_num = 0.10)
even_fwer_10 <- getRun(even_bh_fwer, run_num = 0.10)

# uneven
uneven_10 <- getRun(uneven_bh, run_num = 0.10)
uneven_sig_10 <- getRun(uneven_bh_sig, run_num = 0.10)
uneven_fwer_10 <- getRun(uneven_bh_fwer, run_num = 0.10)

###GEN

#even
even_gen_10 <- getRun(even_gen_bh, run_num = 0.10)
even_gen_sig_10 <- getRun(even_gen_bh_sig, run_num = 0.10)
even_gen_fwer_10 <- getRun(even_gen_bh_fwer, run_num = 0.10)

# uneven
uneven_gen_10 <- getRun(uneven_gen_bh, run_num = 0.10)
uneven_gen_sig_10 <- getRun(uneven_gen_bh_sig, run_num = 0.10)
uneven_gen_fwer_10 <- getRun(uneven_gen_bh_fwer, run_num = 0.10)

###GEN SAM

#even
even_gen_sam_type_10 <- getRun(even_gen_sam_type_bh, run_num = 0.10)
even_gen_sam_type_sig_10 <- getRun(even_gen_sam_type_bh_sig, run_num = 0.10)
even_gen_sam_type_fwer_10 <- getRun(even_gen_sam_type_bh_fwer, run_num = 0.10)

# uneven
uneven_gen_sam_type_10 <- getRun(uneven_gen_sam_type_bh, run_num = 0.10)
uneven_gen_sam_type_sig_10 <- getRun(uneven_gen_sam_type_bh_sig, run_num = 0.10)
uneven_gen_sam_type_fwer_10 <- getRun(uneven_gen_sam_type_bh_fwer, run_num = 0.10)


###GEN SEN

#even
even_gen_sen_type_10 <- getRun(even_gen_sen_type_bh, run_num = 0.10)
even_gen_sen_type_sig_10 <- getRun(even_gen_sen_type_bh_sig, run_num = 0.10)
even_gen_sen_type_fwer_10 <- getRun(even_gen_sen_type_bh_fwer, run_num = 0.10)

# uneven
uneven_gen_sen_type_10 <- getRun(uneven_gen_sen_type_bh, run_num = 0.10)
uneven_gen_sen_type_sig_10 <- getRun(uneven_gen_sen_type_bh_sig, run_num = 0.10)
uneven_gen_sen_type_fwer_10 <- getRun(uneven_gen_sen_type_bh_fwer, run_num = 0.10)

##Type

#even
even_gen_type_10 <- getRun(even_gen_type_bh, run_num = 0.10)
even_gen_type_sig_10 <- getRun(even_gen_type_bh_sig, run_num = 0.10)
even_gen_type_fwer_10 <- getRun(even_gen_type_bh_fwer, run_num = 0.10)

# uneven
uneven_gen_type_10 <- getRun(uneven_gen_type_bh, run_num = 0.10)
uneven_gen_type_sig_10 <- getRun(uneven_gen_type_bh_sig, run_num = 0.10)
uneven_gen_type_fwer_10 <- getRun(uneven_gen_type_bh_fwer, run_num = 0.10)
##########
# quan 15
##########
###NOT GEN
# even
even_15 <- getRun(even_bh, run_num = 0.15)
even_sig_15 <- getRun(even_bh_sig, run_num = 0.15)
even_fwer_15 <- getRun(even_bh_fwer, run_num = 0.15)

# uneven
uneven_15 <- getRun(uneven_bh, run_num = 0.15)
uneven_sig_15 <- getRun(uneven_bh_sig, run_num = 0.15)
uneven_fwer_15 <- getRun(uneven_bh_fwer, run_num = 0.15)

###GEN

#even
even_gen_15 <- getRun(even_gen_bh, run_num = 0.15)
even_gen_sig_15 <- getRun(even_gen_bh_sig, run_num = 0.15)
even_gen_fwer_15 <- getRun(even_gen_bh_fwer, run_num = 0.15)

# uneven
uneven_gen_15 <- getRun(uneven_gen_bh, run_num = 0.15)
uneven_gen_sig_15 <- getRun(uneven_gen_bh_sig, run_num = 0.15)
uneven_gen_fwer_15 <- getRun(uneven_gen_bh_fwer, run_num = 0.15)

###GEN SAM

#even
even_gen_sam_type_15 <- getRun(even_gen_sam_type_bh, run_num = 0.15)
even_gen_sam_type_sig_15 <- getRun(even_gen_sam_type_bh_sig, run_num = 0.15)
even_gen_sam_type_fwer_15 <- getRun(even_gen_sam_type_bh_fwer, run_num = 0.15)

# uneven
uneven_gen_sam_type_15 <- getRun(uneven_gen_sam_type_bh, run_num = 0.15)
uneven_gen_sam_type_sig_15 <- getRun(uneven_gen_sam_type_bh_sig, run_num = 0.15)
uneven_gen_sam_type_fwer_15 <- getRun(uneven_gen_sam_type_bh_fwer, run_num = 0.15)


###GEN SEN

#even
even_gen_sen_type_15 <- getRun(even_gen_sen_type_bh, run_num = 0.15)
even_gen_sen_type_sig_15 <- getRun(even_gen_sen_type_bh_sig, run_num = 0.15)
even_gen_sen_type_fwer_15 <- getRun(even_gen_sen_type_bh_fwer, run_num = 0.15)

# uneven
uneven_gen_sen_type_15 <- getRun(uneven_gen_sen_type_bh, run_num = 0.15)
uneven_gen_sen_type_sig_15 <- getRun(uneven_gen_sen_type_bh_sig, run_num = 0.15)
uneven_gen_sen_type_fwer_15 <- getRun(uneven_gen_sen_type_bh_fwer, run_num = 0.15)

##Type

#even
even_gen_type_15 <- getRun(even_gen_type_bh, run_num = 0.15)
even_gen_type_sig_15 <- getRun(even_gen_type_bh_sig, run_num = 0.15)
even_gen_type_fwer_15 <- getRun(even_gen_type_bh_fwer, run_num = 0.15)

# uneven
uneven_gen_type_15 <- getRun(uneven_gen_type_bh, run_num = 0.15)
uneven_gen_type_sig_15 <- getRun(uneven_gen_type_bh_sig, run_num = 0.15)
uneven_gen_type_fwer_15 <- getRun(uneven_gen_type_bh_fwer, run_num = 0.15)



##########
# quan 20
##########
###NOT GEN
# even
even_20 <- getRun(even_bh, run_num = 0.20)
even_sig_20 <- getRun(even_bh_sig, run_num = 0.20)
even_fwer_20 <- getRun(even_bh_fwer, run_num = 0.20)

# uneven
uneven_20 <- getRun(uneven_bh, run_num = 0.20)
uneven_sig_20 <- getRun(uneven_bh_sig, run_num = 0.20)
uneven_fwer_20 <- getRun(uneven_bh_fwer, run_num = 0.20)

###GEN

#even
even_gen_20 <- getRun(even_gen_bh, run_num = 0.20)
even_gen_sig_20 <- getRun(even_gen_bh_sig, run_num = 0.20)
even_gen_fwer_20 <- getRun(even_gen_bh_fwer, run_num = 0.20)

# uneven
uneven_gen_20 <- getRun(uneven_gen_bh, run_num = 0.20)
uneven_gen_sig_20 <- getRun(uneven_gen_bh_sig, run_num = 0.20)
uneven_gen_fwer_20 <- getRun(uneven_gen_bh_fwer, run_num = 0.20)

###GEN SAM

#even
even_gen_sam_type_20 <- getRun(even_gen_sam_type_bh, run_num = 0.20)
even_gen_sam_type_sig_20 <- getRun(even_gen_sam_type_bh_sig, run_num = 0.20)
even_gen_sam_type_fwer_20 <- getRun(even_gen_sam_type_bh_fwer, run_num = 0.20)

# uneven
uneven_gen_sam_type_20 <- getRun(uneven_gen_sam_type_bh, run_num = 0.20)
uneven_gen_sam_type_sig_20 <- getRun(uneven_gen_sam_type_bh_sig, run_num = 0.20)
uneven_gen_sam_type_fwer_20 <- getRun(uneven_gen_sam_type_bh_fwer, run_num = 0.20)


###GEN SEN

#even
even_gen_sen_type_20 <- getRun(even_gen_sen_type_bh, run_num = 0.20)
even_gen_sen_type_sig_20 <- getRun(even_gen_sen_type_bh_sig, run_num = 0.20)
even_gen_sen_type_fwer_20 <- getRun(even_gen_sen_type_bh_fwer, run_num = 0.20)

# uneven
uneven_gen_sen_type_20 <- getRun(uneven_gen_sen_type_bh, run_num = 0.20)
uneven_gen_sen_type_sig_20 <- getRun(uneven_gen_sen_type_bh_sig, run_num = 0.20)
uneven_gen_sen_type_fwer_20 <- getRun(uneven_gen_sen_type_bh_fwer, run_num = 0.20)

##Type

#even
even_gen_type_20 <- getRun(even_gen_type_bh, run_num = 0.20)
even_gen_type_sig_20 <- getRun(even_gen_type_bh_sig, run_num = 0.20)
even_gen_type_fwer_20 <- getRun(even_gen_type_bh_fwer, run_num = 0.20)

# uneven
uneven_gen_type_20 <- getRun(uneven_gen_type_bh, run_num = 0.20)
uneven_gen_type_sig_20 <- getRun(uneven_gen_type_bh_sig, run_num = 0.20)
uneven_gen_type_fwer_20 <- getRun(uneven_gen_type_bh_fwer, run_num = 0.20)

##########
# quan 25
##########
###NOT GEN
# even
even_25 <- getRun(even_bh, run_num = 0.25)
even_sig_25 <- getRun(even_bh_sig, run_num = 0.25)
even_fwer_25 <- getRun(even_bh_fwer, run_num = 0.25)

# uneven
uneven_25 <- getRun(uneven_bh, run_num = 0.25)
uneven_sig_25 <- getRun(uneven_bh_sig, run_num = 0.25)
uneven_fwer_25 <- getRun(uneven_bh_fwer, run_num = 0.25)

###GEN

#even
even_gen_25 <- getRun(even_gen_bh, run_num = 0.25)
even_gen_sig_25 <- getRun(even_gen_bh_sig, run_num = 0.25)
even_gen_fwer_25 <- getRun(even_gen_bh_fwer, run_num = 0.25)

# uneven
uneven_gen_25 <- getRun(uneven_gen_bh, run_num = 0.25)
uneven_gen_sig_25 <- getRun(uneven_gen_bh_sig, run_num = 0.25)
uneven_gen_fwer_25 <- getRun(uneven_gen_bh_fwer, run_num = 0.25)

###GEN SAM

#even
even_gen_sam_type_25 <- getRun(even_gen_sam_type_bh, run_num = 0.25)
even_gen_sam_type_sig_25 <- getRun(even_gen_sam_type_bh_sig, run_num = 0.25)
even_gen_sam_type_fwer_25 <- getRun(even_gen_sam_type_bh_fwer, run_num = 0.25)

# uneven
uneven_gen_sam_type_25 <- getRun(uneven_gen_sam_type_bh, run_num = 0.25)
uneven_gen_sam_type_sig_25 <- getRun(uneven_gen_sam_type_bh_sig, run_num = 0.25)
uneven_gen_sam_type_fwer_25 <- getRun(uneven_gen_sam_type_bh_fwer, run_num = 0.25)


###GEN SEN

#even
even_gen_sen_type_25 <- getRun(even_gen_sen_type_bh, run_num = 0.25)
even_gen_sen_type_sig_25 <- getRun(even_gen_sen_type_bh_sig, run_num = 0.25)
even_gen_sen_type_fwer_25 <- getRun(even_gen_sen_type_bh_fwer, run_num = 0.25)

# uneven
uneven_gen_sen_type_25 <- getRun(uneven_gen_sen_type_bh, run_num = 0.25)
uneven_gen_sen_type_sig_25 <- getRun(uneven_gen_sen_type_bh_sig, run_num = 0.25)
uneven_gen_sen_type_fwer_25 <- getRun(uneven_gen_sen_type_bh_fwer, run_num = 0.25)

##Type

#even
even_gen_type_25 <- getRun(even_gen_type_bh, run_num = 0.25)
even_gen_type_sig_25 <- getRun(even_gen_type_bh_sig, run_num = 0.25)
even_gen_type_fwer_25 <- getRun(even_gen_type_bh_fwer, run_num = 0.25)

# uneven
uneven_gen_type_25 <- getRun(uneven_gen_type_bh, run_num = 0.25)
uneven_gen_type_sig_25 <- getRun(uneven_gen_type_bh_sig, run_num = 0.25)
uneven_gen_type_fwer_25 <- getRun(uneven_gen_type_bh_fwer, run_num = 0.25)



##############
##########
# remove unneccssary objects
##########
rm(list=ls(pattern="bh"))
rm(bumpHunterBalanced, getProbe, getRun)

# save.image(paste0(model_data, '/bh_feat.RData'))
load(paste0(model_data, '/bh_feat.RData'))

#LOAD THIS AND REMOVE ONES THAT ACTUALLY ARE EMPTY
##########
# remove any object filled with NAs
##########
rm(even_gen_sam_type_sig_20, even_gen_sam_type_sig_25,
   even_gen_sen_type_sig_20, even_gen_sen_type_sig_25,
   even_gen_sig_25, even_gen_type_sig_20, even_gen_type_sig_25,
   even_sig_25, uneven_gen_sam_type_sig_20, uneven_gen_sam_type_sig_25,
   uneven_gen_sen_type_sig_20, uneven_gen_sen_type_sig_25,
   uneven_gen_sig_20, uneven_gen_sig_25,
   uneven_gen_type_sig_25, uneven_sig_20, 
   uneven_sig_25)

rm(uneven_gen_sam_type, uneven_gen_sen_type, uneven_gen_type,
   even_gen_sam_type, even_gen_sen_type, even_gen_type)

rm(cg_locations)

save.image(paste0(model_data, '/bh_feat.RData'))

