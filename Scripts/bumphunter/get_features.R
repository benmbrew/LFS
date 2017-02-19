####### Script will map regions to probes and save the final bh features for both cancer and lfs
# this is 7th step in pipeline

##########
# initialize libraries
##########
library(minfi)
library(dplyr)
library(GenomicRanges)
library(sqldf)


##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
model_data <- paste0(data_folder, '/model_data')

##########
# load bump_hunter_lfs
##########
# save new bh 
raw_bh <- readRDS(paste0(model_data, '/raw_bh.rda'))
quan_bh <- readRDS(paste0(model_data, '/quan_bh.rda'))
swan_bh <- readRDS(paste0(model_data, '/swan_bh.rda'))
funnorm_bh <- readRDS(paste0(model_data, '/funnorm_bh.rda'))
cg_locations <- read.csv(paste0(model_data, '/cg_locations.csv'))

# load(paste0(model_data, '/beta_p53_bh.RData'))

##########
# load bump_hunter_cancer data
##########
# load(paste0(model_data, '/beta_cancer_bh.RData'))

##########
# change column names of cg_locations for merging
##########
cg_locations$X <- NULL
colnames(cg_locations) <- paste0(colnames(cg_locations), '_', 'rgSet')


# ##########
# # get intersections
# ##########
# 
# getIntersection <- function(vec1 = NULL, 
#                             vec2 = NULL, 
#                             vec3 = NULL, 
#                             vec4 = NULL,
#                             num_vecs) 
# {
#   if(num_vecs == 2) {
#     probe_list <- Reduce(intersect, list(vec1,vec2))
#     
#   } 
#   
#   if(num_vecs == 3) {
#     probe_list <- Reduce(intersect, list(vec1,vec2,vec3))
#     
#   } 
#   
#   if (num_vecs == 4) {
#     probe_list <- Reduce(intersect, list(vec1,vec2,vec3,vec4))
#     
#   }
#   combined_probes <- paste(probe_list, collapse = '|')
#   dat <- cg_locations[cg_locations$probe_rgSet %in% probe_list,]
#   stopifnot(length(probe_list) == nrow(dat))
#   return(dat)
# }
# 
# ##########
# # get union
# ##########
# getUnion <- function(vec1 = NULL, 
#                      vec2 = NULL, 
#                      vec3 = NULL, 
#                      vec4 = NULL,
#                      num_vecs) 
# {
#   if(num_vecs == 2) {
#     probe_list <- Reduce(union, list(vec1,vec2))
#     
#   } 
#   
#   if(num_vecs == 3) {
#     probe_list <- Reduce(union, list(vec1,vec2,vec3))
#     
#   } 
#   
#   if (num_vecs == 4) {
#     probe_list <- Reduce(union, list(vec1,vec2,vec3,vec4))
#     
#   }
#   
#   combined_probes <- paste(probe_list, collapse = '|')
#   dat <- cg_locations[cg_locations$probe_rgSet %in% probe_list,]
#   stopifnot(length(probe_list) == nrow(dat))
#   return(dat)
# }

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
  
  # # loop through chromosomes and match then loop through locations and match
  # for (i in unique(data$chr)) {
  # 
  #   sub_data <- data[data$chr == i,] # some region length >0
  #   sub_rg <- cg_locations[cg_locations$seqnames_rgSet == i, ] # no region length >0
  # 
  #   for (j in 1:nrow(sub_data)) {
  # 
  #     chr_data <- sub_data[j,]
  # 
  #     if (any(sub_rg$start_rgSet >= chr_data$start & sub_rg$start_rgSet <= chr_data$end)) {
  # 
  #       results[[j]] <- cbind(sub_rg[sub_rg$start_rgSet >= chr_data$start & sub_rg$start_rgSet <= chr_data$end,], chr_data)
  # 
  #     }
  # 
  #   }
  #   
  #   results_data[[i]] <- do.call('rbind', results)
  #   
  #   print(i)
  # }
  # 
  # combine into data frame 
  # dat_cg <- as.data.frame(do.call('rbind', results_data))
  # 
  # # clean up rownames and create column for cg site 
  # dat_cg$probe_rgSet <- as.character(dat_cg$probe_rgSet)
  # rownames(dat_cg) <- NULL
  
  # # remove chromosome info from probe colum
  # dat_cg$probe_rgSet <- substr(dat_cg$probe_rgSet, nchar(dat_cg$probe_rgSet) - 9, nchar(dat_cg$probe))
  # 
  # # for those that have an extra 1 at the end and missing a c, paste a c back on and remove 1 
  # dat_cg$probe_rgSet <- ifelse(!grepl('c', dat_cg$probe_rgSet), paste0('c', dat_cg$probe_rgSet), dat_cg$probe_rgSet)
  # dat_cg$probe_rgSet <- ifelse(nchar(dat_cg$probe_rgSet) == 11, substr(dat_cg$probe_rgSet, 1, 10), dat_cg$probe_rgSet)
  # 
  # # remove duplicates from dat_cg
  # dat_cg <- dat_cg[!duplicated(dat_cg$probe_rgSet),]
  # 
  # # keep only necessary columns
  # dat_cg <- dat_cg[, c('chr' , 'start_rgSet','end_rgSet', 'probe_rgSet', 'p.value', 'fwer')]
  # 
  # # rename
  # colnames(dat_cg) <- c('chr', 'start', 'end', 'probe', 'p.value', 'fwer')
  # 
  # dat_cg_sig <- dat_cg[dat_cg$p.value < 0.05,]
  
  # keep only necessary columns
  result <- result[, c('chr' , 'start_rgSet','end_rgSet', 'probe_rgSet', 'p.value', 'fwer', 'run')]
  
  # rename
  colnames(result) <- c('chr', 'start', 'end', 'probe', 'p.value', 'fwer', 'run')
  result_sig <- result[result$p.value < 0.05,]
  
  return(list(result, result_sig))
  
}

##########
# apply function to cancer bh for both sig and all on cancer
##########
#raw
raw <- getProbe(raw_bh)
raw_bh <- raw[[1]]
raw_bh_sig <- raw[[2]]


#quan
quan <- getProbe(quan_bh)
quan_bh <- quan[[1]]
quan_bh_sig <- quan[[2]]

#raw
swan <- getProbe(swan_bh)
swan_bh <- swan[[1]]
swan_bh_sig <- swan[[2]]

#funnorm
funnorm <- getProbe(funnorm_bh)
funnorm_bh <- funnorm[[1]]
funnorm_bh_sig <- funnorm[[2]]

##########
# get intersection/union for each preprocessing method
##########

getSet <- function(data, set)
  
{
  # get threshold to loop through
  DELTA_BETA_THRESH = c(0.07, 0.08, 0.09,0.10, 0.11, 0.12, 0.13, 0.14, 0.15) # DNAm difference threshold
  
  # get list 
  feature_list <- list()
  
  if (set == 'intersection') {
    
    for (i in 1:length(DELTA_BETA_THRESH)) {
      result <- data[data$run == DELTA_BETA_THRESH[i],]
      result <- result[!duplicated(result$probe),]
      feature_list[[i]] <- as.character(result$probe)
    }
    
    feat <- Reduce(intersect, list(feature_list[[1]], feature_list[[2]], 
                                   feature_list[[3]], feature_list[[4]],
                                   feature_list[[5]], feature_list[[6]],
                                   feature_list[[7]], feature_list[[8]],
                                                      feature_list[[9]]))
    
  }
  
  if (set == 'union') {
    
    for (i in 1:length(DELTA_BETA_THRESH)) {
      result <- data[data$run == DELTA_BETA_THRESH[i],]
      result <- result[!duplicated(result$probe),]
      feature_list[[i]] <- as.character(result$probe)
    
    }
    feat <- Reduce(union, list(feature_list[[1]], feature_list[[2]], 
                               feature_list[[3]], feature_list[[4]],
                               feature_list[[5]], feature_list[[6]],
                               feature_list[[7]], feature_list[[8]],
                               feature_list[[9]]))
  }
  
  return(feat)
}

##########
# raw
##########
#intersection
raw_int_feat <- getSet(raw_bh, set = 'intersection')
raw_int_sig_feat <- getSet(raw_bh_sig, set = 'intersection')

#union
raw_union_feat <- getSet(raw_bh, set = 'union')
raw_union_sig_feat <- getSet(raw_bh_sig, set = 'union')

##########
# quan
##########
#intersection
quan_int_feat <- getSet(quan_bh, set = 'intersection')
quan_int_sig_feat <- getSet(quan_bh_sig, set = 'intersection')

#union
quan_union_feat <- getSet(quan_bh, set = 'union')
quan_union_sig_feat <- getSet(quan_bh_sig, set = 'union')

##########
# swan
##########
#intersection
swan_int_feat <- getSet(swan_bh, set = 'intersection')
swan_int_sig_feat <- getSet(swan_bh_sig, set = 'intersection')

#union
swan_union_feat <- getSet(swan_bh, set = 'union')
swan_union_sig_feat <- getSet(swan_bh_sig, set = 'union')

##########
# funnorm
##########
#intersection
funnorm_int_feat <- getSet(funnorm_bh, set = 'intersection')
funnorm_int_sig_feat <- getSet(funnorm_bh_sig, set = 'intersection')

#union
funnorm_union_feat <- getSet(funnorm_bh, set = 'union')
funnorm_union_sig_feat <- getSet(funnorm_bh_sig, set = 'union')

##########
# get each run
##########

getRun <- function(data, run_num)
{
  data <- data[data$run == run_num,]
  data <- data[!duplicated(data$probe),]
  data_feat <- as.character(data$probe)
  return(data_feat)
}

##########
# raw
##########
# 0.07
raw_bh_07 <- getRun(raw_bh, run_num = 0.07)
raw_bh_sig_07 <- getRun(raw_bh_sig, run_num = 0.07)

# 0.08
raw_bh_08 <- getRun(raw_bh, run_num = 0.08)
raw_bh_sig_08 <- getRun(raw_bh_sig, run_num = 0.08)

# 0.09
raw_bh_09 <- getRun(raw_bh, run_num = 0.09)
raw_bh_sig_09 <- getRun(raw_bh_sig, run_num = 0.09)

# 0.10
raw_bh_10 <- getRun(raw_bh, run_num = 0.10)
raw_bh_sig_10 <- getRun(raw_bh_sig, run_num = 0.10)

# 0.11
raw_bh_11 <- getRun(raw_bh, run_num = 0.11)
raw_bh_sig_11 <- getRun(raw_bh_sig, run_num = 0.11)

# 0.12
raw_bh_12 <- getRun(raw_bh, run_num = 0.12)
raw_bh_sig_12 <- getRun(raw_bh_sig, run_num = 0.12)

# 0.13
raw_bh_13 <- getRun(raw_bh, run_num = 0.13)
raw_bh_sig_13 <- getRun(raw_bh_sig, run_num = 0.13)

# 0.14
raw_bh_14 <- getRun(raw_bh, run_num = 0.14)
raw_bh_sig_14 <- getRun(raw_bh_sig, run_num = 0.14)

# 0.15
raw_bh_15 <- getRun(raw_bh, run_num = 0.15)
raw_bh_sig_15 <- getRun(raw_bh_sig, run_num = 0.15)

##########
# quan
##########
# 0.07
quan_bh_07 <- getRun(quan_bh, run_num = 0.07)
quan_bh_sig_07 <- getRun(quan_bh_sig, run_num = 0.07)

# 0.08
quan_bh_08 <- getRun(quan_bh, run_num = 0.08)
quan_bh_sig_08 <- getRun(quan_bh_sig, run_num = 0.08)

# 0.09
quan_bh_09 <- getRun(quan_bh, run_num = 0.09)
quan_bh_sig_09 <- getRun(quan_bh_sig, run_num = 0.09)

# 0.10
quan_bh_10 <- getRun(quan_bh, run_num = 0.10)
quan_bh_sig_10 <- getRun(quan_bh_sig, run_num = 0.10)

# 0.11
quan_bh_11 <- getRun(quan_bh, run_num = 0.11)
quan_bh_sig_11 <- getRun(quan_bh_sig, run_num = 0.11)

# 0.12
quan_bh_12 <- getRun(quan_bh, run_num = 0.12)
quan_bh_sig_12 <- getRun(quan_bh_sig, run_num = 0.12)

# 0.13
quan_bh_13 <- getRun(quan_bh, run_num = 0.13)
quan_bh_sig_13 <- getRun(quan_bh_sig, run_num = 0.13)

# 0.14
quan_bh_14 <- getRun(quan_bh, run_num = 0.14)
quan_bh_sig_14 <- getRun(quan_bh_sig, run_num = 0.14)

# 0.15
quan_bh_15 <- getRun(quan_bh, run_num = 0.15)
quan_bh_sig_15 <- getRun(quan_bh_sig, run_num = 0.15)


##########
# swan
##########
# 0.07
swan_bh_07 <- getRun(swan_bh, run_num = 0.07)
swan_bh_sig_07 <- getRun(swan_bh_sig, run_num = 0.07)

# 0.08
swan_bh_08 <- getRun(swan_bh, run_num = 0.08)
swan_bh_sig_08 <- getRun(swan_bh_sig, run_num = 0.08)

# 0.09
swan_bh_09 <- getRun(swan_bh, run_num = 0.09)
swan_bh_sig_09 <- getRun(swan_bh_sig, run_num = 0.09)

# 0.10
swan_bh_10 <- getRun(swan_bh, run_num = 0.10)
swan_bh_sig_10 <- getRun(swan_bh_sig, run_num = 0.10)

# 0.11
swan_bh_11 <- getRun(swan_bh, run_num = 0.11)
swan_bh_sig_11 <- getRun(swan_bh_sig, run_num = 0.11)

# 0.12
swan_bh_12 <- getRun(swan_bh, run_num = 0.12)
swan_bh_sig_12 <- getRun(swan_bh_sig, run_num = 0.12)

# 0.13
swan_bh_13 <- getRun(swan_bh, run_num = 0.13)
swan_bh_sig_13 <- getRun(swan_bh_sig, run_num = 0.13)

# 0.14
swan_bh_14 <- getRun(swan_bh, run_num = 0.14)
swan_bh_sig_14 <- getRun(swan_bh_sig, run_num = 0.14)

# 0.15
swan_bh_15 <- getRun(swan_bh, run_num = 0.15)
swan_bh_sig_15 <- getRun(swan_bh_sig, run_num = 0.15)


##########
# funnorm
##########
# 0.07
funnorm_bh_07 <- getRun(funnorm_bh, run_num = 0.07)
funnorm_bh_sig_07 <- getRun(funnorm_bh_sig, run_num = 0.07)

# 0.08
funnorm_bh_08 <- getRun(funnorm_bh, run_num = 0.08)
funnorm_bh_sig_08 <- getRun(funnorm_bh_sig, run_num = 0.08)

# 0.09
funnorm_bh_09 <- getRun(funnorm_bh, run_num = 0.09)
funnorm_bh_sig_09 <- getRun(funnorm_bh_sig, run_num = 0.09)

# 0.10
funnorm_bh_10 <- getRun(funnorm_bh, run_num = 0.10)
funnorm_bh_sig_10 <- getRun(funnorm_bh_sig, run_num = 0.10)

# 0.11
funnorm_bh_11 <- getRun(funnorm_bh, run_num = 0.11)
funnorm_bh_sig_11 <- getRun(funnorm_bh_sig, run_num = 0.11)

# 0.12
funnorm_bh_12 <- getRun(funnorm_bh, run_num = 0.12)
funnorm_bh_sig_12 <- getRun(funnorm_bh_sig, run_num = 0.12)

# 0.13
funnorm_bh_13 <- getRun(funnorm_bh, run_num = 0.13)
funnorm_bh_sig_13 <- getRun(funnorm_bh_sig, run_num = 0.13)

# 0.14
funnorm_bh_14 <- getRun(funnorm_bh, run_num = 0.14)
funnorm_bh_sig_14 <- getRun(funnorm_bh_sig, run_num = 0.14)

# 0.15
funnorm_bh_15 <- getRun(funnorm_bh, run_num = 0.15)
funnorm_bh_sig_15 <- getRun(funnorm_bh_sig, run_num = 0.15)


##########
# remove unneccssary objects
##########
rm(cg_locations, 
   raw, quan, swan, funnorm,
   funnorm_bh, funnorm_bh_sig,
   quan_bh, quan_bh_sig,
   swan_bh, swan_bh_sig,
   raw_bh, raw_bh_sig)


# save.image(paste0(model_data, 'bh_feat_new.RData'))
# load('/home/benbrew/Desktop/get_feat_temp.RData')

# # beta raw all
# beta_raw_bal_counts_cancer_features <- getProbe(beta_raw_bal_counts_cancer)[[1]]
# beta_raw_bal_cancer_features <- getProbe(beta_raw_bal_cancer)[[1]]
# beta_raw_unbal_cancer_features <- getProbe(beta_raw_unbal_cancer)[[1]]
# 
# # beta swan all
# beta_swan_bal_counts_cancer_features <- getProbe(beta_swan_bal_counts_cancer)[[1]]
# beta_swan_bal_cancer_features <- getProbe(beta_swan_bal_cancer)[[1]]
# beta_swan_unbal_cancer_features <- getProbe(beta_swan_unbal_cancer)[[1]]
# 
# # beta quan all
# beta_quan_bal_counts_cancer_features <- getProbe(beta_quan_bal_counts_cancer)[[1]]
# beta_quan_bal_cancer_features <- getProbe(beta_quan_bal_cancer)[[1]]
# beta_quan_unbal_cancer_features <- getProbe(beta_quan_unbal_cancer)[[1]]
# 
# # beta funnorm all
# beta_funnorm_bal_counts_cancer_features <- getProbe(beta_funnorm_bal_counts_cancer)[[1]]
# beta_funnorm_bal_cancer_features <- getProbe(beta_funnorm_bal_cancer)[[1]]
# beta_funnorm_unbal_cancer_features <- getProbe(beta_funnorm_unbal_cancer)[[1]]
# 
# # beta raw sig
# beta_raw_bal_counts_cancer_sig_features <- getProbe(beta_raw_bal_counts_cancer)[[2]]
# beta_raw_bal_cancer_sig_features <- getProbe(beta_raw_bal_cancer)[[2]]
# beta_raw_unbal_cancer_sig_features <- getProbe(beta_raw_unbal_cancer)[[2]]
# 
# # beta swan sig
# beta_swan_bal_counts_cancer_sig_features <- getProbe(beta_swan_bal_counts_cancer)[[2]]
# beta_swan_bal_cancer_sig_features <- getProbe(beta_swan_bal_cancer)[[2]]
# beta_swan_unbal_cancer_sig_features <- getProbe(beta_swan_unbal_cancer)[[2]]
# 
# # beta quan sig
# beta_quan_bal_counts_cancer_sig_features <- getProbe(beta_quan_bal_counts_cancer)[[2]]
# beta_quan_bal_cancer_sig_features <- getProbe(beta_quan_bal_cancer)[[2]]
# beta_quan_unbal_cancer_sig_features <- getProbe(beta_quan_unbal_cancer)[[2]]
# 
# # beta funnorm sig
# beta_funnorm_bal_counts_cancer_sig_features <- getProbe(beta_funnorm_bal_counts_cancer)[[2]]
# beta_funnorm_bal_cancer_sig_features <- getProbe(beta_funnorm_bal_cancer)[[2]]
# beta_funnorm_unbal_cancer_sig_features <- getProbe(beta_funnorm_unbal_cancer)[[2]]
# 
# ##########
# # apply function to p53 bh for both sig and all on p53
# ##########
# 
# # beta raw all
# beta_raw_bal_counts_p53_features <- getProbe(beta_raw_bal_counts_p53)[[1]]
# beta_raw_bal_p53_features <- getProbe(beta_raw_bal_p53)[[1]]
# beta_raw_unbal_p53_features <- getProbe(beta_raw_unbal_p53)[[1]]
# 
# # beta swan all
# beta_swan_bal_counts_p53_features <- getProbe(beta_swan_bal_counts_p53)[[1]]
# beta_swan_bal_p53_features <- getProbe(beta_swan_bal_p53)[[1]]
# beta_swan_unbal_p53_features <- getProbe(beta_swan_unbal_p53)[[1]]
# 
# # beta quan all
# beta_quan_bal_counts_p53_features <- getProbe(beta_quan_bal_counts_p53)[[1]]
# beta_quan_bal_p53_features <- getProbe(beta_quan_bal_p53)[[1]]
# beta_quan_unbal_p53_features <- getProbe(beta_quan_unbal_p53)[[1]]
# 
# # beta funnorm all
# beta_funnorm_bal_counts_p53_features <- getProbe(beta_funnorm_bal_counts_p53)[[1]]
# beta_funnorm_bal_p53_features <- getProbe(beta_funnorm_bal_p53)[[1]]
# beta_funnorm_unbal_p53_features <- getProbe(beta_funnorm_unbal_p53)[[1]]
# 
# # beta raw sig
# beta_raw_bal_counts_p53_sig_features <- getProbe(beta_raw_bal_counts_p53)[[2]]
# beta_raw_bal_p53_sig_features <- getProbe(beta_raw_bal_p53)[[2]]
# beta_raw_unbal_p53_sig_features <- getProbe(beta_raw_unbal_p53)[[2]]
# 
# # beta swan sig
# beta_swan_bal_counts_p53_sig_features <- getProbe(beta_swan_bal_counts_p53)[[2]]
# beta_swan_bal_p53_sig_features <- getProbe(beta_swan_bal_p53)[[2]]
# beta_swan_unbal_p53_sig_features <- getProbe(beta_swan_unbal_p53)[[2]]
# 
# # beta quan sig
# beta_quan_bal_counts_p53_sig_features <- getProbe(beta_quan_bal_counts_p53)[[2]]
# beta_quan_bal_p53_sig_features <- getProbe(beta_quan_bal_p53)[[2]]
# beta_quan_unbal_p53_sig_features <- getProbe(beta_quan_unbal_p53)[[2]]
# 
# # beta funnorm sig
# beta_funnorm_bal_counts_p53_sig_features <- getProbe(beta_funnorm_bal_counts_p53)[[2]]
# beta_funnorm_bal_p53_sig_features <- getProbe(beta_funnorm_bal_p53)[[2]]
# beta_funnorm_unbal_p53_sig_features <- getProbe(beta_funnorm_unbal_p53)[[2]]
# # 
# # save.image('/home/benbrew/Desktop/temp_bh.RData')
# # load('/home/benbrew/Desktop/temp_bh.RData')
# # 

###########################################################################################################################
# this part of script will get intersection and unions of various sets of probe features

# INTERSECTION

##########
# cancer intersection 
##########

# beta_raw_cancer_intersection
beta_raw_cancer_intersection_features  <- getIntersection(beta_raw_bal_counts_cancer_features$probe,
                                                          beta_raw_bal_cancer_features$probe,
                                                          beta_raw_unbal_cancer_features$probe,
                                                          num_vecs = 3)

# beta_swan_cancer_intersection
beta_swan_cancer_intersection_features  <- getIntersection(beta_swan_bal_counts_cancer_features$probe,
                                                           beta_swan_bal_cancer_features$probe,
                                                           beta_swan_unbal_cancer_features$probe,
                                                           num_vecs = 3)

# beta_quan_cancer_intersection
beta_quan_cancer_intersection_features  <- getIntersection(beta_quan_bal_counts_cancer_features$probe,
                                                           beta_quan_bal_cancer_features$probe,
                                                           beta_quan_unbal_cancer_features$probe,
                                                           num_vecs = 3)

# beta_funnorm_cancer_intersection
beta_funnorm_cancer_intersection_features  <- getIntersection(beta_funnorm_bal_counts_cancer_features$probe,
                                                              beta_funnorm_bal_cancer_features$probe,
                                                              beta_funnorm_unbal_cancer_features$probe,
                                                              num_vecs = 3)

# total cancer intersection
beta_cancer_intersection_features <- getIntersection(beta_raw_cancer_intersection_features$probe_rgSet,
                                                     beta_swan_cancer_intersection_features$probe_rgSet, 
                                                     beta_quan_cancer_intersection_features$probe_rgSet,
                                                     beta_funnorm_cancer_intersection_features$probe_rgSet,
                                                     num_vecs = 4)

# beta_bal_counts_cancer_intersection
beta_bal_counts_cancer_intersection_features  <- getIntersection(beta_raw_bal_counts_cancer_features$probe,
                                                                 beta_swan_bal_counts_cancer_features$probe,
                                                                 beta_quan_bal_counts_cancer_features$probe,
                                                                 beta_funnorm_bal_counts_cancer_features$probe,
                                                                 num_vecs = 4)
##########
# cancer intersection sig
##########

# beta_raw_cancer_intersection
beta_raw_cancer_intersection_sig_features  <- getIntersection(beta_raw_bal_counts_cancer_sig_features$probe,
                                                              beta_raw_bal_cancer_sig_features$probe,
                                                              beta_raw_unbal_cancer_sig_features$probe,
                                                              num_vecs = 3)

# beta_swan_cancer_intersection
beta_swan_cancer_intersection_sig_features  <- getIntersection(beta_swan_bal_counts_cancer_sig_features$probe,
                                                               beta_swan_bal_cancer_sig_features$probe,
                                                               beta_swan_unbal_cancer_sig_features$probe,
                                                               num_vecs = 3)

# beta_quan_cancer_intersection here error
beta_quan_cancer_intersection_sig_features  <- getIntersection(beta_quan_bal_counts_cancer_sig_features$probe,
                                                               beta_quan_bal_cancer_sig_features$probe,
                                                               beta_quan_unbal_cancer_sig_features$probe,
                                                               num_vecs = 3)

# beta_funnorm_cancer_intersection
beta_funnorm_cancer_intersection_sig_features  <- getIntersection(beta_funnorm_bal_counts_cancer_sig_features$probe,
                                                                  beta_funnorm_bal_cancer_sig_features$probe,
                                                                  beta_funnorm_unbal_cancer_sig_features$probe,
                                                                  num_vecs = 3)

# total cancer intersection
beta_cancer_intersection_sig_features <- getIntersection(beta_raw_cancer_intersection_sig_features$probe_rgSet,
                                                         beta_swan_cancer_intersection_sig_features$probe_rgSet, 
                                                         beta_quan_cancer_intersection_sig_features$probe_rgSet,
                                                         beta_funnorm_cancer_intersection_sig_features,
                                                         num_vecs = 4)

# beta_bal_counts_cancer_intersection
beta_bal_counts_cancer_intersection_sig_features  <- getIntersection(beta_raw_bal_counts_cancer_sig_features$probe,
                                                                     beta_swan_bal_counts_cancer_sig_features$probe,
                                                                     beta_quan_bal_counts_cancer_sig_features$probe,
                                                                     beta_funnorm_bal_counts_cancer_sig_features$probe,
                                                                     num_vecs = 4)

##########
# p53 intersection 
##########

# beta_raw_p53_intersection
beta_raw_p53_intersection_features  <- getIntersection(beta_raw_bal_counts_p53_features$probe,
                                                       beta_raw_bal_p53_features$probe,
                                                       beta_raw_unbal_p53_features$probe,
                                                       num_vecs = 3)

# beta_swan_p53_intersection
beta_swan_p53_intersection_features  <- getIntersection(beta_swan_bal_counts_p53_features$probe,
                                                        beta_swan_bal_p53_features$probe,
                                                        beta_swan_unbal_p53_features$probe,
                                                        num_vecs = 3)

# beta_quan_p53_intersection
beta_quan_p53_intersection_features  <- getIntersection(beta_quan_bal_counts_p53_features$probe,
                                                        beta_quan_bal_p53_features$probe,
                                                        beta_quan_unbal_p53_features$probe,
                                                        num_vecs = 3)

# beta_funnorm_p53_intersection
beta_funnorm_p53_intersection_features  <- getIntersection(beta_funnorm_bal_counts_p53_features$probe,
                                                           beta_funnorm_bal_p53_features$probe,
                                                           beta_funnorm_unbal_p53_features$probe,
                                                           num_vecs = 3)

# total p53 intersection
beta_p53_intersection_features <- getIntersection(beta_raw_p53_intersection_features$probe_rgSet,
                                                  beta_swan_p53_intersection_features$probe_rgSet, 
                                                  beta_quan_p53_intersection_features$probe_rgSet,
                                                  beta_funnorm_p53_intersection_features$probe_rgSet,
                                                  num_vecs = 4)

# beta_bal_counts_p53_intersection
beta_bal_counts_p53_intersection_features  <- getIntersection(beta_raw_bal_counts_p53_features$probe,
                                                              beta_swan_bal_counts_p53_features$probe,
                                                              beta_quan_bal_counts_p53_features$probe,
                                                              beta_funnorm_bal_counts_p53_features$probe,
                                                              num_vecs = 4)

##########
# p53 intersection sig
##########

# beta_raw_p53_intersection
beta_raw_p53_intersection_sig_features  <- getIntersection(beta_raw_bal_counts_p53_sig_features$probe,
                                                           beta_raw_bal_p53_sig_features$probe,
                                                           beta_raw_unbal_p53_sig_features$probe,
                                                           num_vecs = 3)

# beta_swan_p53_intersection
beta_swan_p53_intersection_sig_features  <- getIntersection(beta_swan_bal_counts_p53_sig_features$probe,
                                                            beta_swan_bal_p53_sig_features$probe,
                                                            beta_swan_unbal_p53_sig_features$probe,
                                                            num_vecs = 3)

# beta_quan_p53_intersection
beta_quan_p53_intersection_sig_features  <- getIntersection(beta_quan_bal_counts_p53_sig_features$probe,
                                                            beta_quan_bal_p53_sig_features$probe,
                                                            beta_quan_unbal_p53_sig_features$probe,
                                                            num_vecs = 3)

# beta_funnorm_p53_intersection
beta_funnorm_p53_intersection_sig_features  <- getIntersection(beta_funnorm_bal_counts_p53_sig_features$probe,
                                                               beta_funnorm_bal_p53_sig_features$probe,
                                                               beta_funnorm_unbal_p53_sig_features$probe,
                                                               num_vecs = 3)

# total p53 intersection
beta_p53_intersection_sig_features <- getIntersection(beta_raw_p53_intersection_sig_features$probe_rgSet,
                                                      beta_swan_p53_intersection_sig_features$probe_rgSet, 
                                                      beta_quan_p53_intersection_sig_features$probe_rgSet,
                                                      beta_funnorm_p53_intersection_sig_features$probe_rgSet,
                                                      num_vecs = 4)

# beta_bal_counts_p53_intersection
beta_bal_counts_p53_intersection_sig_features  <- getIntersection(beta_raw_bal_counts_p53_sig_features$probe,
                                                                  beta_swan_bal_counts_p53_sig_features$probe,
                                                                  beta_quan_bal_counts_p53_sig_features$probe,
                                                                  beta_funnorm_bal_counts_p53_sig_features$probe,
                                                                  num_vecs = 4)

#############################################################################################################################
# UNION

##########
# cancer union 
##########

# beta_raw_cancer_union
beta_raw_cancer_union_features  <- getUnion(beta_raw_bal_counts_cancer_features$probe,
                                            beta_raw_bal_cancer_features$probe,
                                            beta_raw_unbal_cancer_features$probe,
                                            num_vecs = 3)

# beta_swan_cancer_union
beta_swan_cancer_union_features  <- getUnion(beta_swan_bal_counts_cancer_features$probe,
                                             beta_swan_bal_cancer_features$probe,
                                             beta_swan_unbal_cancer_features$probe,
                                             num_vecs = 3)

# beta_quan_cancer_union
beta_quan_cancer_union_features  <- getUnion(beta_quan_bal_counts_cancer_features$probe,
                                             beta_quan_bal_cancer_features$probe,
                                             beta_quan_unbal_cancer_features$probe,
                                             num_vecs = 3)

# beta_funnorm_cancer_union
beta_funnorm_cancer_union_features  <- getUnion(beta_funnorm_bal_counts_cancer_features$probe,
                                                beta_funnorm_bal_cancer_features$probe,
                                                beta_funnorm_unbal_cancer_features$probe,
                                                num_vecs = 3)

# total cancer union
beta_cancer_union_features <- getUnion(beta_raw_cancer_union_features$probe_rgSet,
                                       beta_swan_cancer_union_features$probe_rgSet, 
                                       beta_quan_cancer_union_features$probe_rgSet,
                                       beta_funnorm_cancer_union_features$probe_rgSet,
                                       num_vecs = 4)

# beta_bal_counts_cancer_union
beta_bal_counts_cancer_union_features  <- getUnion(beta_raw_bal_counts_cancer_features$probe,
                                                   beta_swan_bal_counts_cancer_features$probe,
                                                   beta_quan_bal_counts_cancer_features$probe,
                                                   beta_funnorm_bal_counts_cancer_features$probe,
                                                   num_vecs = 4)

##########
# cancer union sig
##########

# beta_raw_cancer_union
beta_raw_cancer_union_sig_features  <- getUnion(beta_raw_bal_counts_cancer_sig_features$probe,
                                                beta_raw_bal_cancer_sig_features$probe,
                                                beta_raw_unbal_cancer_sig_features$probe,
                                                num_vecs = 3)

# beta_swan_cancer_union
beta_swan_cancer_union_sig_features  <- getUnion(beta_swan_bal_counts_cancer_sig_features$probe,
                                                 beta_swan_bal_cancer_sig_features$probe,
                                                 beta_swan_unbal_cancer_sig_features$probe,
                                                 num_vecs = 3)

# beta_quan_cancer_union
beta_quan_cancer_union_sig_features  <- getUnion(beta_quan_bal_counts_cancer_sig_features$probe,
                                                 beta_quan_bal_cancer_sig_features$probe,
                                                 beta_quan_unbal_cancer_sig_features$probe,
                                                 num_vecs = 3)

# beta_funnorm_cancer_union
beta_funnorm_cancer_union_sig_features  <- getUnion(beta_funnorm_bal_counts_cancer_sig_features$probe,
                                                    beta_funnorm_bal_cancer_sig_features$probe,
                                                    beta_funnorm_unbal_cancer_sig_features$probe,
                                                    num_vecs = 3)

# total cancer union
beta_cancer_union_sig_features <- getUnion(beta_raw_cancer_union_sig_features$probe_rgSet,
                                           beta_swan_cancer_union_sig_features$probe_rgSet, 
                                           beta_quan_cancer_union_sig_features$probe_rgSet,
                                           beta_funnorm_cancer_union_sig_features$probe_rgSet,
                                           num_vecs = 4)

# beta_bal_counts_cancer_union
beta_bal_counts_cancer_union_sig_features  <- getUnion(beta_raw_bal_counts_cancer_sig_features$probe,
                                                       beta_swan_bal_counts_cancer_sig_features$probe,
                                                       beta_quan_bal_counts_cancer_sig_features$probe,
                                                       beta_funnorm_bal_counts_cancer_sig_features$probe,
                                                       num_vecs = 4)

##########
# p53 union 
##########

# beta_raw_p53_union
beta_raw_p53_union_features  <- getUnion(beta_raw_bal_counts_p53_features$probe,
                                         beta_raw_bal_p53_features$probe,
                                         beta_raw_unbal_p53_features$probe,
                                         num_vecs = 3)

# beta_swan_p53_union
beta_swan_p53_union_features  <- getUnion(beta_swan_bal_counts_p53_features$probe,
                                          beta_swan_bal_p53_features$probe,
                                          beta_swan_unbal_p53_features$probe,
                                          num_vecs = 3)

# beta_quan_p53_union
beta_quan_p53_union_features  <- getUnion(beta_quan_bal_counts_p53_features$probe,
                                          beta_quan_bal_p53_features$probe,
                                          beta_quan_unbal_p53_features$probe,
                                          num_vecs = 3)

# beta_funnorm_p53_union
beta_funnorm_p53_union_features  <- getUnion(beta_funnorm_bal_counts_p53_features$probe,
                                             beta_funnorm_bal_p53_features$probe,
                                             beta_funnorm_unbal_p53_features$probe,
                                             num_vecs = 3)

# total p53 union
beta_p53_union_features <- getUnion(beta_raw_p53_union_features$probe_rgSet,
                                    beta_swan_p53_union_features$probe_rgSet, 
                                    beta_quan_p53_union_features$probe_rgSet,
                                    beta_funnorm_p53_union_features$probe_rgSet,
                                    num_vecs = 4)

# beta_bal_counts_p53_union
beta_bal_counts_p53_union_features  <- getUnion(beta_raw_bal_counts_p53_features$probe,
                                                beta_swan_bal_counts_p53_features$probe,
                                                beta_quan_bal_counts_p53_features$probe,
                                                beta_funnorm_bal_counts_p53_features$probe,
                                                num_vecs = 4)

##########
# p53 union sig
##########

# beta_raw_p53_union
beta_raw_p53_union_sig_features  <- getUnion(beta_raw_bal_counts_p53_sig_features$probe,
                                             beta_raw_bal_p53_sig_features$probe,
                                             beta_raw_unbal_p53_sig_features$probe,
                                             num_vecs = 3)

# beta_swan_p53_union
beta_swan_p53_union_sig_features  <- getUnion(beta_swan_bal_counts_p53_sig_features$probe,
                                              beta_swan_bal_p53_sig_features$probe,
                                              beta_swan_unbal_p53_sig_features$probe,
                                              num_vecs = 3)

# beta_quan_p53_union
beta_quan_p53_union_sig_features  <- getUnion(beta_quan_bal_counts_p53_sig_features$probe,
                                              beta_quan_bal_p53_sig_features$probe,
                                              beta_quan_unbal_p53_sig_features$probe,
                                              num_vecs = 3)

# beta_funnorm_p53_union
beta_funnorm_p53_union_sig_features  <- getUnion(beta_funnorm_bal_counts_p53_sig_features$probe,
                                                 beta_funnorm_bal_p53_sig_features$probe,
                                                 beta_funnorm_unbal_p53_sig_features$probe,
                                                 num_vecs = 3)

# total p53 union
beta_p53_union_sig_features <- getUnion(beta_raw_p53_union_sig_features$probe_rgSet,
                                        beta_swan_p53_union_sig_features$probe_rgSet, 
                                        beta_quan_p53_union_sig_features$probe_rgSet,
                                        beta_funnorm_p53_union_sig_features$probe_rgSet,
                                        num_vecs = 4)

# beta_bal_counts_p53_union
beta_bal_counts_p53_union_sig_features  <- getUnion(beta_raw_bal_counts_p53_sig_features$probe,
                                                    beta_swan_bal_counts_p53_sig_features$probe,
                                                    beta_quan_bal_counts_p53_sig_features$probe,
                                                    beta_funnorm_bal_counts_p53_sig_features$probe,
                                                    num_vecs = 4)

##########
# remove unneeded objects
##########
# keep only the objects that have features
rm(list = ls()[!grepl("features", ls())])

# save feaures
save.image(paste0(model_data, '/bh_features.RData'))


