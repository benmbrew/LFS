####### Script will map regions to probes and save the final bh features for both cancer and lfs
# this is 6th step in pipeline

##########
# initialize libraries
##########
library(minfi)
library(dplyr)
library(GenomicRanges)

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
load(paste0(model_data, '/beta_p53_bh.RData'))

##########
# load bump_hunter_cancer data
##########
load(paste0(model_data, '/beta_cancer_bh.RData'))

##########
# change column names of cg_locations for merging
##########
cg_locations$X <- NULL
colnames(cg_locations) <- paste0(colnames(cg_locations), '_', 'rgSet')


##########
# get intersections
##########

getIntersection <- function(vec1 = NULL, 
                            vec2 = NULL, 
                            vec3 = NULL, 
                            vec4 = NULL,
                            num_vecs) 
{
  if(num_vecs == 2) {
    probe_list <- Reduce(intersect, list(vec1,vec2))
    
  } 
  
  if(num_vecs == 3) {
    probe_list <- Reduce(intersect, list(vec1,vec2,vec3))
    
  } 
  
  if (num_vecs == 4) {
    probe_list <- Reduce(intersect, list(vec1,vec2,vec3,vec4))
    
  }
  combined_probes <- paste(probe_list, collapse = '|')
  dat <- cg_locations[grepl(combined_probes, cg_locations$probe_rgSet),]
  stopifnot(length(probe_list) == nrow(dat))
  return(dat)
}

##########
# get union
##########
getUnion <- function(vec1 = NULL, 
                     vec2 = NULL, 
                     vec3 = NULL, 
                     vec4 = NULL,
                     num_vecs) 
{
  
  if(num_vecs == 2) {
    probe_list <- Reduce(union, list(vec1,vec2))
    
  } 
  
  if(num_vecs == 3) {
    probe_list <- Reduce(union, list(vec1,vec2,vec3))
    
  } 
  
  if (num_vecs == 4) {
    probe_list <- Reduce(union, list(vec1,vec2,vec3,vec4))
    
  }
  combined_probes <- paste(probe_list, collapse = '|')
  dat <- cg_locations[grepl(combined_probes, cg_locations$probe_rgSet),]
  stopifnot(length(probe_list) == nrow(dat))
  return(dat)
}

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
  
  # loop through chromosomes and match then loop through locations and match
  for (i in unique(data$chr)) {
    
    sub_data <- data[data$chr == i,] # some region length >0
    sub_rg <- cg_locations[cg_locations$seqnames_rgSet == i, ] # no region length >0
    
    for (j in 1:nrow(sub_data)) {
      
      chr_data <- sub_data[j,]
      
      if (any(sub_rg$start_rgSet >= chr_data$start & sub_rg$start_rgSet <= chr_data$end)) {
        
        results[[j]] <- cbind(sub_rg[sub_rg$start_rgSet >= chr_data$start & sub_rg$start_rgSet <= chr_data$end,], chr_data)
        
      }
      
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
  
  dat_cg_sig <- dat_cg[dat_cg$p.value < 0.05,]
  
  return(list(dat_cg, dat_cg_sig))
  
}

##########
# apply function to cancer bh for both sig and all on cancer
##########

# beta raw all
beta_raw_bal_counts_cancer_features <- getProbe(beta_raw_bal_counts_cancer)[[1]]
beta_raw_bal_cancer_features <- getProbe(beta_raw_bal_cancer)[[1]]
beta_raw_unbal_cancer_features <- getProbe(beta_raw_unbal_cancer)[[1]]

# beta swan all
beta_swan_bal_counts_cancer_features <- getProbe(beta_swan_bal_counts_cancer)[[1]]
beta_swan_bal_cancer_features <- getProbe(beta_swan_bal_cancer)[[1]]
beta_swan_unbal_cancer_features <- getProbe(beta_swan_unbal_cancer)[[1]]

# beta quan all
beta_quan_bal_counts_cancer_features <- getProbe(beta_quan_bal_counts_cancer)[[1]]
beta_quan_bal_cancer_features <- getProbe(beta_quan_bal_cancer)[[1]]
beta_quan_unbal_cancer_features <- getProbe(beta_quan_unbal_cancer)[[1]]

# beta funnorm all
beta_funnorm_bal_counts_cancer_features <- getProbe(beta_funnorm_bal_counts_cancer)[[1]]
beta_funnorm_bal_cancer_features <- getProbe(beta_funnorm_bal_cancer)[[1]]
beta_funnorm_unbal_cancer_features <- getProbe(beta_funnorm_unbal_cancer)[[1]]

# beta raw sig
beta_raw_bal_counts_cancer_sig_features <- getProbe(beta_raw_bal_counts_cancer)[[2]]
beta_raw_bal_cancer_sig_features <- getProbe(beta_raw_bal_cancer)[[2]]
beta_raw_unbal_cancer_sig_features <- getProbe(beta_raw_unbal_cancer)[[2]]

# beta swan sig
beta_swan_bal_counts_cancer_sig_features <- getProbe(beta_swan_bal_counts_cancer)[[2]]
beta_swan_bal_cancer_sig_features <- getProbe(beta_swan_bal_cancer)[[2]]
beta_swan_unbal_cancer_sig_features <- getProbe(beta_swan_unbal_cancer)[[2]]

# beta quan sig
beta_quan_bal_counts_cancer_features <- getProbe(beta_quan_bal_counts_cancer)[[2]]
beta_quan_bal_cancer_features <- getProbe(beta_quan_bal_cancer)[[2]]
beta_quan_unbal_cancer_features <- getProbe(beta_quan_unbal_cancer)[[2]]

# beta funnorm sig
beta_funnorm_bal_counts_cancer_sig_features <- getProbe(beta_funnorm_bal_counts_cancer)[[2]]
beta_funnorm_bal_cancer_sig_features <- getProbe(beta_funnorm_bal_cancer)[[2]]
beta_funnorm_unbal_cancer_sig_features <- getProbe(beta_funnorm_unbal_cancer)[[2]]

##########
# apply function to p53 bh for both sig and all on p53
##########

# beta raw all
beta_raw_bal_counts_p53_features <- getProbe(beta_raw_bal_counts_p53)[[1]]
beta_raw_bal_p53_features <- getProbe(beta_raw_bal_p53)[[1]]
beta_raw_unbal_p53_features <- getProbe(beta_raw_unbal_p53)[[1]]

# beta swan all
beta_swan_bal_counts_p53_features <- getProbe(beta_swan_bal_counts_p53)[[1]]
beta_swan_bal_p53_features <- getProbe(beta_swan_bal_p53)[[1]]
beta_swan_unbal_p53_features <- getProbe(beta_swan_unbal_p53)[[1]]

# beta quan all
beta_quan_bal_counts_p53_features <- getProbe(beta_quan_bal_counts_p53)[[1]]
beta_quan_bal_p53_features <- getProbe(beta_quan_bal_p53)[[1]]
beta_quan_unbal_p53_features <- getProbe(beta_quan_unbal_p53)[[1]]

# beta funnorm all
beta_funnorm_bal_counts_p53_features <- getProbe(beta_funnorm_bal_counts_p53)[[1]]
beta_funnorm_bal_p53_features <- getProbe(beta_funnorm_bal_p53)[[1]]
beta_funnorm_unbal_p53_features <- getProbe(beta_funnorm_unbal_p53)[[1]]

# beta raw sig
beta_raw_bal_counts_p53_sig_features <- getProbe(beta_raw_bal_counts_p53)[[2]]
beta_raw_bal_p53_sig_features <- getProbe(beta_raw_bal_p53)[[2]]
beta_raw_unbal_p53_sig_features <- getProbe(beta_raw_unbal_p53)[[2]]

# beta swan sig
beta_swan_bal_counts_p53_sig_features <- getProbe(beta_swan_bal_counts_p53)[[2]]
beta_swan_bal_p53_sig_features <- getProbe(beta_swan_bal_p53)[[2]]
beta_swan_unbal_p53_sig_features <- getProbe(beta_swan_unbal_p53)[[2]]

# beta quan sig
beta_quan_bal_counts_p53_features <- getProbe(beta_quan_bal_counts_p53)[[2]]
beta_quan_bal_p53_features <- getProbe(beta_quan_bal_p53)[[2]]
beta_quan_unbal_p53_features <- getProbe(beta_quan_unbal_p53)[[2]]

# beta funnorm sig
beta_funnorm_bal_counts_p53_sig_features <- getProbe(beta_funnorm_bal_counts_p53)[[2]]
beta_funnorm_bal_p53_sig_features <- getProbe(beta_funnorm_bal_p53)[[2]]
beta_funnorm_unbal_p53_sig_features <- getProbe(beta_funnorm_unbal_p53)[[2]]

save.image('/home/benbrew/Desktop/temp_bh.RData')
load('/home/benbrew/Desktop/temp_bh.RData')


###########################################################################################################################
# this part of script will get intersection and unions of various sets of probe features

# INTERSECTION

##########
# cancer intersection 
##########

# beta_raw_cancer_intersection
raw_cancer_intersection  <- getIntersection(beta_raw_bal_counts_cancer_features$probe,
                                            beta_raw_bal_cancer_features$probe,
                                            beta_raw_unbal_cancer_features$probe,
                                            num_vecs = 3)

# beta_swan_cancer_intersection
swan_cancer_intersection  <- getIntersection(beta_swan_bal_counts_cancer_features$probe,
                                             beta_swan_bal_cancer_features$probe,
                                             beta_swan_unbal_cancer_features$probe,
                                             num_vecs = 3)

# beta_quan_cancer_intersection
quan_cancer_intersection  <- getIntersection(beta_quan_bal_counts_cancer_features$probe,
                                      beta_quan_bal_cancer_features$probe,
                                      beta_quan_unbal_cancer_features$probe,
                                      num_vecs = 3)

# beta_funnorm_cancer_intersection
funnorm_cancer_intersection  <- getIntersection(beta_funnorm_bal_counts_cancer_features$probe,
                                                beta_funnorm_bal_cancer_features$probe,
                                                beta_funnorm_unbal_cancer_features$probe,
                                                num_vecs = 3)

# total cancer intersection
cancer_intersection <- getIntersection(beta_cancer_intersection,
                                  swan_cancer_intersection, 
                                  quan_cancer_intersection,
                                  funnorm_cancer_intersection,
                                  num_vecs = 4)

# beta_bal_counts_cancer_intersection
bal_counts_cancer_intersection  <- getIntersection(beta_raw_bal_counts_cancer_features$probe,
                                                   beta_swan_bal_counts_cancer_features$probe,
                                                   beta_quan_bal_counts_cancer_features$probe,
                                                   beta_funnorm_bal_counts_cancer_features$probe,
                                                   num_vecs = 4)

##########
# cancer intersection sig
##########

# beta_raw_cancer_intersection
raw_cancer_intersection_sig  <- getIntersection(beta_raw_bal_counts_cancer_sig_features$probe,
                                            beta_raw_bal_cancer_sig_features$probe,
                                            beta_raw_unbal_cancer_sig_features$probe,
                                            num_vecs = 3)

# beta_swan_cancer_intersection
swan_cancer_intersection_sig  <- getIntersection(beta_swan_bal_counts_cancer_sig_features$probe,
                                             beta_swan_bal_cancer_sig_features$probe,
                                             beta_swan_unbal_cancer_sig_features$probe,
                                             num_vecs = 3)

# beta_quan_cancer_intersection
quan_cancer_intersection_sig  <- getIntersection(beta_quan_bal_counts_cancer_sig_features$probe,
                                      beta_quan_bal_cancer_sig_features$probe,
                                      beta_quan_unbal_cancer_sig_features$probe,
                                      num_vecs = 3)

# beta_funnorm_cancer_intersection
funnorm_cancer_intersection_sig  <- getIntersection(beta_funnorm_bal_counts_cancer_sig_features$probe,
                                                beta_funnorm_bal_cancer_sig_features$probe,
                                                beta_funnorm_unbal_cancer_sig_features$probe,
                                                num_vecs = 3)

# total cancer intersection
cancer_intersection_sig <- getIntersection(beta_cancer_intersection_sig,
                                  swan_cancer_intersection_sig, 
                                  quan_cancer_intersection_sig,
                                  funnorm_cancer_intersection_sig,
                                  num_vecs = 4)

# beta_bal_counts_cancer_intersection
bal_counts_cancer_intersection_sig  <- getIntersection(beta_raw_bal_counts_cancer_sig_features$probe,
                                                   beta_swan_bal_counts_cancer_sig_features$probe,
                                                   beta_quan_bal_counts_cancer_sig_features$probe,
                                                   beta_funnorm_bal_counts_cancer_sig_features$probe,
                                                   num_vecs = 4)



##########
# p53 intersection 
##########

# beta_raw_p53_intersection
raw_p53_intersection  <- getIntersection(beta_raw_bal_counts_p53_features$probe,
                                            beta_raw_bal_p53_features$probe,
                                            beta_raw_unbal_p53_features$probe,
                                            num_vecs = 3)

# beta_swan_p53_intersection
swan_p53_intersection  <- getIntersection(beta_swan_bal_counts_p53_features$probe,
                                             beta_swan_bal_p53_features$probe,
                                             beta_swan_unbal_p53_features$probe,
                                             num_vecs = 3)

# beta_quan_p53_intersection
quan_p53_intersection  <- getIntersection(beta_quan_bal_counts_p53_features$probe,
                                             beta_quan_bal_p53_features$probe,
                                             beta_quan_unbal_p53_features$probe,
                                             num_vecs = 3)

# beta_funnorm_p53_intersection
funnorm_p53_intersection  <- getIntersection(beta_funnorm_bal_counts_p53_features$probe,
                                                beta_funnorm_bal_p53_features$probe,
                                                beta_funnorm_unbal_p53_features$probe,
                                                num_vecs = 3)

# total p53 intersection
p53_intersection <- getIntersection(beta_p53_intersection,
                                  swan_p53_intersection, 
                                  quan_p53_intersection,
                                  funnorm_p53_intersection,
                                  num_vecs = 4)

# beta_bal_counts_p53_intersection
bal_counts_p53_intersection  <- getIntersection(beta_raw_bal_counts_p53_features$probe,
                                                   beta_swan_bal_counts_p53_features$probe,
                                                   beta_quan_bal_counts_p53_features$probe,
                                                   beta_funnorm_bal_counts_p53_features$probe,
                                                   num_vecs = 4)

##########
# p53 intersection sig
##########

# beta_raw_p53_intersection
raw_p53_intersection_sig  <- getIntersection(beta_raw_bal_counts_p53_sig_features$probe,
                                                beta_raw_bal_p53_sig_features$probe,
                                                beta_raw_unbal_p53_sig_features$probe,
                                                num_vecs = 3)

# beta_swan_p53_intersection
swan_p53_intersection_sig  <- getIntersection(beta_swan_bal_counts_p53_sig_features$probe,
                                                 beta_swan_bal_p53_sig_features$probe,
                                                 beta_swan_unbal_p53_sig_features$probe,
                                                 num_vecs = 3)

# beta_quan_p53_intersection
quan_p53_intersection_sig  <- getIntersection(beta_quan_bal_counts_p53_sig_features$probe,
                                                 beta_quan_bal_p53_sig_features$probe,
                                                 beta_quan_unbal_p53_sig_features$probe,
                                                 num_vecs = 3)

# beta_funnorm_p53_intersection
funnorm_p53_intersection_sig  <- getIntersection(beta_funnorm_bal_counts_p53_sig_features$probe,
                                                    beta_funnorm_bal_p53_sig_features$probe,
                                                    beta_funnorm_unbal_p53_sig_features$probe,
                                                    num_vecs = 3)

# total p53 intersection
p53_intersection_sig <- getIntersection(beta_p53_intersection_sig,
                                      swan_p53_intersection_sig, 
                                      quan_p53_intersection_sig,
                                      funnorm_p53_intersection_sig,
                                      num_vecs = 4)

# beta_bal_counts_p53_intersection
bal_counts_p53_intersection_sig  <- getIntersection(beta_raw_bal_counts_p53_sig_features$probe,
                                                       beta_swan_bal_counts_p53_sig_features$probe,
                                                       beta_quan_bal_counts_p53_sig_features$probe,
                                                       beta_funnorm_bal_counts_p53_sig_features$probe,
                                                       num_vecs = 4)

# UNION

##########
# cancer union 
##########

# beta_raw_cancer_union
raw_cancer_union  <- getUnion(beta_raw_bal_counts_cancer_features$probe,
                                            beta_raw_bal_cancer_features$probe,
                                            beta_raw_unbal_cancer_features$probe,
                                            num_vecs = 3)

# beta_swan_cancer_union
swan_cancer_union  <- getUnion(beta_swan_bal_counts_cancer_features$probe,
                                             beta_swan_bal_cancer_features$probe,
                                             beta_swan_unbal_cancer_features$probe,
                                             num_vecs = 3)

# beta_quan_cancer_union
quan_cancer_union  <- getUnion(beta_quan_bal_counts_cancer_features$probe,
                                             beta_quan_bal_cancer_features$probe,
                                             beta_quan_unbal_cancer_features$probe,
                                             num_vecs = 3)

# beta_funnorm_cancer_union
funnorm_cancer_union  <- getUnion(beta_funnorm_bal_counts_cancer_features$probe,
                                                beta_funnorm_bal_cancer_features$probe,
                                                beta_funnorm_unbal_cancer_features$probe,
                                                num_vecs = 3)

# total cancer union
cancer_union <- getUnion(beta_cancer_union,
                                  swan_cancer_union, 
                                  quan_cancer_union,
                                  funnorm_cancer_union,
                                  num_vecs = 4)

# beta_bal_counts_cancer_union
bal_counts_cancer_union  <- getUnion(beta_raw_bal_counts_cancer_features$probe,
                                                   beta_swan_bal_counts_cancer_features$probe,
                                                   beta_quan_bal_counts_cancer_features$probe,
                                                   beta_funnorm_bal_counts_cancer_features$probe,
                                                   num_vecs = 4)

##########
# cancer union sig
##########

# beta_raw_cancer_union
raw_cancer_union_sig  <- getUnion(beta_raw_bal_counts_cancer_sig_features$probe,
                                                beta_raw_bal_cancer_sig_features$probe,
                                                beta_raw_unbal_cancer_sig_features$probe,
                                                num_vecs = 3)

# beta_swan_cancer_union
swan_cancer_union_sig  <- getUnion(beta_swan_bal_counts_cancer_sig_features$probe,
                                                 beta_swan_bal_cancer_sig_features$probe,
                                                 beta_swan_unbal_cancer_sig_features$probe,
                                                 num_vecs = 3)

# beta_quan_cancer_union
quan_cancer_union_sig  <- getUnion(beta_quan_bal_counts_cancer_sig_features$probe,
                                                 beta_quan_bal_cancer_sig_features$probe,
                                                 beta_quan_unbal_cancer_sig_features$probe,
                                                 num_vecs = 3)

# beta_funnorm_cancer_union
funnorm_cancer_union_sig  <- getUnion(beta_funnorm_bal_counts_cancer_sig_features$probe,
                                                    beta_funnorm_bal_cancer_sig_features$probe,
                                                    beta_funnorm_unbal_cancer_sig_features$probe,
                                                    num_vecs = 3)

# total cancer union
cancer_union_sig <- getUnion(beta_cancer_union_sig,
                             swan_cancer_union_sig, 
                             quan_cancer_union_sig,
                             funnorm_cancer_union_sig,
                             num_vecs = 4)

# beta_bal_counts_cancer_union
bal_counts_cancer_union_sig  <- getUnion(beta_raw_bal_counts_cancer_sig_features$probe,
                                         beta_swan_bal_counts_cancer_sig_features$probe,
                                         beta_quan_bal_counts_cancer_sig_features$probe,
                                         beta_funnorm_bal_counts_cancer_sig_features$probe,
                                                       num_vecs = 4)



##########
# p53 union 
##########

# beta_raw_p53_union
raw_p53_union  <- getUnion(beta_raw_bal_counts_p53_features$probe,
                                         beta_raw_bal_p53_features$probe,
                                         beta_raw_unbal_p53_features$probe,
                                         num_vecs = 3)

# beta_swan_p53_union
swan_p53_union  <- getUnion(beta_swan_bal_counts_p53_features$probe,
                                          beta_swan_bal_p53_features$probe,
                                          beta_swan_unbal_p53_features$probe,
                                          num_vecs = 3)

# beta_quan_p53_union
quan_p53_union  <- getUnion(beta_quan_bal_counts_p53_features$probe,
                                          beta_quan_bal_p53_features$probe,
                                          beta_quan_unbal_p53_features$probe,
                                          num_vecs = 3)

# beta_funnorm_p53_union
funnorm_p53_union  <- getUnion(beta_funnorm_bal_counts_p53_features$probe,
                                             beta_funnorm_bal_p53_features$probe,
                                             beta_funnorm_unbal_p53_features$probe,
                                             num_vecs = 3)

# total p53 union
p53_union <- getUnion(beta_p53_union,
                               swan_p53_union, 
                               quan_p53_union,
                               funnorm_p53_union,
                               num_vecs = 4)

# beta_bal_counts_p53_union
bal_counts_p53_union  <- getUnion(beta_raw_bal_counts_p53_features$probe,
                                                beta_swan_bal_counts_p53_features$probe,
                                                beta_quan_bal_counts_p53_features$probe,
                                                beta_funnorm_bal_counts_p53_features$probe,
                                                num_vecs = 4)

##########
# p53 union sig
##########

# beta_raw_p53_union
raw_p53_union_sig  <- getUnion(beta_raw_bal_counts_p53_sig_features$probe,
                                             beta_raw_bal_p53_sig_features$probe,
                                             beta_raw_unbal_p53_sig_features$probe,
                                             num_vecs = 3)

# beta_swan_p53_union
swan_p53_union_sig  <- getUnion(beta_swan_bal_counts_p53_sig_features$probe,
                                              beta_swan_bal_p53_sig_features$probe,
                                              beta_swan_unbal_p53_sig_features$probe,
                                              num_vecs = 3)

# beta_quan_p53_union
quan_p53_union_sig  <- getUnion(beta_quan_bal_counts_p53_sig_features$probe,
                                              beta_quan_bal_p53_sig_features$probe,
                                              beta_quan_unbal_p53_sig_features$probe,
                                              num_vecs = 3)

# beta_funnorm_p53_union
funnorm_p53_union_sig  <- getUnion(beta_funnorm_bal_counts_p53_sig_features$probe,
                                                 beta_funnorm_bal_p53_sig_features$probe,
                                                 beta_funnorm_unbal_p53_sig_features$probe,
                                                 num_vecs = 3)

# total p53 union
p53_union_sig <- getUnion(beta_p53_union_sig,
                                   swan_p53_union_sig, 
                                   quan_p53_union_sig,
                                   funnorm_p53_union_sig,
                                   num_vecs = 4)

# beta_bal_counts_p53_union
bal_counts_p53_union_sig  <- getUnion(beta_raw_bal_counts_p53_sig_features$probe,
                                                    beta_swan_bal_counts_p53_sig_features$probe,
                                                    beta_quan_bal_counts_p53_sig_features$probe,
                                                    beta_funnorm_bal_counts_p53_sig_features$probe,
                                                    num_vecs = 4)


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

