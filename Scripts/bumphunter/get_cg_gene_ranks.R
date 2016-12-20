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
bumphunter_data <- paste0(data_folder, '/bumphunter_data')
methyl_data <- paste0(data_folder, '/methyl_data')
model_data <- paste0(data_folder, '/model_data')
imputed_data <- paste0(data_folder, '/imputed_data')
idat_data <- paste0(methyl_data, '/raw_files')
bumphunter_data <- paste0(data_folder, '/bumphunter_data')

#########
# load bump_hunter_lfs
#########
load(paste0(idat_data, '/beta_p53_bh.RData'))

#########
# load bump_hunter_cancer data
#########
load(paste0(idat_data, '/beta_cancer_bh.RData'))

##########
# get probes for bumphunter results
##########

getRgSet <- function() {
  
  #idat files
  idatFiles <- list.files("GSE68777/idat", pattern = "idat.gz$", full = TRUE)
  sapply(idatFiles, gunzip, overwrite = TRUE)
  
  # read into rgSet
  rgSet <- read.450k.exp("GSE68777/idat")
  
  # preprocess quantil
  rgSet <- preprocessQuantile(rgSet)
  
  # get rangers 
  rgSet <- granges(rgSet)
  rgSet <- as.data.frame(rgSet)
  return(rgSet)
}

rgSet <- getRgSet()

###########
# create function that grabs probe site and gene name for results from bumphunter
###########

results <- list()
results_data <- list()

getProbe <- function(data) {
  
  # first make seqnames in rgSet and chr in tab character vectors for identification
  rgSet$seqnames <- as.character(rgSet$seqnames)
  data$chr <- as.character(data$chr)

  # loop through chromosomes and match then loop through locations and match
  for (i in unique(data$chr)) {
    
    sub_data <- data[data$chr == i,]
    sub_rg <- rgSet[rgSet$seqnames == i, ]
    
    for (j in 1:nrow(sub_data)) {
      
      chr_data <- sub_data[j,]
      
      if (any(sub_rg$start >= chr_data$start & sub_rg$start <= chr_data$end)){
        
        results[[j]] <- cbind(sub_rg[sub_rg$start >= chr_data$start & sub_rg$start <= chr_data$end,], chr_data)
        
      }
      
      print(j)
      
    }
    
    results_data[[i]] <- do.call('rbind', results)
    
    print(i)
  }
  
  # combine into data frame 
  dat_cg <- as.data.frame(do.call('rbind', results_data))
  
  # clean up rownames and create column for cg site 
  dat_cg$probe <- rownames(dat_cg)
  rownames(dat_cg) <- NULL
  
  # remove chromosome info from probe colum
  dat_cg$probe <- substr(dat_cg$probe, nchar(dat_cg$probe) - 9, nchar(dat_cg$probe))
  
  # for those that have an extra 1 at the end and missing a c, paste a c back on and remove 1 
  dat_cg$probe <- ifelse(!grepl('c', dat_cg$probe), paste0('c', dat_cg$probe), dat_cg$probe)
  dat_cg$probe <- ifelse(nchar(dat_cg$probe) == 11, substr(dat_cg$probe, 1, 10), dat_cg$probe)
  
  # remove duplicates from dat_cg
  dat_cg <- dat_cg[!duplicated(dat_cg$probe),]
  
  # keep only necessary columns
  dat_cg <- dat_cg[, c('chr' , 'start','end', 'probe', 'p.value', 'fwer')]
  
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
   cg_locations, rgSet)

save.image(paste0(model_data, '/bh_features.RData'))

