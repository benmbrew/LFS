###########################################################################################################
# this script will get the cg probes and genes that bumphunter selected ordered by pvalue and FWER

#######################################################################################
# This script will run bumphunter on methylation probes 
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

# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
bumphunter_data <- paste0(data_folder, '/bumphunter_data')

# load data
load(paste0(bumphunter_data, '/bh_regions.RData'))

################################################
# get probes for bumphunter results
###############################################

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

###################################
# create function that grabs probe site and gene name for results from bumphunter
###################################
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
  
  # ####### Now get nearest gene
  # # First load the 450k data from hg19 (bioconductor, library is FDb.InfiniumMethylation.hg19)
  # hm450 <- get450k()
  # 
  # # Get probe names from our methylation data  
  # probe_names <- as.character(dat_cg$probe)
  # 
  # # remove probes that have less than 10 characters.
  # probe_names <- probe_names[nchar(probe_names) == 10]
  # 
  # # get probes from hm450
  # probes <- hm450[probe_names]
  # #get the nearest gene to each probe location.
  # probe_info <- getNearestGene(probes)
  # probe_info <- cbind(probe = rownames(probe_info), probe_info)
  # rownames(probe_info) <- NULL
  # 
  # # innerjoin dat_cg with probe_info by probe 
  # dat_cg <- inner_join(probe_info, dat_cg, by = 'probe')
  # 
  # # keep only necessary columns
  # dat_cg <- dat_cg[, c('chr', 'probe', 'p.value', 'fwer')]
  # 
  # # order by p.value and save file 
  # dat_cg <- dat_cg[order(dat_cg$p.value, decreasing = F),]
  # 
  # 
  return(dat_cg)
  
}

# apply function to four data bh (balanced)
bh_probe_knn_cancer_features <- getProbe(probe_knn_cancer_bh)
bh_probe_lsa_cancer_features <- getProbe(probe_lsa_cancer_bh)

bh_probe_knn_global_features <- getProbe(probe_knn_global_bh)
bh_probe_lsa_global_features <- getProbe(probe_lsa_global_bh)

# apply function to four data bh (unbalanced)
bh_probe_knn_cancer_unbal_features <- getProbe(probe_knn_cancer_unbal_bh)
bh_probe_lsa_cancer_unbal_features <- getProbe(probe_lsa_cancer_unbal_bh)

bh_probe_knn_global_unbal_features <- getProbe(probe_knn_global_unbal_bh)
bh_probe_lsa_global_unbal_features <- getProbe(probe_lsa_global_unbal_bh)

rm(probe_knn_cancer_bh, probe_lsa_cancer_bh, probe_knn_global_bh, probe_lsa_global_bh, 
   probe_knn_cancer_unbal_bh, probe_lsa_cancer_unbal_bh, probe_knn_global_unbal_bh, probe_lsa_global_unbal_bh,
   cg_locations, rgSet)

# get union of all of them 
bh_union_bal <- rbind(bh_probe_knn_cancer_features,
                  bh_probe_lsa_cancer_features,
                  bh_probe_knn_global_features,
                  bh_probe_lsa_global_features)

bh_union_unbal <- rbind(bh_probe_knn_cancer_unbal_features,
                        bh_probe_lsa_cancer_unbal_features,
                        bh_probe_knn_global_unbal_features,
                        bh_probe_lsa_global_unbal_features)
                  

# remove duplicates to get union 
bh_union_bal <- bh_union_bal[!duplicated(bh_union_bal[,2]),]
bh_union_unbal <- bh_union_unbal[!duplicated(bh_union_unbal[,2]),]


# get intersection of all of them 
intersection_bal <- paste(Reduce(intersect, list(bh_probe_knn_cancer_features$probe, bh_probe_lsa_cancer_features$probe,
                                             bh_probe_knn_global_features$probe, bh_probe_lsa_global_features$probe)),collapse = '|' )

# subset bh_union by these probes to get intersection
bh_intersection_bal <- bh_union_bal[grepl(intersection_bal, bh_union_bal$probe),]

# get intersection of all of them 
intersection_unbal <- paste(Reduce(intersect, list(bh_probe_knn_cancer_unbal_features$probe, bh_probe_lsa_cancer_unbal_features$probe,
                                                 bh_probe_knn_global_unbal_features$probe, bh_probe_lsa_global_unbal_features$probe)),collapse = '|' )

# subset bh_union by these probes to get intersection
bh_intersection_unbal <- bh_union_unbal[grepl(intersection_unbal, bh_union_unbal$probe),]

# save data to modeldata
save.image(paste0(model_data, '/bh_features.RData'))
