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


# read in different bumphunter results # 
bh_cancer <- read.csv(paste0(data_folder, '/bh_cancer.csv'))
bh_global <- read.csv(paste0(data_folder, '/bh_global.csv'))
bh_cancer_sub <- read.csv(paste0(data_folder, '/bh_cancer_sub.csv'))
bh_global_sub <- read.csv(paste0(data_folder, '/bh_global_sub.csv'))
bh_cancer_bal <- read.csv(paste0(data_folder, '/bh_cancer_balanced.csv'))
bh_global_bal <- read.csv(paste0(data_folder, '/bh_global_balanced.csv'))
bh_cancer_ind <-  read.csv(paste0(data_folder, '/bh_cancer_indicator.csv'))

# read in methylation data with probes
# methylation <- read.csv(paste0(data_folder, '/methyl_knn.csv'))
# methylation$X <- NULL
# load(paste0(data_folder, '/methyl_knn.RData'))
################################################
# get probes for bumphunter results
###############################################

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
    sub_rg <- rgSet[rgSet$seqnames ==i, ]
    
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
  
  ####### Now get nearest gene
  # First load the 450k data from hg19 (bioconductor, library is FDb.InfiniumMethylation.hg19)
  hm450 <- get450k()
  
  # Get probe names from our methylation data  
  probe_names <- as.character(dat_cg$probe)
  
  # remove probes that have less than 10 characters.
  probe_names <- probe_names[nchar(probe_names) == 10]
  
  # get probes from hm450
  probes <- hm450[probe_names]
  #get the nearest gene to each probe location.
  probe_info <- getNearestGene(probes)
  probe_info <- cbind(probe = rownames(probe_info), probe_info)
  rownames(probe_info) <- NULL
  
  # innerjoin dat_cg with probe_info by probe 
  dat_cg <- inner_join(probe_info, dat_cg, by = 'probe')
  
  # keep only necessary columns
  dat_cg <- dat_cg[, c('chr', 'probe', 'nearestGeneSymbol', 'p.value', 'fwer')]
  
  # order by p.value and save file 
  dat_cg <- dat_cg[order(dat_cg$p.value, decreasing = F),]
  
  
  return(dat_cg)
  
}

# cancer and global
bh_cancer_full <- getProbe(bh_cancer)
write.csv(bh_cancer_full, paste0(data_folder, '/bh_cancer_full.csv'))
bh_global_full <- getProbe(bh_global)
write.csv(bh_global_full, paste0(data_folder, '/bh_global_full.csv'))

# subselect cancer and global
bh_cancer_sub_full <- getProbe(bh_cancer_sub)
write.csv(bh_cancer_sub_full, paste0(data_folder, '/bh_cancer_sub_full.csv'))
bh_global_sub_full <- getProbe(bh_global_sub)
write.csv(bh_global_sub_full, paste0(data_folder, '/bh_global_sub_full.csv'))

# balanced cancer and global
bh_cancer_bal_full <- getProbe(bh_cancer_bal)
write.csv(bh_cancer_bal_full, paste0(data_folder, '/bh_cancer_bal_full.csv'))
bh_global_bal_full <- getProbe(bh_global_bal)
write.csv(bh_global_bal_full, paste0(data_folder, '/bh_global_bal_full.csv'))

# cancer as indicator instead of p53
bh_cancer_ind_full <- getProbe(bh_cancer_ind)
write.csv(bh_cancer_ind_full, paste0(data_folder, '/bh_cancer_ind_full.csv'))


# get union of all of them 
bh_union <- rbind(bh_cancer_full, 
                  bh_global_full, 
                  bh_cancer_sub_full, 
                  bh_global_sub_full, 
                  bh_cancer_bal_full, 
                  bh_global_bal_full)

# remove duplicates to get union 
bh_union <- bh_union[!duplicated(bh_union[,2]),]
write.csv(bh_union, paste0(data_folder, '/bh_union.csv'))


# get intersection of all of them 
intersection <- paste(Reduce(intersect, list(bh_cancer_full$probe, bh_global_full$probe, bh_cancer_sub_full$probe, bh_global_sub_full$probe,
                       bh_cancer_bal_full$probe, bh_global_bal_full$probe)),collapse = '|' )

# subset bh_union by these probes to get intersection
bh_intersection <- bh_union[grepl(intersection, bh_union$probe),]
write.csv(bh_intersection, paste0(data_folder, '/bh_intersection.csv'))


