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


# read in tab - results from bumphunter
tab <- read.csv(paste0(data_folder, '/bump_hunter_results.csv'))

# read in methylation data with probes
methylation <- read.csv(paste0(data_folder, '/methyl_knn.csv'))
methylation$X <- NULL

# save.image('/home/benbrew/Desktop/methy.RData')

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

# for each chromosome in rgSet find a matching chromosome in tab AND a start and finish that is within the range of tab star
# and finish

# first make seqnames in rgSet and chr in tab character vectors for identification
rgSet$seqnames <- as.character(rgSet$seqnames)
tab$chr <- as.character(tab$chr)

# loop through rgSet and tab and return a data frame with cg sites and their bumphunter significance (p.value and fwer)
results <- list()
results_data <- list()
# loop through chromosomes and match then loop through locations and match
for (i in unique(tab$chr)) {
  
  sub_tab <- tab[tab$chr == i,]
  sub_rg <- rgSet[rgSet$seqnames ==i, ]
  
  for (j in 1:nrow(sub_tab)) {
    
    chr_tab <- sub_tab[j,]
    
    if (any(sub_rg$start >= chr_tab$start & sub_rg$start <= chr_tab$end)){
     
       results[[j]] <- cbind(sub_rg[sub_rg$start >= chr_tab$start & sub_rg$start <= chr_tab$end,], chr_tab)
      
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

# remove duplicates from dat_cg
dat_cg <- dat_cg[!duplicated(dat_cg$probe),]

#################################################
# get nearest gene for dat_cg probes 
#################################################

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

# write.csv(dat_cg, paste0(data_folder, '/probe_gene_dmr.csv'))
