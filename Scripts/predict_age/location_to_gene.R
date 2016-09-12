#######################################################################################
# Map chromosome and regions to nearest gene and save gene names

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
test <- paste0(project_folder, '/Scripts/analyze')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')


# read in tab - results from bumphunter
tab <- read.csv(paste0(data_folder, '/bump_hunter_results.csv'))

# read in methylation data with probes
methylation <- read.csv(paste0(methyl_data, '/methylation.csv'), header = TRUE, check.names = FALSE)

# # load minfi data
# getGEOSuppFiles("GSE68777")
# untar("GSE68777/GSE68777_RAW.tar", exdir = "GSE68777/idat")
# head(list.files("GSE68777/idat", pattern = "idat"))

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

# creat an empty vector to store results
common_loci <- list()

# loop through chromosomes and match then loop through locations and match
for (i in unique(rgSet$seqnames)) {
  
  sub_dat <- rgSet[rgSet$seqnames == i,]
  
  for (j in 1:nrow(tab)) {
    
    sub_tab <- tab[j,]
    
    if(any(sub_dat$seqnames %in% sub_tab$chr & sub_dat$start >= sub_tab$start & sub_dat$start <= sub_tab$end)){
      
      common_chr <- sub_dat[sub_dat$seqnames %in% sub_tab$chr,]
      common_loci[[j]] <- common_chr[common_chr$start >= sub_tab$start & common_chr$start <= sub_tab$end,]
    }
    print(j)
  }
  print(i)
}

# row bind the list and remove repeats
cg_site <- do.call(rbind, common_loci)

# put rownames in first column
cg_site$Probe <- rownames(cg_site)
rownames(cg_site) <- NULL

# innerjoin cg_site and methylation by probe
methyl_dmr <- inner_join(cg_site, methylation, by = 'Probe')

# remove first 5 columns and save data 
methyl_dmr <- methyl_dmr[, 6:ncol(methyl_dmr)]
# write.csv(methyl_dmr, file = paste0(data_folder, '/methyl_dmr.csv'))

# map methyl_dmr probes to nearest gene




