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
data_folder <- paste0(project_folder, '/Data')
clin_data <- paste0(data_folder, '/clin_data')


# read in tab - results from bumphunter
tab <- read.csv(paste0(data_folder, '/bump_hunter_results.csv'))

# read in methylation data with probes
methylation <- read.csv(paste0(data_folder, '/methyl_knn.csv'))
methylation$X <- NULL

# read in clinical data 
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = F)

# save.image('/home/benbrew/Desktop/big_methyl.csv')
# load('/home/benbrew/Desktop/big_methyl.csv')

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
    
    if (any(sub_dat$seqnames %in% sub_tab$chr & sub_dat$start >= sub_tab$start & sub_dat$start <= sub_tab$end)) {
      
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

# transpose methylation 
methylation <- as.data.frame(t(methylation), stringsAsFactors = F)
colnames(methylation) <- methylation[1,]
methylation <- methylation[-1,]
methylation$Probe <- rownames(methylation)

# innerjoin cg_site and methylation by probe
methyl_dmr <- inner_join(cg_site, methylation, by = 'Probe')

# remove first 5 columns and save data 
methyl_dmr <- methyl_dmr[, 6:ncol(methyl_dmr)]

# transpose matrix and save 
methyl_dmr <- as.data.frame(t(methyl_dmr), stringsAsFactors = F)
colnames(methyl_dmr) <- methyl_dmr[1,]
methyl_dmr <- methyl_dmr[-1,]
methyl_dmr$id <- rownames(methyl_dmr)
rownames(methyl_dmr) <- NULL

# standardize ids
methyl_dmr <- methyl_dmr[!grepl('A|B|_', methyl_dmr$id),]
clin$id <- clin$blood_dna_malkin_lab_
clin <- clin[!grepl('A|B|_', clin$blood_dna_malkin_lab_),]

methyl_dmr <- inner_join(clin, methyl_dmr, by = 'id')

# write.csv(methyl_dmr, file = paste0(data_folder, '/methyl_dmr.csv'))

#######################################################################
# map methyl_dmr probes to nearest gene
#######################################################################

# First load the 450k data from hg19 (bioconductor, library is FDb.InfiniumMethylation.hg19)
hm450 <- get450k()

# Get probe names from our methylation data  
probe_names <- as.character(methylation$Probe)

# remove probes that have less than 10 characters.
probe_names <- probe_names[nchar(probe_names) == 10]

# get probes from hm450
probes <- hm450[probe_names]

#get the nearest gene to each probe location.
probe_info <- getNearestGene(probes)
probe_info <- cbind(probe = rownames(probe_info), probe_info)
rownames(probe_info) <- NULL

# join probe_info with methylation. This keeps all of the probes that we could match in hm450 and drops the others.
methyl_gene_dmr <- inner_join(cg_site, methylation, by = c('Probe'= 'Probe'))
methyl_gene_dmr <- inner_join(probe_info, methyl_gene_dmr, by = c('probe'= 'Probe'))

# Get rid of extra variables.
methyl_gene_dmr$probe <- NULL
methyl_gene_dmr$queryHits <- NULL
methyl_gene_dmr$subjectHits <- NULL
methyl_gene_dmr$distance<- NULL
methyl_gene_dmr$seqnames <- NULL
methyl_gene_dmr$start <- NULL
methyl_gene_dmr$end <- NULL
methyl_gene_dmr$width <- NULL
methyl_gene_dmr$strand <- NULL



# transpose matrix and save 
methyl_gene_dmr <- as.data.frame(t(methyl_gene_dmr), stringsAsFactors = F)
colnames(methyl_gene_dmr) <- methyl_gene_dmr[1,]
methyl_gene_dmr <- methyl_gene_dmr[-1,]
methyl_gene_dmr$id <- rownames(methyl_gene_dmr)
rownames(methyl_gene_dmr) <- NULL

methyl_gene_dmr <- inner_join(clin, methyl_gene_dmr, by = 'id')

# save data 
# write.csv(methyl_gene_dmr, paste0(data_folder, '/methyl_gene_dmr.csv'))




