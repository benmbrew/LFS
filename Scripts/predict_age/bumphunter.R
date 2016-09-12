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

# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/analyze')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')

data_name <- '/Chr'

# results_file <- "file.txt"
# constructPath <- function(intermediate_folders, parent=results_folder,
#                           file=results_file) {
#   paste(parent, intermediate_folders, file, sep="/")
# }
# 
# incomplete_file <- constructPath("Incomplete")
# imputed_file <- constructPath("Imputed")
# jvmGBLimit <- 8

source(paste0(project_folder, '/Code/Functions/lsaImputation.R'))
source(paste0(project_folder, '/Code/Functions/knnImputation.R'))


################################################################
# read in and clean raw methylation. 
################################################################

# set empy vectors to fill with methylation
num_sets <- 23 
methylation_raw <-  vector('list', num_sets)

# read in raw methylation files. 
for(i in (1:num_sets)){
  methylation_raw[[i]] <- read.delim(paste(data_folder, data_name, i, '.txt', sep = ''), header = TRUE)
  sub_dat <- methylation_raw[[i]]
  column_split <- strsplit(as.character(sub_dat$X), '#')
  last_digits <- lapply(column_split, function(x) x[length(x)])
  sub_ids <- unlist(last_digits)
  sub_ids <- gsub('RD-', '', sub_ids)
  methylation_raw[[i]][,1] <- sub_ids
  if(i > 1){
    methylation_raw[[i]] <- methylation_raw[[i]][, -1]
  }
  print(i)
}

# concatenate methylation_raw
methylation <- do.call('cbind', methylation_raw)

# make id column for join later
names(methylation)[1] <- 'id'


# clean probe names 
cleanProbe <- function(data){
  
  col_names <- colnames(data)
  data <- data[,!grepl('ch|_|X', colnames(data))]
  data$id <- gsub('-', '_', data$id)
  rownames(data) <- NULL
  
  
  return(data)
  
}

# apply function and get new methylation. also return index for removed methylation
methylation_new <- cleanProbe(methylation)

###################################################################
# impute on methylation_new
###################################################################

# remove duplicated ids to put in rownames of methylation
methylation_new <- methylation_new[!duplicated(methylation_new$id),]
rownames(methylation_new) <- methylation_new$id
methylation_new$id <- NULL

#impute on methyl lsa
# methylation_new <- lsaImputation(incomplete_data = methylation_new, sample_rows = TRUE)
# impute on methyl with knn 
methylation_new <- knnImputation(methylation_new, sample_rows = TRUE)

# join rownames and methyl_impute and then erase rownames
methylation_new <- cbind(id = rownames(methylation_new), methylation_new)
rownames(methylation_new) <- NULL

# make methylation a data frame 
methylation_new <- as.data.frame(methylation_new)

###################################################################
# read in clinical data and merge with methylation_new
###################################################################

# read in clinical
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = TRUE)

# standardize ids
methylation_new <- methylation_new[!grepl('A|B|_', methylation_new$id),]
clin$id <- clin$blood_dna_malkin_lab_
clin <- clin[!grepl('A|B|_', clin$blood_dna_malkin_lab_),]

# inner join methylation and clinical
full_data <- inner_join(methylation_new, clin, by = 'id')

# get mut/wt vector
full_data <- full_data[!is.na(full_data$p53_germline),]
full_data$p53_germline <- factor(full_data$p53_germline, levels = c('WT', 'Mut'))

# get p53 and put into design matrix with intercept 1
p53_vecotr <- full_data$p53_germline
designMatrix <- cbind(rep(1, nrow(full_data)), p53_vecotr)
designMatrix <- as.matrix(designMatrix)

##########################################################
# Get genetic locations
##########################################################

# transpose methylation to join with cg_locations to get genetic location vector.
full_data <- as.data.frame(t(full_data), stringsAsFactors = F)

# make ids column names and remove first row
colnames(full_data) <- full_data[1,]
full_data <- full_data[-1,]

# make probe a column in methyl
full_data$probe <- rownames(full_data)
rownames(full_data) <- NULL

# remove clinical rows at bottom of data set
full_data <- full_data[1:212350,]

# read in cg_locations
cg_locations <- read.csv(paste0(data_folder, '/cg_locations.csv'))

# rename X column as probe for merge
names(cg_locations)[1] <- 'probe'
# inner join methyl and cg_locations by probe
methyl_cg <- inner_join(cg_locations, full_data, by = 'probe')

# get chr and pos vector 
chr <- methyl_cg$seqnames
pos <- methyl_cg$start

# create beta matrix
rownames(methyl_cg) <- methyl_cg$probe
beta <- methyl_cg[, 7:ncol(methyl_cg)]

# make beta numeric 
for (i in 1:ncol(beta)) {
  beta[,i] <- as.numeric(beta[,i])
  print(i)
} 

beta <- as.matrix(beta)


#########################################################################
# Run bumphunter
#########################################################################

# check dimensions 
stopifnot(dim(beta)[2] == dim(designMatrix)[1])
stopifnot(dim(beta)[1] == length(chr))
stopifnot(dim(beta)[1] == length(pos))

# set paramenters 
DELTA_BETA_THRESH = 0.10 # DNAm difference threshold
NUM_BOOTSTRAPS = 3     # number of randomizations

tab <- bumphunter(beta, 
                  designMatrix, 
                  chr = chr, 
                  pos = pos,
                  nullMethod = "bootstrap",
                  cutoff = DELTA_BETA_THRESH,
                  B = NUM_BOOTSTRAPS,
                  type = "Beta")

bump_hunter_results <- tab$table

# write.csv(bump_hunter_results, paste0(data_folder, '/bump_hunter_results.csv'))

# read in bumphunter results 
tab <- read.csv(paste0(data_folder, '/bump_hunter_results.csv'))

# get nearest gene to differentiated methylated regions and use that in model 


# ## Build GRanges with sequence lengths
# regions <- GRanges(seqnames = tab$chr, 
#                    IRanges(start = tab$start, end = tab$end),
#                    strand = '*', value = tab$value, area = tab$area, 
#                    cluster = tab$cluster, L = tab$L, clusterL = tab$clusterL)
# 
# ## Assign chr lengths
# data(hg19Ideogram, package = 'biovizBase')
# seqlengths(regions) <- seqlengths(hg19Ideogram)[names(seqlengths(regions))]
# ## Explore the regions
# regions
