#######################################################################################
# This script will run bumphunter on methylation probes 
library(bumphunter)
library(dplyr)
library(FDb.InfiniumMethylation.hg19)
library(impute)


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



# read in raw methylation
# Read in raw mehtylation data
num_sets <- 23 
methylation_raw <-  vector('list', num_sets)
chrome_vector <- vector('list', num_sets)

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
  chrome_vector[[i]] <- rep(paste0(data_name, i), ncol(methylation_raw[[i]]))
  print(i)
}

methylation <- do.call('cbind', methylation_raw)
# methylation <- as.data.frame(t(methylation), stringsAsFactors = FALSE)
# methylation <- cbind(x = rownames(methylation), methylation)
# colnames(methylation) <- methylation[1,]
# methylation <- methylation[-1,]
# names(methylation)[1] <- 'Probe'
# rownames(methylation) <- NULL  

names(methylation)[1] <- 'id'

temp <- unlist(chrome_vector)
temp.1 <- strsplit(temp, '/')
temp.2 <- lapply(temp.1, function(x) x[length(x)])
chr <- do.call('rbind', temp.2)
chr <- chr[-1]


cleanProbe <- function(data){
  
  col_names <- colnames(data)
  clean_index <- !grepl('ch|_|X', col_names)
  data <- data[,!grepl('ch|_|X', colnames(data))]
  data$id <- gsub('-', '_', data$id)
  rownames(data) <- NULL
  
  
  return(list(data, clean_index))
  
}

methylation_new <- cleanProbe(methylation)[[1]]
index <- cleanProbe(methylation)[[2]]
chr <- chr[index]
chr <- chr[-1]


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

#######################################################
# Transpose methylation and merge with clinical ID 
# Read in data (clinical or clinical_two)
clin <- read.csv(paste0(clin_data, '/clinical_two.csv'), stringsAsFactors = TRUE)

methylation_new <- methylation_new[!grepl('A|B|_', methylation_new$id),]

methylation_new$id <- as.factor(methylation_new$id)



clin$id <- as.factor(clin$blood_dna_malkin_lab_)

# merge clin and methylation
full_data <- inner_join(methylation_new, clin, by = 'id')

# get mut/wt vector
full_data <- full_data[!is.na(full_data$p53_germline),]
full_data$p53_germline <- factor(full_data$p53_germline, levels = c('WT', 'Mut'))


data <- as.data.frame(t(full_data), stringsAsFactors = F)
colnames(data) <- data[1,]
data <- data[-1,]
x <- cbind(rep(1, nrow(full_data)), full_data[212355])
x$p53_germline <- as.numeric(x$p53_germline)
x <- as.matrix(x)
pos <- seq(1, length(chr)*100, by = 10)[1:length(chr)]
pos <- as.matrix(pos)

betas <- data[1:length(chr),]

betas <- apply(betas, 2, as.numeric)

# save data 
DELTA_BETA_THRESH = 0.10 # DNAm difference threshold
NUM_BOOTSTRAPS = 3     # number of randomizations

# save.image(file = '/home/benbrew/Desktop/bumpdata.RData')
# load(file = '/home/benbrew/Desktop/bumpdata.RData')
# y = probes, x = cases and controls, chr = chromosome, pos = genomic locations , cl = cluster
tab <- bumphunter(as.matrix(betas), 
                  x, 
                  chr = chr, 
                  pos = pos,
                  nullMethod = "bootstrap",
                  cutoff = DELTA_BETA_THRESH,
                  B = NUM_BOOTSTRAPS,
                  type = "Beta")


# Yes, I’m using the R bumphunter package for DMRs and minfi package for general data processing. 
# Once you load the IDAT data into an object of the minfi's GenomicRatioSet class (called e.g. grSet) , 
# the main lines that do the bump analysis would look something like this: 

# ANDREI's example
DELTA_BETA_THRESH = 0.10 # DNAm difference threshold
NUM_BOOTSTRAPS = 100     # number of randomizations

betaMatrix <- getBeta(grSet)
annotation <- getAnnotation(grSet)
designMatrix <-  model.matrix(~ Sample_Group + Age + Gender, data=pData(grSet))

dmrs <- bumphunter(
  betaMatrix,
  design = designMatrix, 
  chr = annotation$chr, 
  pos = annotation$pos,
  nullMethod = "bootstrap",
  cutoff = DELTA_BETA_THRESH,
  B = NUM_BOOTSTRAPS,
  type = "Beta"
)


head(dmrs$table)

# Then just examine the results in the dmrs$table, post-process and filter them as you wish, etc.
# Let me know if this look clear.

# Example: We first generate an example of typical genomic locations

pos <- list(pos1=seq(1,1000,35),
            pos2=seq(2001,3000,35),
            pos3=seq(1,1000,50))
chr <- rep(paste0("chr",c(1,1,2)), times = sapply(pos,length))
pos <- unlist(pos, use.names=FALSE)

cl <- clusterMaker(chr, pos, maxGap = 300)
table(cl)

Indexes <- split(seq_along(cl), cl)
beta1 <- rep(0, length(pos))
for(i in seq(along=Indexes)){
  ind <- Indexes[[i]]
  x <- pos[ind]
  z <- scale(x, median(x), max(x)/12)
  beta1[ind] <- i*(-1)^(i+1)*pmax(1-abs(z)^3,0)^3 ##multiply by i to vary size
}

# Bumphunter is a more complicated function. In addition to regionFinder and clusterMaker
# it also implements a statistical model as well as permutation testing to assess uncertainty.
# We start by creating a simulated data set of 10 cases and 10 controls (recall that beta1 was defined above).


beta0 <- 3*sin(2*pi*pos/720)
X <- cbind(rep(1,20), rep(c(0,1), each=10))
error <- matrix(rnorm(20*length(beta1), 0, 1), ncol=20)
y <- t(X[,1])%x%beta0 + t(X[,2])%x%beta1 + error
# Now we can run bumphunter
# y = probes, x = cases and controls, chr = chromosome, pos = genomic locations , cl = cluster
tab <- bumphunter(y, X, chr, pos, cl, cutoff=.5)
