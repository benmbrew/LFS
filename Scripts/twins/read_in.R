#Twins data

##########
# initialize libraries
##########
library(dplyr)
library(sva)
library(minfi)
library(methyAnalysis)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Scripts/twins/data')
cg_folder <- paste0(project_folder, '/Scripts/twins')

##########
# load data batch data
##########
rgSet <- read.metharray.exp(paste0(data_folder))

# get cg_locations
cg_locations <- read.csv(paste0(cg_folder, '/cg_locations.csv'))

# rename to probe
names(cg_locations)[1] <- 'probe'

##########
# normalize data
##########

getNorm <- function(data, type) {
  
  if (type == 'raw') {
    ratioSet <- preprocessRaw(data)
  }
  
  if (type == 'quan') {
    ratioSet <- preprocessQuantile(data)
  }
  
  
  gset <- mapToGenome(ratioSet) 
  beta <- getBeta(gset)
}

raw <- getNorm(rgSet, type = 'raw')

# get complete cases
raw <- raw[complete.cases(raw),]


##########
# transpose data and join
# R03C01 has cancer,  RO2C01 is unnaffected
##########

raw <- as.data.frame(t(raw), stringsAsFactors = F)


##########
# explore methylation with methyAnalysis
##########
source("https://bioconductor.org/biocLite.R")
biocLite("methyAnalysis")












##########
# run bumphunter to get DMR
##########

# make type column
raw$type <- c('cancer', 'no_cancer')

# get indicator and put into design matrix with intercept 1
indicator_vector <- as.factor(raw$type)
designMatrix <- cbind(rep(1, nrow(raw)), indicator_vector)
designMatrix <- as.matrix(designMatrix)

raw <- as.data.frame(t(raw), stringsAsFactors = F)

# make probe a column in methyl
raw$probe <- rownames(raw)
rownames(raw) <- NULL

# inner join methyl and cg_locations by probe
methyl_cg <- inner_join(raw, cg_locations, by = 'probe')

# get chr and pos vector 
chr <- methyl_cg$seqnames
pos <- methyl_cg$start

# create beta matrix
beta <- methyl_cg[, 1:2]

# make beta numeric 
for (i in 1:ncol(beta)) {
  beta[,i] <- as.numeric(beta[,i])
  print(i)
} 

beta <- as.matrix(beta)

##########
# Run bumphunter
##########

# check dimensions 
stopifnot(dim(beta)[2] == dim(designMatrix)[1])
stopifnot(dim(beta)[1] == length(chr))
stopifnot(dim(beta)[1] == length(pos))

# set paramenters 
DELTA_BETA_THRESH = .10 # DNAm difference threshold
NUM_BOOTSTRAPS = 4   # number of randomizations

# create tab list
tab <- bumphunter(beta, 
                  designMatrix, 
                  chr = chr, 
                  pos = pos,
                  nullMethod = "bootstrap",
                  cutoff = DELTA_BETA_THRESH,
                  B = NUM_BOOTSTRAPS,
                  type = "Beta")

bump_hunter_results <- tab$table

# order by p value 
bump_hunter_results <- bump_hunter_results[order(bump_hunter_results$p.value),]
