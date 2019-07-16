#Twins data
# ##########
# # explore methylation with methyAnalysis
# ##########
# source("https://bioconductor.org/biocLite.R")
# biocLite("methyAnalysis")

##########
# initialize libraries
##########
library(dplyr)
library(sva)
library(minfi)
library(methyAnalysis)
library(sqldf)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Scripts/twins/data')
data_folder_2 <- paste0(project_folder, '/Data')
model_data <- paste0(data_folder_2, '/model_data')
results_folder <- paste0(project_folder, '/Results')
cg_folder <- paste0(project_folder, '/Scripts/twins')

# read cases
quan_cases_sam <- readRDS(paste0(model_data, '/quan_cases_sam.rda'))

model_feat <- colnames(quan_cases_sam)[10:ncol(quan_cases_sam)]

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

colnames(raw) <- c('no_cancer', 'cancer')
##########
# run bumphunter to get DMR
##########
raw <- as.data.frame(t(raw), stringsAsFactors = F)

# subset raw by model_features

# get intersection
raw_feat <- colnames(raw)
int_feat <- intersect(raw_feat, model_feat)
raw <- raw[, int_feat]

# make type column
raw$type <- rownames(raw)

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
NUM_BOOTSTRAPS = 2   # number of randomizations

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

colnames(cg_locations) <- paste0(colnames(cg_locations), '_', 'rgSet')


# use sql data table to merger validators with model_data based on age range
result = sqldf("select * from cg_locations
                 inner join bump_hunter_results
                 on cg_locations.start_rgSet between bump_hunter_results.start and bump_hunter_results.end")

# keep only necessary columns
result <- result[, c('chr' , 'start_rgSet','end_rgSet', 'probe_rgSet', 'p.value', 'fwer')]

# rename cols
colnames(result) <- c('chr', 'start', 'end', 'probe', 'p.value', 'fwer')

# get sig results
result_sig <- result[result$p.value < 0.05,]

# get fwer results
result_fwer <- result[result$fwer == 0,]


##########
# look at full_feat and bh results
##########
# read in full data
full_feat <- readRDS('/home/benbrew/Desktop/full_feat.rda')

# get features greater than 5
feat_5 <- full_feat

# subset for reg, 48, 60, 72
feat_reg <- feat_5[feat_5$V3 == 'reg',]
feat_48 <- feat_5[feat_5$V3 == '48',]
feat_60 <- feat_5[feat_5$V3 == '60',]
feat_72 <- feat_5[feat_5$V3 == '72',]

# merge with 
twins_reg <- inner_join(feat_reg, result, by = 'probe')
twins_48 <- inner_join(feat_48, result, by = 'probe')
twins_60 <- inner_join(feat_60, result, by = 'probe')
twins_72 <- inner_join(feat_72, result, by = 'probe')


##########
# look at mod_features probes in twins.
##########

raw_reg <- inner_join(feat_reg, raw, by = 'probe')
raw_48 <- inner_join(feat_48, raw, by = 'probe')
raw_60 <- inner_join(feat_60, raw, by = 'probe')
raw_72 <- inner_join(feat_72, raw, by = 'probe')

##########
# get distance 
##########
getDist <- function(data) {
  
  ##########
  # calculate avg distance between beta values in raw
  ##########
  raw$no_cancer <- as.numeric(raw$no_cancer)
  raw$cancer <- as.numeric(raw$cancer)
  
  dis_vec_raw <- abs(raw$no_cancer - raw$cancer)
  
  ##########
  # calculate avg distance between beta values in raw_reg
  ##########
  raw_reg$no_cancer <- as.numeric(raw_reg$no_cancer)
  raw_reg$cancer <- as.numeric(raw_reg$cancer)
  
  dis_vec_reg <- abs(raw_reg$no_cancer - raw_reg$cancer)
  
  ###########
  # imbalanced t test
  ###########
  t.test(dis_vec_reg, dis_vec_raw, alternative = 'greater', paired = F)
  
}

# should use all available probes chosen
# should exclude those probes in raw 


