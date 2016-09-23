#####################################################################################################
# This script will load bumphunter full data and run different variations of bumhunter on cancer status and save results 
library(minfi)
library(bumphunter)
library(dplyr)
library(FDb.InfiniumMethylation.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(GenomicRanges)
library(biovizBase)
library(MatchIt)

# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
clin_data <- paste0(data_folder, '/clin_data')

# read in cg_locations 
# read in cg_locations
cg_locations <- read.csv(paste0(data_folder, '/cg_locations.csv'))

# rename X column as probe for merge
names(cg_locations)[1] <- 'probe'

# load bh data
load(paste0(data_folder, '/full_data_bh.RData'))

# make function that selects the appropriate WT population and runs bumphunter
bumpHunter <- function(selection) {
  
  # create factor variable to indicate cancer
  full_data <- full_data[!is.na(full_data$cancer_diagnosis_diagnoses),]
  full_data$cancer <- ifelse(full_data$cancer_diagnosis_diagnoses != 'Unaffected', 'yes', 'no')
  
  # get clinical data 
  bump_clin <- full_data[,212352:212380]
  
  # get p53 and put into design matrix with intercept 1
  cancer_vector <- as.factor(full_data$cancer)
  designMatrix <- cbind(rep(1, nrow(full_data)), cancer_vector)
  designMatrix <- as.matrix(designMatrix)
  
  ######################
  # Get genetic locations
  ######################
  
  # transpose methylation to join with cg_locations to get genetic location vector.
  full_data <- as.data.frame(t(full_data), stringsAsFactors = F)
  
  # make ids column names and remove first row
  colnames(full_data) <- full_data[1,]
  full_data <- full_data[-1,]
  
  # make probe a column in methyl
  full_data$probe <- rownames(full_data)
  rownames(full_data) <- NULL
  
  # remove clinical rows at bottom of data set
  bump_data <- full_data[1:212351,]
  
  # inner join methyl and cg_locations by probe
  methyl_cg <- inner_join(cg_locations, bump_data, by = 'probe')
  
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
  
  
  ######################
  # Run bumphunter
  ######################
  
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
  
  write.csv(bump_hunter_results, paste0(data_folder, '/bh_cancer_indicator.csv'))
 
  # analyze the distribution of cancer and age in the bumphunter data 
  
  # get counts for mut and WT 
  wt_count <- nrow(bump_clin[which(bump_clin$p53_germline == 'WT'),])
  mut_count <- nrow(bump_clin[which(bump_clin$p53_germline == 'Mut'),])
  
  
  # difference in age between 
  wt_age_summary <- summary(bump_clin$age_diagnosis[which(bump_clin$p53_germline == 'WT')])
  mut_age_summary <- summary(bump_clin$age_diagnosis[which(bump_clin$p53_germline == 'Mut')])
  
  
  # difference in cancer cancer 
  wt_cancer_summary <- summary(bump_clin$cancer_diagnosis_diagnoses[which(bump_clin$p53_germline == 'WT')])
  mut_cancer_summary <- summary(bump_clin$cancer_diagnosis_diagnoses[which(bump_clin$p53_germline == 'Mut')])
  
  
  return(list(bump_hunter_results, wt_count, mut_count, wt_age_summary, mut_age_summary, wt_cancer_summary, mut_cancer_summary))
  
}

cancer <- bumpHunter(selection = NULL)

