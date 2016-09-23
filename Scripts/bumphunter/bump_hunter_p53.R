#####################################################################################################
# This script will load bumphunter full data and run different variations of bumhunter on p53 status and save results 
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
  
  if (selection == 'cancer') {
    full_data <- full_data[!is.na(full_data$p53_germline),]
    full_data <- full_data[full_data$cancer_diagnosis_diagnoses != 'Unaffected',]
  }
  if (selection == 'global') {
    full_data <- full_data[!is.na(full_data$p53_germline),]
  }
  
  # get clinical data 
  bump_clin <- full_data[,212352:212379]
  
  
  # get p53 and put into design matrix with intercept 1
  p53_vector <- full_data$p53_germline
  designMatrix <- cbind(rep(1, nrow(full_data)), p53_vector)
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
  
  if (selection == 'cancer') {
    write.csv(bump_hunter_results, paste0(data_folder, '/bh_cancer.csv'))
  }
  if (selection == 'global') {
    write.csv(bump_hunter_results, paste0(data_folder, '/bh_global.csv')) 
  }

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

cancer <- bumpHunter(selection = 'cancer')
global <- bumpHunter(selection = 'global')

# look at summary stats of global and cancer

# Cancer 
# counts for wt
cancer[[2]]
# counts for mut
cancer[[3]]
# summary for age wt
cancer[[4]]
#summary for age mut
cancer[[5]]
#summary for cancer wt
cancer[[6]]
#summary for cancer mut
cancer[[7]]

## Global
# counts for wt
global[[2]]
# counts for mut
global[[3]]
# summary for age wt
global[[4]]
#summary for age mut
global[[5]]
#summary for cancer wt
global[[6]]
#summary for cancer mut
global[[7]]


####################################
# function that matches mut and wt counts for global and cancer 
####################################

# make function that selects the appropriate WT population and runs bumphunter
bumpHunterMatch <- function(selection,
                            iterations) {
  
  if (selection == 'cancer') {
    full_data <- full_data[!is.na(full_data$p53_germline),]
    full_data <- full_data[full_data$cancer_diagnosis_diagnoses != 'Unaffected',]
  }
  if (selection == 'global') {
    full_data <- full_data[!is.na(full_data$p53_germline),]
  }
  
  # create lists to store results of loop
  bump_hunter_results <- list()
  wt_age_summary <- list()
  mut_age_summary <- list()
  wt_cancer_summary <- list()
  mut_cancer_summary <- list()
 
  # loop that runs bumphunter a number of times, each time selecting a different number of Mutants 
  # that match the number of WT
  for ( i in 1:iterations) {
    
    set.seed(i)
    # subset full data so that counts of Mut match counts of WT
    num_of_wt <- nrow(full_data[full_data$p53_germline == 'WT',])
    wt_index <- which(full_data$p53_germline == 'WT')
    mut_index <- which(full_data$p53_germline == 'Mut')
    
    # sample mutants 14 times 
    equal_index <- sample(mut_index, num_of_wt)
    
    data_index <- sort(append(equal_index, wt_index))
    
    data <- full_data[data_index,]
    
    # get clinical data 
    bump_clin <- data[,212352:212379]
    # analyze the distribution of cancer and age in the bumphunter data 
    
    # get counts for mut and WT 
    wt_count <- nrow(bump_clin[which(bump_clin$p53_germline == 'WT'),])
    mut_count <- nrow(bump_clin[which(bump_clin$p53_germline == 'Mut'),])
    
    
    # difference in age between 
    wt_age_summary[[i]] <- summary(bump_clin$age_diagnosis[which(bump_clin$p53_germline == 'WT')])
    mut_age_summary[[i]] <- summary(bump_clin$age_diagnosis[which(bump_clin$p53_germline == 'Mut')])
    
    
    # difference in cancer cancer 
    wt_cancer_summary[[i]] <- summary(bump_clin$cancer_diagnosis_diagnoses[which(bump_clin$p53_germline == 'WT')])
    mut_cancer_summary[[i]] <- summary(bump_clin$cancer_diagnosis_diagnoses[which(bump_clin$p53_germline == 'Mut')])
    
    
    stopifnot(nrow(data[data$p53_germline == 'Mut',]) == nrow(data[data$p53_germline == 'WT',]))
    # get p53 and put into design matrix with intercept 1
    p53_vector <- data$p53_germline
    designMatrix <- cbind(rep(1, nrow(data)), p53_vector)
    designMatrix <- as.matrix(designMatrix)
    
    ######################
    # Get genetic locations
    ######################
    
    # transpose methylation to join with cg_locations to get genetic location vector.
    data <- as.data.frame(t(data), stringsAsFactors = F)
    
    # make ids column names and remove first row
    colnames(data) <- data[1,]
    data <- data[-1,]
    
    # make probe a column in methyl
    data$probe <- rownames(data)
    rownames(data) <- NULL
    
    # remove clinical rows at bottom of data set
    bump_data <- data[1:212351,]
    
    # inner join methyl and cg_locations by probe
    methyl_cg <- inner_join(cg_locations, bump_data, by = 'probe')
    
    # get chr and pos vector 
    chr <- methyl_cg$seqnames
    pos <- methyl_cg$start
    
    # create beta matrix
    rownames(methyl_cg) <- methyl_cg$probe
    beta <- methyl_cg[, 7:ncol(methyl_cg)]
    
    # make beta numeric 
    for (j in 1:ncol(beta)) {
      beta[,j] <- as.numeric(beta[,j])
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
    
    bump_hunter_results[[i]] <- tab$table
    
    print(i)
  }
  
  # combine bumphunter results from all iterations and save
  
  if (selection == 'cancer') {
    bump_hunter_results <- do.call('rbind', bump_hunter_results)
    bump_hunter_results <- bump_hunter_results[!duplicated(bump_hunter_results[1:2]),]
    write.csv(bump_hunter_results, paste0(data_folder, '/bh_cancer_sub.csv'))
  }
  if (selection == 'global') {
    bump_hunter_results <- do.call('rbind', bump_hunter_results)
    bump_hunter_results <- bump_hunter_results[!duplicated(bump_hunter_results[1:2]),]
    write.csv(bump_hunter_results, paste0(data_folder, '/bh_global_sub.csv')) 
  }
  
  
  return(list(bump_hunter_results, wt_count, mut_count, wt_age_summary, mut_age_summary, wt_cancer_summary, mut_cancer_summary))
  
}

cancer <- bumpHunterMatch(selection = 'cancer', iterations = 10)
global <- bumpHunterMatch(selection = 'global', iterations = 10)


# look at summary stats of global and cancer

# Cancer 
# counts for wt
cancer[[2]]
# counts for mut
cancer[[3]]
# summary for age wt
cancer[[4]]
#summary for age mut
age_cancer <- as.data.frame(do.call('rbind', cancer[[5]]))
mean(age_cancer$Min.)
mean(age_cancer$Median)
mean(age_cancer$Mean)
mean(age_cancer$Max)

#summary for cancer wt
cancer[[6]]
#summary for cancer mut
diagnosis <- as.data.frame(do.call('rbind', cancer[[7]]))
mean(diagnosis$ACC)


# global 
# counts for wt
global[[2]]
# counts for mut
global[[3]]
# summary for age wt
global[[4]]
#summary for age mut
age_global <- as.data.frame(do.call('rbind', global[[5]]))
mean(age_global$Min.)
mean(age_global$Median)
mean(age_global$Mean)
mean(age_global$Max)

#summary for global wt
global[[6]]
#summary for global mut
diagnosis <- as.data.frame(do.call('rbind', global[[7]]))
mean(diagnosis$ACC)



####################################
# function that matches mut and wt counts as well as creates relatively similar samples for each category
####################################
# load(paste0(data_folder, '/full_data_bh.RData'))

# make function that selects the appropriate WT population and runs bumphunter
bumpHunterBalanced <- function(selection) {
  
  if (selection == 'cancer') {
    full_data <- full_data[!is.na(full_data$p53_germline),]
    full_data <- full_data[full_data$cancer_diagnosis_diagnoses != 'Unaffected',]
    # first subset by ACC and Other (the biggest difference is between "Other" cancers - so subset by just ACC and Unaffected)
    full_data <- full_data[grepl('ACC', full_data$cancer_diagnosis_diagnoses),]
    full_data <- full_data[full_data$age_diagnosis <= 180,]
  }
  if (selection == 'global') {
    full_data <- full_data[!is.na(full_data$p53_germline),]
    # first subset by ACC and Other (the biggest difference is between "Other" cancers - so subset by just ACC and Unaffected)
    full_data <- full_data[grepl('ACC|Unaffected', full_data$cancer_diagnosis_diagnoses),]
    full_data <- full_data[full_data$age_diagnosis <= 180 | is.na(full_data$age_diagnosis),]
    # full_data[full_data$cancer_diagnosis_diagnoses != 'Unaffected',] <- full_data$cancer_diagnosis_diagnoses != 'Unaffected'[full_data$age_diagnosis <= 180,]
  }
  
  # Make a mut and WT population have "similar" age
  #########################################
  
  # remove rows that are all NA 
  all_na_ind <- apply(full_data, 1, function(x) all(is.na(x)))
  full_data <- full_data[ !all_na_ind, ]
  
  # get clinical data 
  bump_clin <- full_data[,212352:212379]
  
  # get p53 and put into design matrix with intercept 1
  p53_vector <- full_data$p53_germline
  designMatrix <- cbind(rep(1, nrow(full_data)), p53_vector)
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
  
  if (selection == 'cancer') {
    write.csv(bump_hunter_results, paste0(data_folder, '/bh_cancer_balanced.csv'))
  }
  if (selection == 'global') {
    write.csv(bump_hunter_results, paste0(data_folder, '/bh_global_balanced.csv')) 
  }
  
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

cancer <- bumpHunterBalanced(selection = 'cancer')
global <- bumpHunterBalanced(selection = 'global')

# look at summary stats of global and cancer

# Cancer 
# counts for wt
cancer[[2]]
# counts for mut
cancer[[3]]
# summary for age wt
cancer[[4]]
#summary for age mut
cancer[[5]]
#summary for cancer wt
cancer[[6]]
#summary for cancer mut
cancer[[7]]

## Global
# counts for wt
global[[2]]
# counts for mut
global[[3]]
# summary for age wt
global[[4]]
#summary for age mut
global[[5]]
#summary for cancer wt
global[[6]]
#summary for cancer mut
global[[7]]

# 
# bh_cancer_sub
# bh_global_sub
# bh_cancer_balanced
# bh_global_balanced
# bh_cancer
# bh_global
