#####################################################################################################
# This script will load bumphunter full data and run different variations of bumhunter
# using lfs patients and controls are those without cancer

##########
# initialize libraries
##########
library(minfi)
library(bumphunter)
library(dplyr)
library(MatchIt)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
clin_data <- paste0(data_folder, '/clin_data')
methyl_data <- paste0(data_folder, '/methyl_data')
model_data <- paste0(data_folder, '/model_data')
imputed_data <- paste0(data_folder, '/imputed_data')
idat_data <- paste0(methyl_data, '/raw_files')
bumphunter_data <- paste0(data_folder, '/bumphunter_data')

#########
# Load model data which contains all variations of methylation data
#########

# load controls and rename
load(paste0(idat_data, '/imputed_idat_betas_final_control.RData'))
beta_raw_controls <- beta_raw
beta_swan_controls <- beta_swan
beta_quan_controls <- beta_quan
beta_funnorm_controls <- beta_funnorm
# remove old named controls
rm(beta_raw, beta_swan, beta_quan, beta_funnorm)

# load original
load(paste0(idat_data, '/imputed_idat_betas_final.RData'))


##########
# function to find overlapping probes and subset control and regular by those probes
##########
findOverlap <- function(dat, dat_controls) {
  control_probes <- colnames(dat_controls)[5:ncol(dat_controls)]
  orig_probes <- colnames(dat)[6:ncol(dat)]
  overlaps <- intersect(control_probes, orig_probes)
  dat_controls <- dat_controls[,  c('id', 'p53_germline', 'cancer_diagnosis_diagnoses', 'age_sample_collection', overlaps)]
  dat_controls$age_diagnosis <- NA
  dat <- dat[,  c('id', 'p53_germline', 'age_diagnosis', 'cancer_diagnosis_diagnoses', 'age_sample_collection', overlaps)]
  dat_full <- rbind(dat, dat_controls)
  return(dat_full)
}

# get data
beta_raw_full <- findOverlap(beta_raw, beta_raw_controls)
beta_swan_full <- findOverlap(beta_swan, beta_swan_controls)
beta_quan_full <- findOverlap(beta_quan, beta_quan_controls)
beta_funnorm_full <- findOverlap(beta_funnorm, beta_funnorm_controls)

# remove smaller object
rm(beta_raw, beta_raw_controls, beta_swan, beta_swan_controls, beta_quan, beta_quan_controls, 
   beta_funnorm, beta_funnorm_controls)

##########
# function that takes LFS patients and run balanced and unbalanced bumphunter on cancer and controls
##########

bumpHunterBalanced <- function(dat,
                               balanced,
                               even_counts) {
  
  # get p53 mutatnts
  dat <- dat[dat$p53_germline == 'Mut',]
  
  # create column for controls and cases
  dat$indicator <- ifelse(dat$cancer_diagnosis_diagnoses != 'Unaffected', 'cases', 'controls')
  
  if (balanced) {
    # # balance age 
    # hist(dat$age_sample_collection[dat$indicator == 'cases'])
    # hist(dat$age_sample_collection[dat$indicator == 'controls'])
    # remove a few from ranges 100-200, 300-400
    # randomly remove Mut that have a less than 50 month age of diganosis to have balanced classes
    remove_index <- which(dat$indicator == 'controls' & (dat$age_sample_collection >= 100 & dat$age_sample_collection <= 200) |
                            (dat$age_sample_collection >= 300 & dat$age_sample_collection <= 400))
    
    remove_index <- sample(remove_index, 15, replace = F )
    dat <- dat[-remove_index,]
    
    if (even_counts) {
      
 # balance numbers
      # nrow(dat[dat$indicator == 'cases',]) # 80
      # nrow(dat[dat$indicator == 'controls',]) #50
      remove_index <- which(dat$indicator == 'cases')
      remove_index <- sample(nrow(dat), 30, replace = F )
      dat <- dat[-remove_index,]
      
    } 
    
    
    
  } 

  
  ##########
  # get clinical dat 
  ##########
  bump_clin <- dat[,1:5]
  
  ##########
  # get indicator and put into design matrix with intercept 1
  ##########
  indicator_vector <- as.factor(dat$indicator)
  designMatrix <- cbind(rep(1, nrow(dat)), indicator_vector)
  designMatrix <- as.matrix(designMatrix)
  
  ##########
  # Get genetic locations
  ##########
  dat$p53_germline <- dat$age_diagnosis <- dat$cancer_diagnosis_diagnoses <-
    dat$age_sample_collection <- dat$id <- NULL
  # transpose methylation to join with cg_locations to get genetic location vector.
  dat <- as.data.frame(t(dat), stringsAsFactors = F)
  
  # make probe a column in methyl
  dat$probe <- rownames(dat)
  rownames(dat) <- NULL
  
  # inner join methyl and cg_locations by probe
  methyl_cg <- inner_join(dat, cg_locations, by = 'probe')
  
  # get chr and pos vector 
  chr <- methyl_cg$seqnames
  pos <- methyl_cg$start
  
  # create beta matrix
  beta <- methyl_cg[, 1:(ncol(methyl_cg) - 6)]
  
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
  
  return(bump_hunter_results)
  
}


#########
# Now Apply to Original IDAT data
#########

# beta raw
beta_raw_unbal <- bumpHunterBalanced(beta_raw, unbalanced = F)
beta_raw_bal_ <- bumpHunterBalanced(beta_raw, unbalanced = F)

beta_raw_global_bal <- bumpHunterBalanced(beta_raw, unbalanced = F)[[1]]
beta_raw_cancer_unbal <- bumpHunterBalanced(beta_raw, unbalanced = T)[[1]]
beta_raw_global_unbal <- bumpHunterBalanced(beta_raw, unbalanced = T)[[1]]


# save.image(paste0(idat_data, '/imputed_idat_betas_bh.RData'))

