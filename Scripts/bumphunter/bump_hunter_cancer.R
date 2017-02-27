####### Script will run bumphunter on p53 Mut and use non cancer (family members) as controls
# this is 7th  step in pipeline

##########
# initialize libraries
##########
library(minfi)
library(bumphunter)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
model_data <- paste0(data_folder, '/model_data')


#########
# Load model data batch, non batch, raw, quan, cases, controls
#########
# load image
load(paste0(model_data, '/model_data.RData'))

# ge cg_locations
cg_locations <- read.csv(paste0(model_data, '/cg_locations.csv'))

##########
# function that takes LFS patients and run balanced and unbalanced bumphunter on cancer and controls
# type is the indicator
##########

# will return one unbalanced, one balanced by age, and balanced by age and counts
bumpHunterBalanced <- function(dat_cases,
                               dat_controls,
                               bal_age) {
  
  
  # combine data
  dat <- rbind(dat_cases, dat_controls)
  
  if (bal_age) {
    # # balance age 
    hist(dat$age_sample_collection[dat$type == 'cases'])
    hist(dat$age_sample_collection[dat$type == 'controls'])
    
    # remove a few from ranges 100-200, 300-400
    # randomly remove controls that have a less than 50 month age of diganosis to have balanced classes
    remove_index <- which(dat$type == 'controls' & ((dat$age_sample_collection >= 100 & dat$age_sample_collection <= 200) |
                            (dat$age_sample_collection >= 300 & dat$age_sample_collection <= 400)))
    
    remove_index <- sample(remove_index, 8, replace = F )
    dat <- dat[-remove_index,]
    
  } 

  ##########
  # get clinical dat 
  ##########
  bump_clin <- dat[,1:7]
  
  ##########
  # get indicator and put into design matrix with intercept 1
  ##########
  indicator_vector <- as.factor(dat$type)
  designMatrix <- cbind(rep(1, nrow(dat)), indicator_vector)
  designMatrix <- as.matrix(designMatrix)
  
  ##########
  # Get genetic locations
  ##########
  dat$p53_germline <- dat$age_diagnosis <- dat$cancer_diagnosis_diagnoses <-
    dat$age_sample_collection <- dat$id <- dat$type <- dat$gender <-  NULL
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
  beta <- methyl_cg[, 1:(ncol(methyl_cg) - 7)]
  
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
  DELTA_BETA_THRESH = c(0.10, 0.20, 0.30, 0.40, 0.50) # DNAm difference threshold
  NUM_BOOTSTRAPS = 4   # number of randomizations
  
  # create tab list
  tab <- list()
  bump_hunter_results <- list()
  for(i in 1:length(DELTA_BETA_THRESH)) {
    tab[[i]] <- bumphunter(beta, 
                           designMatrix, 
                           chr = chr, 
                           pos = pos,
                           nullMethod = "bootstrap",
                           cutoff = DELTA_BETA_THRESH[i],
                           B = NUM_BOOTSTRAPS,
                           type = "Beta")
    
    bump_hunter_results[[i]] <- tab[[i]]$table
    bump_hunter_results[[i]]$run <- DELTA_BETA_THRESH[i]
  }
  
  bh_results <- do.call(rbind, bump_hunter_results)
  
  return(bh_results)
  
}


#########
# Now Apply to Original IDAT data
#########

##########
# raw
##########
# RAW no batch
raw_even <- bumpHunterBalanced(raw_cases, raw_controls, bal_age = T)
raw_uneven <- bumpHunterBalanced(raw_cases, raw_controls, bal_age = F)

# RAW batch
raw_batch_even <- bumpHunterBalanced(raw_cases_batch, raw_controls_batch, bal_age = T)
raw_batch_uneven <- bumpHunterBalanced(raw_cases_batch, raw_controls_batch, bal_age = F)

##########
# quantile
##########
# quan no batch
quan_even <- bumpHunterBalanced(quan_cases, quan_controls, bal_age = T)
quan_uneven <- bumpHunterBalanced(quan_cases, quan_controls, bal_age = F)

# quan batch
quan_batch_even <- bumpHunterBalanced(quan_cases_batch, quan_controls_batch, bal_age = T)
quan_batch_uneven <- bumpHunterBalanced(quan_cases_batch, quan_controls_batch, bal_age = F)

##########
# quality ontrol of bumphunter results (feature intersection - are higher thresholds a subset of lower ones?)
##########

##########
# remove non model objects and save
##########
rm(cg_locations, raw_cases, raw_cases_batch,
   quan_cases, quan_cases_batch, raw_controls, 
   raw_controls_batch, quan_controls, quan_controls_batch)

###########
# save image of bh_features
###########
save.image(paste0(model_data, '/modal_feat.RData'))

