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
# Load data
#########

# read cases
quan_cases <- readRDS(paste0(model_data, '/quan_cases.rda'))

# read contorls
quan_controls <- readRDS(paste0(model_data, '/quan_controls.rda'))

# read batch corrected data for gender
quan_cases_gen <- readRDS(paste0(model_data, '/quan_cases_gen.rda'))

# read controls
quan_controls_gen <- readRDS(paste0(model_data, '/quan_controls_gen.rda'))

# read batch corrected data for sentrix id and SAM
quan_cases_sen <- readRDS(paste0(model_data, '/quan_cases_sen.rda'))

quan_cases_sam <- readRDS(paste0(model_data, '/quan_cases_sam.rda'))

# read batch corrected data for sentrix id and SAM and gender!
quan_cases_sen_gen <- readRDS(paste0(model_data, '/quan_cases_sen_gen.rda'))

quan_cases_sam_gen <- readRDS(paste0(model_data, '/quan_cases_sam_gen.rda'))

# ge cg_locations
cg_locations <- read.csv(paste0(model_data, '/cg_locations.csv'))

##########
# function that takes LFS patients and run balanced and unbalanced bumphunter on cancer and controls
# type is the indicator
##########
# dat_cases <- quan_cases
# dat_controls <- quan_controls
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
    
    remove_index <- sample(remove_index, 10, replace = F )
    dat <- dat[-remove_index,]
    
  } 

  ##########
  # get clinical dat 
  ##########
  bump_clin <- dat[,1:9]
  
  ##########
  # get indicator and put into design matrix with intercept 1
  ##########
  indicator_vector <- as.factor(dat$type)
  designMatrix <- cbind(rep(1, nrow(dat)), indicator_vector)
  designMatrix <- as.matrix(designMatrix)
  
  ##########
  # Get genetic locations
  ##########
  dat$p53_germline <- dat$age_diagnosis <- dat$cancer_diagnosis_diagnoses <- dat$ids <- dat$batch <- 
    dat$age_sample_collection <- dat$id <- dat$type <- dat$gender <-  dat$sentrix_id <-  NULL
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
  DELTA_BETA_THRESH = c(0.10, 0.20, 0.30) # DNAm difference threshold
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


##########
# quantile
##########
# quan no batch
quan_even <- bumpHunterBalanced(quan_cases, quan_controls, bal_age = T)
quan_uneven <- bumpHunterBalanced(quan_cases, quan_controls, bal_age = F)

# quan gender
quan_even_gen <- bumpHunterBalanced(quan_cases_gen, quan_controls_gen, bal_age = T)
quan_uneven_gen <- bumpHunterBalanced(quan_cases_gen, quan_controls_gen, bal_age = F)

# quan gender sentrix
quan_even_gen_sen <- bumpHunterBalanced(quan_cases_sen_gen, quan_controls_gen, bal_age = T)
quan_uneven_gen_sen <- bumpHunterBalanced(quan_cases_sen_gen, quan_controls_gen, bal_age = F)

# quan gender sam
quan_even_gen_sam <- bumpHunterBalanced(quan_cases_sam_gen, quan_controls_gen, bal_age = T)
quan_uneven_gen_sam <- bumpHunterBalanced(quan_cases_sam_gen, quan_controls_gen, bal_age = F)


##########
# remove non model objects and save
##########
rm(quan_cases_gen, quan_controls_gen,
   quan_cases, quan_controls,
   quan_cases_sen, quan_cases_sen_gen,
   quan_cases_sam_gen, quan_cases_sam)

###########
# save image of bh_features
###########
save.image(paste0(model_data, '/modal_feat.RData'))

