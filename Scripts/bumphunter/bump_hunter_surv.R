####### Script will run bumphunter on (1) LFS patients with cancer (cases) and (2) LFS patients wtihout cancer (controls)
# this is 5th  step in pipeline

##########
# initialize libraries
##########
library(minfi)
library(bumphunter)
library(dplyr)

##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
model_data <- paste0(data_folder, '/model_data')
clin_data <- paste0(data_folder, '/clin_data')

##########
# load cases 
##########
# full
quan_cases_full <- readRDS(paste0(model_data, '/quan_cases_full.rda'))

funnorm_cases_full <- readRDS(paste0(model_data, '/funnorm_cases_full.rda'))

raw_cases_full <- readRDS(paste0(model_data, '/raw_cases_full.rda'))


# sub
quan_cases_sub <- readRDS(paste0(model_data, '/quan_cases_sub.rda'))

funnorm_cases_sub <- readRDS(paste0(model_data, '/funnorm_cases_sub.rda'))

raw_cases_sub <- readRDS(paste0(model_data, '/raw_cases_sub.rda'))


##########
# load controls 
##########
# full
quan_controls_full <- readRDS(paste0(model_data, '/quan_controls_full.rda'))

funnorm_controls_full <- readRDS(paste0(model_data, '/funnorm_controls_full.rda'))

raw_controls_full <- readRDS(paste0(model_data, '/raw_controls_full.rda'))


# sub
quan_controls_sub <- readRDS(paste0(model_data, '/quan_controls_sub.rda'))

funnorm_controls_sub <- readRDS(paste0(model_data, '/funnorm_controls_sub.rda'))

raw_controls_sub <- readRDS(paste0(model_data, '/raw_controls_sub.rda'))


# ge cg_locations
cg_locations <- read.csv(paste0(model_data, '/cg_locations.csv'))


##########
# remove samples to get balanced age
##########

getBalAge <- function(data_controls, full)
{
  # # balance age
  # hist(quan_cases_full$age_sample_collection)
  # hist(data_controls$age_sample_collection)

  # remove a few from ranges 100-200, 300-400
  # randomly remove controls that have a less than 50 month age of diganosis to have balanced classes
  remove_index <- which((data_controls$age_sample_collection >= 100 & data_controls$age_sample_collection <= 250) |
                       (data_controls$age_sample_collection >= 300 & data_controls$age_sample_collection <= 400))
  
  if(full) {
    remove_index <- sample(remove_index, 13, replace = F)
    
  } else {
    remove_index <- sample(remove_index, 10, replace = F)
    
  }
  
  data_controls <- data_controls[-remove_index,]
  
  return(data_controls)
  
}

##########
# get balanced ages for each control group
##########

# full
quan_controls_full_bal <- getBalAge(quan_controls_full, full = T)
funnorm_controls_full_bal <- getBalAge(funnorm_controls_full, full = T)

# sub
quan_controls_sub_bal <- getBalAge(quan_controls_sub, full = F)
funnorm_controls_sub_bal <- getBalAge(funnorm_controls_sub, full = F)


##########
# function that takes LFS patients and run balanced and unbalanced bumphunter on cancer and controls
##########
# HERE
# check histograms  
dat_cases <- quan_cases_full
dat_controls <- quan_controls_full_bal
# will return one unbalanced, one balanced by age, and balanced by age and counts
bumpHunterBalanced <- function(dat_cases,
                               dat_controls) {
  
  
  # combine data
  dat <- rbind(dat_cases, dat_controls)


  ##########
  # get clinical dat 
  ##########
  bump_clin <- dat[,1:9]
  
  # recode type
  dat$type <- ifelse(grepl('Unaffected', dat$cancer_diagnosis_diagnoses), 'controls', 'cases')
  
  ##########
  # get indicator and put into design matrix with intercept 1
  #########
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
  DELTA_BETA_THRESH = c(0.10, 0.15, 0.20, 0.25) # DNAm difference threshold
  NUM_BOOTSTRAPS = 3   # number of randomizations
  
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
# quantile even
##########
# quan no batch
even <- bumpHunterBalanced(quan_cases, quan_controls_bal)

# quan gender
even_gen <- bumpHunterBalanced(quan_cases_gen, quan_controls_gen_bal)

# quan
even_type <- bumpHunterBalanced(quan_cases_gen, quan_controls_type_bal)

# quan gender sentrix
even_gen_sen_type <- bumpHunterBalanced(quan_cases_sen_gen, quan_controls_type_bal)

# quan gender sam
even_gen_sam_type <- bumpHunterBalanced(quan_cases_sam_gen, quan_controls_type_bal)


##########
# quantile uneven
##########
# quan no batch
uneven <- bumpHunterBalanced(quan_cases, quan_controls)

# quan gender
uneven_gen <- bumpHunterBalanced(quan_cases_gen, quan_controls_gen)

# quan
uneven_type <- bumpHunterBalanced(quan_cases_gen, quan_controls_type)

# quan gender sentrix
uneven_gen_sen_type <- bumpHunterBalanced(quan_cases_sen_gen, quan_controls_type)

# quan gender sam
uneven_gen_sam_type <- bumpHunterBalanced(quan_cases_sam_gen, quan_controls_type)


##########
# remove non model objects and save
##########
rm(quan_cases_gen, quan_controls_gen,
   quan_cases, quan_controls,
   quan_cases_sen, quan_cases_sen_gen,
   quan_cases_sam_gen, quan_cases_sam,
   quan_controls, quan_controls_bal, 
   quan_controls_gen, quan_controls_gen_bal,
   quan_controls_type, quan_controls_type_bal,
   cg_locations)

###########
# save image of bh_features
###########
save.image(paste0(model_data, '/modal_feat.RData'))

