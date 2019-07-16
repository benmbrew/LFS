####### Script will run bumphunter on (1) WT patients with no cancer and (2) LFS with no cancer

##########
# initialize libraries
##########
library(minfi)
library(bumphunter)
library(dplyr)

##########
# Initialize folders
##########
home_folder <- '/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
model_data <- paste0(data_folder, '/model_data')
clin_data <- paste0(data_folder, '/clin_data')

##########
# load controls 
##########

#load
controls <- readRDS(paste0(model_data, '/controls.rda'))

controls_full <- readRDS(paste0(model_data, '/controls_full.rda'))

controls_wt <- readRDS(paste0(model_data, '/controls_wt.rda'))

# ge cg_locations
cg_locations <- read.csv(paste0(model_data, '/cg_locations.csv'))

##########
# KS test for difference in age of sample collection
##########

# testKS <- function(x, y)
# {
#   y <- y[!is.na(y)]
#   x <- x[!is.na(x)]
#   
#   # Do x and y come from the same distribution?
#   ks.test(jitter(x), jitter(y), alternative = 'two.sided')
#   
# }
# 
# testKS(controls_full$age_sample_collection, controls_wt$age_sample_collection)
# testKS(controls$age_sample_collection, controls_wt$age_sample_collection)
# 
# testKS(controls_full_bal$age_sample_collection, controls_wt$age_sample_collection)
# testKS(controls_bal$age_sample_collection, controls_wt$age_sample_collection)



##########
# remove samples to get balanced age
##########
# data_wt <- controls_wt
# data_controls <- controls
getBalAge <- function(data_wt, data_controls, full)
{
  
  # # # # # balance age
  hist(data_wt$age_sample_collection)
  hist(data_controls$age_sample_collection)
  # 
  if (full){
    
    remove_index_2 <- which(data_controls$age_sample_collection > 200 & data_controls$age_sample_collection <= 300) 
    remove_index_2 <- sample(remove_index_2, 3, replace = F)
    
    data_controls <- data_controls[-remove_index_2,]
    
    # remove a few from ranges 100-200, 300-400
    # randomly remove controls that have a less than 50 month age of diganosis to have balanced classes
    remove_index <- which(data_controls$age_sample_collection >= 100 & data_controls$age_sample_collection <= 200) 
    remove_index <- sample(remove_index, 8, replace = F)
    
    data_controls <- data_controls[-remove_index,]
    
    
  } else {
    
    remove_index_2 <- which(data_controls$age_sample_collection > 200 & data_controls$age_sample_collection <= 300) 
    remove_index_2 <- sample(remove_index_2, 2, replace = F)
    
    data_controls <- data_controls[-remove_index_2,]
    
    # remove a few from ranges 100-200, 300-400
    # randomly remove controls that have a less than 50 month age of diganosis to have balanced classes
    remove_index <- which(data_controls$age_sample_collection >= 100 & data_controls$age_sample_collection <= 200) 
    remove_index <- sample(remove_index, 8, replace = F)
    
    data_controls <- data_controls[-remove_index,]
    
  }
  
  
  return(data_controls)
  
}



##########
# get balanced ages for each control group
##########

# full
controls_full_bal <- getBalAge(controls_wt, controls_full, full = T)

# sub
controls_bal <- getBalAge(controls_wt, controls, full = F)


##########
# function that takes LFS patients and run balanced and unbalanced bumphunter on cancer and controls
##########
# HERE
# check histograms  
# dat_controls_wt <- controls_wt
# dat_controls <- controls_full_bal
# will return one unbalanced, one balanced by age, and balanced by age and counts
bumpHunterBalanced <- function(dat_controls_wt,
                               dat_controls) {
  
  
  # combine data
  dat <- rbind(dat_controls_wt, dat_controls)
  
  
  ##########
  # get clinical dat 
  ##########
  bump_clin <- dat[,1:9]
  
  # recode type
  dat$type <- dat$p53_germline
  
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
  DELTA_BETA_THRESH = .10 # DNAm difference threshold
  NUM_BOOTSTRAPS = 2   # number of randomizations
  
  # create tab list
  tab <- list()
  bump_hunter_results <- list()
  for (i in 1:length(DELTA_BETA_THRESH)) {
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
# controls full and controls 
##########

# all cases
wt_bal_full <-  bumpHunterBalanced(controls_wt, controls_full_bal)

# sub cases
wt_bal <- bumpHunterBalanced(controls_wt, controls_bal)


###########
# rempve cases and controls
###########
rm(list=ls(pattern="cases"))
rm(list=ls(pattern="controls"))


###########
# save image of bh_features
###########
save.image(paste0(model_data, '/modal_feat_pred_10.RData'))


