####### Script will run bumphunter on p53 Mut and use non cancer (family members) as controls
# this is 6th (part A) step in pipeline

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
# dat_cases <- raw_cases_batch
# dat_controls <- raw_controls_batch
# will return one unbalanced, one balanced by age, and balanced by age and counts
bumpHunterBalanced <- function(dat_cases,
                               dat_controls,
                               balanced,
                               even_counts) {
  
  
  # combine data
  dat <- rbind(dat_cases, dat_controls)
  
  if (balanced) {
    # # balance age 
    hist(dat$age_sample_collection[dat$type == 'cases'])
    hist(dat$age_sample_collection[dat$type == 'controls'])
    
    # remove a few from ranges 100-200, 300-400
    # randomly remove controls that have a less than 50 month age of diganosis to have balanced classes
    remove_index <- which(dat$type == 'controls' & (dat$age_sample_collection >= 100 & dat$age_sample_collection <= 200) |
                            (dat$age_sample_collection >= 300 & dat$age_sample_collection <= 400))
    
    remove_index <- sample(remove_index, 15, replace = F )
    dat <- dat[-remove_index,]
    
    if (even_counts) {
      
      # balance numbers
      nrow(dat[dat$type == 'cases',]) # 70
      nrow(dat[dat$type == 'controls',]) #26
      remove_index <- which(dat$type == 'cases')
      remove_index <- sample(remove_index, 41, replace = F )
      dat <- dat[-remove_index,]
      
    } 
  } 

  ##########
  # get clinical dat 
  ##########
  bump_clin <- dat[,1:8]
  
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
    dat$age_sample_collection <- dat$id <- dat$type <- dat$batch <-  NULL
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
  DELTA_BETA_THRESH = c(0.10, 0.15, 0.20) # DNAm difference threshold
  NUM_BOOTSTRAPS = 3    # number of randomizations
  
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

# bal, even
raw_bh <- bumpHunterBalanced(raw_cases_batch, raw_controls_batch, balanced = T, even_counts = T)
quan_bh <- bumpHunterBalanced(quan_cases_batch, quan_controls_batch, balanced = T, even_counts = T)
swan_bh <- bumpHunterBalanced(swan_cases_batch, swan_controls_batch, balanced = T, even_counts = T)
funnorm_bh <- bumpHunterBalanced(funnorm_cases_batch, funnorm_controls_batch, balanced = T, even_counts = T)

# unbal, uneven
raw_unbal_bh <- bumpHunterBalanced(raw_cases_batch, raw_controls_batch, balanced = F, even_counts = F)
quan_unbal_bh <- bumpHunterBalanced(quan_cases_batch, quan_controls_batch, balanced = F, even_counts = F)
swan_unbal_bh <- bumpHunterBalanced(swan_cases_batch, swan_controls_batch, balanced = F, even_counts = F)
funnorm_unbal_bh <- bumpHunterBalanced(funnorm_cases_batch, funnorm_controls_batch, balanced = F, even_counts = F)


# save new bh 
saveRDS(raw_bh, paste0(model_data, '/raw_bh.rda'))
saveRDS(quan_bh, paste0(model_data, '/quan_bh.rda'))
saveRDS(swan_bh, paste0(model_data, '/swan_bh.rda'))
saveRDS(funnorm_bh, paste0(model_data, '/funnorm_bh.rda'))

# save new bh  unbal
saveRDS(raw_unbal_bh, paste0(model_data, '/raw_unbal_bh.rda'))
saveRDS(quan_unbal_bh, paste0(model_data, '/quan_unbal_bh.rda'))
saveRDS(swan_unbal_bh, paste0(model_data, '/swan_unbal_bh.rda'))
saveRDS(funnorm_unbal_bh, paste0(model_data, '/funnorm_unbal_bh.rda'))





# # beta raw
# beta_raw_bal_cancer <- bumpHunterBalanced(beta_raw_full, balanced = T, even_counts = F)
# beta_raw_bal_counts_cancer <- bumpHunterBalanced(beta_raw_full, balanced = T, even_counts = T)
# beta_raw_unbal_cancer <- bumpHunterBalanced(beta_raw_full, balanced = F, even_counts = F)
# 
# # beta swan
# beta_swan_bal_cancer <- bumpHunterBalanced(beta_swan_full, balanced = T, even_counts = F)
# beta_swan_bal_counts_cancer <- bumpHunterBalanced(beta_swan_full, balanced = T, even_counts = T)
# beta_swan_unbal_cancer <- bumpHunterBalanced(beta_swan_full, balanced = F, even_counts = F)
# 
# # beta quan
# beta_quan_bal_cancer <- bumpHunterBalanced(beta_quan_full, balanced = T, even_counts = F)
# beta_quan_bal_counts_cancer <- bumpHunterBalanced(beta_quan_full, balanced = T, even_counts = T)
# beta_quan_unbal_cancer <- bumpHunterBalanced(beta_quan_full, balanced = F, even_counts = F)
# 
# # beta funnorm
# beta_funnorm_bal_cancer <- bumpHunterBalanced(beta_funnorm_full, balanced = T, even_counts = F)
# beta_funnorm_bal_counts_cancer <- bumpHunterBalanced(beta_funnorm_full, balanced = T, even_counts = T)
# beta_funnorm_unbal_cancer <- bumpHunterBalanced(beta_funnorm_full, balanced = F, even_counts = F)
# 
# # remove larger objects
# rm(beta_raw_full, beta_swan_full, beta_quan_full, beta_funnorm_full, clin)
# 
# load(paste0(model_data, '/beta_cancer_bh.RData'))



# save.image(paste0(model_data, '/beta_cancer_bh.RData'))

