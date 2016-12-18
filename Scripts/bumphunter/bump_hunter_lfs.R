#####################################################################################################
# This script will load bumphunter full data and run different variations of bumhunter on p53 status and save results 
library(minfi)
library(bumphunter)
library(dplyr)
library(MatchIt)

# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
clin_data <- paste0(data_folder, '/clin_data')
methyl_data <- paste0(data_folder, '/methyl_data')
model_data <- paste0(data_folder, '/model_data')
imputed_data <- paste0(data_folder, '/imputed_data')
idat_data <- paste0(methyl_data, '/raw_files')
bumphunter_data <- paste0(data_folder, '/bumphunter_data')


##########
# load imputed methylation data and cg_locations csv
##########
load(paste0(model_data, '/model_data_cases.RData'))
load(paste0(model_data, '/model_data_controls.RData'))

cg_locations <- read.csv(cg_locations, paste0(model_data, '/cg_locations.csv'))


#########
# function that looks within cancer and sets WT as controls
#########

bumpHunterBalanced <- function(dat,
                               balanced,
                               even_counts) {
  
  # subset data so its just cancer patietns
  dat <- dat[dat$cancer_diagnosis_diagnoses != 'Unaffected',]
  
  
  if (balanced) {
    
    hist(dat$age_sample_collection[dat$p53_germline == 'Mut'])
    hist(dat$age_sample_collection[dat$p53_germline == 'WT'])
    
    # randomly remove 3 WT that have age of sample collection between 100 and 200
    remove_index <- which(dat$p53_germline == 'WT' & (dat$age_sample_collection >= 100 & 
                                                        dat$age_sample_collection <= 200))
    remove_index <- sample(remove_index, 3, replace = F )
    dat <- dat[-remove_index,]
    
    if (even_counts) {
      
      # balance numbers
      # nrow(dat[dat$p53_germline == 'Mut',]) # 80
      # nrow(dat[dat$p53_germline == 'WT',]) #21
      remove_index <- which(dat$p53_germline == 'Mut')
      remove_index <- sample(remove_index, 59, replace = F )
      dat <- dat[-remove_index,]
      
    } 
    
  } 
  
  # get clinical dat 
  bump_clin <- dat[,1:5]
  
  # get p53 and put into design matrix with intercept 1
  p53_vector <- dat$p53_germline
  designMatrix <- cbind(rep(1, nrow(dat)), p53_vector)
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


##########
# Now Apply to Original IDAT data
##########

# beta raw
beta_raw_bal <- bumpHunterBalanced(beta_raw, balanced = T, even_counts = F)
beta_raw_bal_counts <- bumpHunterBalanced(beta_raw, balanced = T, even_counts = T)
beta_raw_unbal <- bumpHunterBalanced(beta_raw, balanced = F, even_counts = F)

# beta swan
beta_swan_bal <- bumpHunterBalanced(beta_swan, balanced = T, even_counts = F)
beta_swan_bal_counts <- bumpHunterBalanced(beta_swan, balanced = T, even_counts = T)
beta_swan_unbal <- bumpHunterBalanced(beta_swan, balanced = F, even_counts = F)

# beta quan
beta_quan_bal <- bumpHunterBalanced(beta_quan, balanced = T, even_counts = F)
beta_quan_bal_counts <- bumpHunterBalanced(beta_quan, balanced = T, even_counts = T)
beta_quan_unbal <- bumpHunterBalanced(beta_quan, balanced = F, even_counts = F)

# beta funnorm
beta_funnorm_bal <- bumpHunterBalanced(beta_funnorm, balanced = T, even_counts = F)
beta_funnorm_bal_counts <- bumpHunterBalanced(beta_funnorm, balanced = T, even_counts = T)
beta_funnorm_unbal <- bumpHunterBalanced(beta_funnorm, balanced = F, even_counts = F)

# remove unneeded objects
rm(beta_raw, beta_swan, beta_quan, beta_funnorm, cg_locations)

# save data
save.image(paste0(idat_data, '/imputed_bh_lfs.RData'))

