
# source all_functions.R to load libraries and my functions
source('all_functions.R')

##########
# get base directory for 4 batch
##########
path_to_cases_tor <- '../../Data/methyl_data/cases_toronto'
path_to_cases_mon <- '../../Data/methyl_data/cases_montreal'
path_to_controls <- '../../Data/methyl_data/controls'
path_to_valid <- '../../Data/methyl_data/validation'

##########
# read in meth array - Data/methyl_data/cases_toronto, cases_montreal, controls, validation
##########

# cases 
rgCasesT <- read.metharray.exp(path_to_cases_tor, recursive = T)
rgCasesM <- read.metharray.exp(path_to_cases_mon, recursive = T)

# combine cases arrays 
rgCases <- combineArrays(rgCasesT, rgCasesM)
rm(rgCasesT, rgCasesM)

# controls
rgControls <- read.metharray.exp(path_to_controls, recursive = T)

rgValid <- read.metharray.exp(path_to_valid, recursive = T)

##########
# load genomic methyl set (from controls) - you need genetic locations by probe from this object
##########
ratio_set <- readRDS('../../Data/model_data/raw_ratio_set.rda')

# get granges object
g_ranges <- as.data.frame(getLocations(ratio_set))

# get probes from rownames
g_ranges$probe <- rownames(g_ranges)

# remove ch and duplicatee
g_ranges <- g_ranges[!duplicated(g_ranges$start),]
g_ranges <- g_ranges[!grepl('ch', g_ranges$probe),]

##########
# read in clinical data
##########
clin <- read.csv('../../Data/clin_data/clinical_two.csv', stringsAsFactors = F)

# clean clinical ids
clin$ids <-  gsub('A|B|_|-', '', clin$blood_dna_malkin_lab_)


##########
# cases 
##########

# cases batch1
id_map_tor <- read.csv(paste0(path_to_cases_tor, '/SampleSheet.csv'), stringsAsFactors = F)

#cases batch2
id_map_mon <- read.csv(paste0(path_to_cases_mon, '/SampleSheet.csv'), stringsAsFactors = F)
id_map_mon$Project <- NULL

# combine id_map and id_map_other
id_map_cases <- rbind(id_map_tor, id_map_mon)
rm(id_map_tor, id_map_mon)

# clean id map
id_map_cases <- cleanIdMap(id_map_cases)


##########
# Controls batch1
##########
id_map_con <- read.csv(paste0(path_to_controls, '/SampleSheet.csv'), stringsAsFactors = F)

# clean idmap
id_map_con <- cleanIdMap(id_map_con)

##########
# valid
##########
id_map_val <- read.csv(paste0(path_to_valid, '/SampleSheet.csv'))

# homogenize valid map data with cases and controls
id_map_val <- id_map_val[, c('Sample.ID', 'Sample.Group', 'Sentrix.Barcode', 'Sample.Section',
                             'Project', 'Pool_ID', 'Sample_Well')]

# sub out '.' for '_'
colnames(id_map_val) <- gsub('.', '_', colnames(id_map_val), fixed = T)

# change 'Sample_ID' to 'Sample_Name' and 'Sentrix_Barcode' to 'Sentrix_ID'
colnames(id_map_val)[1] <- 'Sample_Name'
colnames(id_map_val)[3] <- 'Sentrix_ID'
colnames(id_map_val)[4] <- 'Sentrix_Position'
colnames(id_map_val)[5] <- 'Sample_Plate'

# clean idmap
id_map_val <- cleanIdMap(id_map_val)

##########
# remove outliers (previously determined) from rgset before normalization
##########
rgControls <- remove_outliers(rgSet = rgControls,
                              id_map = id_map_con,
                              method = 'doesnt_matter',
                              type = 'controls')

rgValid <- remove_outliers(rgSet = rgValid,
                           id_map = id_map_val,
                           method = 'doesnt_matter',
                           type = 'valid')


##########
# subset data - remove controls probes on each data set only if raw preprocessing
##########

# save.image('~/Desktop/temp_450_850.RData')
load('~/Desktop/all_new.RData')

full_pipeline <- function(rgCases, 
                          rgControls, 
                          rgValid, 
                          method,
                          survival,
                          random_forest,
                          rf_surv_fac,
                          rf_surv_con,
                          age_cutoff,
                          cg_gene_regions,
                          remove_age_cgs_lit,
                          remove_age_cgs_lm,
                          gender,
                          tech,
                          base_change,
                          exon_intron,
                          control_for_family,
                          k_folds,
                          beta_thresh,
                          controls) {
  
  
  
  # load in gene cpgs
  gene_probes <- read_csv('../../Data/all_gene_cpg_loc.csv')
  
  gene_region <- paste(cg_gene_regions, collapse = '|')
  # get probe gene region
  gene_probes <- gene_probes[grepl(gene_region, gene_probes$focal_gene_regions),]
  
  gene_probes <- as.character(gene_probes$focal_CpGs[!duplicated(gene_probes$focal_CpGs)])
  
  # dont control for gender in model if using funnorm
  # control for gender if you use raw or noob
  if (method == 'funnorm') {
    keep_gender <- 
      keep_controls <- T
      keep_snps <- F
  } else if (method == 'noob') {
    keep_gender <- F
      keep_controls <- T
      keep_snps <- F
  } else {
    keep_gender <- 
    keep_controls <-
    keep_snps <- F
  }
  
  
  # cases
  rg_cases <- subset_rg_set(rg_set = rgCases,
                            keep_gender = keep_gender,
                            keep_controls = keep_controls,
                            keep_snps = keep_snps,
                            get_island = NULL,
                            get_chr = NULL,
                            get_type = NULL,
                            gene_probes = gene_probes)
  
  # controls
  rg_controls <- subset_rg_set(rg_set = rgControls,
                               keep_gender = keep_gender,
                               keep_controls = keep_controls,
                               keep_snps = keep_snps,
                               get_island = NULL,
                               get_chr = NULL,
                               get_type = NULL,
                               gene_probes = gene_probes)
  
  # valid
  rg_valid <- subset_rg_set(rg_set = rgValid,
                            keep_gender = keep_gender,
                            keep_controls = keep_controls,
                            keep_snps = keep_snps,
                            get_island = NULL,
                            get_chr = NULL,
                            get_type = NULL,
                            gene_probes = gene_probes)
  
  
  # list to store cv results
  temp_cases <- list()
  temp_controls <- list()
  surv_results <- list()
  
  # preprocess controls and valid
  beta_cases <-  preprocessMethod(rg_cases, preprocess = method)
  beta_controls <- preprocessMethod(rg_controls, preprocess = method)
  beta_valid <- preprocessMethod(rg_valid, preprocess = method)
  
  # get controls
  beta_cases <- process_rg_set_single(beta_data = beta_cases, 
                                      id_map = id_map_cases, 
                                      clin = clin)
  # get controls
  beta_controls_mod <- process_rg_set_single(beta_data = beta_controls, 
                                             id_map = id_map_con, 
                                             clin = clin)
  
  # get valid
  beta_valid_mod <- process_rg_set_single(beta_data = beta_valid, 
                                          id_map = id_map_val, 
                                          clin = clin)
  
  if(method == 'raw'){
    # remove NAs
    beta_cases <- removeNA(beta_cases, probe_start = 12)
    beta_controls_mod <- removeNA(beta_controls_mod, probe_start = 12)
    beta_valid_mod <- removeNA(beta_valid_mod, probe_start = 12)
    # remove infinite values 
    beta_cases <- removeInf(beta_cases, probe_start = 12)
    beta_controls_mod <- removeInf(beta_controls_mod, probe_start = 12)
    beta_valid_mod <- removeInf(beta_valid_mod, probe_start = 12)
  }
  

  # get intersecting name
  intersect_names <- Reduce(intersect, list(colnames(beta_cases)[12:ncol(beta_cases)],
                                            colnames(beta_controls_mod)[12:ncol(beta_controls_mod)],
                                            colnames(beta_valid_mod)[12:ncol(beta_valid_mod)]))
  
  # get the probes that are associated with genes
  intersect_names <- intersect(intersect_names, 
                               gene_probes)
  
  # read in probes associated with age
  if (remove_age_cgs_lit) {
    age_cgs <- readRDS('../../Data/age_probes.rda')
    intersect_names <- intersect_names[!intersect_names %in% age_cgs]
  }
  
  
  # remove age
  if (remove_age_cgs_lm) {
    age_cgs_lm <- readRDS('../../Data/age_cgs_lm.rda')
    intersect_names <- intersect_names[!intersect_names %in% age_cgs_lm]
  }
  
  # get clin names
  clin_names <- colnames(beta_cases)[1:11]
  
  
  # train data 
  beta_cases <- beta_cases[, c(clin_names,
                               intersect_names)]
  
  # controls data 
  beta_controls_mod <- beta_controls_mod[, c(clin_names,
                                             intersect_names)]
  
  # train data 
  beta_valid_mod <- beta_valid_mod[, c(clin_names,
                                       intersect_names)]
  
  # get old controls (mut and no cancer from cases)
  beta_controls_mod_old <- rbind(subset(beta_cases, p53_germline == 'Mut' & 
                                          cancer_diagnosis_diagnoses == 'Unaffected'))
  
  # get model data in cases for training and test
  beta_cases <- getModData(beta_cases)
  
  # get rid of cancer samples in controls 
  beta_controls_mod <- beta_controls_mod[grepl('Unaffected', beta_controls_mod$cancer_diagnosis_diagnoses),]
  
  #subset valid - get ids from train and test
  case_ids <- beta_cases$ids
  beta_valid_mod <- beta_valid_mod[!beta_valid_mod$ids %in% case_ids,]
  
  # remove NAs 
  beta_cases <-beta_cases[complete.cases(beta_cases),]

  # combine beta cases and beta valid
  beta_cases_full <- rbind(beta_cases,
                           beta_valid_mod)
  
  # combine beta_controls and beta_controls_old
  beta_controls_full <- rbind(beta_controls_mod, 
                              beta_controls_mod_old)
  beta_controls_full <- beta_controls_full[!duplicated(beta_controls_full$ids),]
  
  # add an indicator for 450 and 850
  beta_cases_full$tech <- ifelse(grepl('^57|97', beta_cases_full$sentrix_id), 'a', 'b')
  beta_controls_full$tech <- ifelse(grepl('^57|97', beta_controls_full$sentrix_id), 'a', 'b')
  
  # recode gdna.base.change for to get base letter after mutation
  beta_cases_full$gdna.base.change <- gsub(' b', 'none', substr(beta_cases_full$gdna.base.change, 3,4))
  beta_controls_full$gdna.base.change <- gsub(' b', 'none', substr(beta_controls_full$gdna.base.change, 3,4))
  
  
  # recode exon intron for eith exon or intron
  beta_cases_full$gdna.exon.intron <- ifelse(grepl('Exon', beta_cases_full$gdna.exon.intron), 'exon',
                                             ifelse(grepl('Intron', beta_cases_full$gdna.exon.intron), 'intron', 'not_clear'))
  beta_controls_full$gdna.exon.intron <- ifelse(grepl('Exon', beta_controls_full$gdna.exon.intron), 'exon',
                                             ifelse(grepl('Intron', beta_controls_full$gdna.exon.intron), 'intron', 'not_clear'))
    
  # get a base change variable for each data set
  beta_cases_full <- cbind(as.data.frame(class.ind(beta_cases_full$gdna.base.change)), beta_cases_full)
  beta_controls_full <- cbind(as.data.frame(class.ind(beta_controls_full$gdna.base.change)), beta_controls_full)
  
  # get a exon_intron variable for each data set 
  beta_cases_full <- cbind(as.data.frame(class.ind(beta_cases_full$gdna.exon.intron)), beta_cases_full)
  beta_controls_full <- cbind(as.data.frame(class.ind(beta_controls_full$gdna.exon.intron)), beta_controls_full)
  
  # get gender variable for each data set
  beta_cases_full <- cbind(as.data.frame(class.ind(beta_cases_full$gender)), beta_cases_full)
  beta_controls_full <- cbind(as.data.frame(class.ind(beta_controls_full$gender)), beta_controls_full)

  # get tech variable for each data set
  beta_cases_full <- cbind(as.data.frame(class.ind(beta_cases_full$tech)), beta_cases_full)
  beta_controls_full <- cbind(as.data.frame(class.ind(beta_controls_full$tech)), beta_controls_full)
  
  # remove na in both
  beta_cases_full <- beta_cases_full[!is.na(beta_cases_full$age_sample_collection),]
  beta_controls_full <- beta_controls_full[!is.na(beta_controls_full$age_sample_collection),]
  
  
  if(control_for_family) {
    
    # group by family name for both data sets (subset 1:12 because it will go quicker)
    # and get counts for each family
    temp_family_cases <- beta_cases_full[, 1:100] %>%
      group_by(family_name) %>%
      summarise(counts_cases = n()) 
    
    temp_family_controls <- beta_controls_full[, 1:100] %>%
      group_by(family_name) %>%
      summarise(counts_controls = n()) 
    
    # join by faimly name and get the number of time families in cases have more counts than families in controls
    # cases has more 
    temp_all <- inner_join(temp_family_cases, temp_family_controls, by = 'family_name')
    temp_all$more_cases <- ifelse(temp_all$counts_cases > temp_all$counts_controls, TRUE, FALSE)
    
    remove_from_controls <- paste0('^', temp_all$family_name[which(temp_all$more_cases)], '$',collapse  = '|')
    remove_from_cases <- paste0('^', temp_all$family_name[which(!temp_all$more_cases)], '$',collapse = '|')
    
    beta_cases_full <- beta_cases_full[!grepl(remove_from_cases, beta_cases_full$family_name), ]
    beta_controls_full <- beta_controls_full[!grepl(remove_from_controls, beta_controls_full$family_name), ]
    
    # remove cases that are over a certain age
    # beta_cases_full <- beta_cases_full[beta_cases_full$age_sample_collection < 400,]
  
  }
  
  if (survival) {
    
    full_data <- rbind(beta_cases_full,
                       beta_controls_full)
    
    full_data <- full_data[full_data$age_sample_collection < 400,]
    
    fold_vec <- sample(1:k_folds, nrow(full_data), replace = T)
    
  } else {
    # get vector of random folds
    fold_vec <- sample(1:k_folds, nrow(beta_cases_full), replace = T)
    
  }
 
 
  # combine 
  for(i in 1:k_folds) {
    
    # get train and test index
    train_index <- !grepl(i, fold_vec)
    test_index <- !train_index
    
    # split into training and test 
    if(survival){
      full_train_cases <- full_data[train_index, ]
      full_test_cases <- full_data[test_index, ]
      
      # function to predict with all test, controls, controls old, and valid
      surv_results[[i]]  <- run_coxreg(training_dat = full_train_cases,
                                       test_dat = full_test_cases,
                                       age_cutoff = age_cutoff,
                                       random_forest = random_forest,
                                       rf_surv_fac = rf_surv_fac,
                                       rf_surv_con = rf_surv_con,
                                       gender = gender,
                                       tech = tech,
                                       base_change = base_change,
                                       exon_intron = exon_intron,
                                       intersect_names = intersect_names)
      

    } else {
      
      beta_train_cases <- beta_cases_full[train_index, ]
      beta_test_cases <- beta_cases_full[test_index, ]
      
      
      # use cases training and controls to get bumphunter features
      bh_feats <- bump_hunter(dat_1 = beta_train_cases, 
                              dat_2 = beta_controls_full, 
                              bump = 'cancer', 
                              boot_num = 5, 
                              thresh = beta_thresh,
                              g_ranges = g_ranges)
      
      
      # get feature list
      colnames(bh_feats)[1] <- 'chr'
      remove_features <- inner_join(bh_feats, g_ranges)$probe
      
      # take remove features out of colnames 
      bh_features <- intersect_names[!intersect_names %in% remove_features]
      
      # function to predict with all test, controls, controls old, and valid
      mod_result  <- run_enet_450_850(training_dat = beta_train_cases,
                                      controls_dat = beta_controls_full,
                                      test_dat = beta_test_cases,
                                      age_cutoff = age_cutoff,
                                      gender = gender,
                                      tech = tech,
                                      base_change = base_change,
                                      exon_intron = exon_intron,
                                      bh_features = bh_features)
      
    
      
      
      temp_cases[[i]] <- mod_result[[1]]
      temp_controls[[i]] <- mod_result[[2]]
      
      print(i)
    }
    
  }
  
  if(survival) {
    surv_final <- do.call(rbind, surv_results)
    return(surv_final)
    
  }
  full_cases <- do.call(rbind, temp_cases)
  full_controls <- do.call(rbind, temp_controls)

  return(list(full_cases, full_controls))
}


load('~/Desktop/all_new.RData')

source('all_functions.R')

##########
# fixed variables
##########
#1stExon   3'UTR   5'UTR    Body TSS1500  TSS200
method = 'noob'
age_cutoff = 72
cg_gene_regions <- c("Body")
random_forest = F
rf_surv_fac = F
rf_surv_con = F
survival = F
remove_age_cgs = F
remove_age_lit = T
gender = T
tech = T
base_change = F
exon_intron = F
control_for_family = T
k_folds = 5
beta_thresh = 0.1

# run full pipeline
full_results <- full_pipeline(rgCases = rgCases,
                              rgControls = rgControls,
                              rgValid = rgValid,
                              method = method,
                              survival = survival,
                              random_forest = random_forest,
                              rf_surv_fac = rf_surv_fac,
                              rf_surv_con = rf_surv_con,
                              age_cutoff = age_cutoff,
                              cg_gene_regions = cg_gene_regions,
                              remove_age_cgs_lit = remove_age_lit,
                              remove_age_cgs = remove_age_cgs,
                              gender = gender,
                              tech = tech,
                              base_change = base_change,
                              exon_intron = exon_intron,
                              control_for_family = control_for_family,
                              k_folds = k_folds,
                              beta_thresh = beta_thresh,
                              controls = control_type)

full_results <-full_results[order(full_results$test_pred, decreasing = TRUE), ]


# random forest fac
controls <- full_results %>% 
  filter(cancer_diagnosis_diagnoses == 'Unaffected')

controls <- controls[order(controls$pred_y, decreasing = T),]
# get the individuals 
 #save results
# saveRDS(full_results, paste0('../../Data/results_data/noob_survival_72.rda'))

# saveRDS(full_results, paste0('../../Data/results_data/',age_cutoff,'_',method,'_', cg_gene_regions,'_',survival ,'_', remove_age_lit ,'_', remove_age_cgs ,'_',gender, '_', tech, '_', base_change,'_',exon_intron, '_', control_for_family,'.rda'))






# get results from list 
temp_cases <- full_results[[1]]

temp_cases <- temp_cases[ , c('test_pred', 'test_label', 
                              'age_diagnosis' ,  'age_sample_collection')]
temp_cases$test_pred_label <- ifelse(temp_cases$test_pred > .5, 1, 0)

temp_cases$pred_is <- ifelse(temp_cases$test_pred_label == temp_cases$test_label, 
                             'good',
                             'bad')


temp_controls <- full_results[[2]]

temp_controls <- temp_controls[ , c('controls_age_pred', 'controls_age_label', 
                                    'age_sample_collection')]
temp_controls$controls_pred_label <- ifelse(temp_controls$controls_age_pred > .5, 1, 0)

temp_controls$pred_is <- ifelse(temp_controls$controls_pred_label == temp_controls$controls_age_label, 
                                'good',
                                'bad')

# get person identfier 
temp_controls$p_id <- rep.int(seq(1, 44, 1), 5)

# group by fold and get mean 
temp_controls_pred <- temp_controls %>%
  group_by(p_id) %>%
  summarise(mean_pred = mean(controls_age_pred, na.rm =T)) %>%
  cbind(temp_controls[1:44,])

# remove original prediction 
temp_controls_pred$controls_age_pred <- NULL
rm(temp_controls,temp_results)


##########
# examin cases with prediction objects (ROC, TRP, etc)
##########

temp_pred_cases <- prediction(temp_cases$test_pred, temp_cases$test_label)
temp_roc <- performance(temp_pred_cases, measure = 'tpr', x.measure = 'tnr')
cutoffs <- data.frame(cut=temp_roc@alpha.values[[1]], tpr=temp_roc@x.values[[1]], 
                      tnr=temp_roc@y.values[[1]])

plot(cutoffs$cut, cutoffs$tpr,type="l",col="red")
par(new=TRUE)
plot(cutoffs$cut, cutoffs$tnr,type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="")

temp_pr <- performance(temp_pred_cases, measure = 'prec', x.measure = 'rec')
temp_lc <- performance(temp_pred_cases, measure="lift")

plot(temp_roc)
plot(temp_pr)
plot(temp_lc)
caret::confusionMatrix(round(temp_cases$test_pred, 0.1), temp_cases$test_label)

