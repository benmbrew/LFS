

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


# preprocess controls and valid
data_cases <-  preprocessMethod(rg_cases, preprocess = method, methyl_type = methyl_type)
data_controls <- preprocessMethod(rg_controls, preprocess = method, methyl_type = methyl_type)
data_valid <- preprocessMethod(rg_valid, preprocess = method, methyl_type = methyl_type)

# get controls
data_cases <- process_rg_set_single(beta_data = data_cases, 
                                    id_map = id_map_cases, 
                                    clin = clin)
# get controls
data_controls_mod <- process_rg_set_single(beta_data = data_controls, 
                                           id_map = id_map_con, 
                                           clin = clin)

# get valid
data_valid_mod <- process_rg_set_single(beta_data = data_valid, 
                                        id_map = id_map_val, 
                                        clin = clin)

############

if(method == 'raw'){
  # remove NAs
  data_cases <- removeNA(data_cases, probe_start = 12)
  data_controls_mod <- removeNA(data_controls_mod, probe_start = 12)
  data_valid_mod <- removeNA(data_valid_mod, probe_start = 12)
  # remove infinite values 
  data_cases <- removeInf(data_cases, probe_start = 12)
  data_controls_mod <- removeInf(data_controls_mod, probe_start = 12)
  data_valid_mod <- removeInf(data_valid_mod, probe_start = 12)
}


# get intersecting name
intersect_names <- Reduce(intersect, list(colnames(data_cases)[12:ncol(data_cases)],
                                          colnames(data_controls_mod)[12:ncol(data_controls_mod)],
                                          colnames(data_valid_mod)[12:ncol(data_valid_mod)]))

# get the probes that are associated with genes
intersect_names <- intersect(intersect_names, 
                             gene_probes)

# get clin names
clin_names <- colnames(data_cases)[1:11]


# train data 
data_cases <- data_cases[, c(clin_names,
                             intersect_names)]

# controls data 
data_controls_mod <- data_controls_mod[, c(clin_names,
                                           intersect_names)]

# train data 
data_valid_mod <- data_valid_mod[, c(clin_names,
                                     intersect_names)]

# get old controls (mut and no cancer from cases)
# data_controls_mod_old <- rbind(subset(data_cases, p53_germline == 'Mut' & 
#                                         cancer_diagnosis_diagnoses == 'Unaffected'))

# get model data in cases for training and test
data_cases <- getModData(data_cases)

# get rid of cancer samples in controls 
data_controls_mod <- data_controls_mod[grepl('Unaffected', data_controls_mod$cancer_diagnosis_diagnoses),]

#subset valid - get ids from train and test
case_ids <- data_cases$tm_donor_
data_valid_mod <- data_valid_mod[!data_valid_mod$tm_donor_ %in% case_ids,]

# remove NAs 
data_cases <-data_cases[complete.cases(data_cases),]


# add an indicator for 450 and 850
data_cases$tech <- ifelse(grepl('^57|97', data_cases$sentrix_id), 'a', 'b')
data_controls_mod$tech <- ifelse(grepl('^57|97', data_controls_mod$sentrix_id), 'a', 'b')
# data_controls_mod_old$tech <- ifelse(grepl('^57|97', data_controls_mod_old$sentrix_id), 'a', 'b')
data_valid_mod$tech <- ifelse(grepl('^57|97', data_valid_mod$sentrix_id), 'a', 'b')

# get gender variable for each data set
data_cases <- cbind(as.data.frame(class.ind(data_cases$gender)), data_cases)
data_controls_mod <- cbind(as.data.frame(class.ind(data_controls_mod$gender)), data_controls_mod)
# data_controls_mod_old <- cbind(as.data.frame(class.ind(data_controls_mod_old$gender)), data_controls_mod_old)
data_valid_mod <- cbind(as.data.frame(class.ind(data_valid_mod$gender)), data_valid_mod)

# get tech variable for each data set
data_cases <- cbind(as.data.frame(class.ind(data_cases$tech)), data_cases)
data_controls_mod <- cbind(as.data.frame(class.ind(data_controls_mod$tech)), data_controls_mod)
# data_controls_mod_old <- cbind(as.data.frame(class.ind(data_controls_mod_old$tech)), data_controls_mod_old)
data_valid_mod <- cbind(as.data.frame(class.ind(data_valid_mod$tech)), data_valid_mod)

# remove tech variable at end of data frame
data_cases$tech <- NULL
data_controls_mod$tech <- NULL
# data_controls_mod_old$tech <- NULL
data_valid_mod$tech <- NULL

# remove na in both
data_cases <- data_cases[!is.na(data_cases$age_sample_collection),]
data_controls_mod <- data_controls_mod[!is.na(data_controls_mod$age_sample_collection),]
# data_controls_mod_old <- data_controls_mod_old[!is.na(data_controls_mod_old$age_sample_collection),]
data_valid_mod <- data_valid_mod[!is.na(data_valid_mod$age_sample_collection),]


# remove duplicates from each one
data_cases <- data_cases[!duplicated(data_cases$tm_donor_),]
data_controls_mod <- data_controls_mod[!duplicated(data_controls_mod$tm_donor_),]
# data_controls_mod_old <- data_controls_mod_old[!duplicated(data_controls_mod_old$tm_donor_),]
data_valid_mod <- data_valid_mod[!duplicated(data_valid_mod$tm_donor_),]

# save.image(paste0(data_dir,paste0(data_used,'_',methyl_type, '_processed_temp', '.RData')))
load(paste0(data_dir,paste0(data_used,'_',methyl_type, '_processed_temp', '.RData')))

# get cancer indicator and ratio for family members
get_family_list <- get_family_cancer_old(data_cases, data_controls_mod, data_valid_mod)

# extract from list 
data_cases <- get_family_list[[1]]
data_controls_mod <- get_family_list[[2]]
data_valid_mod <- get_family_list[[3]]

# remove duplicates
data_cases <- data_cases[!duplicated(data_cases$tm_donor_),]
#subset valid - get ids from train and test
case_ids <- data_cases$tm_donor_
data_valid_mod <- data_valid_mod[!data_valid_mod$tm_donor_ %in% case_ids,]


# add a and b variable 
data_cases$b <- 0
data_controls_mod$a <- 0
data_valid_mod$a <- 0



# remove unneeded objects
rm(list=ls(pattern="^rg"))
rm(list=ls(pattern="^id"))
rm(data_valid, data_controls, ratio_set, gene_probes, get_family_list,
   age_cgs, case_ids, clin_names, data_dir, gene_region,
   intersect_names, keep_controls, keep_gender,
   keep_snps, remove_age_cgs_lm, remove_age_cgs_lit)


# save data
save.image(paste0(data_dir,paste0(data_used,'_',methyl_type, '_processed', '.RData')))

