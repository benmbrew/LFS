#1stExon   3'UTR   5'UTR    Body TSS1500  TSS200

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
clin <- read.csv('../../Data/clin_data/new_clin.csv', stringsAsFactors = F)

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
# rgControls <- remove_outliers(rgSet = rgControls,
#                               id_map = id_map_con,
#                               method = 'doesnt_matter',
#                               type = 'controls')
# 
# rgValid <- remove_outliers(rgSet = rgValid,
#                            id_map = id_map_val,
#                            method = 'doesnt_matter',
#                            type = 'valid')

# load in gene cpgs
gene_probes <- read_csv('../../Data/all_gene_cpg_loc.csv')

gene_region <- paste('Body', collapse = '|')
# get probe gene region
gene_probes <- gene_probes[grepl(gene_region, gene_probes$focal_gene_regions),]

gene_probes <- as.character(gene_probes$focal_CpGs[!duplicated(gene_probes$focal_CpGs)])

# dont control for gender in model if using funnorm
# control for gender if you use raw or noob
if (method == 'funnorm') {
  keep_gender <- T
    keep_controls <- T
  keep_snps <- T
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
# HERE
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
data_controls_mod_old <- rbind(subset(data_cases, p53_germline == 'MUT' &
                                        cancer_diagnosis_diagnoses == 'Unaffected'))

#HERE
# # get the data here - combine when done running 
# full_data <- rbind(data_cases,
#                    data_controls_mod,
#                    data_valid_mod)
# 
# saveRDS(full_data, paste0(data_dir, 'all_data.rda'))

# data_cases_all <- data_cases
# saveRDS(data_cases_all, paste0(data_dir, 'cases_all.rda'))
# get model data in cases for training and test
data_cases <- getModData(data_cases)

# get rid of cancer samples in controls 
data_controls_mod <- data_controls_mod[grepl('Unaffected', data_controls_mod$cancer_diagnosis_diagnoses),]

# remove NAs 
data_cases <-data_cases[complete.cases(data_cases),]

# combine beta cases and beta valid
data_cases_full <- rbind(data_cases,
                         data_valid_mod)

# combine data_controls and data_controls_old
data_controls_full <- rbind(data_controls_mod, 
                            data_controls_mod_old)

# data_controls_full <- data_controls_full[!duplicated(data_controls_full$ids),]

# add an indicator for 450 and 850
data_cases_full$tech <- ifelse(grepl('^57|97', data_cases_full$sentrix_id), 'a', 'b')
data_controls_full$tech <- ifelse(grepl('^57|97', data_controls_full$sentrix_id), 'a', 'b')

# get gender variable for each data set
data_cases_full <- cbind(as.data.frame(class.ind(data_cases_full$gender)), data_cases_full)

data_controls_full <- cbind(as.data.frame(class.ind(data_controls_full$gender)), data_controls_full)

# get tech variable for each data set
data_cases_full <- cbind(as.data.frame(class.ind(data_cases_full$tech)), data_cases_full)
data_controls_full <- cbind(as.data.frame(class.ind(data_controls_full$tech)), data_controls_full)
# remove tech variable at end of data frame
data_cases_full$tech <- NULL
data_controls_full$tech <- NULL

# remove na in both
data_cases_full <- data_cases_full[!is.na(data_cases_full$age_sample_collection),]

data_controls_full <- data_controls_full[!is.na(data_controls_full$age_sample_collection),]

# remove duplicates from each one
temp <- data_cases_full[, 1:30]
tm_dup <- temp$tm_donor_[duplicated(temp$tm_donor_)]
dups <- temp %>% filter(tm_donor_ %in% tm_dup) %>% arrange(tm_donor_)


full_data <- rbind(data_cases_full,
                   data_controls_full)

temp_image <- clin[clin$image_status == 'yes',]

shared_ids <- full_data$ids[full_data$ids %in% temp_image$blood_dna_malkin_lab]

full_data1 <- full_data[full_data$ids %in% shared_ids, ]
length(which(duplicated(full_data1$ids)))

data_cases_full_last <- data_cases_full[!duplicated(data_cases_full$tm_donor_, fromLast = TRUE),]
data_cases_full_first <- data_cases_full[!duplicated(data_cases_full$tm_donor_, fromLast = FALSE),]

data_controls_full <- data_controls_full[!duplicated(data_controls_full$tm_donor_),]

save.image(paste0(data_dir,paste0(data_used,'_',methyl_type, '_processed_temp_outlier_more', '.RData')))

# create indicator for number of members in family with cancer and ratio of people in family with cancer
# each individual's cancer status will be unknown, so leave them out.
get_family_list_first <- get_family_cancer(data_cases_full_first, data_controls_full)
get_family_list_last <- get_family_cancer(data_cases_full_last, data_controls_full)

data_cases_full_first <- get_family_list_first[[1]] #cases
data_cases_full_last <- get_family_list_last[[1]] #cases
data_controls_full <- get_family_list_first[[2]] #controls



# remove unneeded objects
rm(list=ls(pattern="^rg"))
rm(list=ls(pattern="^id"))
rm(data_valid, data_controls, ratio_set, get_family_list,
   data_controls_mod, data_controls_mod_old,
   data_valid_mod, data_cases, gene_probes,
   age_cgs, case_ids, clin_names, data_dir, gene_region,
   intersect_names, keep_controls, keep_gender,
   keep_snps, remove_age_cgs_lit, remove_age_cgs_lit)


# savae data
save.image(paste0(data_dir,paste0(data_used,'_',methyl_type, '_processed_beta_first_last_gen', '.RData')))
############

