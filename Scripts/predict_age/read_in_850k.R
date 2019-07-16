# This script will read in 850k methylation data and preprocess it

# get path to validation data 
path_to_controls <- '../../Data/methyl_data/controls'
path_to_valid <- '../../Data/methyl_data/validation'

# set preprocessing method and methylation data type 
method <- 'quan'
methyl_type <- 'beta'

# source all functions
source('all_functions.R')

# read in 450k data to get probe names in order to subset 850k data by those probe names
cases <- readRDS(paste0('../../Data/', method,'/cases_450_beta.rda'))

# get probe names
features_450 <- colnames(cases)[12:ncol(cases)]
rm(cases)
##########
# read in 850k controls (LFS no cancer) and cases (LFS cancer)
##########

# controls
rgControls <- read.metharray.exp(path_to_controls, recursive = T)

# valid
rgValid <- read.metharray.exp(path_to_valid, recursive = T)

# read in clinical data
clin <- read.csv('../../Data/clin_data/new_clin.csv', stringsAsFactors = F)

# clean clinical ids
clin$ids <-  gsub('A|B|_|-', '', clin$blood_dna_malkin_lab)

# remove "ids" column because it will be created later.
clin$ids <- NULL

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

gene_region <- paste('Body', collapse = '|')
# get probe gene region
gene_probes <- gene_probes[grepl(gene_region, gene_probes$focal_gene_regions),]

gene_probes <- as.character(gene_probes$focal_CpGs[!duplicated(gene_probes$focal_CpGs)])


# controls
rg_controls <- subset_rg_set(rg_set = rgControls,
                            keep_gender = TRUE,
                            keep_controls = TRUE,
                            keep_snps = FALSE,
                            get_island = NULL,
                            get_chr = NULL,
                            get_type = NULL,
                            gene_probes = gene_probes)

# validation
rg_val <- subset_rg_set(rg_set = rgValid,
                             keep_gender = TRUE,
                             keep_controls = TRUE,
                             keep_snps = FALSE,
                             get_island = NULL,
                             get_chr = NULL,
                             get_type = NULL,
                             gene_probes = gene_probes)


load('~/Desktop/temp_850.RData')
method = 'illumina'

rm(rgControls, rgValid)

# preprocess controls and valid
data_controls <-  preprocessMethod(rg_controls, preprocess = method, methyl_type = 'beta')
data_valid <-  preprocessMethod(rg_val, preprocess = method, methyl_type = 'beta')


# subset by feature_450
data_controls <- data_controls[rownames(data_controls) %in% features_450,]
data_valid <- data_valid[rownames(data_valid) %in% features_450,]

# remove large objects to speed up computation
rm(g_ranges, ratio_set, rgControls, rgValid, rg_controls, rg_val)

# get controls and valid mapped to clinical data 
data_controls <- process_rg_set_single(beta_data = data_controls, 
                                    id_map = id_map_con, 
                                    clin = clin)

data_valid <- process_rg_set_single(beta_data = data_valid, 
                                       id_map = id_map_val, 
                                       clin = clin)
# save data
saveRDS(data_controls, paste0('../../Data/', method,'/controls_850_beta.rda'))
saveRDS(data_valid, paste0('../../Data/', method,'/cases_850_beta.rda'))

