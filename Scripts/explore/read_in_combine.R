#1stExon   3'UTR   5'UTR    Body TSS1500  TSS200

##########
# get base directory for 4 batch
##########
path_to_cases_tor <- '../../Data/methyl_data/cases_toronto'
path_to_cases_mon <- '../../Data/methyl_data/cases_montreal'
path_to_controls <- '../../Data/methyl_data/controls'


# get functions
source('all_functions.R')

load('~/Desktop/temp_450.RData')
load('~/Desktop/temp_850.RData')

rm(rgControls, rgValid)

# set preprocessing method
method <- 'noob'
methyl_type <- 'beta'


##########
# load genomic methyl set (from controls) - you need genetic locations by probe from this object
##########
ratio_set <- readRDS('../../Data/g_ranges.rda')

# # get granges object
# g_ranges <- as.data.frame(getLocations(ratio_set))
# 
# # get probes from rownames
# g_ranges$probe <- rownames(g_ranges)
# 
# # remove ch and duplicatee
# g_ranges <- g_ranges[!duplicated(g_ranges$start),]
# g_ranges <- g_ranges[!grepl('ch', g_ranges$probe),]

##########
# read in clinical data
##########
clin <- read.csv('../../Data/clin_data/new_clin.csv', stringsAsFactors = F)

# clean clinical ids
clin$ids <-  gsub('A|B|_|-', '', clin$blood_dna_malkin_lab)


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



# load in gene cpgs
gene_probes <- read_csv('../../Data/all_gene_cpg_loc.csv')

gene_region <- paste('Body', collapse = '|')
# get probe gene region
gene_probes <- gene_probes[grepl(gene_region, gene_probes$focal_gene_regions),]

gene_probes <- as.character(gene_probes$focal_CpGs[!duplicated(gene_probes$focal_CpGs)])

#cases
rg_cases <- subset_rg_set(rg_set = rgCases,
                          keep_gender = TRUE,
                          keep_controls = TRUE,
                          keep_snps = FALSE,
                          get_island = NULL,
                          get_chr = NULL,
                          get_type = NULL,
                          gene_probes = gene_probes)
rm(rgCases)







# combine 450 to 850 controls (casting to 850k technology)
rg_cases_450_to_850 <- convertArray(rg_cases,
                          outType = c("IlluminaHumanMethylationEPIC"),
                          verbose = TRUE)

rg_controls_850_to_450 <- convertArray(rg_controls,
                                    outType = c("IlluminaHumanMethylation450k"),
                                    verbose = TRUE)

rg_valid_850_to_450 <- convertArray(rg_val,
                                       outType = c("IlluminaHumanMethylation450k"),
                                       verbose = TRUE)

# get overlapping probes
# get overlapping probes 
int_names <- intersect(rg_controls@NAMES, rg_cases@NAMES)

rg_cases <- rg_cases[rownames(rg_cases) %in% int_names,]
rg_controls <- rg_controls[rownames(rg_controls) %in% int_names,]
rg_val <- rg_val[rownames(rg_val) %in% int_names,]


# combine cases and controls
rg_cases_controls <- combineArrays(rg_cases, rg_controls)
rg_cases_controls_converted <- combineArrays(rg_cases_450_to_850, rg_controls)


# normalize 


# get controls
data_cases1 <- process_rg_set_single(beta_data = data_cases, 
                                     id_map = id_map_cases, 
                                     clin = clin)
# get controls
data_controls_mod <- process_rg_set_single(beta_data = data_controls, 
                                           id_map = id_map_con, 
                                           clin = clin)


