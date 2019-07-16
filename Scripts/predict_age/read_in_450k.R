#1stExon   3'UTR   5'UTR    Body TSS1500  TSS200

##########
# get base directory for 4 batch
##########
path_to_cases_tor <- '../../Data/methyl_data/cases_toronto'
path_to_cases_mon <- '../../Data/methyl_data/cases_montreal'

# set preprocessing method
method <- 'illumina'
methyl_type <- 'beta'

# get functions
source('all_functions.R')

##########
# read in meth arrayd - Data/methyl_data/cases_toronto, cases_montreal, controls, validation
##########

# cases 
rgCasesT <- read.metharray.exp(path_to_cases_tor, recursive = T)
rgCasesM <- read.metharray.exp(path_to_cases_mon, recursive = T)

# combine cases arrays 
rgCases <- combineArrays(rgCasesT, rgCasesM)
rm(rgCasesT, rgCasesM)

save.image('~/Desktop/temp_450.RData')

##########
# load genomic methyl set (from controls) - you need genetic locations by probe from this object
# ##########
# ratio_set <- readRDS('../../Data/model_data/raw_ratio_set.rda')
# 
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

# remove "ids" column because it will be created later.
clin$ids <- NULL


temp <- clin[grepl('4749|4324|3714', clin$tm_donor),]

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

# load in gene cpgs
gene_probes <- read_csv('../../Data/all_gene_cpg_loc.csv')

gene_region <- paste('Body', collapse = '|')
# get probe gene region
gene_probes <- gene_probes[grepl(gene_region, gene_probes$focal_gene_regions),]

gene_probes <- as.character(gene_probes$focal_CpGs[!duplicated(gene_probes$focal_CpGs)])

# remove outliers
rgCases <- remove_outliers(rgSet = rgCases,
                              id_map = id_map_cases,
                              method = 'doesnt_matter',
                              type = 'cases')
# cases
rg_cases <- subset_rg_set(rg_set = rgCases,
                          keep_gender = TRUE,
                          keep_controls = TRUE,
                          keep_snps = FALSE,
                          get_island = NULL,
                          get_chr = NULL,
                          get_type = NULL,
                          gene_probes = gene_probes)
rm(rgCases)


# preprocess controls and valid
data_cases <-  preprocessMethod(rg_cases, preprocess = method, methyl_type = 'beta')

# temp <- GenomicRatioSet(data_cases)
# get cases mapped to clinical data 
data_cases <- process_rg_set_single(beta_data = data_cases, 
                                    id_map = id_map_cases, 
                                    clin = clin)

# get cases 450 and controls 450
data_cases_450 <- data_cases[!data_cases$cancer_diagnosis_diagnoses %in% 'Unaffected',]
data_controls_450 <- data_cases[data_cases$cancer_diagnosis_diagnoses %in% 'Unaffected',]

# get WT 
data_wt_cases_450 <- data_cases_450[data_cases_450$p53_germline == 'WT',]
data_wt_controls_450 <- data_controls_450[data_controls_450$p53_germline == 'WT',]


# save all four data sets: 1) LFS cancer, LFS controls, WT cancer, WT controls
saveRDS(data_cases_450, paste0('../../Data/', method,'/cases_450_beta.rda'))
saveRDS(data_controls_450, paste0('../../Data/', method,'/controls_450_beta.rda'))
saveRDS(data_wt_cases_450, paste0('../../Data/', method,'/cases_wt_450_beta.rda'))
saveRDS(data_controls_450,paste0('../../Data/', method,'/controls_wt_450_beta.rda'))
