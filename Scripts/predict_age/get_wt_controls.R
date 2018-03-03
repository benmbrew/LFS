# get WT controls from 450k data

# validation or combined data 
data_used <- 'new'

# get cg regions
cg_gene_regions = 'Body'

# set preprocessing method
method <- 'noob'

# set type of data, beta or m
methyl_type <- 'm'

# source all_functions.R to load libraries and my functions
source('all_functions.R')

# set data directory
data_dir <- '../../Data/'

##########
# get base directory for 4 batch
##########
path_to_cases_tor <- '../../Data/methyl_data/cases_toronto'
path_to_cases_mon <- '../../Data/methyl_data/cases_montreal'

##########
# read in meth array - Data/methyl_data/cases_toronto, cases_montreal, controls, validation
##########

# cases 
rgCasesT <- read.metharray.exp(path_to_cases_tor, recursive = T)
rgCasesM <- read.metharray.exp(path_to_cases_mon, recursive = T)

# combine cases arrays 
rgCases <- combineArrays(rgCasesT, rgCasesM)
rm(rgCasesT, rgCasesM)


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


# load in gene cpgs
gene_probes <- read_csv('../../Data/all_gene_cpg_loc.csv')

gene_region <- paste('Body', collapse = '|')
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

# preprocess controls and valid
data_cases <-  preprocessMethod(rg_cases, preprocess = method, methyl_type = methyl_type)

# get controls
data_cases <- process_rg_set_single(beta_data = data_cases, 
                                    id_map = id_map_cases, 
                                    clin = clin)

# if raw data remove NA and Inf data
if(method == 'raw'){
  # remove NAs
  data_cases <- removeNA(data_cases, probe_start = 12)
  
  # remove infinite values 
  data_cases <- removeInf(data_cases, probe_start = 12)

}

# get wt 
data_wt <- data_cases[data_cases$p53_germline =='WT',]

# remove duplicates 
data_wt <- data_wt[!duplicated(data_wt$tm_donor_),]

# add an indicator for 450 and 850
data_wt$tech <- ifelse(grepl('^57|97', data_wt$sentrix_id), 'a', 'b')

# get gender variable for each data set
data_wt <- cbind(as.data.frame(class.ind(data_wt$gender)), data_wt)

# get tech variable for each data set
data_wt <- cbind(as.data.frame(class.ind(data_wt$tech)), data_wt)
# remove tech variable at end of data frame
data_wt$tech <- NULL

# savae data
saveRDS(data_wt, paste0(data_dir,paste0(data_used,'_',methyl_type, '_wild_type', '.rda')))
############

