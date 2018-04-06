

# validation or combined data 
data_used <- 'new'

# get cg regions
cg_gene_regions = 'Body'

# set preprocessing method
method <- 'noob'

# set type of data, beta or m
methyl_type <- 'm'

# set data directory
data_dir <- '../../Data/'

# get data
if(paste0(data_used,'_',methyl_type, '_final', '.RData') %in% dir(data_dir)) {
  load(paste0(data_dir,paste0(data_used,'_',methyl_type, '_final', '.RData')))
} else {
  source(paste0('get_combat.R'))
}

# read in wt data
data_wt <- readRDS(paste0(data_dir,paste0(data_used,'_',methyl_type, '_wild_type', '.rda')))

# source all_functions.R to load libraries and my functions
source('all_functions.R')
##########
# run model
##########

# get controls_lfs from full_Data
controls_lfs <- full_data[full_data$cancer_diagnosis_diagnoses == 'Unaffected',]
controls_lfs <- controls_lfs[controls_lfs$a == 1,]


# use cases training and controls to get bumphunter features
bh_feats <- bump_hunter_surv(dat = controls_lfs, 
                             wild_type = data_wt,
                             boot_num = 5, 
                             thresh = 0.1,
                             g_ranges = g_ranges)

# get feature list
colnames(bh_feats)[1] <- 'chr'
these_features <- inner_join(bh_feats, g_ranges)$probe

# save different lfs features lengths
saveRDS(these_features, '../../Data/lfs_features_01.rda')

