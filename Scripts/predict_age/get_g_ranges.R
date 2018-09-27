# This script will read in raw methylation data from 450k and get a ratio set

#1stExon   3'UTR   5'UTR    Body TSS1500  TSS200

##########
# get base directory for 4 batch
##########
path_to_cases_tor <- '../../Data/methyl_data/cases_toronto'
path_to_cases_mon <- '../../Data/methyl_data/cases_montreal'

# set preprocessing method
method <- 'noob'
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

# map to genome to get ratio set
ratio_set <- mapToGenome(rgCases)

# get genomic locations
g_ranges <- as.data.frame(getLocations(ratio_set))

# save g_ranges
saveRDS(g_ranges, '../../Data/g_ranges.rda')
