
# this script will force read in all methylation data and get column intersection and remove outliers before normalization

##########
# load libraries
##########
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(biovizBase)
library(GEOquery)
library(IlluminaHumanMethylation450kmanifest)
library(preprocessCore)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

##########
# initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
idat_data <- paste0(methyl_data, '/raw_files')
idat_data_con <- paste0(methyl_data, '/controls')
idat_data_valid <- paste0(methyl_data, '/validation/idat_files')
model_data <- paste0(data_folder, '/model_data')
feat_data <- paste0(data_folder, '/feat_data')


# save.image('~/Desktop/temp_raw_annotation.RData')
# 
# load('~/Desktop/temp_cases_raw.RData')

##########
# source all_functions.R script
##########
source(paste0(project_folder, '/Scripts/predict_age/all_functions.R'))

##########
# fixed variables
##########
method = 'raw'
##########

##########
# read in idat 
##########

# read in cases, controls, and valids
rg_set_cases <- read.metharray.exp(idat_data)
rg_set_con <- read.metharray.exp(idat_data_con)
rg_set_val <- read.metharray.exp(idat_data_valid)

# get annotation
temp_cases <- getAnnotation(rg_set_cases)
temp_con <- getAnnotation(rg_set_con)
temp_val <- getAnnotation(rg_set_val)

# tet chr, Name, Type,NextBase, Color, Relation_to_island, Methyl27_Locit
temp_cases <- temp_cases[, c('chr', 'Name', 'Type', 'NextBase', 'Color', 'Relation_to_Island', 'Methyl27_Loci')]
temp_con <- temp_con[, c('chr', 'Name', 'Type', 'NextBase', 'Color', 'Relation_to_Island', 'Methyl27_Loci', 'Methyl450_Loci')]
temp_val <- temp_val[, c('chr', 'Name', 'Type', 'NextBase', 'Color', 'Relation_to_Island', 'Methyl27_Loci', 'Methyl450_Loci')]

##########
# explore
##########

# get chr17 probes only for 450k and 850k 
cases_17 <- rownames(temp_cases)[temp_cases$chr == 'chr17']
con_17 <- rownames(temp_con)[temp_con$chr == 'chr17']
val_17 <- rownames(temp_val)[temp_val$chr == 'chr17']

# get intersection
chr_17_int <- Reduce(intersect, list(cases_17, 
                                     con_17,
                                     val_17))

# get chr17 island 
cases_island_17 <- rownames(temp_cases)[temp_cases$chr == 'chr17' & 
                                         temp_cases$Relation_to_Island == 'Island']

con_island_17 <- rownames(temp_con)[temp_con$chr == 'chr17' & 
                                         temp_con$Relation_to_Island == 'Island']

val_island_17 <- rownames(temp_val)[temp_val$chr == 'chr17' &
                                       temp_val$Relation_to_Island == 'Island']

# get intersetion
chr_17_island_int <- Reduce(intersect, list(cases_island_17, 
                                            con_island_17,
                                            val_island_17))


# how many choromosomes represented in 3 data types 27, 450, 850
temp_con <- as.data.frame(temp_con)

case_chrome <- 
  temp_con %>%
  group_by(chr) %>%
  summarise(counts = n(),
            methyl_27 = sum(Methyl27_Loci == T),
            methyl_450 = sum(Methyl450_Loci == T))



# get probes that are in islands 
temp_cases <- temp_cases[temp_cases$Methyl27_Loci == T,]
temp_con <- temp_con[temp_con$Methyl27_Loci == T,]
temp_val <- temp_val[temp_val$Methyl27_Loci == T,]

case_island <- temp_cases$Name[temp_cases$chr == 'chr17']
con_island <- temp_con$Name[temp_cases$chr == 'chr17']
val_island <- temp_val$Name[temp_cases$chr == 'chr17']

# get intersection
island_probes <- Reduce(intersect, list(case_island,
                                        con_island,
                                        val_island))


# get probes that are in islands 
case_island <- temp_cases$Name[temp_cases$Relation_to_Island == 'Island']
con_island <- temp_con$Name[temp_con$Relation_to_Island == 'Island']
val_island <- temp_val$Name[temp_val$Relation_to_Island == 'Island']

# get intersection
island_probes <- Reduce(intersect, list(case_island,
                                        con_island,
                                        val_island))

# get probes from 27k analysis
case_27k <- temp_cases$Name[temp_cases$Methyl27_Loci == TRUE]
con_27k <- temp_con$Name[temp_con$Methyl27_Loci == TRUE]
val_27k <- temp_val$Name[temp_val$Methyl27_Loci == TRUE]

loci27k_probes <- Reduce(intersect, list(case_27k,
                                         con_27k,
                                         val_27k))




case_27k <- temp_cases[temp_cases$Methyl27_Loci == TRUE,]
con_27k <- temp_con[temp_con$Methyl27_Loci == TRUE,]
val_27k <- temp_val[temp_val$Methyl27_Loci == TRUE,]


# get probes that are in islands 
case_island <- case_27k$Name[case_27k$Relation_to_Island == 'Island']
con_island <- con_27k$Name[con_27k$Relation_to_Island == 'Island']
val_island <- val_27k$Name[val_27k$Relation_to_Island == 'Island']

# get intersection
island_probes <- Reduce(intersect, list(case_island,
                                        con_island,
                                        val_island))


# save features
saveRDS(island_probes, paste0(feat_data, paste0('/', method, '_', 'island_probes.rda')))
saveRDS(loci27k_probes, paste0(feat_data, paste0('/', method, '_', 'loci27k_probes.rda')))
saveRDS(island_probes, paste0(feat_data, paste0('/', method, '_', 'loci27k_island_probes.rda')))
saveRDS(chr_17_int, paste0(feat_data, paste0('/', method, '_', 'chr_17_int.rda')))
saveRDS(chr_17_island_int, paste0(feat_data, paste0('/', method, '_', 'chr_17_island_int.rda')))





