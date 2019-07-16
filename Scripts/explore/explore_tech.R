
# get functions
source('../predict_age/all_functions.R')
library(RColorBrewer)
path_to_cases_tor <- '../../Data/methyl_data/cases_toronto'
path_to_cases_mon <- '../../Data/methyl_data/cases_montreal'
path_to_controls <- '../../Data/methyl_data/controls'
path_to_valid <- '../../Data/methyl_data/validation'

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

# get ids from id_maps
get_ids <- function(data) {
  column_split <- strsplit(as.character(data$sample_name), '#')
  last_digits <- lapply(column_split, function(x) x[length(x)])
  sub_ids <- unlist(last_digits)
  sub_ids <- gsub('RD-', '', sub_ids)
  data$ids <- sub_ids
  data$identifier <- NULL
  return(data)
}
id_map_cases <- get_ids(id_map_cases)
id_map_con <- get_ids(id_map_con)
id_map_val<- get_ids(id_map_val)


##### -------- beta
method = 'noob'
# read in 850k
con_850 <- readRDS(paste0('../../Data/', method,'/controls_850_beta.rda'))
cases_850 <- readRDS(paste0('../../Data/', method,'/cases_850_beta.rda'))

# read in 450k
cases_450 <-readRDS(paste0('../../Data/', method,'/cases_450_beta.rda'))
con_450 <- readRDS(paste0('../../Data/', method,'/controls_450_beta.rda'))
cases_wt_450 <- readRDS(paste0('../../Data/', method,'/cases_wt_450_beta.rda'))
con_wt_450 <- readRDS(paste0('../../Data/', method,'/controls_wt_450_beta.rda'))


# read in sample sheets to get array data (sample well, etc)

# ---------------homogenize column names by subsetting by 850k
features_850 <- colnames(cases_850)
cases_450 <- cases_450[, features_850]
cases_wt_450 <- cases_wt_450[, features_850]
con_450 <- con_450[, features_850]
con_wt_450 <- con_wt_450[, features_850]


# combine with id maps
cases_450 <- inner_join(id_map_cases, cases_450, by = 'ids')
cases_wt_450 <- inner_join(id_map_cases, cases_wt_450, by = 'ids')

cases_850 <- inner_join(id_map_val, cases_850, by = 'ids')
con_450 <- inner_join(id_map_con, con_450, by = 'ids')
con_850 <- inner_join(id_map_con, con_850, by = 'ids')
con_wt_450 <- inner_join(id_map_con, con_wt_450, by = 'ids')

# read structure data
restructure_data <- function(temp_dat){
  clin_dat <- temp_dat[, !grepl('^cg', names(temp_dat))]
  cpg_dat <- temp_dat[, grepl('^cg', names(temp_dat))]
  final_dat <- as.data.frame(cbind(clin_dat, cpg_dat))
  return(final_dat)
}

cases_450 <- restructure_data(cases_450)
cases_wt_450 <- restructure_data(cases_wt_450)
cases_850 <- restructure_data(cases_850)
con_450 <- restructure_data(con_450)
con_wt_450 <- restructure_data(con_wt_450)
con_850 <- restructure_data(con_850)

# remove sentrix id.x 
names(cases_450) <- gsub('.x', '',names(cases_450), fixed = TRUE )
names(cases_wt_450) <- gsub('.x', '',names(cases_wt_450), fixed = TRUE )
names(cases_850) <- gsub('.x', '',names(cases_850), fixed = TRUE )

names(con_450) <- gsub('.x', '',names(con_450), fixed = TRUE )
names(con_wt_450) <- gsub('.x', '',names(con_wt_450), fixed = TRUE )
names(con_850) <- gsub('.x', '',names(con_850), fixed = TRUE )

cases_450$sentrix_id.y <- NULL
cases_wt_450$sentrix_id.y <- NULL
cases_850$sentrix_id.y <- NULL

con_450$sentrix_id.y <- NULL
con_wt_450$sentrix_id.y <- NULL
con_850$sentrix_id.y <- NULL

load('~/Desktop/tempall.RData')

# run pca on all data sets, plotting different potential batch effects

##########
# CASES individual
#########


cases_450$sentrix_id <- as.character(cases_450$sentrix_id)
cases_wt_450$sentrix_id <- as.character(cases_wt_450$sentrix_id)
cases_850$sentrix_id <- as.character(cases_850$sentrix_id)
cases_450$pool_id <- ifelse(is.na(cases_450$pool_id), 'other', cases_450$pool_id)
cases_wt_450$pool_id <- ifelse(is.na(cases_wt_450$pool_id), 'other', cases_wt_450$pool_id)

# recode sample_well, by sepearing letter from numbers into two separate columns
cases_450$sample_well_letter <- substr(cases_450$sample_well, 1, 1)
cases_450$sample_well_num <- substr(cases_450$sample_well, 2, 3)

cases_wt_450$sample_well_letter <- substr(cases_wt_450$sample_well, 1, 1)
cases_wt_450$sample_well_num <- substr(cases_wt_450$sample_well, 2, 3)

cases_850$sample_well_letter <- substr(cases_850$sample_well, 1, 1)
cases_850$sample_well_num <- substr(cases_850$sample_well, 2, 3)

# recode sentrix position = RO1C0, R01C, RO1, C01, 
cases_450$sentrix_ro1 <- ifelse(grepl('^R01C0', cases_450$sentrix_position), 'R01C0',
                                   ifelse(grepl('^R02C0', cases_450$sentrix_position), 'R02C0',
                                          ifelse(grepl('^R03C0', cases_450$sentrix_position), 'R03C0',
                                                 ifelse(grepl('^R04C0', cases_450$sentrix_position), 'R04C0',
                                                        ifelse(grepl('^R05C0', cases_450$sentrix_position), 'R05C0', 'R06C0')))))

# recode sentrix position = RO1C0, R01C, RO1, C01, 
cases_450$sentrix_co1 <- ifelse(grepl('C01', cases_450$sentrix_position), 'C01','C02')

# recode sentrix position = RO1C0, R01C, RO1, C01, 
cases_wt_450$sentrix_ro1 <- ifelse(grepl('^R01C0', cases_wt_450$sentrix_position), 'R01C0',
                                ifelse(grepl('^R02C0', cases_wt_450$sentrix_position), 'R02C0',
                                       ifelse(grepl('^R03C0', cases_wt_450$sentrix_position), 'R03C0',
                                              ifelse(grepl('^R04C0', cases_wt_450$sentrix_position), 'R04C0',
                                                     ifelse(grepl('^R05C0', cases_wt_450$sentrix_position), 'R05C0', 'R06C0')))))

# recode sentrix position = RO1C0, R01C, RO1, C01, 
cases_wt_450$sentrix_co1 <- ifelse(grepl('C01', cases_wt_450$sentrix_position), 'C01','C02')

# recode sentrix position = RO1C0, R01C, RO1, C01, 
cases_850$sentrix_r3 <- ifelse(grepl('R01|R02|R03', cases_850$sentrix_position), 'R123', 
                                ifelse(grepl('R04|R05|R06', cases_850$sentrix_position), 'R356', 'R78'))
     
# recode sentrix position = RO1C0, R01C, RO1, C01, 
cases_850$sentrix_r2 <- ifelse(grepl('R01|R02|R03|R04', cases_850$sentrix_position), 'R1234','R5678') 

# restructure
cases_450 <- restructure_data(cases_450)
cases_wt_450 <- restructure_data(cases_wt_450)
cases_850 <- restructure_data(cases_850)

# cases_450 = sample_well_letter, sample_well_num, sample_plate, sample_group, pool_id (edit), sentrix_id, sentrix_position
# cases_wt_450 = sample_well_letter, sample_well_num,, sample_plate, sample_group, pool_id (edit), sentrix_id, sentrix_position
# cases_850 = sample_well_letter, sample_well_num,, sample_plate,  sentrix_id, sentrix_position
# get_pca(pca_data = cases_450, 
#         column_name = 'sample_well_letter', 
#         show_variance = TRUE, 
#         pc_x = 1, 
#         pc_y = 2, 
#         main_title = 'sample_well_letter')

# cases_450$sentrix_id <- ifelse(grepl('^576029', cases_450$sentrix_id), '576029', 
#                                ifelse(grepl('^576066', cases_450$sentrix_id), '576066', 
#                                       ifelse(grepl('^577', cases_450$sentrix_id), '577', '972')))
# 
# cases_450$sentrix_id <- ifelse(grepl('^576', cases_450$sentrix_id), '576', 
#                                ifelse(grepl('^577', cases_450$sentrix_id), '577', '972')) 


cases_850$sentrix_id <- ifelse(grepl('^200', cases_850$sentrix_id), '200', 
                               ifelse(grepl('^20100', cases_850$sentrix_id), '20100', 
                                      ifelse(grepl('^20104', cases_850$sentrix_id), '20104', cases_850$sentrix_id))) 




get_pca(pca_data = cases_850, 
        column_name = 'sentrix_r2', 
        show_variance = FALSE, 
        pc_x = 1, 
        pc_y = 2, 
        main_title = 'sentrix_r2')

get_pca(pca_data = cases_850, 
        column_name = 'sentrix_r3', 
        show_variance = FALSE, 
        pc_x = 1, 
        pc_y = 2, 
        main_title = 'sentrix_r3')

get_pca(pca_data = cases_450, 
        column_name = 'sentrix_co1', 
        show_variance = FALSE, 
        pc_x = 1, 
        pc_y = 2, 
        main_title = 'sentrix_co1')

# get_pca(pca_data = cases_wt_450, 
#         column_name = 'sample_well_letter', 
#         show_variance = TRUE, 
#         pc_x = 1, 
#         pc_y = 2, 
#         main_title = 'sample_well_letter')

get_pca(pca_data = cases_wt_450, 
        column_name = 'sample_well_letter', 
        show_variance = FALSE, 
        pc_x = 1, 
        pc_y = 2, 
        main_title = 'sample_well_letter')

# get_pca(pca_data = cases_850, 
#         column_name = 'sample_well_letter', 
#         show_variance = TRUE, 
#         pc_x = 1, 
#         pc_y = 2, 
#         main_title = 'sample_well_letter')

get_pca(pca_data = cases_850, 
        column_name = 'sentrix_id', 
        show_variance = FALSE, 
        pc_x = 2, 
        pc_y = 3, 
        main_title = 'sample_plate')

##########
# CASES combine 450 to see p53 effect
#########
all_cases_450 <- rbind(cases_450, cases_wt_450)

get_pca(pca_data = all_cases_450, 
        column_name = 'p53_germline', 
        show_variance = FALSE, 
        pc_x = 2, 
        pc_y = 3, 
        main_title = 'p53_germline')



#####
# CONTROLS
con_850$sentrix_id <- as.character(con_850$sentrix_id)


