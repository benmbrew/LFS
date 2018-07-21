# ids_450 <- con_450$ids
# ids_850 <- con_850$ids
# re <- intersect(ids_450, ids_850)
# 
# con_450_m <- con_450_m[!con_450_m$ids %in% re, ]
# con_850_m <- con_850_m[!con_850_m$ids %in% re, ]
# 
# con_450_m <- con_450_m[con_450_m$p53_germline == 'MUT',]
# con_850_m <- con_850_m[con_850_m$p53_germline == 'MUT',]
# 
# con_450_m <- con_450_m[!is.na(con_450_m$age_sample_collection),]
# con_850_m <- con_850_m[!is.na(con_850_m$age_sample_collection),]
# 
# cases_450_m <- cases_450_m[cases_450_m$p53_germline == 'MUT',]
# cases_450_m <- cases_450_m[!is.na(cases_450_m$age_sample_collection),]
# ids_450 <- cases_450$ids
# cases_850_m <- cases_850_m[!cases_850_m$ids %in% ids_450,]
# cases_850_m <- cases_850_m[!duplicated(cases_850_m$ids),]
# 
# 
# saveRDS(con_450_m, '~/Desktop/for_valli/controls_m_450k.rda')
# saveRDS(con_850_m, '~/Desktop/for_valli/controls_m_850k.rda')
# saveRDS(cases_450_m, '~/Desktop/for_valli/cases_m_450k.rda')
# saveRDS(cases_850_m, '~/Desktop/for_valli/cases_m_850k.rda')
# 
# 
# cases_450 <- cases_450[cases_450$p53_germline == 'MUT',]
# cases_450 <- cases_450[!is.na(cases_450$age_sample_collection),]
# ids_450 <- cases_450$ids
# cases_850 <- cases_850[!cases_850$ids %in% ids_450,]
# cases_850 <- cases_850[!duplicated(cases_850$ids),]
# 
# ids_450 <- con_450$ids
# ids_850 <- con_850$ids
# re <- intersect(ids_450, ids_850)
# 
# con_450 <- con_450[!con_450$ids %in% re, ]
# con_850 <- con_850[!con_850$ids %in% re, ]
# 
# con_450 <- con_450[con_450$p53_germline == 'MUT',]
# con_850 <- con_850[con_850$p53_germline == 'MUT',]
# 
# con_450 <- con_450[!is.na(con_450$age_sample_collection),]
# con_850 <- con_850[!is.na(con_850$age_sample_collection),]
# 
# saveRDS(cases_450, '~/Desktop/for_valli/cases_bmiq_beta_450k.rda')
# saveRDS(cases_850, '~/Desktop/for_valli/cases_bmiq_beta_850k.rda')
# saveRDS(con_450, '~/Desktop/for_valli/controls_bmiq_beta_450k.rda')
# saveRDS(con_850, '~/Desktop/for_valli/controls_bmiq_beta_850k.rda')


# get functions
source('all_functions.R')

##### -------- beta

# read in 850k
con_850 <- readRDS('../../Data/controls_850.rda')
cases_850 <- readRDS('../../Data/cases_850.rda')

# read in 450k 
cases_450 <-readRDS('../../Data/cases_450.rda')
con_450 <- readRDS('../../Data/controls_450.rda')
cases_wt_450 <- readRDS('../../Data/cases_wt_450.rda')
con_wt_450 <- readRDS('../../Data/controls_wt_450.rda')

##### -------- m values (log2(meth/unmeth))

# read in 850k
con_850_m <- readRDS('../../Data/controls_850_m.rda')
cases_850_m <- readRDS('../../Data/cases_850_m.rda')

# read in 450k 
cases_450_m <-readRDS('../../Data/cases_450_m.rda')
con_450_m <- readRDS('../../Data/controls_450_m.rda')
cases_wt_450_m  <- readRDS('../../Data/cases_wt_450_m.rda')
con_wt_450_m <- readRDS('../../Data/controls_wt_450_m.rda')


# ---------------homogenize column names by subsetting by 850k
features_850 <- colnames(cases_850)
cases_450 <- cases_450[, features_850]
cases_wt_450 <- cases_wt_450[, features_850]
con_450 <- con_450[, features_850]
con_wt_450 <- con_wt_450[, features_850]

# save wild type
saveRDS(con_wt_450, '../../Data/con_wt_450_beta.rda')
saveRDS(cases_wt_450, '../../Data/cases_wt_450_beta.rda')

rm(con_wt_450, cases_wt_450)

# homogenize column names by subsetting by 850k
features_850_m <- colnames(cases_850)
cases_450_m <- cases_450_m[, features_850]
cases_wt_450_m <- cases_wt_450_m[, features_850]
con_450_m <- con_450_m[, features_850]
con_wt_450_m <- con_wt_450_m[, features_850]

# save wild type
saveRDS(con_wt_450_m, '../../Data/con_wt_450_m.rda')
saveRDS(cases_wt_450_m, '../../Data/cases_wt_450_m.rda')

rm(con_wt_450_m, cases_wt_450_m)

###----------------------------------------------------------------------------------
# get pc of all data sets with lambda chart

# run PCA
get_pca(pca_data = cases_450, 
        age_cutoff = 72,
        column_name = 'gender',
        show_variance = FALSE,
        pc_x = 1,
        pc_y = 2,
        main_title = 'PC 1 and 2')

get_pca(pca_data = cases_450, 
        age_cutoff = 72,
        column_name = 'gender',
        show_variance = TRUE,
        pc_x = 1,
        pc_y = 2,
        main_title = 'PC 1 and 2')



###----------------------------------------------------------------------------------
# combine cases, controls, valid and pca

# add indicator for 850 or 450
cases_450$gdna.base.change <- '450k'
con_450$gdna.base.change <- '450k'
con_850$gdna.base.change <- '850k'
cases_850$gdna.base.change <- '850k'

# rename columns
names(cases_450)[11] <- 'tech'
names(con_450)[11] <- 'tech'
names(cases_850)[11] <- 'tech'
names(con_850)[11] <- 'tech'

# add indicator for cases and controls
cases_450$gdna.exon.intron <- 'cases'
con_450$gdna.exon.intron <- 'controls'
con_850$gdna.exon.intron <- 'controls'
cases_850$gdna.exon.intron <- 'cases'

# rename columns
names(cases_450)[10] <- 'cancer_status'
names(con_450)[10] <- 'cancer_status'
names(cases_850)[10] <- 'cancer_status'
names(con_850)[10] <- 'cancer_status'

# add indicator for 850 or 450
cases_450_m$gdna.base.change <- '450k'
con_450_m$gdna.base.change <- '450k'
con_850_m$gdna.base.change <- '850k'
cases_850_m$gdna.base.change <- '850k'

# rename columns
names(cases_450_m)[11] <- 'tech'
names(con_450_m)[11] <- 'tech'
names(cases_850_m)[11] <- 'tech'
names(con_850_m)[11] <- 'tech'

# add indicator for cases and controls
cases_450_m$gdna.exon.intron <- 'cases'
con_450_m$gdna.exon.intron <- 'controls'
con_850_m$gdna.exon.intron <- 'controls'
cases_850_m$gdna.exon.intron <- 'cases'

# rename columns
names(cases_450_m)[10] <- 'cancer_status'
names(con_450_m)[10] <- 'cancer_status'
names(cases_850_m)[10] <- 'cancer_status'
names(con_850_m)[10] <- 'cancer_status'

# create cases only and controls only
# cases
all_cases_beta <- rbind(cases_450,
                        cases_850)
all_cases_m <- rbind(cases_450_m,
                        cases_850_m)

# cases
all_con_beta <- rbind(con_450,
                        con_850)
all_con_m <- rbind(con_450_m,
                     con_850_m)

# remove outliers '2564', '3010'
all_con_m <- all_con_m[all_con_m$ids != '2564',]
all_con_m <- all_con_m[all_con_m$ids != '3010',]

all_con_beta<- all_con_beta[all_con_beta$ids != '2564',]
all_con_beta<- all_con_beta[all_con_beta$ids != '3010',]


# remove datasets
rm(cases_450,
   con_450, 
   cases_850,
   con_850,
   cases_450_m,
   con_450_m, 
   cases_850_m,
   con_850_m)

# remove duplicates from datasets 

# Cases - for tm donor make sure you keep 450k over 850k. 
all_cases_beta <- all_cases_beta[!duplicated(all_cases_beta$tm_donor),]
all_cases_m <- all_cases_m[!duplicated(all_cases_m$tm_donor),]

# Controls - for tm donor make sure you keep 450k over 850k. 
all_con_beta <- all_con_beta[!duplicated(all_con_beta$tm_donor, fromLast = TRUE),]
all_con_m <- all_con_m[!duplicated(all_con_m$tm_donor, fromLast = TRUE),]

# plot pc
# run PCA
get_pca(pca_data = all_cases_beta, 
        age_cutoff = 72,
        column_name = 'tech',
        show_variance = FALSE,
        pc_x = 1,
        pc_y = 2,
        main_title = 'PC 1 and 2 cases beta values')

get_pca(pca_data = all_cases_m, 
        age_cutoff = 72,
        column_name = 'tech',
        show_variance = FALSE,
        pc_x = 1,
        pc_y = 2,
        main_title = 'PC 1 and 2 cases m values')

# run PCA
get_pca(pca_data = all_con_beta, 
        age_cutoff = 72,
        column_name = 'tech',
        show_variance = FALSE,
        pc_x = 1,
        pc_y = 2,
        main_title = 'PC 1 and 2 controls beta values')

get_pca(pca_data = all_con_m, 
        age_cutoff = 72,
        column_name = 'tech',
        show_variance = FALSE,
        pc_x = 1,
        pc_y = 2,
        main_title = 'PC 1 and 2 controls m values')




###----------------------------------------------------------------------------------
# run combat to fix techology - just run batch correction on cases 450k using cases 850k, if its truly
# removing the tech bach effec due to arrya differnece, than 450k should be homogenized to 850k controls, so # we can test the data

# run combat on cases (beta and m) to obtain a batch corrected 450k cases.
all_cases_beta_combat <- run_combat(all_cases_beta)
all_cases_m_combat <- run_combat(all_cases_m)

# run combat on controls (beta and m) to obtain a batch corrected 850k controls.
all_con_beta_combat <- run_combat(all_con_beta)
all_con_m_combat <- run_combat(all_con_m)

# save beta data
# cases
saveRDS(all_cases_beta, '../../Data/all_cases_beta.rda')
saveRDS(all_cases_beta_combat, '../../Data/all_cases_beta_combat.rda')

# controls
saveRDS(all_con_beta, '../../Data/all_con_beta.rda')
saveRDS(all_con_beta_combat, '../../Data/all_con_beta_combat.rda')

# save m data
# cases
saveRDS(all_cases_m, '../../Data/all_cases_m.rda')
saveRDS(all_cases_m_combat, '../../Data/all_cases_m_combat.rda')

# controls
saveRDS(all_con_m, '../../Data/all_con_m.rda')
saveRDS(all_con_m_combat, '../../Data/all_con_m_combat.rda')


# check PCs
get_pca(pca_data = all_cases_m, 
        age_cutoff = 72,
        column_name = 'tech',
        show_variance = FALSE,
        pc_x = 1,
        pc_y = 2,
        main_title = 'PC 1 and 2 cases beta values without combat')
###----------------------------------------------------------------------------------
# remove first PC from analysis

