
# get functions
source('all_functions.R')

##### -------- beta

# read in 850k
con_850 <- readRDS('../../Data/controls_850_beta.rda')
cases_850 <- readRDS('../../Data/cases_850.rda')

# read in 450k
cases_450 <-readRDS('../../Data/cases_450_beta.rda')
con_450 <- readRDS('../../Data/controls_450_beta.rda')
cases_wt_450 <- readRDS('../../Data/cases_wt_450_beta.rda')
con_wt_450 <- readRDS('../../Data/controls_wt_450_beta.rda')
# 
# ##### -------- m values (log2(meth/unmeth))

# read in 850k
con_850_m <- readRDS('../../Data/controls_850.rda')
cases_850_m <- readRDS('../../Data/cases_850.rda')


# read in 450k
cases_450_m <-readRDS('../../Data/cases_450.rda')
con_450_m <- readRDS('../../Data/controls_450.rda')
cases_wt_450_m  <- readRDS('../../Data/cases_wt_450.rda')
con_wt_450_m <- readRDS('../../Data/controls_wt_450.rda')
# 

# ---------------homogenize column names by subsetting by 850k
features_850 <- colnames(cases_850)
cases_450 <- cases_450[, features_850]
cases_wt_450 <- cases_wt_450[, features_850]
con_450 <- con_450[, features_850]
con_wt_450 <- con_wt_450[, features_850]


# homogenize column names by subsetting by 850k
features_850_m <- colnames(cases_850)
cases_450_m <- cases_450_m[, features_850]
cases_wt_450_m <- cases_wt_450_m[, features_850]
con_450_m <- con_450_m[, features_850]
con_wt_450_m <- con_wt_450_m[, features_850]



# save indivudual (not combined data)

# # cases
# saveRDS(cases_450_m, '../../Data/cases_450_mod_m.rda')
# saveRDS(cases_450, '../../Data/cases_450_mod_beta.rda')
# 
# saveRDS(cases_850_m, '../../Data/cases_850_mod_m.rda')
# saveRDS(cases_850, '../../Data/cases_850_mod_beta.rda')
# 
# #

# ################# temp family
# temp <- cases_450_m
# temp <- temp[temp$p53_germline != 'WT',]
# temp <- temp[!is.na(temp$family_name),]
# temp <- temp[, 1:20] %>% group_by(family_name) %>% summarise(counts = n())
# temp$family_name <- gsub('family_', '', temp$family_name)
# temp$family_name <- as.character(seq(1, 71, 1))
# # temp$family_name <- as.factor(temp$family_name)
# # temp$counts <- as.numeric(temp$counts)
# # temp <- as.data.frame(temp)
# ggplot(temp, aes(reorder(family_name, -counts),
#                  counts)) +
#   geom_bar(stat = 'identity', fill = 'black', color = 'darkgrey', alpha = 0.6) +
#   xlab('Family Id') +
#   ylab('# of Family Members') +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=45, hjust=1, size = 5))
# 
# 
# ###----------------------------------------------------------------------------------
# # get pc of all data sets with lambda chart
# 
# # run PCA
# get_pca(pca_data = cases_450_m, 
#         age_cutoff = 72,
#         column_name = 'gender',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'PC 1 and 2')
# 
# get_pca(pca_data = cases_450_m, 
#         age_cutoff = 72,
#         column_name = 'gender',
#         show_variance = TRUE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'PC 1 and 2')
# 
# 

###----------------------------------------------------------------------------------
# combine cases, controls, valid and pca

# add indicator for 850 or 450
cases_450$gdna.base.change <- '450k'
con_450$gdna.base.change <- '450k'
con_wt_450$gdna.base.change <- '450k'
con_850$gdna.base.change <- '850k'
cases_850$gdna.base.change <- '850k'

# rename columns
names(cases_450)[11] <- 'tech'
names(con_450)[11] <- 'tech'
names(con_wt_450)[11] <- 'tech'
names(cases_850)[11] <- 'tech'
names(con_850)[11] <- 'tech'

# add indicator for cases and controls
cases_450$gdna.exon.intron <- 'cases'
con_450$gdna.exon.intron <- 'controls'
con_wt_450$gdna.exon.intron <- 'controls'
con_850$gdna.exon.intron <- 'controls'
cases_850$gdna.exon.intron <- 'cases'

# rename columns
names(cases_450)[10] <- 'cancer_status'
names(con_450)[10] <- 'cancer_status'
names(con_wt_450)[10] <- 'cancer_status'
names(cases_850)[10] <- 'cancer_status'
names(con_850)[10] <- 'cancer_status'

# # add indicator for 850 or 450
cases_450_m$gdna.base.change <- '450k'
con_450_m$gdna.base.change <- '450k'
con_wt_450_m$gdna.base.change <- '450k'
con_850_m$gdna.base.change <- '850k'
cases_850_m$gdna.base.change <- '850k'

# rename columns
names(cases_450_m)[11] <- 'tech'
names(con_450_m)[11] <- 'tech'
names(con_wt_450_m)[11] <- 'tech'
names(cases_850_m)[11] <- 'tech'
names(con_850_m)[11] <- 'tech'

# add indicator for cases and controls
cases_450_m$gdna.exon.intron <- 'cases'
con_450_m$gdna.exon.intron <- 'controls'
con_wt_450_m$gdna.exon.intron <- 'controls'
con_850_m$gdna.exon.intron <- 'controls'
cases_850_m$gdna.exon.intron <- 'cases'

# rename columns
names(cases_450_m)[10] <- 'cancer_status'
names(con_450_m)[10] <- 'cancer_status'
names(con_wt_450_m)[10] <- 'cancer_status'
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

# wt controls 
all_con_beta_wt <- rbind(con_wt_450,
                         con_850)

all_con_m_wt <- rbind(con_wt_450_m,
                      con_850_m)

# remove outliers '2564', '3010'
all_con_m <- all_con_m[all_con_m$ids != '2564',]
all_con_m <- all_con_m[all_con_m$ids != '3010',]

all_con_beta<- all_con_beta[all_con_beta$ids != '2564',]
all_con_beta<- all_con_beta[all_con_beta$ids != '3010',]

all_con_beta_wt <- all_con_beta_wt[all_con_beta_wt$ids != '2564',]
all_con_beta_wt <- all_con_beta_wt[all_con_beta_wt$ids != '3010',]

all_con_m_wt <- all_con_beta_wt[all_con_beta_wt$ids != '2564',]
all_con_m_wt <- all_con_m_wt[all_con_m_wt$ids != '3010',]

con_wt_450_m <- con_wt_450_m[con_wt_450_m$ids != '2564',]
con_wt_450_m <- con_wt_450_m[con_wt_450_m$ids != '3010',]

con_wt_450 <- con_wt_450[con_wt_450$ids != '2564',]
con_wt_450 <- con_wt_450[con_wt_450$ids != '3010',]

saveRDS(con_wt_450_m, '../../Data/con_wt_450_m.rda')
saveRDS(cases_wt_450_m, '../../Data/cases_wt_450_m.rda')
# #
rm(con_wt_450_m, cases_wt_450_m)

# # save wild type
saveRDS(con_wt_450, '../../Data/con_wt_450_beta.rda')
saveRDS(cases_wt_450, '../../Data/cases_wt_450_beta.rda')

rm(con_wt_450, cases_wt_450)


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

# remove case controls from control dataset 
all_con_beta <- all_con_beta[grepl('Unaffected', all_con_beta$cancer_diagnosis_diagnoses),]
all_con_m <- all_con_m[grepl('Unaffected', all_con_m$cancer_diagnosis_diagnoses),]

# Controls - for tm donor make sure you keep 450k over 850k. 
all_con_beta_wt <- all_con_beta_wt[!duplicated(all_con_beta_wt$tm_donor, fromLast = TRUE),]
all_con_m_wt <- all_con_m_wt[!duplicated(all_con_m_wt$tm_donor, fromLast = TRUE),]

# remove case controls from control dataset 
all_con_beta_wt <- all_con_beta_wt[grepl('Unaffected', all_con_beta_wt$cancer_diagnosis_diagnoses),]
all_con_m_wt <- all_con_m_wt[grepl('Unaffected', all_con_m_wt$cancer_diagnosis_diagnoses),]



# remove WT 
remove_wild_type <- function(m_or_beta_values){
  m_or_beta_values <- m_or_beta_values[m_or_beta_values$p53_germline == 'MUT',]
  return(m_or_beta_values)
}
all_cases_beta <- remove_wild_type(all_cases_beta)
all_cases_m <- remove_wild_type(all_cases_m)

all_con_beta <- remove_wild_type(all_con_beta)
all_con_m <- remove_wild_type(all_con_m)

# # plot pc
# # run PCA
# get_pca(pca_data = all_cases_m, 
#         age_cutoff = 72,
#         column_name = 'tech',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'PC 1 and 2')
# 
# # get_pca(pca_data = all_cases_m, 
# #         age_cutoff = 72,
# #         column_name = 'tech',
# #         show_variance = FALSE,
# #         pc_x = 1,
# #         pc_y = 2,
# #         main_title = 'PC 1 and 2 cases m values')
# 
# # run PCA
# get_pca(pca_data = all_con_m, 
#         age_cutoff = 72,
#         column_name = 'tech',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'PC 1 and 2 controls')
# 
# # get_pca(pca_data = all_con_m, 
# #         age_cutoff = 72,
# #         column_name = 'tech',
# #         show_variance = FALSE,
# #         pc_x = 1,
# #         pc_y = 2,
# #         main_title = 'PC 1 and 2 controls m values')


# # # save wild type

###----------------------------------------------------------------------------------
# run combat to fix techology - just run batch correction on cases 450k using cases 850k, if its truly
# removing the tech bach effec due to arrya differnece, than 450k should be homogenized to 850k controls, so # we can test the data

# run combat on cases (beta and m) to obtain a batch corrected 450k cases.
all_cases_beta_combat <- run_combat(all_cases_beta)
all_cases_m_combat <- run_combat(all_cases_m)

# run combat on controls (beta and m) to obtain a batch corrected 850k controls.
all_con_beta_combat <- run_combat(all_con_beta)
all_con_m_combat <- run_combat(all_con_m)

# run combat on controls for wt and 850 
all_con_beta_wt_combat <- run_combat(all_con_beta_wt)
all_con_m_wt_combat <- run_combat(all_con_m_wt)

# subset each data set by age diagnosis if cases or age of sample collection if controls- all 8 datasets
all_cases_beta <- all_cases_beta[!is.na(all_cases_beta$age_diagnosis),]
all_cases_beta_combat <- all_cases_beta_combat[!is.na(all_cases_beta_combat$age_diagnosis),]

all_cases_m <- all_cases_m[!is.na(all_cases_m$age_diagnosis),]
all_cases_m_combat <- all_cases_m_combat[!is.na(all_cases_m_combat$age_diagnosis),]

all_con_beta <- all_con_beta[!is.na(all_con_beta$age_sample_collection),]
all_con_beta_combat <- all_con_beta_combat[!is.na(all_con_beta_combat$age_sample_collection),]

all_con_m <- all_con_m[!is.na(all_con_m$age_sample_collection),]
all_con_m_combat <- all_con_m_combat[!is.na(all_con_m_combat$age_sample_collection),]

all_con_beta_wt <- all_con_beta_wt[!is.na(all_con_beta_wt$age_sample_collection),]
all_con_beta_wt_combat <- all_con_beta_wt_combat[!is.na(all_con_beta_wt_combat$age_sample_collection),]

all_con_m_wt <- all_con_m_wt[!is.na(all_con_m_wt$age_sample_collection),]
all_con_m_wt_combat <- all_con_m_wt_combat[!is.na(all_con_m_wt_combat$age_sample_collection),]



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
 
# controls wt
saveRDS(all_con_beta_wt, '../../Data/all_con_beta_wt.rda')
saveRDS(all_con_beta_wt_combat, '../../Data/all_con_beta_wt_combat.rda')
saveRDS(all_con_m_wt, '../../Data/all_con_m_wt.rda')
saveRDS(all_con_m_wt_combat, '../../Data/all_con_m_wt_combat.rda')

# # check PCs
# get_pca(pca_data = all_cases_beta, 
#         age_cutoff = 72,
#         column_name = 'tech',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'PC 1 and 2 cases beta no combat')
# 
# 
# # check PCs
# get_pca(pca_data = all_cases_beta_combat, 
#         age_cutoff = 72,
#         column_name = 'tech',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'PC 1 and 2 cases beta combat')
# 
# # check PCs
# get_pca(pca_data = all_cases_m,
#         age_cutoff = 72,
#         column_name = 'tech',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'PC 1 and 2 cases')
# 
# 
# # check PCs
# # all_cases_m_combat$tech <- ifelse(all_cases_m_combat$tech == 'batch_1', '450k', '850k')
# get_pca(pca_data = all_cases_m_combat,
#         age_cutoff = 72,
#         column_name = 'tech',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'PC 1 and 2')
# 
# 
# # check PCs
# get_pca(pca_data = all_con_beta, 
#         age_cutoff = 72,
#         column_name = 'tech',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'PC 1 and 2 con beta no combat')
# 
# 
# 
# # check PCs
# get_pca(pca_data = all_con_beta_combat, 
#         age_cutoff = 72,
#         column_name = 'tech',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'PC 1 and 2 con beta combat')
# 
# # # check PCs
# get_pca(pca_data = all_con_m,
#         age_cutoff = 72,
#         column_name = 'tech',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'PC 1 and 2 con m no combat')
# 
# 
# 
# # check PCs
# get_pca(pca_data = all_con_m_combat,
#         age_cutoff = 72,
#         column_name = 'tech',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'PC 1 and 2 con m combat')
# 
