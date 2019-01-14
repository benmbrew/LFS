
# get functions
source('all_functions.R')

method = 'swan'
##### -------- beta

# read in 850k
con_850 <- readRDS(paste0('../../Data/', method,'/controls_850_beta.rda'))
cases_850 <- readRDS(paste0('../../Data/', method,'/cases_850_beta.rda'))

# read in 450k
cases_450 <- readRDS(paste0('../../Data/', method,'/cases_450_beta.rda'))
con_450 <- readRDS(paste0('../../Data/', method,'/controls_450_beta.rda'))
cases_wt_450 <- readRDS(paste0('../../Data/', method,'/cases_wt_450_beta.rda'))
con_wt_450 <- readRDS(paste0('../../Data/', method,'/controls_wt_450_beta.rda'))

# 
# # ##### -------- m values (log2(meth/unmeth))
# 
# # read in 850k
# con_850_m <- readRDS('../../Data/controls_850.rda')
# cases_850_m <- readRDS('../../Data/cases_850.rda')
# 
# 
# # read in 450k
# cases_450_m <-readRDS('../../Data/cases_450.rda')
# con_450_m <- readRDS('../../Data/controls_450.rda')
# cases_wt_450_m  <- readRDS('../../Data/cases_wt_450.rda')
# con_wt_450_m <- readRDS('../../Data/controls_wt_450.rda')
#  
# 
# # ---------------homogenize column names by subsetting by 850k
# features_850 <- colnames(cases_850)
# cases_450 <- cases_450[, features_850]
# cases_wt_450 <- cases_wt_450[, features_850]
# con_450 <- con_450[, features_850]
# con_wt_450 <- con_wt_450[, features_850]
# 
# 
# # homogenize column names by subsetting by 850k
# features_850_m <- colnames(cases_850)
# cases_450_m <- cases_450_m[, features_850]
# cases_wt_450_m <- cases_wt_450_m[, features_850]
# con_450_m <- con_450_m[, features_850]
# con_wt_450_m <- con_wt_450_m[, features_850]

###----------------------------------------------------------------------------------
# combine cases, controls, valid and pca

# add indicator for 850 or 450
cases_450$gdna.base.change <- '450k'
cases_wt_450$gdna.base.change <- '450k'

con_450$gdna.base.change <- '450k'
con_wt_450$gdna.base.change <- '450k'
con_850$gdna.base.change <- '850k'
cases_850$gdna.base.change <- '850k'

# rename columns
names(cases_450)[11] <- 'tech'
names(cases_wt_450)[11] <- 'tech'

names(con_450)[11] <- 'tech'
names(con_wt_450)[11] <- 'tech'
names(cases_850)[11] <- 'tech'
names(con_850)[11] <- 'tech'

# add indicator for cases and controls
cases_450$gdna.exon.intron <- 'cases'
cases_wt_450$gdna.exon.intron <- 'cases'

con_450$gdna.exon.intron <- 'controls'
con_wt_450$gdna.exon.intron <- 'controls'
con_850$gdna.exon.intron <- 'controls'
cases_850$gdna.exon.intron <- 'cases'

# rename columns
names(cases_450)[10] <- 'cancer_status'
names(cases_wt_450)[10] <- 'cancer_status'

names(con_450)[10] <- 'cancer_status'
names(con_wt_450)[10] <- 'cancer_status'
names(cases_850)[10] <- 'cancer_status'
names(con_850)[10] <- 'cancer_status'

# # # add indicator for 850 or 450
# cases_450_m$gdna.base.change <- '450k'
# con_450_m$gdna.base.change <- '450k'
# con_wt_450_m$gdna.base.change <- '450k'
# con_850_m$gdna.base.change <- '850k'
# cases_850_m$gdna.base.change <- '850k'
# 
# # rename columns
# names(cases_450_m)[11] <- 'tech'
# names(con_450_m)[11] <- 'tech'
# names(con_wt_450_m)[11] <- 'tech'
# names(cases_850_m)[11] <- 'tech'
# names(con_850_m)[11] <- 'tech'
# 
# # add indicator for cases and controls
# cases_450_m$gdna.exon.intron <- 'cases'
# con_450_m$gdna.exon.intron <- 'controls'
# con_wt_450_m$gdna.exon.intron <- 'controls'
# con_850_m$gdna.exon.intron <- 'controls'
# cases_850_m$gdna.exon.intron <- 'cases'
# 
# # rename columns
# names(cases_450_m)[10] <- 'cancer_status'
# names(con_450_m)[10] <- 'cancer_status'
# names(con_wt_450_m)[10] <- 'cancer_status'
# names(cases_850_m)[10] <- 'cancer_status'
# names(con_850_m)[10] <- 'cancer_status'


# get shared validation set
shared_valid_beta <- get_shared_data(cases_450, 
                                     cases_850, 
                                     cases_or_controls = 'cases')

# get shared controls set
shared_controls_beta <- get_shared_data(con_450, 
                                        con_850, 
                                        cases_or_controls = 'controls')

# # get shared validation set
# shared_valid_m <- get_shared_data(cases_450_m, 
#                                 cases_850_m, 
#                                 cases_or_controls = 'cases')
# 
# # get shared controls set
# shared_controls_m <- get_shared_data(con_450_m, 
#                                    con_850_m, 
#                                    cases_or_controls = 'controls')

# save shared cases and controls
saveRDS(shared_valid_beta, paste0('../../Data/', method,'/shared_valid_beta.rda'))
saveRDS(shared_controls_beta, paste0('../../Data/', method,'/shared_controls_beta.rda'))

# saveRDS(shared_valid_m, 'all_data/shared_valid_m.rda')
# saveRDS(shared_controls_m, 'all_data/shared_controls_m.rda')

# cases
cases_450 <- clean_dat(cases_450, tech = '450k', cases_or_controls = 'cases', mut_or_wt = 'mut')
# cases_450_m <- clean_dat(cases_450_m, tech = '450k', cases_or_controls = 'cases', mut_or_wt = 'mut')
cases_wt_450 <- clean_dat(cases_wt_450, tech = '450k',cases_or_controls = 'cases', mut_or_wt = 'wt')
# cases_wt_450_m <- clean_dat(cases_wt_450_m, tech = '450k',cases_or_controls = 'cases', mut_or_wt = 'wt')

# get tm donor from 450
tm_donor_450 <- unique(cases_450$tm_donor)

# cases 850
cases_850 <- clean_dat(cases_850, tech = '850k',cases_or_controls = 'cases', mut_or_wt = 'mut')
# cases_850_m <- clean_dat(cases_850_m, tech = '850k',cases_or_controls = 'cases', mut_or_wt = 'mut')

# con
con_450 <- clean_dat(con_450, tech = '450k', cases_or_controls = 'con', mut_or_wt = 'mut')
# con_450_m <- clean_dat(con_450_m, tech = '450k', cases_or_controls = 'con', mut_or_wt = 'mut')
con_wt_450 <- clean_dat(con_wt_450, tech = '450k',cases_or_controls = 'con', mut_or_wt = 'wt')
# con_wt_450_m <- clean_dat(con_wt_450_m, tech = '450k',cases_or_controls = 'con', mut_or_wt = 'wt')

# con 850
con_850 <- clean_dat(con_850, tech = '850k',cases_or_controls = 'con', mut_or_wt = 'mut')
# con_850_m <- clean_dat(con_850_m, tech = '850k',cases_or_controls = 'con', mut_or_wt = 'mut')


# run PCA
get_pca(pca_data = cases_450,
        column_name = 'gender',
        show_variance = FALSE,
        pc_x = 1,
        pc_y = 2,
        main_title = 'cases_450')

# get_pca(pca_data = cases_450_m,
#         column_name = 'tech',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'case_450_m')

get_pca(pca_data = cases_850,
        column_name = 'gender',
        show_variance = FALSE,
        pc_x = 1,
        pc_y = 2,
        main_title = 'cases_850')

# get_pca(pca_data = cases_850_m,
#         column_name = 'tech',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'cases_850_m')

# run PCA
get_pca(pca_data = con_450,
        column_name = 'gender',
        show_variance = FALSE,
        pc_x = 1,
        pc_y = 2,
        main_title = 'con_450')

# get_pca(pca_data = con_450_m,
#         column_name = 'tech',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_450_m')

get_pca(pca_data = con_850,
        column_name = 'gender',
        show_variance = FALSE,
        pc_x = 1,
        pc_y = 2,
        main_title = 'con_850')

# get_pca(pca_data = con_850_m,
#         column_name = 'tech',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_850_m')

# remove outliers for cases_850: 3740, 4122
cases_850 <- cases_850[!grepl('3740|4122', cases_850$tm_donor),]
# cases_850_m <- cases_850_m[!grepl('3740|4122', cases_850_m$tm_donor),]

# remove outliers for cases_450: 4089
cases_450 <- cases_450[!grepl('4089', cases_450$tm_donor),]
# cases_450_m <- cases_450_m[!grepl('4089', cases_450_m$tm_donor),]

# remove outliers for con_450: 3847
con_450 <- con_450[!grepl('3847', con_450$tm_donor),]
# con_450_m <- con_450_m[!grepl('3847', con_450_m$tm_donor),]

# get overlapping names
cg_names <- names(cases_850)[12:ncol(cases_850)]
clin_names <- names(cases_850)[1:11]


# now subset cases_450, cases_wt_450, con_450, con_wt_450
cases_450 <- cases_450[, c(clin_names, cg_names)]
cases_wt_450 <- cases_wt_450[, c(clin_names, cg_names)]

con_450 <- con_450[, c(clin_names, cg_names)]
con_wt_450 <- con_wt_450[, c(clin_names, cg_names)]


# combine data
all_cases_beta <- rbind(cases_450,
                        cases_850)
# all_cases_m <- rbind(cases_450_m,
#                         cases_850_m)

# cases
all_con_beta <- rbind(con_450,
                        con_850)
# all_con_m <- rbind(con_450_m,
#                      con_850_m)

# wt controls
all_con_beta_wt <- rbind(con_wt_450,
                         con_850)

# all_con_m_wt <- rbind(con_wt_450_m,
#                       con_850_m)

rm(cases_450, cases_450_m, cases_850, cases_850_m,
   con_450, con_450_m, con_850, con_850_m,
   cases_wt_450, cases_wt_450_m, con_wt_450, con_wt_450_m)

###----------------------------------------------------------------------------------
# run combat to fix techology - just run batch correction on cases 450k using cases 850k, if its truly
# removing the tech bach effec due to arrya differnece, than 450k should be homogenized to 850k controls, so # we can test the data

# run combat on cases (beta and m) to obtain a batch corrected 450k cases.
all_cases_beta_combat <- run_combat(all_cases_beta, type = 'tech')
all_cases_beta_combat_sen <- run_combat(all_cases_beta_combat, type = 'sentrix_id')
all_cases_beta_combat_gen <- run_combat(all_cases_beta_combat, type = 'gender')

# all_cases_m_combat <- run_combat(all_cases_m)

# run combat on controls (beta and m) to obtain a batch corrected 850k controls.
all_con_beta_combat <- run_combat(all_con_beta, type = 'tech')
all_con_beta_combat_sen <- run_combat(all_con_beta_combat, type = 'sentrix_id')
all_con_beta_combat_gen <- run_combat(all_con_beta_combat, type = 'gender')



# all_con_m_combat <- run_combat(all_con_m)

# run combat on controls for wt and 850 
all_con_beta_wt_combat <- run_combat(all_con_beta_wt, type = 'tech')
all_con_beta_wt_combat_sen <- run_combat(all_con_beta_wt_combat, type = 'sentrix_id')
all_con_beta_wt_combat_gen <- run_combat(all_con_beta_wt_combat, type = 'gender')


# all_con_m_wt_combat <- run_combat(all_con_m_wt)


# save beta data
# cases
saveRDS(all_cases_beta, paste0('../../Data/', method,'/all_cases_beta.rda'))
saveRDS(all_cases_beta_combat, paste0('../../Data/', method,'/all_cases_beta_combat.rda'))
saveRDS(all_cases_beta_combat_sen, paste0('../../Data/', method,'/all_cases_beta_combat_sen.rda'))
saveRDS(all_cases_beta_combat_gen, paste0('../../Data/', method,'/all_cases_beta_combat_gen.rda'))


# controls
saveRDS(all_con_beta, paste0('../../Data/', method,'/all_con_beta.rda'))
saveRDS(all_con_beta_combat, paste0('../../Data/', method,'/all_con_beta_combat.rda'))
saveRDS(all_con_beta_combat_sen, paste0('../../Data/', method,'/all_con_beta_combat_sen.rda'))
saveRDS(all_con_beta_combat_gen, paste0('../../Data/', method,'/all_con_beta_combat_gen.rda'))


# controls wt
saveRDS(all_con_beta_wt, paste0('../../Data/', method,'/all_con_beta_wt.rda'))
saveRDS(all_con_beta_wt_combat, paste0('../../Data/', method,'/all_con_beta_wt_combat.rda'))
saveRDS(all_con_beta_wt_combat_sen, paste0('../../Data/', method,'/all_con_beta_wt_combat_sen.rda'))
saveRDS(all_con_beta_wt_combat_gen, paste0('../../Data/', method,'/all_con_beta_wt_combat_gen.rda'))




# save m data
# # cases
# saveRDS(all_cases_m, 'all_data/all_cases_m.rda')
# saveRDS(all_cases_m_combat, 'all_data/all_cases_m_combat.rda')
# 
# # controls
# saveRDS(all_con_m, 'all_data/all_con_m.rda')
# saveRDS(all_con_m_combat, 'all_data/all_con_m_combat.rda')
#  

# saveRDS(all_con_m_wt, 'all_data/all_con_m_wt.rda')
# saveRDS(all_con_m_wt_combat, 'all_data/all_con_m_wt_combat.rda')

# check PCs
get_pca(pca_data = all_cases_beta_combat_gen,
        column_name = 'gender',
        show_variance = FALSE,
        pc_x = 1,
        pc_y = 2,
        main_title = 'beta cancer no combat')


# check PCs
get_pca(pca_data = all_cases_beta_combat,
        column_name = 'tech',
        show_variance = FALSE,
        pc_x = 1,
        pc_y = 2,
        main_title = 'beta cancer with combat')

# # check PCs
# get_pca(pca_data = all_cases_m,
#         column_name = 'tech',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'm cancer no combat')


# check PCs
# # all_cases_m_combat$tech <- ifelse(all_cases_m_combat$tech == 'batch_1', '450k', '850k')
# get_pca(pca_data = all_cases_m_combat,
#         column_name = 'tech',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'm cancer with combat')


# check PCs
get_pca(pca_data = all_con_beta,
        column_name = 'tech',
        show_variance = FALSE,
        pc_x = 1,
        pc_y = 2,
        main_title = 'beta no cancer no combat')

# check PCs
get_pca(pca_data = all_con_beta_combat,
        column_name = 'tech',
        show_variance = FALSE,
        pc_x = 1,
        pc_y = 2,
        main_title = 'beta no cancer with combat')

# # # check PCs
# get_pca(pca_data = all_con_m,
#         column_name = 'tech',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'm no cancer no combat')
# 
# 
# 
# # check PCs
# get_pca(pca_data = all_con_m_combat,
#         column_name = 'tech',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'm no cancer with combat')

