source('all_functions.R')

# remove WT 
remove_wild_type <- function(m_or_beta_values){
  m_or_beta_values <- m_or_beta_values[m_or_beta_values$p53_germline == 'MUT',]
  return(m_or_beta_values)
}


# set fixed variables
method = 'swan'
combat = 'combat_sen'

# save beta data
# cases

if(combat == 'normal'){
  all_cases <- readRDS(paste0('../../Data/', method,'/all_cases_beta.rda'))
  all_con <- readRDS(paste0('../../Data/', method,'/all_con_beta.rda'))
  all_con_wt <- readRDS(paste0('../../Data/', method,'/all_con_beta_wt.rda'))
  
}

if(combat == 'combat_1'){
  all_cases <- readRDS(paste0('../../Data/', method,'/all_cases_beta_combat.rda'))
  all_con <- readRDS(paste0('../../Data/', method,'/all_con_beta_combat.rda'))
  all_con_wt <- readRDS(paste0('../../Data/', method,'/all_con_beta_wt.rda'))
  
  all_cases$tech <- ifelse(all_cases$tech == 'batch_1', '450k', '850k')
  all_con$tech <- ifelse(all_con$tech == 'batch_1', '450k', '850k')
  all_con_wt$tech <- ifelse(all_con_wt$tech == 'batch_1', '450k', '850k')
}

if(combat == 'combat_gen'){
  all_cases <- readRDS(paste0('../../Data/', method,'/all_cases_beta_combat_gen.rda'))
  all_con <- readRDS(paste0('../../Data/', method,'/all_con_beta_combat_gen.rda'))
  all_con_wt <- readRDS(paste0('../../Data/', method,'/all_con_beta_wt_combat_gen.rda'))
  
  all_cases$tech <- ifelse(all_cases$tech == 'batch_1', '450k', '850k')
  all_con$tech <- ifelse(all_con$tech == 'batch_1', '450k', '850k')
  all_con_wt$tech <- ifelse(all_con_wt$tech == 'batch_1', '450k', '850k')
}

if(combat == 'combat_sen'){
  all_cases <- readRDS(paste0('../../Data/', method,'/all_cases_beta_combat_sen.rda'))
  all_con <- readRDS(paste0('../../Data/', method,'/all_con_beta_combat_sen.rda'))
  all_con_wt <- readRDS(paste0('../../Data/', method,'/all_con_beta_wt_combat_sen.rda'))
  
  all_cases$tech <- ifelse(all_cases$tech == 'batch_1', '450k', '850k')
  all_con$tech <- ifelse(all_con$tech == 'batch_1', '450k', '850k')
  all_con_wt$tech <- ifelse(all_con_wt$tech == 'batch_1', '450k', '850k')
}


# get controls 
all_con_wt <- all_con_wt[!duplicated(all_con_wt$tm_donor),]
con_wt <- all_con_wt[all_con_wt$p53_germline == 'WT',]
rm(all_con_wt)

#  remove duplcates
get_data <- function(cases, controls) {
  
  # split up by tech
  cases_450 <- cases[cases$tech == '450k',]
  cases_850 <- cases[cases$tech == '850k',]
  
  cases_450 <- cases_450[!duplicated(cases_450$tm_donor),]
  cases_850 <- cases_850[!duplicated(cases_850$tm_donor),]
  cases_450 <- remove_wild_type(cases_450)
  cases_850 <- remove_wild_type(cases_850)
  cases_450 <- cases_450[!is.na(cases_450$age_diagnosis),]
  cases_450 <- cases_450[!is.na(cases_450$age_sample_collection),]
  cases_850 <- cases_850[!is.na(cases_850$age_diagnosis),]
  cases_850 <- cases_850[!is.na(cases_850$age_sample_collection),]
  
  # controls
  con_850 <- controls[controls$tech == '850k',]
  con_450 <- controls[controls$tech == '450k',]
  
  
  con_450 <- remove_wild_type(con_450)
  con_850 <- remove_wild_type(con_850)
  con_450 <- con_450[!is.na(con_450$age_sample_collection),]
  con_850 <- con_850[!is.na(con_850$age_sample_collection),]
  
  return(list(cases_450, cases_850, con_450, con_850))
}

data_list <- get_data(all_cases, all_con)

# unlist data 
cases_450 <- data_list[[1]]
cases_850 <- data_list[[2]]
con_450 <- data_list[[3]]
con_850 <- data_list[[4]]

# convert age to months
cases_450$age_diagnosis <- 
  round(cases_450$age_diagnosis*12, 2)
cases_450$age_sample_collection <- 
  round(cases_450$age_sample_collection*12, 2)

# controls
con_850$age_diagnosis <- 
  round(con_850$age_diagnosis*12, 2)
con_850$age_sample_collection <- 
  round(con_850$age_sample_collection*12, 2)

# valie
cases_850$age_diagnosis <- 
  round(cases_850$age_diagnosis*12, 2)
cases_850$age_sample_collection <- 
  round(cases_850$age_sample_collection*12, 2)

# conmut and conmw
con_450$age_sample_collection <- 
  round(con_450$age_sample_collection*12, 2)
con_wt$age_sample_collection <- 
  round(con_wt$age_sample_collection*12, 2)

# apply to wt controls
con_wt <- con_wt[!is.na(con_wt$age_sample_collection),]
con_wt <- con_wt[!is.na(con_wt$gender),]
con_450 <- con_450[!is.na(con_450$age_sample_collection),]
con_450 <- con_450[!is.na(con_450$gender),]


rm(data_list)
# first create a function to check how the variation in age
# is accounted for by the PCs

plot_pc <- function(temp_dat, 
                    column_name, 
                    show_variance,
                    age,
                    pc_x, 
                    pc_y, 
                    main_title){
  
  # create age fac for plotting
  # create age fac for plotting
  temp_dat$age_fac <- ifelse(temp_dat$age_sample_collection > 0 &
                                temp_dat$age_sample_collection <= 12, 'first_age', 
                              ifelse(temp_dat$age_sample_collection > 12 &
                                       temp_dat$age_sample_collection <= 36, 'second_age',
                                     ifelse(temp_dat$age_sample_collection > 36 &
                                              temp_dat$age_sample_collection <= 140, 'third_age',
                                            ifelse(temp_dat$age_sample_collection > 140 &
                                                     temp_dat$age_sample_collection <= 211, 'fourth_age',
                                                   ifelse(temp_dat$age_sample_collection > 211 &
                                                            temp_dat$age_sample_collection <= 348, 'fifth_age', 'sixth_age')))))
  
  # run pca plot and return                             
  get_pca(pca_data = temp_dat, 
          column_name = column_name, 
          show_variance = show_variance, 
          pc_x = pc_x, 
          pc_y = pc_y, 
          main_title = main_title)
  
}
# pdf(paste0('~/Desktop/before_pc', '_', method, '_', combat, '.pdf'))
# 
# plot dataset PCs colored by age for PC 1, 2, 3
# plot_pc(temp_dat = cases_450,
#         column_name = 'age_fac',
#         show_variance = TRUE,
#         age= 180,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'cases_450 before pc')
# 
# plot_pc(temp_dat = cases_850,
#         column_name = 'age_fac',
#         show_variance = TRUE,
#         age= 180,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'cases_850 before pc')
# 
# plot_pc(temp_dat = con_450,
#         column_name = 'age_fac',
#         show_variance = TRUE,
#         age= 180,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_450 before pc')
# 
# plot_pc(temp_dat = con_wt,
#         column_name = 'age_fac',
#         show_variance = TRUE,
#         age= 180,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_wt before pc')
# 
# plot_pc(temp_dat = con_850,
#         column_name = 'age_fac',
#         show_variance = TRUE,
#         age= 180,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_850 before pc')
# 
# plot_pc(temp_dat = cases_450,
#         column_name = 'age_fac',
#         show_variance = FALSE,
#         age= 180,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'cases_450 before pc')
# 
# plot_pc(temp_dat = cases_850,
#         column_name = 'age_fac',
#         show_variance = FALSE,
#         age= 180,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'cases_850 before pc')
# 
# plot_pc(temp_dat = con_450,
#         column_name = 'age_fac',
#         show_variance = FALSE,
#         age= 180,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_450 before pc')
# 
# plot_pc(temp_dat = con_850,
#         column_name = 'age_fac',
#         show_variance = FALSE,
#         age= 180,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_850 before pc')
# 
# plot_pc(temp_dat = con_wt,
#         column_name = 'age_fac',
#         show_variance = FALSE,
#         age= 180,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_wt before pc')
# 


# run combat on age
cases_450 <- run_combat(cases_450, type = 'age2')
cases_850 <- run_combat(cases_850, type = 'age2')
con_450 <- run_combat(con_450, type = 'age2')
con_850 <- run_combat(con_850, type = 'age2')
con_wt <- run_combat(con_wt, type = 'age2')

# plot_pc(temp_dat = cases_450,
#         column_name = 'age_fac',
#         show_variance = TRUE,
#         age= 180,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'cases_450 after combat')
# 
# plot_pc(temp_dat = cases_850,
#         column_name = 'age_fac',
#         show_variance = TRUE,
#         age= 180,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'cases_850 after combat')
# 
# plot_pc(temp_dat = con_450,
#         column_name = 'age_fac',
#         show_variance = TRUE,
#         age= 180,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_450 after combat')
# 
# plot_pc(temp_dat = con_wt,
#         column_name = 'age_fac',
#         show_variance = TRUE,
#         age= 180,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_wt after combat')
# 
# plot_pc(temp_dat = con_850,
#         column_name = 'age_fac',
#         show_variance = TRUE,
#         age= 180,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_850 after combat')
# 
# plot_pc(temp_dat = cases_450,
#         column_name = 'age_fac',
#         show_variance = FALSE,
#         age= 180,
#         pc_x = 1,
#         pc_y = 3,
#         main_title = 'cases_450 after combat')
# 
# plot_pc(temp_dat = cases_850,
#         column_name = 'age_fac',
#         show_variance = FALSE,
#         age= 180,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'cases_850 after combat')
# 
# plot_pc(temp_dat = con_450,
#         column_name = 'age_fac',
#         show_variance = FALSE,
#         age= 180,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_450 after combat')
# 
# plot_pc(temp_dat = con_850,
#         column_name = 'age_fac',
#         show_variance = FALSE,
#         age= 180,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_850 after combat')
# 
# plot_pc(temp_dat = con_wt,
#         column_name = 'age_fac',
#         show_variance = FALSE,
#         age= 180,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_wt after combat')
# 



saveRDS(cases_450, paste0('../../Data/', method,'/cases_450_age_combat_2_', combat,'.rda'))
saveRDS(cases_850, paste0('../../Data/', method,'/cases_850_age_combat_2_', combat,'.rda'))
saveRDS(con_450, paste0('../../Data/', method,'/con_450_age_combat_2_', combat,'.rda'))
saveRDS(con_850, paste0('../../Data/', method,'/con_850_age_combat_2_', combat,'.rda'))
saveRDS(con_wt, paste0('../../Data/', method,'/con_wt_age_combat_2_', combat,'.rda'))

