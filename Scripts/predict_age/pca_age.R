# this script will analyze the PCs that are associated with age and remove it.
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/
# https://cran.r-project.org/web/packages/superpc/superpc.pdf
# get functions
# get functions
# prepare pc data
source('all_functions.R')

# remove WT 
remove_wild_type <- function(m_or_beta_values){
  m_or_beta_values <- m_or_beta_values[m_or_beta_values$p53_germline == 'MUT',]
  return(m_or_beta_values)
}


# set fixed variables
method = 'swan'
combat = 'combat_sen'
remove_leading_pcs = 'first'

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

rm(data_list)
# first create a function to check how the variation in age
# is accounted for by the PCs

plot_pc <- function(temp_dat, 
                    column_name, 
                    show_variance,
                    pc_x, 
                    pc_y, 
                    main_title){
  
  # create age fac for plotting
  temp_dat$age_fac <- ifelse(temp_dat$age_sample_collection > 0 &
                               temp_dat$age_sample_collection <= 10, 'first_age',
                             ifelse(temp_dat$age_sample_collection > 10 &
                                      temp_dat$age_sample_collection <= 30, 'second_age', 'third_age'))
  
  # run pca plot and return                             
  get_pca(pca_data = temp_dat, 
          column_name = column_name, 
          show_variance = show_variance, 
          pc_x = pc_x, 
          pc_y = pc_y, 
          main_title = main_title)
  
}
pdf(paste0('~/Desktop/before_pc', '_', method, '_', combat, '.pdf'))

# plot dataset PCs colored by age for PC 1, 2, 3
plot_pc(temp_dat = cases_450,
        column_name = 'age_fac',
        show_variance = TRUE,
        pc_x = 1,
        pc_y = 2,
        main_title = 'cases_450 before pc')

plot_pc(temp_dat = cases_850,
        column_name = 'age_fac',
        show_variance = TRUE,
        pc_x = 1,
        pc_y = 2,
        main_title = 'cases_850 before pc')

plot_pc(temp_dat = con_450,
        column_name = 'age_fac',
        show_variance = TRUE,
        pc_x = 1,
        pc_y = 2,
        main_title = 'con_450 before pc')

plot_pc(temp_dat = con_850,
        column_name = 'age_fac',
        show_variance = TRUE,
        pc_x = 1,
        pc_y = 2,
        main_title = 'con_850 before pc')

plot_pc(temp_dat = cases_450,
        column_name = 'age_fac',
        show_variance = FALSE,
        pc_x = 1,
        pc_y = 2,
        main_title = 'cases_450 before pc')

plot_pc(temp_dat = cases_850,
        column_name = 'age_fac',
        show_variance = FALSE,
        pc_x = 1,
        pc_y = 2,
        main_title = 'cases_850 before pc')

plot_pc(temp_dat = con_450,
        column_name = 'age_fac',
        show_variance = FALSE,
        pc_x = 1,
        pc_y = 2,
        main_title = 'con_450 before pc')

plot_pc(temp_dat = con_850,
        column_name = 'age_fac',
        show_variance = FALSE,
        pc_x = 1,
        pc_y = 2,
        main_title = 'con_850 before pc')

dev.off()


# creat function that takes a data set and computes the PCs and removes 
# the one associated with age

remove_pc_age <- function(temp_dat, remove_leading_pcs){
  clin_dat <- temp_dat[,!grepl('^cg', names(temp_dat))]
  feat_matrix <- temp_dat[,grepl('^cg', names(temp_dat))]
  mu = colMeans(feat_matrix)
  
  Xpca <- prcomp(feat_matrix)
  
  if(remove_leading_pcs == 'first'){
    pc_start <- 2
  } else if(remove_leading_pcs == 'first_two'){
    pc_start <- 3
  } else {
    pc_start <- 6
  }
  
  nComp = nrow(temp_dat)
  Xhat = Xpca$x[,pc_start:nComp] %*% t(Xpca$rotation[,pc_start:nComp])
  Xhat = scale(Xhat, center = -mu, scale = FALSE)
  
  final_dat <- as.data.frame(cbind(clin_dat, Xhat))
  
  return(final_dat)
}

rm(all_cases, all_con)

cases_450 <- remove_pc_age(cases_450, remove_leading_pcs = remove_leading_pcs)
cases_850 <- remove_pc_age(cases_850, remove_leading_pcs = remove_leading_pcs)
con_450 <- remove_pc_age(con_450, remove_leading_pcs = remove_leading_pcs)
con_850 <- remove_pc_age(con_850, remove_leading_pcs = remove_leading_pcs)
con_wt <- remove_pc_age(con_wt, remove_leading_pcs = remove_leading_pcs)

pdf(paste0('~/Desktop/after_pc', '_', method, '_',combat, '.pdf'))


# plot dataset PCs colored by age for PC 1, 2, 3
plot_pc(temp_dat = cases_450,
        column_name = 'age_fac',
        show_variance = TRUE,
        pc_x = 1,
        pc_y = 2,
        main_title =  paste0('cases_850', '_removed_', remove_leading_pcs))

plot_pc(temp_dat = cases_850,
        column_name = 'age_fac',
        show_variance = TRUE,
        pc_x = 1,
        pc_y = 2,
        main_title =  paste0('cases_850', '_removed_', remove_leading_pcs))

plot_pc(temp_dat = con_450,
        column_name = 'age_fac',
        show_variance = TRUE,
        pc_x = 1,
        pc_y = 2,
        main_title =  paste0('con_450', '_removed_', remove_leading_pcs))

plot_pc(temp_dat = con_850,
        column_name = 'age_fac',
        show_variance = TRUE,
        pc_x = 1,
        pc_y = 2,
        main_title =  paste0('con_850', '_removed_', remove_leading_pcs))

plot_pc(temp_dat = cases_450,
        column_name = 'age_fac',
        show_variance = FALSE,
        pc_x = 1,
        pc_y = 2,
        main_title =  paste0('cases_450', '_removed_', remove_leading_pcs))

plot_pc(temp_dat = cases_850,
        column_name = 'age_fac',
        show_variance = FALSE,
        pc_x = 1,
        pc_y = 2,
        main_title =  paste0('cases_850', '_removed_', remove_leading_pcs))

plot_pc(temp_dat = con_450,
        column_name = 'age_fac',
        show_variance = FALSE,
        pc_x = 1,
        pc_y = 2,
        main_title =  paste0('con_450', '_removed_', remove_leading_pcs))

plot_pc(temp_dat = con_850,
        column_name = 'age_fac',
        show_variance = FALSE,
        pc_x = 1,
        pc_y = 2,
        main_title = paste0('con_850', '_removed_', remove_leading_pcs))

dev.off()

saveRDS(cases_450, paste0('../../Data/', method,'/cases_450_', combat,'.rda'))
saveRDS(cases_850, paste0('../../Data/', method,'/cases_850_', combat,'.rda'))
saveRDS(con_450, paste0('../../Data/', method,'/con_450_', combat,'.rda'))
saveRDS(con_850, paste0('../../Data/', method,'/con_850_', combat,'.rda'))
saveRDS(con_wt, paste0('../../Data/', method,'/con_wt_', combat,'.rda'))

