# source functions script
source('all_functions.R')

# remove WT 
remove_wild_type <- function(m_or_beta_values){
  m_or_beta_values <- m_or_beta_values[m_or_beta_values$p53_germline == 'MUT',]
  return(m_or_beta_values)
}


# set fixed variables
size = 'used_bh'
method = 'noob'
methyl_type = 'beta'

if(size =='used_bh'){
  cases_450 <- readRDS(paste0('../../Data/', method,'/cases_450_small.rda'))
  con_transform <- readRDS(paste0('../../Data/', method,'/con_transform_small.rda'))
  valid_transform <- readRDS(paste0('../../Data/', method,'/valid_transform_small.rda'))
  con_wt <- readRDS(paste0('../../Data/', method,'/con_wt_small.rda'))
  con_mut <- readRDS(paste0('../../Data/', method,'/con_mut_small.rda'))
  
} else {
  
  cases_450 <- readRDS(paste0('../../Data/', method,'/cases_450.rda'))
  con_transform <- readRDS(paste0('../../Data/', method,'/con_transform.rda'))
  valid_transform <- readRDS(paste0('../../Data/', method,'/valid_transform.rda'))
  con_wt <- readRDS(paste0('../../Data/', method,'/con_wt.rda'))
  con_mut <- readRDS(paste0('../../Data/', method,'/con_mut.rda'))
  
  # add dummy tech variable for data sets with only one, replace family_name
  names(cases_450)[9] <- 'tech'
  names(con_transform)[9] <- 'tech'
  names(valid_transform)[9] <- 'tech'
  
  # fill them with Zero
  cases_450$tech <- '450k'
  con_transform$tech <- '850k'
  valid_transform$tech <- '850k'
  
  # do the same to con_mut and con_wt
  names(con_mut)[9] <- 'tech'
  names(con_wt)[9] <- 'tech'
  
  # fill new variable with right tech indication
  con_mut$tech <- '450k'
  con_wt$tech <- '450k'
  
  # randomly sample from all cgs
  clin_names <- names(cases_450)[1:10]
  r_cgs <- sample(names(cases_450)[11:ncol(cases_450)], 3000)
  cases_450 <- cases_450[c(clin_names, r_cgs)]
  valid_transform <- valid_transform[c(clin_names, r_cgs)]
  con_transform <- con_transform[c(clin_names, r_cgs)]
  con_mut <- con_mut[c(clin_names, r_cgs)]
  con_wt <- con_wt[c(clin_names, r_cgs)]
  
  
}


plot_pc <- function(temp_dat, 
                    column_name, 
                    show_variance,
                    pc_x, 
                    pc_y, 
                    main_title){
  
  # create age fac for plotting
  temp_dat$age_fac <- ifelse(temp_dat$age_sample_collection > 0 &
                               temp_dat$age_sample_collection <= 72, 'first_age',
                             ifelse(temp_dat$age_sample_collection > 72 &
                                      temp_dat$age_sample_collection <= 216, 'second_age', 'third_age'))
  
  # run pca plot and return                             
  get_pca(pca_data = temp_dat, 
          column_name = column_name, 
          show_variance = show_variance, 
          pc_x = pc_x, 
          pc_y = pc_y, 
          main_title = main_title)
  
}

# # plot dataset PCs colored by age for PC 1, 2, 3
# plot_pc(temp_dat = cases_450,
#         column_name = 'age_fac',
#         show_variance = TRUE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'cases_450 before pc')
# 
# plot_pc(temp_dat = cases_450,
#         column_name = 'age_fac',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'cases_450 before pc')
# 
# 
# # plot dataset PCs colored by age for PC 1, 2, 3
# plot_pc(temp_dat = valid_transform,
#         column_name = 'age_fac',
#         show_variance = TRUE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'valid_transform before pc')
# 
# plot_pc(temp_dat = valid_transform,
#         column_name = 'age_fac',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'valid_transform before pc')
# 
# # plot dataset PCs colored by age for PC 1, 2, 3
# plot_pc(temp_dat = con_transform,
#         column_name = 'age_fac',
#         show_variance = TRUE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_transform before pc')
# 
# plot_pc(temp_dat = con_transform,
#         column_name = 'age_fac',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_transform before pc')
# 
# # plot dataset PCs colored by age for PC 1, 2, 3
# plot_pc(temp_dat = con_mut,
#         column_name = 'age_fac',
#         show_variance = TRUE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_mut before pc')
# 
# plot_pc(temp_dat = con_mut,
#         column_name = 'age_fac',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_mut before pc')
# 
# 
# # plot dataset PCs colored by age for PC 1, 2, 3
# plot_pc(temp_dat = con_wt,
#         column_name = 'age_fac',
#         show_variance = TRUE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_wt before pc')
# 
# plot_pc(temp_dat = con_wt,
#         column_name = 'age_fac',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_wt before pc')
# 

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


cases_450 <- remove_pc_age(cases_450, remove_leading_pcs = 'first')
valid_transform <- remove_pc_age(valid_transform, remove_leading_pcs = 'first')
con_transform <- remove_pc_age(con_transform, remove_leading_pcs = 'first')
con_mut <- remove_pc_age(con_mut, remove_leading_pcs = 'first')
con_wt <- remove_pc_age(con_wt, remove_leading_pcs = 'first')




remove_leading_pcs = 'first'
# plot dataset PCs colored by age for PC 1, 2, 3
# plot dataset PCs colored by age for PC 1, 2, 3
# plot_pc(temp_dat = cases_450,
#         column_name = 'age_fac',
#         show_variance = TRUE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'cases_450 after pc')
# 
# plot_pc(temp_dat = cases_450,
#         column_name = 'age_fac',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'cases_450 after pc')
# 
# 
# # plot dataset PCs colored by age for PC 1, 2, 3
# plot_pc(temp_dat = valid_transform,
#         column_name = 'age_fac',
#         show_variance = TRUE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'valid_transform after pc')
# 
# plot_pc(temp_dat = valid_transform,
#         column_name = 'age_fac',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'valid_transform after pc')
# 
# # plot dataset PCs colored by age for PC 1, 2, 3
# plot_pc(temp_dat = con_transform,
#         column_name = 'age_fac',
#         show_variance = TRUE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_transform after pc')
# 
# plot_pc(temp_dat = con_transform,
#         column_name = 'age_fac',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_transform after pc')
# 
# # plot dataset PCs colored by age for PC 1, 2, 3
# plot_pc(temp_dat = con_mut,
#         column_name = 'age_fac',
#         show_variance = TRUE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_mut after pc')
# 
# plot_pc(temp_dat = con_mut,
#         column_name = 'age_fac',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_mut after pc')
# 
# 
# # plot dataset PCs colored by age for PC 1, 2, 3
# plot_pc(temp_dat = con_wt,
#         column_name = 'age_fac',
#         show_variance = TRUE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_wt after pc')
# 
# plot_pc(temp_dat = con_wt,
#         column_name = 'age_fac',
#         show_variance = FALSE,
#         pc_x = 1,
#         pc_y = 2,
#         main_title = 'con_wt after pc')


if(size == 'used_bh'){
  # save trainig  set 
  saveRDS(cases_450, paste0('../../Data/', method,'/cases_450_small_transform_pc.rda'))
  
  # save validation set 
  saveRDS(con_transform, paste0('../../Data/', method,'/con_transform_small_transform_pc.rda'))
  
  # save validation set 
  saveRDS(valid_transform, paste0('../../Data/', method,'/valid_transform_small_transform_pc.rda'))
  
  # save wt controls 
  saveRDS(con_wt, paste0('../../Data/', method,'/con_wt_small_transform_pc.rda'))
  
  # save mut controls 
  saveRDS(con_mut, paste0('../../Data/', method,'/con_mut_small_transform_pc.rda'))
  
} else {
  # save trainig  set 
  saveRDS(cases_450, paste0('../../Data/', method,'/cases_450_transform_pc.rda'))
  
  # save validation set 
  saveRDS(con_transform, paste0('../../Data/', method,'/con_transform_pc.rda'))
  
  # save validation set 
  saveRDS(valid_transform, paste0('../../Data/', method,'/valid_transform_pc.rda'))
  
  # save wt controls 
  saveRDS(con_wt, paste0('../../Data/', method,'/con_wt_transform_pc.rda'))
  
  # save mut controls 
  saveRDS(con_mut, paste0('../../Data/', method,'/con_mut_transform_pc.rda'))
  
  
}



