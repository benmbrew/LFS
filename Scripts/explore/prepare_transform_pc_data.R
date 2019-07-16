# source functions script
source('all_functions.R')

# create fixed objects to model and pipeline inputs and saving  

# set fixed variables
data_type = 'beta'
remove_leading_pcs = 'first'


# read in all data
cases_450 <- readRDS(paste0('transform_data_cv/', 'cases_450_',data_type, '.rda'))
con_transform <- readRDS( paste0('transform_data_cv/', 'con_transform_',data_type, '.rda'))
valid_transform <- readRDS(paste0('transform_data_cv/', 'valid_transform_',data_type,'.rda'))
con_wt <-  readRDS(paste0('transform_data_cv/', 'con_wt_',data_type, '.rda'))
con_mut <- readRDS( paste0('transform_data_cv/', 'con_mut_',data_type,'.rda'))
lfs_bump_probes <- readRDS(paste0('transform_data_cv/', 'lfs_bumps_', data_type, '_', '.rda'))


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
# pdf(paste0('~/Desktop/before_pc', '_', data_type, '_',used_combat, '.pdf'))


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


cases_450 <- remove_pc_age(cases_450, remove_leading_pcs = remove_leading_pcs)
valid_transform <- remove_pc_age(valid_transform, remove_leading_pcs = remove_leading_pcs)
con_450 <- remove_pc_age(con_mut, remove_leading_pcs = remove_leading_pcs)
con_transform <- remove_pc_age(con_transform, remove_leading_pcs = remove_leading_pcs)
con_wt <- remove_pc_age(con_wt, remove_leading_pcs = remove_leading_pcs)




# cases 450
saveRDS(cases_450, paste0('transform_pc_data/cases_450_pc_', 
                          remove_leading_pcs,
                          '_',
                          data_type,  '.rda'))

saveRDS(valid_transform, paste0('transform_pc_data/cases_850_pc_', 
                          remove_leading_pcs,
                          '_',
                          data_type, '.rda'))

# con 450
saveRDS(con_mut, paste0('transform_pc_data/con_450_pc_', 
                        remove_leading_pcs,
                        '_',
                        data_type,  '.rda'))

saveRDS(con_transform, paste0('transform_pc_data/con_850_pc_', 
                        remove_leading_pcs,
                        '_',
                        data_type,'.rda'))

saveRDS(con_wt, paste0('transform_pc_data/con_wt_pc_', 
                       remove_leading_pcs,
                       '_',
                       data_type,  '.rda'))
