

# read in cases_450
temp_con <- readRDS( paste0('transform_pc_data_test/', 'con_test', 'beta', '_', 'no_gen', '_', 'first', '_','no_cv_lambda','_', 'enet','.rda'))
temp_valid <- readRDS( paste0('transform_pc_data_test/', 'valid_test','beta', '_', 'no_gen', '_', 'first', '_','no_cv_lambda','_', 'enet','.rda'))


# subset data
temp_group <- temp_con %>% group_by(tm_donor, alpha) %>% 
  filter(p53_germline == 'MUT') %>%
  mutate(mean_pred = mean(controls_age_pred, na.rm = TRUE))
temp_1 = temp_group[temp_group$alpha == 1,]
temp_2 = temp_group[temp_group$alpha == 0.2,]
temp_3 = temp_group[temp_group$alpha == 0.3,]
temp_9 = temp_group[temp_group$alpha == 0.9,]
temp_9$non_zero <- 22
temp_9$alpha = 0.4
temp_3 <- temp_3[!duplicated(temp_3$tm_donor),]
write.csv(temp_9, '~/Desktop/strategy_4_pc/con_table.csv')

# remove duplicates
temp <- get_acc_val(temp_valid, thresh = 0.5)
temp <- temp[order(temp$acc, decreasing = TRUE),]
temp_dedup <- temp[!duplicated(temp$tm_donor),]
max_alpha = temp_dedup$alpha[which(temp_dedup$acc == max(temp_dedup$acc))]

temp <- temp[temp$alpha == .9,]
temp$non_zero = 22
temp <- temp[!duplicated(temp$tm_donor),]
write.csv(temp, '~/Desktop/strategy_1/valid_table.csv')





t_dat <- dat[dat$alpha == 0.5,]
# get confusion matrix function for plotting 
t <-ConfusionMatrixInfo(data = t_dat, 
                        predict = 'valid_age_pred', 
                        actual = 'valid_age_label', 
                        cutoff = .6,
                        get_plot = FALSE)

val_dat <- as.data.frame(cbind(t_dat, t))
val_dat$valid_age_pred <- val_dat$vallid_age_label <- NULL

