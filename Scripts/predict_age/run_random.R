####### This script will predict age of onset with raw preprocessed data
# this is part of 7th step in pipleline


##########
# Initialize folders
##########
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
data_folder <- paste0(project_folder, '/Data')
scripts_folder <- paste0(project_folder, '/Scripts')

##########
# source model_functions.R and run_models.R
##########
source(paste0(scripts_folder, '/predict_age/model_functions.R'))
source(paste0(scripts_folder, '/predict_age/run_models.R'))


#########################################################################################################################
# Random features - 100, 200, 500, 1000, 2000

###########
# function for random
###########

getRand <- function(data, rand_num, num_feat, data_name) 
{
  
  model_holder = list()
  table_holder = list()
  
  for (rand_num in 1:rand_num) {
    
    model_holder[[rand_num]] <- runModels(data, 
                                          bump_hunter = F,
                                          num_feat = num_feat,
                                          random = T,
                                          seed_num = rand_num)
    
    table_holder[[rand_num]] <- extractResults(model_holder[[rand_num]],
                                               data_name = data_name)
    
  }
  
  full_table <- do.call(rbind, table_holder)
  
  return(full_table)
}


###########
# raw
###########

# raw 100 features
raw_rand_100_table <- getRand(beta_raw, 
                              rand_num = 2, 
                              num_feat = 100, 
                              data_name = "raw_rand_100")


# group by data_name and get mean scores
raw_rand_100_table <- raw_rand_100_table %>%
                      group_by(p53_status, age, type, data) %>%
                      summarise(mean_score = mean(score))

# raw 200 features
raw_rand_200_table <- getRand(beta_raw, 
                              rand_num = 10, 
                              num_feat = 200, 
                              data_name = "raw_rand_200")


# group by data_name and get mean scores
raw_rand_200_table <- raw_rand_200_table %>%
                      group_by(p53_status, age, type, data) %>%
                      summarise(mean_score = mean(score))

# raw 500 features
raw_rand_500_table <- getRand(beta_raw, 
                              rand_num = 10, 
                              num_feat = 500, 
                              data_name = "raw_rand_500")


# group by data_name and get mean scores
raw_rand_500_table <- raw_rand_500_table %>%
                      group_by(p53_status, age, type, data) %>%
                      summarise(mean_score = mean(score))
    

# raw 1000 features
raw_rand_1000_table <- getRand(beta_raw, 
                               rand_num = 10, 
                               num_feat = 1000, 
                               data_name = "raw_rand_1000")


# group by data_name and get mean scores
raw_rand_1000_table <- raw_rand_1000_table %>%
                       group_by(p53_status, age, type, data) %>%
                       summarise(mean_score = mean(score))


###########
# rbind tables and save RDA file
###########
raw_rand_table <- rbind(raw_rand_100_table, raw_rand_200_table, 
                        raw_rand_500_table, raw_rand_1000_table)

# remove data 
rm(raw_rand_100_table, raw_rand_200_table, 
   raw_rand_500_table, raw_rand_1000_table)


#save table 
saveRDS(raw_rand_table, 
        file = paste0(rand_folder, '/raw_rand_table.rda'))

rm(list = ls(pattern = "beta_raw_*"))
rm(list = ls(pattern = "raw_*"))

###########
# swan
###########

# swan 100 features
swan_rand_100_table <- getRand(beta_swan, 
                               rand_num = 10, 
                               num_feat = 100, 
                               data_name = "swan_rand_100")


# group by data_name and get mean scores
swan_rand_100_table <- swan_rand_100_table %>%
                       group_by(p53_status, age, type, data) %>%
                       summarise(mean_score = mean(score))

# swan 200 features
swan_rand_200_table <- getRand(beta_swan, 
                               rand_num = 10, 
                               num_feat = 200, 
                               data_name = "swan_rand_200")


# group by data_name and get mean scores
swan_rand_200_table <- swan_rand_200_table %>%
                       group_by(p53_status, age, type, data) %>%
                       summarise(mean_score = mean(score))

# swan 500 features
swan_rand_500_table <- getRand(beta_swan, 
                               rand_num = 10, 
                               num_feat = 500, 
                               data_name = "swan_rand_500")


# group by data_name and get mean scores
swan_rand_500_table <- swan_rand_500_table %>%
                       group_by(p53_status, age, type, data) %>%
                       summarise(mean_score = mean(score))


# swan 1000 features
swan_rand_1000_table <- getRand(beta_swan, 
                                rand_num = 10, 
                                num_feat = 1000, 
                                data_name = "swan_rand_1000")


# group by data_name and get mean scores
swan_rand_1000_table <- swan_rand_1000_table %>%
                        group_by(p53_status, age, type, data) %>%
                        summarise(mean_score = mean(score))


###########
# rbind tables and save RDA file
###########
swan_rand_table <- rbind(swan_rand_100_table, swan_rand_200_table, 
                        swan_rand_500_table, swan_rand_1000_table)

# remove data 
rm(swan_rand_100_table, swan_rand_200_table, 
   swan_rand_500_table, swan_rand_1000_table)


#save table 
saveRDS(swan_rand_table, 
        file = paste0(rand_folder, '/swan_rand_table.rda'))


rm(list = ls(pattern = "beta_swan_*"))
rm(list = ls(pattern = "swan_*"))

###########
# quan
###########

# quan 100 features
quan_rand_100_table <- getRand(beta_quan, 
                               rand_num = 10, 
                               num_feat = 100, 
                               data_name = "quan_rand_100")


# group by data_name and get mean scores
quan_rand_100_table <- quan_rand_100_table %>%
                       group_by(p53_status, age, type, data) %>%
                       summarise(mean_score = mean(score))

# quan 200 features
quan_rand_200_table <- getRand(beta_quan, 
                               rand_num = 10, 
                               num_feat = 200, 
                               data_name = "quan_rand_200")


# group by data_name and get mean scores
quan_rand_200_table <- quan_rand_200_table %>%
                       group_by(p53_status, age, type, data) %>%
                       summarise(mean_score = mean(score))

# quan 500 features
quan_rand_500_table <- getRand(beta_quan, 
                               rand_num = 10, 
                               num_feat = 500, 
                               data_name = "quan_rand_500")


# group by data_name and get mean scores
quan_rand_500_table <- quan_rand_500_table %>%
                       group_by(p53_status, age, type, data) %>%
                       summarise(mean_score = mean(score))


# quan 1000 features
quan_rand_1000_table <- getRand(beta_quan, 
                                rand_num = 10, 
                                num_feat = 1000, 
                                data_name = "quan_rand_1000")


# group by data_name and get mean scores
quan_rand_1000_table <- quan_rand_1000_table %>%
                        group_by(p53_status, age, type, data) %>%
                        summarise(mean_score = mean(score))


###########
# rbind tables and save RDA file
###########
quan_rand_table <- rbind(quan_rand_100_table, quan_rand_200_table, 
                         quan_rand_500_table, quan_rand_1000_table)

# remove data 
rm(quan_rand_100_table, quan_rand_200_table, 
   quan_rand_500_table, quan_rand_1000_table)


#save table 
saveRDS(quan_rand_table, 
        file = paste0(rand_folder, '/quan_rand_table.rda'))


###########
# funnorm
###########

# funnorm 100 features
funnorm_rand_100_table <- getRand(beta_funnorm, 
                                  rand_num = 10, 
                                  num_feat = 100, 
                                  data_name = "funnorm_rand_100")


# group by data_name and get mean scores
funnorm_rand_100_table <- funnorm_rand_100_table %>%
                          group_by(p53_status, age, type, data) %>%
                          summarise(mean_score = mean(score))

# funnorm 200 features
funnorm_rand_200_table <- getRand(beta_funnorm, 
                                  rand_num = 10, 
                                  num_feat = 200, 
                                  data_name = "funnorm_rand_200")


# group by data_name and get mean scores
funnorm_rand_200_table <- funnorm_rand_200_table %>%
                          group_by(p53_status, age, type, data) %>%
                          summarise(mean_score = mean(score))

# funnorm 500 features
funnorm_rand_500_table <- getRand(beta_funnorm, 
                                  rand_num = 10, 
                                  num_feat = 500, 
                                  data_name = "funnorm_rand_500")


# group by data_name and get mean scores
funnorm_rand_500_table <- funnorm_rand_500_table %>%
                          group_by(p53_status, age, type, data) %>%
                          summarise(mean_score = mean(score))


# funnorm 1000 features
funnorm_rand_1000_table <- getRand(beta_funnorm, 
                                   rand_num = 10, 
                                   num_feat = 1000, 
                                   data_name = "funnorm_rand_1000")


# group by data_name and get mean scores
funnorm_rand_1000_table <- funnorm_rand_1000_table %>%
                           group_by(p53_status, age, type, data) %>%
                           summarise(mean_score = mean(score))


###########
# rbind tables and save RDA file
###########
funnorm_rand_table <- rbind(funnorm_rand_100_table, funnorm_rand_200_table, 
                            funnorm_rand_500_table, funnorm_rand_1000_table)

# remove data 
rm(funnorm_rand_100_table, funnorm_rand_200_table, 
   funnorm_rand_500_table, funnorm_rand_1000_table)


#save table 
saveRDS(funnorm_rand_table, 
        file = paste0(rand_folder, '/funnorm_rand_table.rda'))

