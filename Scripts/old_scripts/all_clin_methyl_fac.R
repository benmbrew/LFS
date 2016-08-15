####################################################################################################3
# This script will predict a multinomal age model ibrary(dplyr)
library(doParallel)
library(randomForest)
library(caret)


registerDoParallel(1)

# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/Analyze')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')

# Read in 3 different data sets 
full_data <- read.csv(paste0(data_folder, '/full_data.csv'), stringsAsFactors = F)
full_data_cor <- read.csv(paste0(data_folder, '/full_data_cor.csv'), stringsAsFactors = F)
full_data_rf <- read.csv(paste0(data_folder, '/full_data_rf.csv'), stringsAsFactors = F)
# Load in clinical data
# remove extra columns
full_data$X <- NULL
full_data_cor$X <- NULL
full_data_rf$X <- NULL

# Create binary variables for age of diagnosis and age of sample collection on both data sets
full_data_cor$age_diagnosis_fac <- ifelse(full_data_cor$age_diagnosis <= 48, 1, 2)

full_data_cor$age_sample_fac<- ifelse(full_data_cor$age_sample_collection <= 48, 1, 2)

full_data_rf$age_diagnosis_fac <- ifelse(full_data_rf$age_diagnosis <= 48, 1, 2)

full_data_rf$age_sample_fac <- ifelse(full_data_rf$age_sample_collection <= 48, 1, 2)

# Create binomial variables for age of diagnosis and age of sample collection
full_data$age_diagnosis_multi <- ifelse(full_data$age_diagnosis <= 30, 1, 
                                        ifelse(full_data$age_diagnosis > 30 & full_data$age_diagnosis <= 50, 2,
                                               ifelse(full_data$age_diagnosis > 50 & full_data$age_diagnosis <= 300, 3, 4)))

full_data$age_sample_multi <- ifelse(full_data$age_sample_collection <= 30, 1, 
                                        ifelse(full_data$age_sample_collection > 30 & full_data$age_sample_collection <= 50, 2,
                                               ifelse(full_data$age_sample_collection > 50 & full_data$age_sample_collection <= 300, 3, 4)))

full_data_rf$age_diagnosis_multi <- ifelse(full_data_rf$age_diagnosis <= 30, 1, 
                                        ifelse(full_data_rf$age_diagnosis > 30 & full_data_rf$age_diagnosis <= 50, 2,
                                               ifelse(full_data_rf$age_diagnosis > 50 & full_data_rf$age_diagnosis <= 300, 3, 4)))

full_data_rf$age_sample_multi <- ifelse(full_data_rf$age_sample_collection <= 30, 1, 
                                     ifelse(full_data_rf$age_sample_collection > 30 & full_data_rf$age_sample_collection <= 50, 2,
                                            ifelse(full_data_rf$age_sample_collection > 50 & full_data_rf$age_sample_collection <= 300, 3, 4)))

summary(as.factor(full_data_rf$age_diagnosis_multi))
summary(as.factor(full_data_rf$age_sample_multi))

# Create multinomal variables for age of diagnosis and age of sample collection
full_data$age_diagnosis_year <- ifelse(full_data$age_diagnosis <= 24, 1, 
                                          ifelse(full_data$age_diagnosis > 24 & full_data$age_diagnosis <= 48, 2,
                                                 ifelse(full_data$age_diagnosis > 48 & full_data$age_diagnosis <= 180, 3,
                                                        ifelse(full_data$age_diagnosis > 180 & full_data$age_diagnosis <= 360, 4, 5))))

full_data$age_sample_year <- ifelse(full_data$age_diagnosis <= 24, 1, 
                                       ifelse(full_data$age_diagnosis > 24 & full_data$age_diagnosis <= 48, 2,
                                              ifelse(full_data$age_diagnosis > 48 & full_data$age_diagnosis <= 180, 3,
                                                     ifelse(full_data$age_diagnosis > 180 & full_data$age_diagnosis <= 360, 4, 5))))

full_data_rf$age_diagnosis_year <- ifelse(full_data_rf$age_diagnosis <= 24, 1, 
                                        ifelse(full_data_rf$age_diagnosis > 24 & full_data_rf$age_diagnosis <= 48, 2,
                                               ifelse(full_data_rf$age_diagnosis > 48 & full_data_rf$age_diagnosis <= 180, 3,
                                                    ifelse(full_data_rf$age_diagnosis > 180 & full_data_rf$age_diagnosis <= 360, 4, 5))))

full_data_rf$age_sample_year <- ifelse(full_data_rf$age_diagnosis <= 24, 1, 
                                          ifelse(full_data_rf$age_diagnosis > 24 & full_data_rf$age_diagnosis <= 48, 2,
                                                 ifelse(full_data_rf$age_diagnosis > 48 & full_data_rf$age_diagnosis <= 180, 3,
                                                        ifelse(full_data_rf$age_diagnosis > 180 & full_data_rf$age_diagnosis <= 360, 4, 5))))

summary(as.factor(full_data_rf$age_diagnosis_year))
summary(as.factor(full_data_rf$age_sample_year))

# load in residual methyl
methyl_resid_cor <- read.csv(paste0(data_folder, '/methyl_resid_cor.csv'), stringsAsFactors = F)
methyl_resid_rf <- read.csv(paste0(data_folder, '/methyl_resid_rf.csv'), stringsAsFactors = F)

methyl_resid_cor$X <- NULL
methyl_resid_rf$X <- NULL

# Random Forest - this is training and testing on clinical data using k fold cross validation
predictAll <- function(data,
                       subset, 
                       selected_features,
                       binary,
                       multi,
                       resid,
                       iterations) {
  
  model <- list()
  predictions <- list()
  test.ground_truth <- list()
  rf.test_acc <- list()
  rf.test_stats <- list()

  if(resid){
    
    
    data <- data[, subset]
    
    # remove rows where age of diagnosis missinf
    data <- data[complete.cases(data),]
    data <- cbind(data, methyl_resid_cor)
    
  } else{
    
    genes <- colnames(data)[27:ncol(data)]
    
    data <- data[, c(subset, genes)]
    
    # remove rows where age of diagnosis missinf
    data <- data[!(is.na(data$age_diagnos)),]
  }
  
  obs <- nrow(data)
  # # convert characters to factors 
  # for ( i in 1:ncol(data)){
  #   
  #   if(typeof(data[,i]) == 'character' || grepl('num', names(data[i]))) {
  #     data[,i] <- as.factor(data[,i])
  #     print(i)
  #   } 
  # }
  # 
  for (i in 1:iterations){
    
    set.seed(i)
    temp_levels <- 1
    loop_count <- i
    
    variable <- 'age_diagnosis_fac'
    
    if (binary){
      level_count <- 2
    } else if (multi) {
      level_count <- 4
    } else {
      level_count <- 5
    }
    # 
    while (length(temp_levels) < level_count) {
      set.seed(loop_count)
      train_index <- sample(nrow(data), nrow(data) *.7, replace = F)
      temp_levels <- levels(as.factor(data[, variable][-train_index]))
      # if(length(temp_levels) < level_count) {
      # loop_count <- i + 1
      # }
      loop_count <- loop_count + 1
    }
    
    rf_y = make.names(as.factor(data[, variable])[train_index])

    # 4) Random Forest 
    
    if (length(levels(rf_y)) == 2) {
      summaryFunc <- twoClassSummary
    } else {
      summaryFunc <- multiClassSummary
    }
    
    # determines how you train the model.
    fitControl <- trainControl( 
      method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
      number = 2, 
      classProbs = TRUE,     
      repeats = 1,
      allowParallel = TRUE,
      summaryFunction = summaryFunc)
    
    mtry <- sqrt(ncol(data))
    tunegrid <- expand.grid(.mtry=mtry)
    
    model[[i]] <- train(x = data[train_index, 3:ncol(data)]
                        , y = rf_y
                        , method = "rf"
                        , trControl = fitControl                   
                        , verbose = FALSE
                        , metric = "logLoss")
    
    predictions[[i]] <- predict(model[[i]] 
                                , newdata = data[-train_index, 3:ncol(data)]
                                , type = "prob")
    
    test.ground_truth[[i]] <- as.factor(data[, variable][-train_index])
    test.ground_truth_1inN <- as.factor(class.ind(data[, variable][-train_index]))
    #print(table(levels(test.ground_truth)[apply(rf.predictions, 1, which.is.max)], test.ground_truth))
    # AUC
    # #create ROCR prediction object
    # temp.predict <- prediction(unlist(predictions[[i]]), test.ground_truth_1inN)  
    # 
    # if(auc){
    #   rf.test_auc <- unlist(slot(performance(temp.predict, "auc"), "y.values"))  
    #   print(paste("RF test AUC:", rf.test_auc))  
    # }
    
    # Accuracy
    rf.test_acc[[i]] <- sum(levels(test.ground_truth[[i]])[apply(predictions[[i]], 1, which.is.max)] == test.ground_truth[[i]]) / dim(predictions[[i]])[1]
    # print(paste("RF test acc:", rf.test_acc))  
    # Compute Confusion Matrix and Statistics
    #confusionMatrix(pred, truth)
    rf.test_stats[[i]] <- confusionMatrix(levels(test.ground_truth[[i]])[apply(predictions[[i]], 1, which.is.max)], test.ground_truth[[i]])
    # print(rf.test_stats)
    
    print(i)
    
  }
  
  return(list(predictions, test.ground_truth, rf.test_acc, model, rf.test_stats, obs))
  
}


#################################################################################################
# Just methylation
##########################
# Use Methylation
rf_methyl_fac <- predictAll(data = full_data_cor,
                        subset = c('age_diagnosis_fac', 'age_sample_fac'),
                        selected_features = NULL,
                        binary = TRUE,
                        resid = T,
                        iterations = 10)

mean(unlist(rf_methyl_fac[[3]]))

##############
# confusion matrix
# look at confustion matrices
temp <- list()
for (i in 1:10){
  temp[[i]] <- rf_methyl_fac[[5]][[i]]$table
}
mat <- unlist(temp)
new_mat <- matrix(, 2, 2)
spot_1<- (mat[1] + mat[5] + mat[9] + mat[13] + mat[17] + mat[21] +mat[25] + mat[29] + mat[33] + mat[37])/10
spot_2 <- (mat[2] + mat[6] + mat[10] + mat[14] + mat[18] + mat[22] +mat[26] + mat[30] + mat[34] + mat[38])/10
spot_3<- (mat[3] + mat[7] + mat[11] + mat[15] + mat[19] + mat[23] +mat[27] + mat[31] + mat[35] + mat[39])/10
spot_4 <- (mat[4] + mat[8] + mat[12] + mat[16] + mat[20] + mat[24] +mat[28] + mat[32] + mat[36] + mat[40])/10

new_mat[1,1] <- spot_1
new_mat[2,1] <- spot_2
new_mat[1,2] <- spot_3
new_mat[2,2] <- spot_4

#########################

rf_methyl_multi <- predictAll(data = full_data_rf,
                            subset = c('age_diagnosis_multi'),
                            selected_features = NULL,
                            binary = F,
                            multi = T,
                            iterations = 10)

mean(unlist(rf_methyl_multi[[3]]))

# confusion matrices 
temp <- list()
for (i in 1:10){
  temp[[i]] <- rf_methyl_multi[[5]][[i]]$table
}

mat <- unlist(temp)
new_mat <- matrix(, 4, 4)
spot_1<- (mat[1] + mat[17] + mat[33] + mat[49] + mat[65] + mat[81] +mat[97] + mat[113] + mat[129] + 
            mat[145])/10
spot_2 <- (mat[2] + mat[18] + mat[34] + mat[50] + mat[66] + mat[82] +mat[98] + mat[114] + mat[130] + 
            mat[146])/10
spot_3<- (mat[3] + mat[19] + mat[35] + mat[51] + mat[67] + mat[83] +mat[99] + mat[115] + mat[131] + 
            mat[147])/10
spot_4<- (mat[4] + mat[20] + mat[36] + mat[52] + mat[68] + mat[84] +mat[100] + mat[116] + mat[132] + 
            mat[148])/10
spot_5<- (mat[5] + mat[21] + mat[37] + mat[53] + mat[69] + mat[85] +mat[101] + mat[117] + mat[133] + 
            mat[149])/10
spot_6<- (mat[6] + mat[22] + mat[38] + mat[54] + mat[70] + mat[86] +mat[102] + mat[118] + mat[134] + 
            mat[150])/10
spot_7<- (mat[7] + mat[23] + mat[39] + mat[55] + mat[71] + mat[87] +mat[103] + mat[119] + mat[135] + 
            mat[151])/10

spot_8<- (mat[8] + mat[24] + mat[40] + mat[56] + mat[72] + mat[88] +mat[104] + mat[120] + mat[136] + 
            mat[152])/10
spot_9<- (mat[9] + mat[25] + mat[41] + mat[57] + mat[73] + mat[89] +mat[105] + mat[121] + mat[137] + 
            mat[153])/10

spot_10<- (mat[10] + mat[26] + mat[42] + mat[58] + mat[74] + mat[90] +mat[106] + mat[122] + mat[138] + 
            mat[154])/10

spot_11<- (mat[11] + mat[27] + mat[43] + mat[59] + mat[75] + mat[91] +mat[107] + mat[123] + mat[139] + 
             mat[155])/10

spot_12<- (mat[12] + mat[28] + mat[44] + mat[60] + mat[76] + mat[92] +mat[108] + mat[124] + mat[140] + 
             mat[156])/10
spot_13<- (mat[13] + mat[29] + mat[45] + mat[61] + mat[77] + mat[93] +mat[109] + mat[125] + mat[141] + 
             mat[157])/10
spot_14<- (mat[14] + mat[30] + mat[46] + mat[62] + mat[78] + mat[94] +mat[110] + mat[126] + mat[142] + 
             mat[158])/10
spot_15<- (mat[15] + mat[31] + mat[47] + mat[63] + mat[79] + mat[95] +mat[111] + mat[127] + mat[143] + 
             mat[159])/10
spot_16<- (mat[16] + mat[32] + mat[48] + mat[64] + mat[80] + mat[96] +mat[112] + mat[128] + mat[144] + 
             mat[160])/10

new_mat[1,1] <- spot_1
new_mat[2,1] <- spot_2
new_mat[3,1] <- spot_3
new_mat[4,1] <- spot_4
new_mat[1,2] <- spot_5
new_mat[2,2] <- spot_6
new_mat[3,2] <- spot_7
new_mat[4,2] <- spot_8
new_mat[1,3] <- spot_9
new_mat[2,3] <- spot_10
new_mat[3,3] <- spot_11
new_mat[4,3] <- spot_12
new_mat[1,4] <- spot_13
new_mat[2,4] <- spot_14
new_mat[3,4] <- spot_15
new_mat[4,4] <- spot_16

##########################

rf_methyl_year <- predictAll(data = full_data_rf,
                              subset = c('age_diagnosis_year'),
                              selected_features = NULL,
                              binary = F,
                              multi = F,
                              iterations = 10)

mean(unlist(rf_methyl_year[[3]]))

# confusion matrices 
temp <- list()
for (i in 1:10){
  temp[[i]] <- rf_methyl_year[[5]][[i]]$table
}

mat <- unlist(temp)
new_mat <- matrix(, 5, 5)

spot_1<- (mat[1] + mat[26] + mat[51] + mat[76] + mat[101] + mat[126] +mat[151] + mat[176] + mat[201] + 
            mat[226])/10
spot_2<- (mat[2] + mat[27] + mat[52] + mat[77] + mat[102] + mat[127] +mat[152] + mat[177] + mat[202] + 
            mat[227])/10
spot_3<- (mat[3] + mat[28] + mat[53] + mat[78] + mat[103] + mat[128] +mat[153] + mat[178] + mat[203] + 
            mat[228])/10
spot_4<- (mat[4] + mat[29] + mat[54] + mat[79] + mat[104] + mat[129] +mat[154] + mat[179] + mat[204] + 
            mat[229])/10
spot_5<- (mat[5] + mat[30] + mat[55] + mat[80] + mat[105] + mat[130] +mat[155] + mat[180] + mat[205] + 
            mat[230])/10
spot_6<- (mat[6] + mat[31] + mat[56] + mat[81] + mat[106] + mat[131] +mat[156] + mat[181] + mat[206] + 
            mat[231])/10
spot_7<- (mat[7] + mat[32] + mat[57] + mat[82] + mat[107] + mat[132] +mat[157] + mat[182] + mat[207] + 
            mat[232])/10
spot_8<- (mat[8] + mat[33] + mat[58] + mat[83] + mat[108] + mat[133] +mat[158] + mat[183] + mat[208] + 
            mat[233])/10
spot_9<- (mat[9] + mat[34] + mat[59] + mat[84] + mat[109] + mat[134] +mat[159] + mat[184] + mat[209] + 
            mat[234])/10
spot_10<- (mat[10] + mat[35] + mat[60] + mat[85] + mat[110] + mat[135] +mat[160] + mat[185] + mat[210] + 
            mat[235])/10
spot_11<- (mat[11] + mat[36] + mat[61] + mat[86] + mat[111] + mat[136] +mat[161] + mat[186] + mat[211] + 
             mat[236])/10
spot_12<- (mat[12] + mat[37] + mat[62] + mat[87] + mat[112] + mat[137] +mat[162] + mat[187] + mat[212] + 
             mat[237])/10
spot_13<- (mat[13] + mat[38] + mat[63] + mat[88] + mat[113] + mat[138] +mat[163] + mat[188] + mat[213] + 
             mat[238])/10
spot_14<- (mat[14] + mat[39] + mat[64] + mat[89] + mat[114] + mat[139] +mat[164] + mat[189] + mat[214] + 
             mat[239])/10
spot_15<- (mat[15] + mat[40] + mat[65] + mat[90] + mat[115] + mat[140] +mat[165] + mat[190] + mat[215] + 
             mat[240])/10
spot_16<- (mat[16] + mat[41] + mat[66] + mat[91] + mat[116] + mat[141] +mat[166] + mat[191] + mat[216] + 
             mat[241])/10
spot_17<- (mat[17] + mat[42] + mat[67] + mat[92] + mat[117] + mat[142] +mat[167] + mat[192] + mat[217] + 
            mat[242])/10
spot_18<- (mat[18] + mat[43] + mat[68] + mat[93] + mat[118] + mat[143] +mat[168] + mat[193] + mat[218] + 
             mat[243])/10
spot_19<- (mat[19] + mat[44] + mat[69] + mat[94] + mat[119] + mat[144] +mat[169] + mat[194] + mat[219] + 
             mat[244])/10
spot_20<- (mat[20] + mat[45] + mat[70] + mat[95] + mat[120] + mat[145] +mat[170] + mat[195] + mat[220] + 
             mat[245])/10
spot_21<- (mat[21] + mat[46] + mat[71] + mat[96] + mat[121] + mat[146] +mat[171] + mat[196] + mat[221] + 
             mat[246])/10
spot_22<- (mat[22] + mat[47] + mat[72] + mat[97] + mat[122] + mat[147] +mat[172] + mat[197] + mat[222] + 
             mat[247])/10
spot_23<- (mat[23] + mat[48] + mat[73] + mat[98] + mat[123] + mat[148] +mat[173] + mat[198] + mat[223] + 
             mat[248])/10
spot_24<- (mat[24] + mat[49] + mat[74] + mat[99] + mat[124] + mat[149] +mat[174] + mat[199] + mat[224] + 
             mat[249])/10
spot_25<- (mat[25] + mat[50] + mat[75] + mat[100] + mat[125] + mat[150] +mat[175] + mat[200] + mat[225] + 
             mat[250])/10



new_mat[1,1] <- spot_1
new_mat[2,1] <- spot_2
new_mat[3,1] <- spot_3
new_mat[4,1] <- spot_4
new_mat[5,1] <- spot_5
new_mat[1,2] <- spot_6
new_mat[2,2] <- spot_7
new_mat[3,2] <- spot_8
new_mat[4,2] <- spot_9
new_mat[5,2] <- spot_10
new_mat[1,3] <- spot_11
new_mat[2,3] <- spot_12
new_mat[3,3] <- spot_13
new_mat[4,3] <- spot_14
new_mat[5,3] <- spot_15
new_mat[1,4] <- spot_16
new_mat[2,4] <- spot_17
new_mat[3,4] <- spot_18
new_mat[4,4] <- spot_19
new_mat[5,4] <- spot_20
new_mat[1,5] <- spot_21
new_mat[2,5] <- spot_22
new_mat[3,5] <- spot_23
new_mat[4,5] <- spot_24
new_mat[5,5] <- spot_25






###########################
# unlist and combine 

all <- rbind (
  append('binary_50', mean(unlist(rf_methyl_fac[[3]]))),
  append('binary_50_sample', mean(unlist(rf_methyl_fac_samp[[3]]))) ,
  append('multi_30_50_300', mean(unlist(rf_methyl_multi[[3]]))),
  append('multi_50_50_300_sample', mean(unlist(rf_methyl_multi[[3]])))

)


# rf_methyl_multi_more <- predictAll(data = full_data_rf,
#                               subset = c('age_diagnosis_multi_more'),
#                               selected_features = NULL,
#                               binary = FALSE,
#                               multi = FALSE,
#                               iterations = 20)
# 
# rf_methyl_multi_samp <- predictAll(data = full_data_rf,
#                                    subset = c('age_diagnosis_multi_more_samp'),
#                                    selected_features = NULL,
#                                    binary = FALSE,
#                                    iterations = 20)

# # unlist accuracy results and put into data frame 
# 
# rf_mdm2.nG_clin_reg <- predictAll(data = full_data_rf,
#                                   use_genes = FALSE,
#                                   subset = c('age_diagnosis_fac', 'codon72.npro', 'gdna.exon.intron', 'gdna.base.change', 'gdna.codon',
#                                              'protein.codon.change', 'protein.codon.num', 'splice.delins.snv',
#                                              'mdm2.nG'),
#                                   selected_features = c('codon72.npro', 'gdna.exon.intron', 'gdna.base.change', 'gdna.codon',
#                                                         'protein.codon.change', 'protein.codon.num', 'splice.delins.snv',
#                                                         'mdm2.nG'),
#                                   iterations = 5)
