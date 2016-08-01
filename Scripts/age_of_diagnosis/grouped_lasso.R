### 
# this script will run grouped lasso
library(grplasso)


# Initialize folders
home_folder <- '/home/benbrew/hpf/largeprojects/agoldenb/ben/Projects/'
project_folder <- paste0(home_folder, '/LFS')
test <- paste0(project_folder, '/Scripts/regression_template')
data_folder <- paste0(project_folder, '/Data')
methyl_data <- paste0(data_folder, '/methyl_data')
clin_data <- paste0(data_folder, '/clin_data')
results_folder <- paste0(test, '/Results')


# read in methyl impute raw
full_data <- read.csv(paste0(data_folder, '/full_data.csv'), stringsAsFactors = F)
kmeans <- read.csv(paste0(data_folder, '/kmeans_labels.csv'), stringsAsFactors = F)
hier <- read.csv(paste0(data_folder, '/hier_labels.csv'), stringsAsFactors = F)

# add 1 to each label because grplasso cant deal with label 1. 
kmeans <- kmeans$V1 
hier <- hier$labels 

# remove unnecessary columns 
full_data$X <- NULL

###############################################################################################
# first run on just methylation 

# take subset to do test run 
test_data <- full_data[, c(6, 27:ncol(full_data))]
test_data <- test_data[complete.cases(test_data),]

y <- test_data$age_diagnosis

test_dat <- cbind(1, as.matrix(test_data[,-1]))


## Use a multiplicative grid for the penalty parameter lambda, starting
## at the maximal lambda value
index <- as.numeric(c(NA, kmeans))
colnames(test_dat)[1] <- 'Intercept'

dim(test_data)
length(index)
length(y)

lambda <- lambdamax(test_dat, y, index = index, penscale = sqrt,
                    model = LinReg()) * 0.5^(0:7)

## Fit the solution path on the lambda grid
fit <- grplasso(test_dat,  y, index = index, lambda = lambda, model = LinReg(),
                penscale = sqrt,
                control = grpl.control(update.hess = "lambda", trace = 0))
## Plot coefficient paths
plot(fit)



pred <- predict(fit)
pred.resp <- predict(fit, type = "response")
## The following points should lie on the sigmoid curve
plot(pred, pred.resp)





##########################################################################################3
# Run grouped lasso with model formula, throws error because of stack overflow. 

###############################################################################################
# first run on just methylation 

# take subset to do test run 
full_data <- full_data[, c(6, 27:5000)]
full_data <- full_data[complete.cases(full_data),]

y <- as.numeric(full_data$age_diagnosis)

test_dat <- cbind(1, full_data[,-1])
colnames(test_dat)[1] <- 'Intercept'
test_dat <- as.matrix(test_dat)

## Use a multiplicative grid for the penalty parameter lambda, starting

lambda <- lambdamax(age_diagnosis~., nonpen = ~1, full_data)* 0.5^(0:5)


fit <- grplasso(age_diagnosis~., nonpen = ~1, full_data, lambda = lambda, model = LinReg(),
                penscale = sqrt,
                control = grpl.control(update.hess = "lambda", trace = 0))

## Plot coefficient paths
plot(fit)


pred <- predict(fit)
pred.resp <- predict(fit, type = "response")
## The following points should lie on the sigmoid curve
plot(pred, pred.resp)





































