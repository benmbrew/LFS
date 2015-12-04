
kFolds <- length(folds) # list of data split into 5 

for (test.i in 1:kFolds) {
  test <- folds[[test.i]]
  train.folds <- folds[-test.i]
  for (valid.i in 1:(kFolds - 1)) {
    valid <- train.folds[[valid.i]]
    train <- do.call("rbind", train.folds[-valid.i])
    for (parameter.i in 1:length(parameters)) {
      # train model on "train" with the this parameter
      # predict outcomes on "valid" and save the score
    }
  }
  train.total <- do.call("rbind", train.folds)
  # train model on "train.total" using the parameter value with the best average score on the validation sets
  # predict outcomes on "test" and save the score
}
# average the scores on the test set to get an estimate of your procedure's ability to generalize to new data

Note: you might not want to follow this procedure if you have a lot of data to start with

