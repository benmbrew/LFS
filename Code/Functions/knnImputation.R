# Keeping track of the splitting conditions would be easier if there
# was a single variable which could be used to keep track of the logic

knnImputation <- function(incompleteData, sample_rows) {
  # Transpose the data for KNN
  # KNN requires samples in the columns
  sampleRowsRequired <- FALSE
  transposeData <- sample_rows != sampleRowsRequired
  # incompleteData <- as.matrix(incompleteData)
  if (transposeData) {
    incompleteData <- t(incompleteData)
  }
  
  # make data numeric
  for(i in 1:ncol(incompleteData)) {
    incompleteData[, i] <- as.numeric(incompleteData[, i])
    print(i)
  }
  
  # Concatenate the data by rows because the samples are in columns
  concatenatedIncompleteData <- incompleteData
  
  # Impute the missing data
  set.seed(10)
  concatenatedImputedData <- impute.knn(concatenatedIncompleteData, k = 10)$data
  
  
  # Reorient the data to its original form
  if (transposeData) {
    imputedData <- t(concatenatedImputedData)
  }
  
  return(imputedData)
}