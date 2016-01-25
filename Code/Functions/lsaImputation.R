lsaImputation <- function(incomplete_data, sample_rows) {
  # Variables which have not been initialized:
  # imputedFile, incompleteFile, projectFolder, jvmGBLimit
  
  # Transpose the data for LSA
  # LSA requires samples in rows
  sampleRowsRequired <- TRUE
  transposeData <- sample_rows != sampleRowsRequired
  if (transposeData) {
    incomplete_data <- t(incomplete_data)
  }
  
  # Concatenate the data by columns because the samples are in rows
  # Row and column names are required for the LSA code
  #concatenatedIncompleteData <- do.call(cbind, incompleteData) # cols
  #concatenatedRownames <- rownames(concatenatedIncompleteData) # rows
  concatenatedRownames <- rownames(incomplete_data) # rows
#   colnames(incomplete_data) <-
#     as.character(1:ncol(incomplete_data))
#   rownames(incomplete_data) <-
#     as.character(1:nrow(incomplete_data))
  
#   colnames(concatenatedIncompleteData) <-
#    as.character(1:ncol(concatenatedIncompleteData))
#   rownames(concatenatedIncompleteData) <-
#    as.character(1:nrow(concatenatedIncompleteData))
#   
  # Save the incomplete data to a file
  cat("xxx","\t", file=incomplete_file, sep="")
  write.table(incomplete_data, incomplete_file,
              col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE,
              na="NULL",append=TRUE)
  
  # Impute the missing data
  lsaCodePath <- paste(project_folder, "Code/Functions/LSimpute.jar",
                       sep="/")
  jvmMemLimit <- paste("-Xmx", jvmGBLimit, "G", sep="")
  javaCommand <- paste("java -jar", jvmMemLimit, "-server ",
                       lsaCodePath, incomplete_file, imputed_file, "6",
                       sep=" ")
  system(javaCommand)
  
  # Read the imputed file
  #concatenatedImputedData <-
  imputed_data <- read.table(imputed_file, header=TRUE, sep="\t")
  #concatenatedImputedData <- concatenatedImputedData[, -1]
  imputed_data <- imputed_data[, -1]
  rownames(imputed_data) <- concatenatedRownames # rows
  
  # Remove the files used for imputation
  system(paste("rm", imputed_file, sep=" "))
  system(paste("rm", incomplete_file, sep=" "))
  
#   # Split the concatenated data by columns
#   featureSize <- sapply(incompleteData, ncol) # cols
#   imputedData <- splitConcatenatedData(concatenatedImputedData, 
#                                        featureSize,
#                                        sampleRowsRequired)
  
  # Reorient the data to its original form
  if (transposeData) {
    imputedData <- lapply(imputedData, t)
  }
  
  return(imputed_data)
}