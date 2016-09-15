plot_models_performance <- function(models, num.of.partition, features.name, pdf.name, plot.title) {
  #plot_models_performance(cp2.other_model.all, NUM_OF_PARTITION, "exp+metyl", "Classifiers_Group4_c4_expANDmethyl.pdf", "Group4 c4 using Top1% of exp. & methyl.")
  #plot_models_performance(models.all, NUM_OF_PARTITION, "", "Classifiers_Group4_c4_expANDmethyl.pdf", "Group4 c4 using Top1% of exp. & methyl.")
  # This function takes models = models.methyl, num partitions, features.name ="",  "",
  # pdf.name = paste0(results_folder,"/Classifiers_test_results.pdf"), plot.title = "Classification using methyl. markers"
  require("doParallel")
  #require("vioplot")

  # Boxplot Data ---------------------------------
  # for each parition run this to create matrix with auc values 
  boxplot_matrix_auc <- foreach (temp.ind = 1:num.of.partition, .combine=rbind) %do% {
    c(
      mean(models[[temp.ind]]$elasticNet$test_auc)
      , mean(models[[temp.ind]]$lasso$test_auc)
      , mean(models[[temp.ind]]$ridge$test_auc)
      , mean(models[[temp.ind]]$rf$test_auc)
      , mean(models[[temp.ind]]$dtree$test_auc)
      , mean(models[[temp.ind]]$svmLinear$test_auc)
      , mean(models[[temp.ind]]$svmRadial$test_auc)
    )  
  }
  boxplot_matrix_acc <- foreach (temp.ind = 1:num.of.partition, .combine=rbind) %do% {
    c(
      mean(models[[temp.ind]]$elasticNet$test_acc)
      , mean(models[[temp.ind]]$lasso$test_acc)
      , mean(models[[temp.ind]]$ridge$test_acc)
      , mean(models[[temp.ind]]$rf$test_acc)
      , mean(models[[temp.ind]]$dtree$test_acc)
      , mean(models[[temp.ind]]$svmLinear$test_acc)
      , mean(models[[temp.ind]]$svmRadial$test_acc)
    )  
  }
  temp.boxplot_data_auc <- data.frame(boxplot_matrix_auc)
  temp.boxplot_data_acc <- data.frame(boxplot_matrix_acc)
  #print(dim(temp.boxplot_data))

  temp.names <- c(
    paste("Elastic Net", features.name)
    , paste("Lasso", features.name)
    , paste("Ridge", features.name)
    , paste("Rand.Forest", features.name)
    , paste("Decis.Tree", features.name)
    , paste("SVM lin", features.name)
    , paste("SVM rbf", features.name)
  )
  #temp.order <- order(apply(temp.boxplot_data, 2, median))
  #with(temp.boxplot_data, reorder(TYPE, temp.boxplot_data, median))

  #print(rank(apply(temp.boxplot_data, 2, median), ties.method = c("first")))

  ## download and source SystematicInvestor package for nice table plotting
  # require(RCurl)
  # sit = getURLContent('https://github.com/systematicinvestor/SIT/raw/master/sit.gz', binary=TRUE, followlocation = TRUE, ssl.verifypeer = FALSE)
  #     con = gzcon(rawConnection(sit, 'rb'))
  #     source(con)
  # close(con)
  ## read SystematicInvestor from local file instead of downloading it
  sit <- gzcon(file(paste0(test,"/sit.gz"), "rb"))
  source(sit)

  alltabs <- foreach (temp.ind = 1:num.of.partition, .combine=rbind) %do% {
    list(
      models[[temp.ind]]$elasticNet$test_stats$table
      , models[[temp.ind]]$lasso$test_stats$table
      , models[[temp.ind]]$ridge$test_stats$table
      , models[[temp.ind]]$rf$test_stats$table
      , models[[temp.ind]]$dtree$test_stats$table
      , models[[temp.ind]]$svmLinear$test_stats$table
      , models[[temp.ind]]$svmRadial$test_stats$table
    )  
  }
  
  # print(alltabs)
  # print(models[[1]]$lasso$test_stats$table)
  # print(Reduce(`+`, alltabs[, 1]))
  # A <- Reduce(`+`, alltabs[, 1]) 
  # A <- round(A/rowSums(A), digits=3)
  # print(A)

  pdf(pdf.name)


  ### violin plot
  # library(reshape)
  # library(easyGgplot2)

  # df <- temp.boxplot_data_auc
  # names(df) <- temp.names
  # at.order <- temp.names
  # at.order[rank(apply(df, 2, mean), ties.method = c("first"))] <- temp.names
  # df <- df[, at.order]

  # #print(melt(df))
  # #print(at.order)
  # print(ggplot2.violinplot(data=melt(df),
  #                   xName='variable', yName='value',
  #                   mainTitle=plot.title,
  #                   ytitle="AUC on test set",
  #                   xtitle="classifiers")
  # )

  boxplot(temp.boxplot_data_auc, las=2,
          main=plot.title,
          par(mar = c(10, 5, 4, 2) + 0.1), 
          ylab = "AUC on test set",
          names = temp.names,
          at = rank(apply(temp.boxplot_data_auc, 2, mean), ties.method = c("first"))
  )
  
  boxplot(temp.boxplot_data_acc, las=2, 
          main=plot.title,
          par(mar = c(10, 5, 4, 2) + 0.1), 
          ylab = "accuracy on test set",
          names = temp.names,
          at = rank(apply(temp.boxplot_data_acc, 2, mean), ties.method = c("first"))
  )

  for (i in 1:7){
    A <- Reduce(`+`, alltabs[, i]) # adds up confusion matrix for each parition
    #A <- A/rowSums(A)        # row normalize
    A <- t(t(A)/colSums(A))   # column normalize c
    A <- round(A, digits=3)
    plot.table(A, temp.names[i])
  }

  dev.off()

}

# temp.median_results <- apply(temp.boxplot_data, 2, median)[temp.order]
# names(temp.median_results) <- temp.names[temp.order]
# temp.median_results

## USED for plotting varying number of training data
# temp.median_results <- apply(temp.boxplot_data, 2, median)
# names(temp.median_results) <- temp.names
# temp.median_results
# #varying_training_matrix <- matrix(nrow = 0, ncol = length(temp.median_results))
# load("~/bor_var_cp.RData")
# colnames(varying_training_matrix)
# dim(varying_training_matrix)
# varying_training_matrix <- rbind(varying_training_matrix, temp.median_results)
# dim(varying_training_matrix)
# getwd()
# save(varying_training_matrix, file="bor_var_cp.RData")

# # use when done
# #varying_training_matrix <- varying_training_matrix[c(5:8, 4:1), ]

# stopifnot(dim(varying_training_matrix)[2] == 21)
# plot(0, 0, xlim = c(30, 300), ylim = c(0.5, 0.8), type = "n", 
#      xlab = "# of Cell Lines in Training Data",
#      ylab = "Test AUC")
# cl <- rainbow(7)
# temp.points <- rep(1, 7)
# temp.points <- c(temp.points, rep(2, 7))
# temp.points <- c(temp.points, rep(3, 7))
# for (temp.ind in 1:21) {
#   lines(c(30, 40, 50, 70, 90, 110, 150, 200, 250, 300), varying_training_matrix[, temp.ind], col = cl[temp.ind %% 7], type = 'b', pch=temp.points[temp.ind])
# }
# legend("topleft", legend = colnames(varying_training_matrix), col=cl, pch=temp.points, ncol=3) # optional legend
## END

## Comparing c2p, cp2p, and p2p (produce plots for best classifiers for each)
# best for c2p
# load("~/AeroFS/Drug Prediction Project/Output/output_bortezomib/bor_var_cp.RData")

# best <- list()

# best_ind <- which(max(varying_training_matrix[10, ]) == varying_training_matrix[10, ])
# best_ind
# best$cp <- varying_training_matrix[, best_ind]

# load("~/AeroFS/Drug Prediction Project/Output/output_bortezomib/bor_var_cpp.RData")
# best$cpp <- varying_training_matrix[, best_ind]
# load("~/AeroFS/Drug Prediction Project/Output/output_bortezomib/bor_var_pp.RData")
# best$pp <- varying_training_matrix[, best_ind]
# min(length(best$cpp), length(best$pp))
# load("~/bor_var_cp.RData")
# best$cp <-rep(max(varying_training_matrix), 6)

# plot(0, 0, xlim = c(30, 110), ylim = c(0.5, 0.8), type = "n", 
#      xlab = "# of Patients in Training Data",
#      ylab = "Test AUC")
# cl <- rainbow(3)
# lines(c(30, 40, 50, 70, 90, 110), best$pp, col = cl[1], type = 'b', pch=1)
# lines(c(30, 40, 50, 70, 90, 110), best$cpp[1:6], col = cl[2], type = 'b', pch=2)
# lines(c(30, 40, 50, 70, 90, 110), best$cp, col = cl[3], type = 'o', pch=3, lty=3)
# legend("topleft", legend = c("P2P", "C2P", "C2P"), col=cl, pch=1:3, ncol=3) # optional legend

# ##
# best <- list()
# load("~/AeroFS/Drug Prediction Project/Output/output_bortezomib/bor_var_cpp.RData")
# best$cpp <- varying_training_matrix[, 13]
# load("~/AeroFS/Drug Prediction Project/Output/output_bortezomib/bor_var_pp.RData")
# best$pp <- varying_training_matrix[, 13]
# max(length(best$cpp), length(best$pp))

# load("~/AeroFS/Drug Prediction Project/Output/output_bortezomib/bor_var_cp.RData")
# best$cp <-rep(max(varying_training_matrix[, 13]), 12)

# plot(0, 0, xlim = c(0, 150), ylim = c(0.5, 0.9), type = "n", 
#      xlab = "# of Patients in Training Data",
#      ylab = "Test AUC")
# cl <- rainbow(3)
# lines(c(30, 40, 50, 70, 90, 110), best$pp, col = cl[1], type = 'b', pch=1)
# lines(c(10, 20, 30, 40, 50, 70, 90, 110, 120, 130, 140, 150), best$cpp, col = cl[2], type = 'b', pch=2)
# lines(c(10, 20, 30, 40, 50, 70, 90, 110, 120, 130, 140, 150), best$cp, col = cl[3], type = 'o', pch=3, lty=3)
# legend("topleft", legend = c("P2P", "CP2P", "C2C"), col=cl, pch=1:3, ncol=3) # optional legend


## END



# setwd("/home/cheng/")
# write.table(temp.median_results, file="cp09_median_auc.txt", quote=FALSE, sep="\t")


# # Which patients are being misclassified? ---------------------------------
# cp.misclassified.matrix = matrix(nrow = length(partition$patients), ncol = dim(boxplot_matrix)[2])
# rownames(cp.misclassified.matrix) = names(binaryResponse)
# colnames(cp.misclassified.matrix) = temp.names 

# source('compute_misclassified.R')

# rm(list = ls(pattern="temp*"))

# temp.matrix = compute_misclassified(model = models, ground_truth = ground_truth, partition = partition, response_with_name = binaryResponse)
# cp.misclassified.matrix[, 1:6] = temp.matrix


# temp.heatmap_colors <- rep("red", length(binaryResponse))
# temp.heatmap_colors[binaryResponse == 0] = "blue"
# dev.off()
# heatmap.2(cp.misclassified.matrix, RowSideColors = temp.heatmap_colors, margins = c(10, 5), col=terrain.colors(50))

# # which patients are always misclassified? 
# temp.acc = apply(cp.misclassified.matrix, 1, mean)
# temp.always_misclassified <- names(which(temp.acc < 0.5))
# write.table(temp.always_misclassified, file="cp_misclassified.txt", quote=FALSE)
# length(names(which(temp.acc < 0.5)))
# length(names(which(temp.acc < 0.25)))

