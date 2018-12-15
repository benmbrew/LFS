# Useful functions when working with logistic regression

library(plotly)
library(scatterplot3d)
library(car)
library(rgl)
library(ROCR)
library(caret)
library(pROC)
library(tidyverse)
library(grid)
library(broom)
library(scales)
library(gridExtra)
library(data.table)

get_acc_val <- function(temp_dat, thresh){
  all_alphas <- (1:10)/10
  result_list <- list()
  for(i in 1:length(all_alphas)){
    this_alpha <- all_alphas[i]
    sub_dat <- temp_dat[temp_dat$alpha == this_alpha,]
    sub_dat$pred_label <-as.factor(ifelse(sub_dat$valid_age_pred > thresh, 'positive', 'negative'))
    sub_dat$pred_label <- factor(sub_dat$pred_label, levels = c('positive', 'negative'))
    sub_dat$valid_age_label <- factor(sub_dat$valid_age_label, levels = c('positive', 'negative'))
    sub_dat$acc <- caret::confusionMatrix(sub_dat$pred_label, sub_dat$valid_age_label)$overall[1]
    result_list[[i]] <- sub_dat
    print(i)
  }
  temp <- do.call('rbind', result_list)
  return(temp)
}
 
get_young_labels <- function(temp_dat, thresh, age){
  all_alphas <- (1:10)/10
  result_list <- list()
  for(i in 1:length(all_alphas)){
    this_alpha <- all_alphas[i]
    sub_dat <- temp_dat[temp_dat$alpha == this_alpha,]
    sub_dat <- sub_dat[sub_dat$age_sample_collection < age,]
    sub_dat$pred_label <-as.factor(ifelse(sub_dat$controls_age_pred > thresh, 'positive', 'negative'))
    sub_dat$pred_label <- factor(sub_dat$pred_label, levels = c('positive', 'negative'))
    sub_dat$controls_age_label <- factor(sub_dat$controls_age_label, levels = c('positive', 'negative'))
    sub_dat <- sub_dat[, c('tm_donor','alpha','age_sample_collection','controls_age_pred', 'controls_age_label', 'pred_label')]
    result_list[[i]] <- sub_dat
    print(i)
  }
  temp <- do.call('rbind', result_list)
  return(temp)
}

plot_acc <- function(temp_dat, acc_column, column, bar) {
  
  column_name <- column
  
  
  if(bar){
    temp_dat <- temp_dat[, c(column, acc_column)]
    names(temp_dat) <- c('V1', 'Avg Accuracy')
    
    g1 <- ggplot(temp_dat, 
                 aes(reorder(V1, -`Avg Accuracy`),
                     `Avg Accuracy`)) +
      geom_bar(alpha = 0.6,
               color = 'black',
               fill = 'grey',
               stat = 'identity') +
      geom_smooth(method = 'lm',
                  color = 'red',
                  linetype = 1) +
      labs(title = paste0(column_name, ' and Accuracy'),
           x = column_name,
           y = 'Accuracy') +
      theme(axis.text.x = element_text(angle = 45, hjust = 0.5, size = 5))
  } else {
    temp_dat <- temp_dat[, c(column, acc_column)]
    names(temp_dat) <- c('V1', 'Accuracy')
    
    g1 <- ggplot(temp_dat, 
                 aes(V1, Accuracy)) +
      geom_point(pch = 21,
                 size = 2,
                 alpha = 0.6,
                 color = 'black',
                 fill = 'grey') +
      geom_smooth(method = 'lm',
                  color = 'red',
                  linetype = 1) +
      labs(title = paste0(column_name, ' and Accuracy'),
           x = column_name,
           y = 'Accuracy')
    
  }
  
  return(g1)
  
}


plot_3d_model_means <- function(temp_dat){
  with(temp_dat, {
    s3d <- scatterplot3d(mean_lambda, mean_alpha, mean_acc,        # x y and z axis
                         color="darkgrey", 
                         pch=1,        # filled blue circles
                         type="h",
                         main="Alpha and Lambda choice",
                         xlab="Enet lambda",
                         ylab="Enet alpha",
                         zlab="Model accuracy")
    s3d.coords <- s3d$xyz.convert(mean_lambda, mean_alpha, mean_acc) # convert 3D coords to 2D projection
    my.lm <- lm(temp_dat$mean_acc ~ temp_dat$mean_lambda + temp_dat$mean_alpha )
    s3d$plane3d(my.lm)
    s3d$points3d(mean_lambda, mean_alpha, mean_acc,
                 col = adjustcolor("black", alpha.f = 0.8), type = 'p', pch = 16)
  })
}


plot_3d_model <- function(temp_dat, type){
  with(temp_dat, {
    s3d <- scatterplot3d(lambda, alpha, accuracy,        # x y and z axis
                         color=adjustcolor('grey', alpha.f = 0.2), 
                         pch=1,        # filled blue circles
                         type=type,
                         main="Alpha and Lambda choice",
                         xlab="Enet lambda",
                         ylab="Enet alpha",
                         zlab="Model accuracy")
    s3d.coords <- s3d$xyz.convert(lambda, alpha, accuracy) # convert 3D coords to 2D projection
    my.lm <- lm(temp_dat$accuracy ~ temp_dat$lambda + temp_dat$alpha )
    s3d$plane3d(my.lm)
    s3d$points3d(lambda, alpha, accuracy,
                 col = adjustcolor("black", alpha.f = 0.1), type = 'p', pch = 16)
  })
}


AccuracyCutoffInfo <- function( test, predict, actual )
{
  # change the cutoff value's range as you please 
  cutoff <- seq( .4, .8, by = .05 )
  
  accuracy <- lapply( cutoff, function(c)
  {
    # use the confusionMatrix from the caret package
    cm_test  <- confusionMatrix( as.numeric( test[[predict]]  > c ), test[[actual]]  )
    
    dt <- data.table( cutoff = c,
                      test   = cm_test$overall[["Accuracy"]] )
    return(dt)
  }) %>% rbindlist()
  
  # visualize the accuracy of the train and test set for different cutoff value 
  # accuracy in percentage.
  accuracy_long <- gather( accuracy, "data", "accuracy", -1 )
  
  plot <- ggplot( accuracy_long, aes( cutoff, accuracy) ) + 
    geom_line( size = 1 ) + geom_point( size = 3 ) +
    scale_y_continuous( label = percent ) +
    ggtitle( "Cutoff" )
  
  return( list( data = accuracy, plot = plot ) )
}


# ------------------------------------------------------------------------------------------
# [ConfusionMatrixInfo] : 
# Obtain the confusion matrix plot and data.table for a given
# dataset that already consists the predicted score and actual outcome.
# @data    : your data.table or data.frame type data that consists the column
#            of the predicted score and actual outcome 
# @predict : predicted score's column name
# @actual  : actual results' column name
# @cutoff  : cutoff value for the prediction score 
# return   : 1. data : a data.table consisting of three column
#            		   the first two stores the original value of the prediction and actual outcome from
#			 		   the passed in data frame, the third indicates the type, which is after choosing the 
#			 		   cutoff value, will this row be a true/false positive/ negative 
#            2. plot : plot that visualizes the data.table 



ConfusionMatrixInfo <- function( data, predict, actual, cutoff, get_plot)
{	
  # extract the column ;
  # relevel making 1 appears on the more commonly seen position in 
  # a two by two confusion matrix	
  predict <- data[[predict]]
  actual  <- relevel( as.factor( data[[actual]] ), "positive" )
  
  result <- data.table( actual = actual, predict = predict )
  
  # caculating each pred falls into which category for the confusion matrix
  result[ , type := ifelse( predict >= cutoff & actual == 'positive', "TP",
                            ifelse( predict >= cutoff & actual == 'negative', "FP", 
                                    ifelse( predict <  cutoff & actual == 'positive', "FN", "TN" ) ) ) %>% as.factor() ]
  
  # jittering : can spread the points along the x axis 
  plot <- ggplot( result, aes( actual, predict, color = type ) ) + 
    geom_violin( fill = "white", color = NA ) +
    geom_jitter( shape = 1 ) + 
    geom_hline( yintercept = cutoff, color = "blue", alpha = 0.6 ) + 
    scale_y_continuous( limits = c( 0, 1 ) ) + 
    scale_color_discrete( breaks = c( "TP", "FN", "FP", "TN" ) ) + # ordering of the legend 
    guides( col = guide_legend( nrow = 2 ) ) + # adjust the legend to have two rows  
    ggtitle( sprintf( "Confusion Matrix with Cutoff at %.2f", cutoff ) )
  
  if(get_plot) {
    return(plot)
  } else {
    return(as.data.frame(result))
  }
}


# [ROCInfo] : 
# Pass in the data that already consists the predicted score and actual outcome.
# to obtain the ROC curve 
# @data    : your data.table or data.frame type data that consists the column
#            of the predicted score and actual outcome
# @predict : predicted score's column name
# @actual  : actual results' column name
# @cost.fp : associated cost for a false positive 
# @cost.fn : associated cost for a false negative 
# return   : a list containing  
#			 1. plot        : a side by side roc and cost plot, title showing optimal cutoff value
# 				 	   		  title showing optimal cutoff, total cost, and area under the curve (auc)
# 		     2. cutoff      : optimal cutoff value according to the specified fp/fn cost 
#		     3. totalcost   : total cost according to the specified fp/fn cost
#			 4. auc 		: area under the curve
#		     5. sensitivity : TP / (TP + FN)
#		     6. specificity : TN / (FP + TN)


# data = dat_sample 
# predict = 'mean_preds' 
# actual = 'real_label' 
# cost.fp = cost_fp
# cost.fn = cost_fn
ROCInfo <- function( data, predict, actual, cost.fp, cost.fn )
{
  # calculate the values using the ROCR library
  # true positive, false postive 
  pred <- prediction( data[[predict]], data[[actual]] )
  perf <- performance( pred, "tpr", "fpr" )
  roc_dt <- data.frame( fpr = perf@x.values[[1]], tpr = perf@y.values[[1]] )
  
  # cost with the specified false positive and false negative cost 
  # false postive rate * number of negative instances * false positive cost + 
  # false negative rate * number of positive instances * false negative cost
  cost <- perf@x.values[[1]] * cost.fp * sum( data[[actual]] == 'negative' ) + 
    ( 1 - perf@y.values[[1]] ) * cost.fn * sum( data[[actual]] == 'positive' )
  
  cost_dt <- data.frame( cutoff = pred@cutoffs[[1]], cost = cost )
  
  # optimal cutoff value, and the corresponding true positive and false positive rate
  best_index  <- which.min(cost)
  best_cost   <- cost_dt[ best_index, "cost" ]
  best_tpr    <- roc_dt[ best_index, "tpr" ]
  best_fpr    <- roc_dt[ best_index, "fpr" ]
  best_cutoff <- pred@cutoffs[[1]][ best_index ]
  
  # area under the curve
  auc <- performance( pred, "auc" )@y.values[[1]]
  
  # normalize the cost to assign colors to 1
  normalize <- function(v) ( v - min(v) ) / diff( range(v) )
  
  # create color from a palette to assign to the 100 generated threshold between 0 ~ 1
  # then normalize each cost and assign colors to it, the higher the blacker
  # don't times it by 100, there will be 0 in the vector
  col_ramp <- colorRampPalette( c( "green", "orange", "red", "black" ) )(100)   
  col_by_cost <- col_ramp[ ceiling( normalize(cost) * 99 ) + 1 ]
  
  roc_plot <- ggplot( roc_dt, aes( fpr, tpr ) ) + 
    geom_line( color = rgb( 0, 0, 1, alpha = 0.3 ) ) +
    geom_point( color = col_by_cost, size = 4, alpha = 0.2 ) + 
    geom_segment( aes( x = 0, y = 0, xend = 1, yend = 1 ), alpha = 0.8, color = "royalblue" ) + 
    labs( title = "ROC", x = "False Postive Rate", y = "True Positive Rate" ) +
    geom_hline( yintercept = best_tpr, alpha = 0.8, linetype = "dashed", color = "steelblue4" ) +
    geom_vline( xintercept = best_fpr, alpha = 0.8, linetype = "dashed", color = "steelblue4" )				
  
  cost_plot <- ggplot( cost_dt, aes( cutoff, cost ) ) +
    geom_line( color = "blue", alpha = 0.5 ) +
    geom_point( color = col_by_cost, size = 4, alpha = 0.5 ) +
    ggtitle( "Cost" ) +
    scale_y_continuous( labels = comma ) +
    geom_vline( xintercept = best_cutoff, alpha = 0.8, linetype = "dashed", color = "steelblue4" )	
  
  options(scipen = '999')
  # the main title for the two arranged plot
  sub_title <- sprintf( "Cutoff at %.2f - Total Cost = %a, AUC = %.3f", 
                        best_cutoff, best_cost, auc )
  
  # arranged into a side by side plot
  plot <- arrangeGrob( roc_plot, cost_plot, ncol = 2, 
                       top = textGrob( sub_title, gp = gpar( fontsize = 16, fontface = "bold" ) ) )
  
  return( list( plot 		  = plot, 
                cutoff 	  = best_cutoff, 
                totalcost   = best_cost, 
                auc         = auc,
                sensitivity = best_tpr, 
                specificity = 1 - best_fpr ) )
}
