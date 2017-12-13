##### Supervised PCA
library(superpc)
library(DAAG)
# tutorial here
# http://statweb.stanford.edu/~tibs/superpc/tutorial.html

# could add in clinical variables here too

# Supervised principal components is a generalization of principal components regression. The first (or first few) 
# principal components are the linear combinations of the features that capture the directions of largest variation in a 
# dataset. But these directions may or may not be related to an outcome variable of interest. To find linear combinations 
# that are related to an outcome variable, we compute univariate scores for each gene and then retain only those features 
# whose score exceeds a threshold. A principal components analysis is carried out using only the data from these selected 
# features. Finally, these "supervised principal components" are used in a regression model to predict the outcome. 
# To summarize, the steps are:
  
# 1) Compute (univariate) standard regression coefficients for each feature
# 2) Form a reduced data matrix consisting of only those features whose univariate coefficient exceeds a threshold theta in 
# absolute value (theta is estimated by cross-validation) 
# 3) Compute the first (or first few) principal components of the 
# reduced data matrix 
# 4) Use these principal component(s) in a regression model to predict the outcome

# This idea can be used in standard regression problems with a quantitative outcome, and also in generalized 
# regression problems such as survival analysis. In the latter problem, the regression coefficients in step (1) are 
# obtained from a proportional hazards model. The superpc R package handles these two cases: standard regression and survival data.

# There is one more important point: the features (e.g genes) which important in the prediction are not necessarily 
# the ones that passed the screen in step 2. There are other features that may have as high a correlation with the 
# supervised PC predictor. So we compute an importance score for each feature equal to its correlation with the supervised 
# PC predictor. A reduced predictor is formed by soft-thresholding the importance scores, and using these shrunken scores
# as weights. The soft-thresholding sets the weight of some features to zero, hence throwing them out of the model. 
# The amount of shrinkage is determined by cross-validation. The reduced predictor often performs as well or better than 
# than the supervised PC predictor, and is more interpretable.



# generate some synthetic survival data.
#
# there are 1000 features and 60 samples
#  the outcome is highly correlated with the first principal component of the
#  first 80 features


set.seed(464)


x<-matrix(rnorm(1000*100),ncol=100)
v1<- svd(x[1:80,])$v[,1]

y<-2+5*v1+ .05*rnorm(100)

xtest<-x[, 1:50]
ytest<-(2+5*v1+ .05*rnorm(100))[1:50]

featurenames <- paste("feature",as.character(1:1000),sep="")



# create train and test data objects. censoring.status=1 means the event occurred;
#  censoring.status=0 means censored

data<-list(x=x,y=y, featurenames=featurenames)
data.test<-list(x=xtest,y=ytest,featurenames= featurenames)


# train  the model. This step just computes the  scores for each feature

# Data object with components x- p by n matrix of features, one observation per column; y- n-vector of outcome measurements; 
# censoring.status- n-vector of censoring censoring.status (1= died or event occurred, 0=survived, or event was censored),
# needed for a censored survival outcome type	 Problem type: "survival" for censored survival outcome, or 
# "regression" for simple quantitative outcomes 0.perc	
# Factor for denominator of score statistic, between 0 and 1: the percentile of standard deviation values 
# added to the denominator. Default is 0.5 (the median)

# training the model - score for each feature
train.obj<- superpc.train(data, type="regression")

# note for regression (non-survival) data, we leave the component "censoring.status"
# out of the data object, and call superpc.train with type="regression".
# otherwise the superpc commands are all the same




# get optimal threshold
cv.obj <- superpc.cv(train.obj, data)

#plot the cross-validation curves. From this plot we see that the 1st 
# principal component is significant and the best threshold  is around 0.7

# For each value of theta you do feature selection and then PCA, and then compute the prediction error on the test 
# data. After cross-validation is finished, you will have one value of error for each value of theta. 
# Select the theta that gives the lowest error. After that you can apply this theta to build your prediction model 
# with all available data (and "to reduce the matrix size" if that is what you are after).
superpc.plotcv(cv.obj)


#See pdf version of the cv plot 

# here we have the luxury of  test data, so we can compute the likelihood ratio statistic
# over the test data and plot them. We see that the threshold of 0.7
# works pretty well
#  
# at each threshold (20 values) we compute the LRST from the features at that threshold 
# (threshold increase, number of features decrease). 
# this might be using a label in the test data, but this is pretty much just a consistency test i think.
# 
# lrtest.obj <- superpc.lrtest.curv(train.obj, data, data.test)
# lrtest.obj$lrtest
# 
# superpc.plot.lrtest(lrtest.obj)
# 
# #See pdf version of the lrtest plot


# at the optimal threshold on the training data (0.7), 
# we compute the princicpal components from the features selcted.
fit.cts <- superpc.predict(object = train.obj, 
                           data = data, 
                           newdata = data.test, 
                           threshold = 1.5, 
                           n.components = 3, 
                           prediction.type ="continuous")

# these are the PCAs you should use for predicting on the test set 
fit.cts$v.pred

# Fit predictive model using outcome of supervised principal components, via lm (for regression data)
# this is using chosen PCAs for test data and estimating a linear model with test y regressed on 3 chosen PCAs


temp <- superpc_fit_outcome_cv(y = data.test$y, 
                               score = fit.cts$v.pred)

# Finally, we look for a predictor of survival a small number of
#genes (rather than all 1000 genes). We do this by computing an importance
# score for each equal its correlation with the supervised PC predictor.
# Then we soft threshold the importance scores, and use the shrunken
# scores as gene weights to from a reduced predictor. Cross-validation
# gives us an estimate of the best amount to shrink and an idea of
#how well the shrunken predictor works.


fit.red<- superpc.predict.red(train.obj, data, data.test, threshold=0.7)

fit.redcv<- superpc.predict.red.cv(fit.red, cv.obj,  data,  threshold=0.7)

superpc.plotred.lrtest(fit.redcv)

#See pdf version of this plot


# Finally we list the significant genes, in order of decreasing importance score


# 
# A note on interpretation:
#   
# The signs of the scores (latent factors) v.pred returned by superpc.predict are chosen so that the regression 
# of the outcome on each factor has a positive coefficient. This is just a convention to aid in interpretation.
# 
# For regression data, this means 
# Higher score => higher mean of the outcome 
# How about the direction of effect for each individual feature (gene)? The function superpc.listfeatures reports 
# an importance score equal to the correlation between each feature and the latent factor.
# 
# For regression data,
# 
# Importance score positive means 
# increase in value of feature => higher mean of the outcome
# 
# The reverse are true for Importance score negative.


