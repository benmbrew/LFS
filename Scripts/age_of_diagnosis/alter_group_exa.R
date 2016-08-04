## Use the Logistic Group Lasso on the splice data set
library(grplasso)
data(splice)
## Define a list with the contrasts of the factors
contr <- rep(list("contr.sum"), ncol(splice) - 1)
names(contr) <- names(splice)[-1]

## Fit a logistic model
fit.splice <- grplasso(y ~ ., data = splice, model = LogReg(), lambda = 20,
                       contrasts = contr, center = TRUE, standardize = TRUE)



## Perform the Logistic Group Lasso on a random dataset
set.seed(79)
n <- 50 ## observations
p <- 299 ## variables
## First variable (intercept) not penalized, two groups having 2 degrees
## of freedom each
index <- c(NA, rep.int(c(1,2,3,4,5,6,7), 301/7))
index <- index[-c(2)]
## Create a random design matrix, including the intercept (first column)
x <- cbind(1, matrix(rnorm(p * n), nrow = n))
colnames(x) <- c("Intercept", paste("X", 1:p+1, sep = ""))

# make response vector
y <- rnorm(n)^2

dim(x)
length(index)
length(y)

## Use a multiplicative grid for the penalty parameter lambda, starting
## at the maximal lambda value
lambda <- lambdamax(x, y = y, index = index[-2], penscale = sqrt,
                    model = LinReg()) * 0.5^(0:5)


## Fit the solution path on the lambda grid
fit <- grplasso(x, y = y, index = index, lambda = lambda, model = LinReg(),
                penscale = sqrt,
                control = grpl.control(update.hess = "lambda", trace = 0))
## Plot coefficient paths
plot(fit)
