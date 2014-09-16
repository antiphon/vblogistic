#' test fitter
n<-100
p <- 10
set.seed(1)
X <- matrix(rnorm(n*p), ncol=p)
theta <- rnorm(p)
logit <- function(x) 1/(1+exp(-x))
offset <- rnorm(n)
pr <- logit( X%*%theta +offset)

y <- rbinom(n, 1, pr)


source("R/logistic_regression.R")

source("R/classes.R")
library(Matrix)
f <- vblogit(y, X, verb=T, offset=offset, eps=0.001)
g <- glm(y~-1+X, family=binomial, offset=offset)
