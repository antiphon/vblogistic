#' test fitter
library(vblogistic)
library(devtools)
load_all(".")


n <- 800 + 4*800
p <- 30
np <- 3
set.seed(2)
#' add noise components
X <- matrix(rnorm(n*(p+np)), ncol=p+np)
theta <- c(rep(0, np), rnorm(p))
logit <- function(x) 1/(1+exp(-x))
offset <- rnorm(n, 3.5)
pr <- logit( X%*%theta +offset)

y <- rbinom(n, 1, pr)

library(Matrix)
g<-f <- list()
m <- (p+np)-2
for(i in 1:m) f[[i]]<- vblogit(y, X[,-c(1:i)], verb=T, offset=offset, S0=diag(1, (p+np)-i))
for(i in 1:m) g[[i]]<- glm(y~-1+X[,-c(1:i)], family=binomial, offset=offset)
for(i in 1:m) print(round(c(logLik(f[[i]]), f[[i]]$logLik_ML2, logLik(g[[i]]))))

par(mfrow=c(3,1))
plot(sapply(f, logLik), col=cl<-3-1*(theta[-1]==0))
plot(sapply(g, logLik), col=cl)
plot(theta[-1], g[[1]]$coef)
points(theta[-1], f[[1]]$coeffi, col=3)
abline(0,1)
s<-round(cbind(diag( f[[1]]$S), diag(summary(g[[1]])$cov.unscaled)),2)
print(head(s))
