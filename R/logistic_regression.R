#' Fit logistic regression model using VB approximation
#' 
#' Bayesian fit of logistic regression model. p coefficients, n observations.
#' 
#' @param y binary vector of responses, length n
#' @param X n x p matrix of covariates, including 1-column for intercept
#' @param offset n-vector of offsets (or 1-vector which will be replicated)
#' @param eps convergence criterion, increase in log-likelihood is no more than this
#' @param m0 p-vector of prior means
#' @param S0 p x p prior covariance matrix
#' @param xi0 p-vector of initial 
#' @param verb verbose output, logical
#' @param maxiter upper limit for iterations
#' @param ... ignored.
#' 
#' Computes the posterior distribution of regression coefficients in logistic regression
#' using the method of Jaakkola&Jordan 1996.
#' 
#' @examples
#' ## some data
#' n <- 100
#' p <- 10
#' X <- matrix( rnorm(n*p), ncol=p)
#' theta <- rnorm(p)
#' prob <- 1/(1+exp(-X%*%theta))
#' y <- rbinom(n, 1, prob)
#' 
#' ## See that it works:
#' ## vb:
#' fit_vb <- vblogit(y, X, verb=TRUE)
#' ## glm:
#' fit_glm <- glm(y ~ -1+X, family=binomial)
#' 
#' coefs <- cbind(vb=fit_vb$coef, glm=fit_glm$coef)
#' 
#' summary(fit_vb)
#' 
#' ## compare vb and glm
#' plot(coefs, main="Estimates")
#' abline(0,1)
#' 
#' ## Compare to true coefficients
#' plot(coefs[,1]-theta)
#' points(coefs[,2]-theta, col=3, pch=4)
#' abline(h=0)
#' legend("topright", c("glm","vblogit"), col=c(1,3), pch=c(1,4))
#' 
#' @import Matrix 
#' @export

vblogit<-function(y, X, offset, eps=1e-2, m0, S0, S0i, xi0, verb=FALSE, maxiter=1000, ...) {
  ### Logistic regression using JJ96 idea. Ormeron00 notation.
  ## p(y, w, t) = p(y | w) p(w | t) p(t) 
  ##
  ## Y ~ Bern(logit(Xw + offset))
  ## w  ~ N(m0, S0) iid
  ##
  ## "*0" are fixed priors.
  ##
  cat2 <- if(verb) cat else function(...) NULL
  varnames <- colnames(data.frame(as.matrix(X[1:2,])))
  
  ## Write 
  X <- Matrix(X)
  X <- drop0(X)
  N <- length(y)
  K <- ncol(X)
  #'
  #'
  #' offset
  if(missing('offset')) offset <- 0
  if(length(offset)<N) offset <- rep(offset, N)[1:N]
  #
  #'
  #' Priors and initial estimates.
  if(missing(S0))S0   <- Diagonal(1e5, n=K)
  if(missing(S0i))S0i <- solve(S0)
  if(missing(m0))m0   <- rep(0, K)
  #' Constants:
  OO <- offset%*%t(offset)
  oo2 <- offset*offset
  LE_CONST <- -0.5*t(m0)%*%S0i%*%m0 -0.5*determinant(S0)$mod - 0.5*sum(oo2)
  Sm0 <- S0i%*%m0
  #' start values for xi:
  if(missing(xi0))xi0   <- rep(4, N)
  if(length(xi0)!=N) xi0 <- rep(xi0, N)[1:N]
  
  est <- list(m=m0, S=S0, Si=S0i, xi=xi0)
  #'
  #
  ## helper functions needed:
  lambda <- function(x)  -tanh(x/2)/(4*x)
  gamma <- function(x)  x/2 - log(1+exp(x)) + x*tanh(x/2)/4
  #' Use Matrices
  m <- Matrix(m0)
  S <- Matrix(S0)
  Si <- Matrix(S0i)
  xi <- xi0
  ###
  ## loop
  le <- -Inf
  le_hist <- le
  loop <- TRUE
  iter <- 0
  la <- lambda(xi)
  while(loop){
    old <- le
    # update all
    #est <- update(est)
    L <- Diagonal( x = la )
    #' update post covariance
    Si <- S0i - 2 * t(X)%*%L%*%X
    S <- solve(Si)
    #' update post mean
    m <- S%*%( t(X)%*%( (y-0.5) + 2*L%*%offset ) + Sm0  )
    #' update variational parameters
    xi2 <- diag(  X%*%(S+m%*%t(m))%*%t(X)+ 2*(X%*%m)%*%t(offset) ) + oo2
    xi <- sqrt(xi2)
    la <- lambda(xi)
    ## compute the log evidence
    le <-  as.numeric(  0.5*determinant(S)$mod  + sum( gamma(xi) ) + 
      sum(oo2*la) + 0.5*t(m)%*%Si%*%m + LE_CONST    )
    #' check convergence
    d <- le - old
    if(d < 0) warning("Log-evidence decreasing, try different starting values for xi.")
    loop <- abs(d) > eps & (iter<-iter+1) <= maxiter
    le_hist <- c(le_hist, le)
    cat2("diff:", d, "             \r")
  }
  if(iter == maxiter) warning("Maximum iteration limit reached.")
  cat2("\n")
  ## done. Compile:
  est <- list(m=m, S=S, Si=Si, xi=xi, L=L)
  est$logLik <- le[1]
  est$logp_hist <- le_hist
  
  est$coefficients <- est$m[,1]
  names(est$coefficients) <- varnames
  est$call <- sys.call()
  est$parameters <- list(eps=eps, maxiter=maxiter)
  est$priors <- list(m=m0, S=S0)
  est$converged <- !(maxiter==iter)
  est$iterations <- iter
  class(est) <- "vblogitfit"
  ## return
  
  est
}


