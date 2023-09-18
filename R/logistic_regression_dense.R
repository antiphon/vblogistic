#' vb-logit basic
#' No sparse matrix class
#  
#' @export

vblogit_dense <- function(y, X, offset, eps=1e-2, m0, S0, S0i, xi0, verb=FALSE, maxiter=1000, ...) {
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
  N <- length(y)
  K <- ncol(X)
  #
  #
  # offset
  if(missing('offset')) offset <- 0
  if(length(offset)<N) offset <- rep(offset, N)[1:N]
  #
  #
  # Priors and initial estimates.
  if(missing(S0))S0   <- diag(1e5, K, K)
  if(missing(S0i))S0i <- solve(S0)
  if(missing(m0))m0   <- rep(0, K)
  # Constants:
  oo2 <- offset^2
  LE_CONST <- as.numeric( -0.5*t(m0)%*%S0i%*%m0 - 0.5*determinant(S0)$mod + sum((y-0.5)*offset) ) 
  Sm0 <- S0i%*%m0
  # start values for xi:
  if(missing(xi0))xi0   <- rep(4, N) # something positive
  if(length(xi0)!=N) xi0 <- rep(xi0, N)[1:N]
  
  est <- list(m=m0, S=S0, Si=S0i, xi=xi0)
  #
  #
  ## helper functions needed:
  lambda <- function(x)  -tanh(x/2)/(4*x)
  gamma <- function(x)  x/2 - log(1+exp(x)) + x*tanh(x/2)/4
  ###
  ## loop
  le <- -Inf
  le_hist <- le
  loop <- TRUE
  iter <- 0
  # initials:
  la <- lambda(xi0)
  Si <- S0i - 2 * t(X*la)%*%X
  S <- solve(Si)
  m <- S%*%( t(X)%*%( (y-0.5) + 2*la*offset ) + Sm0  )
  #
  # Main loop:
  while(loop){
    old <- le
    #  update variational parameters
    xi2 <- diag(  X%*%( S+m%*%t(m) )%*%t(X)+ 2*(X%*%m)%*%t(offset) ) + oo2
    xi <- sqrt(xi2)
    la <- lambda(xi)
    #  update post covariance
    Si <- S0i - 2 * t(X*la)%*%X
    S <- solve(Si)
    #  update post mean
    m <- S%*%( t(X)%*%( (y-0.5) + 2*la*offset ) + Sm0  )
    #  compute the log evidence
    le <-  as.numeric( 0.5*determinant(S)$mod + sum( gamma(xi) ) + sum(oo2*la) + 0.5*t(m)%*%Si%*%m + LE_CONST)
    #  check convergence 
    devi <- le - old
    if(devi < 0) warning("Log-evidence decreasing, try different starting values for xi.")
    loop <- abs(devi) > eps & (iter<-iter+1) <= maxiter
    le_hist <- c(le_hist, le)
    cat2("diff:", devi, "             \r")
  }
  if(iter == maxiter) warning("Maximum iteration limit reached.")
  cat2("\n")
  ## done. Compile:
  est <- list(m=m, S=S, Si=Si, xi=xi, lambda_xi=la)
  #  Marginal evidence
  est$logLik <- le
  #  Compute max logLik with the Bernoulli model, this should be what glm gives:
  est$logLik_ML <- as.numeric( t(y)%*%(X%*%m+offset) - sum( log( 1 + exp(X%*%m+offset)) ) )
  #  Max loglik with the approximation
  est$logLik_ML2 <- as.numeric(  t(y)%*%(X%*%m + offset)  + 
                                   t(m)%*%t(X*la)%*%X%*%m - 
                                   0.5*sum(X%*%m) + sum(gamma(xi)) +
                                   2*t(offset*la)%*%X%*%m + 
                                   t(offset*la)%*%offset - 
                                   0.5 * sum(offset)  )
  #  some additional parts, like in glm output
  est$coefficients <- est$m[,1]
  names(est$coefficients) <- varnames
  est$call <- sys.call()
  est$converged <- !(maxiter==iter)
  #  more additional stuff
  est$logp_hist <- le_hist
  est$parameters <- list(eps=eps, maxiter=maxiter)
  est$priors <- list(m=m0, S=S0)
  est$iterations <- iter
  class(est) <- "vblogit"
  ## return
  est
}

