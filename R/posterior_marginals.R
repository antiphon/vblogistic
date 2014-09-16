#' Posterior marginals of the parameters in a vblogit fit
#' 
#' Use the fact that the variables are gaussian.
#' @param x vblogit object
#' @param log If FALSE, exponentiate the values.
#' @param ... ignored
#' 
#' @export
#' @exportMethod marginals

marginals.vblogifit <- function(x, log=TRUE, ...) {
  #' extract the posterior mean and S
  m <- as.numeric(x$m)
  S <- as.matrix(x$S)
  
  varnames <- names(x$coef)
  
  den <- function(i){
    r <-  m[i]+4*c(-1,1)*sqrt(S[i,i]) 
    if(!log) r <- exp(r)
    x <- seq(r[1], r[2], length=100)
    d <- if(log) dnorm else dlnorm
    y <- d(x, m[i], sqrt(S[i,i]))
    cbind(x=x, y=y)
  }
  v <- lapply(1:length(m), den)
  names(v) <- varnames
  class(v) <- "vblogitfit_marginals"
  v
}

