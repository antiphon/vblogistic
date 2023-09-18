#' Posterior marginals of the parameters in a vblogit fit
#' 
#' Use the fact that the variables are gaussian.
#' @param x vblogit object
#' @param which If NULL, create all marginals. Otherwise the indices to create.
#' @param log If FALSE, exponentiate the values.
#' @param ... ignored
#' 
#' @export

marginals.vblogifit <- function(x, which = NULL, log=TRUE, ...) {
  # extract the posterior mean and S
  m <- as.numeric(x$m)
  S <- as.matrix(x$S)
  if(is.null(which)) which <- 1:length(m)
  varnames <- names(x$coef)
  

  den <- function(i){
    r <-  m[i]+4*c(-1,1)*sqrt(S[i,i]) 
    if(!log) r <- exp(r)
    x <- seq(r[1], r[2], length=100)
    d <- if(log) dnorm else dlnorm
    y <- d(x, m[i], sqrt(S[i,i]))
    cbind(x=x, y=y)
  }
  v <- lapply(which, den)
  names(v) <- varnames[which]
  class(v) <- "vblogitfit_marginals"
  v
}

