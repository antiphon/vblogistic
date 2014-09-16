#' vblogit fit print method
#' @param x vblogit object
#' @param ... ignored.
#' @exportMethod print
#' @export

print.vblogitfit <- function(x, ...) {
  cat("Variational Bayes logistic regression fit\n")
  cat("\nCall: ")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\n")
  cat("Log-likelihood:", x$logLik)
  cat("\nConverged:", x$converged)
  cat("\nConvergence threshold:", x$parameters$eps)
  cat("\nIterations / max:", x$iterations, "/", x$parameters$maxiter)
  cat("\n* Caution: the estimates are conditional on convergence.\n")
}

#' vblogit fit summary method
#' @param x vblogit object
#' @param ... ignored.
#' @exportMethod summary
#' @export

summary.vblogitfit <- function(x, ...) {
  cat("Variational Bayes logistic regression fit\n")
  cat("\nCall: ")
  print(x$call)
  cat("\nCoefficients and posterior 95% central regions:\n")
  vna <- names(x$coefficients)
  s <- sqrt(diag(x$S))
  q0 <- qnorm(c(0.025, 0.975))
  m <- as.numeric(x$m)
  df <- data.frame(estimate=m, "low 0.05"=m+s*q0[1], "high 97.5"=m+s*q0[2], "prior mean"=x$priors$m, "prior var"=diag(x$priors$S))
  rownames(df) <- vna
  print(df)
  cat("\n")
  cat("Log-likelihood:", x$logLik)
  cat("\n")
}


#' Coef
#' @param x object from vblogit
#' @param ... ignored.
#' 
#' @exportMethod coef
coef.vblogitfit <- function(x, ...) x$coefficients


#' vblogit fit plot
#' 
#' Plot something.
#' 
#' @exportMethod plot
#' @param x object from vblogit
#' @param ... passed to plot.vblogitfit_marginals.
#' @export
plot.vblogitfit <- function(x, ...){
  y <- marginals.vblogifit(x, ...)
  plot.vblogitfit_marginals(y, ...)
}

#' Marginals plot
#' @param x object from vblogit
#' @param ... ignored.
#' 
#' @exportMethod plot
#' @export
plot.vblogitfit_marginals <- function(x, ...){
  nn <- names(x)
  n <- length(nn)
  ncol <- 2
  nrow <- ceiling(n/ncol)
  par(mfrow=c(nrow, ncol))
  for(i in 1:n) plot(x[[i]], xlab=nn[i], ylab="posterior density", main=nn[i], type="l", col=4)
}
#' Marginals print
#' @param x object from vblogit
#' @param ... ignored.
#' @exportMethod print
#' @export
print.vblogitfit_marginals <- function(x, ...){
  nn <- names(x)
  cat("Marginals for variables:", paste(nn, sep=", "))
}

#' Log-evidence
#' @param x object from vblogit
#' @param ... ignored.
#' @exportMethod logLik
#' @export

logLik.vblogitfit <- function(x, ...) {
  x$logLik
}

