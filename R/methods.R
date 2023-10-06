#' vblogit fit print method
#' @param x vblogit object
#' @param ... ignored.
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
#' Summary method for logistic regression with approximately Normal posterior
#' 
#' @param x vblogit object
#' @param ... ignored.
#' @importFrom Matrix diag
#' @export
summary.vblogitfit <- function(object, ...) {
  x <- object
  vna <- names(x$coefficients)
  s <- sqrt(diag(x$S))
  q0 <- qnorm(c(0.025, 0.975))
  m <- as.numeric(x$m)
  df <- data.frame(estimate=m, sd=s, "low 0.05"=m+s*q0[1], "high 97.5"=m+s*q0[2], "prior mean"=x$priors$m, "prior var"=diag(x$priors$S))
  rownames(df) <- vna
  x$table <- df
  class(x) <- c("sum_vblogitfit", class(x))
  x
}

#' Print for vblogit fit summary
#' @param x vblogit object
#' @param ... ignored.
#' @export

print.sum_vblogitfit <- function(x, ...) {
  cat("Variational Bayes logistic regression fit\n")
  cat("\nCall: ")
  print(x$call)
  cat("\nCoefficients and posterior 95% central regions:\n")
  print(x$table)
  cat("\n")
  cat("Log-likelihood:", x$logLik)
  cat("\n")
  
}

#' Coef
#' @param x object from vblogit
#' @param ... ignored.
#' 
#' @export
coef.vblogitfit <- function(x, ...) x$coefficients


#' vblogit fit plot
#' 
#' Plot something.
#' 
#' @param x object from vblogit
#' @param ... passed to plot.vblogitfit_marginals.
#' @export
plot.vblogitfit <- function(x, ...){
  y <- marginals.vblogifit(x, ...)
  plot.vblogitfit_marginals(y, ...)
}

#' Marginals plot
#' @param x object from vblogit
#' @param ncol default:2
#' @param ... ignored.
#' 
#' @export
plot.vblogitfit_marginals <- function(x, ncol = 2, ...){
  nn <- names(x)
  n <- length(nn)
  nrow <- ceiling(n/ncol)
  par(mfrow=c(nrow, ncol))
  for(i in 1:n) plot(x[[i]], xlab=nn[i], ylab="posterior density", main=nn[i], type="l", col=4)
}
#' Marginals print
#' @param x object from vblogit
#' @param ... ignored.
#' @export
print.vblogitfit_marginals <- function(x, ...){
  nn <- names(x)
  cat("Marginals for variables:", paste(nn, sep=", "))
}

#' Log-evidence
#' @param x object from vblogit
#' @param ... ignored.
#' @export

logLik.vblogitfit <- function(x, ...) {
  x$logLik
}

#' Var-Cov Method
#' 
#' @param object 'vblogitfit' object 
#' @param ... not used 
#' @export
vcov.vblogitfit <- function(object, ...) {
  S <- as.matrix(object$S)
  # dimnames
  n <- names(coef(object))
  rownames(S) <- colnames(S) <- n
  S
}

