#' The generalized normal distribution
#'
#' @name gnorm
#'
#' @description Density, distribution function and random generation for the generalized normal distribution
#'
#' @aliases dgnorm pgnorm qgnorm rgnorm
#'
#' @usage \code{dgnorm(x, mu = 0, alpha = 1, beta = 1, log = FALSE)}
#' \code{pgnorm(q, mu = 0, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)}
#' \code{qgnorm(p, mu = 0, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)}
#' \code{rgnorm(n, mu = 0, alpha = 1, beta = 1)}
#'
#' @param \code{x} vector of quantiles
#' @param \code{p} vector of probabilities
#' @param \code{n} number of observations
#' @param \code{mu} location parameter
#' @param \code{alpha} scale parameter
#' @param \code{beta} shape parameter
#'
#' @source \code{dgnorm}, \code{pgnorm}, \code{qgnorm} and\code{rgnorm} are all parametrized as in Version 1 of the following: https://en.wikipedia.org/wiki/Generalized_normal_distribution
#' Functions and documentation modeled after GammaDist functions
#' @export
dgnorm <- function(x, mu = 0, alpha = 1, beta = 1, log = FALSE) {
  if (alpha <= 0 | beta <= 0) {
    cat("Not defined for negative values of alpha and/or beta.\n")
    return(rep(NaN, length(x)))
  }
  if (!log) {
    return(exp(-(abs(x - mu)/alpha)^beta)*beta/(2*alpha*gamma(1/beta)))
  } else {
    if (!log) {
      return(-(abs(x - mu)/alpha)^beta + log(beta) - (log(2) + log(alpha) + log(gamma(1/beta))))
    }
  }
}
#' @export
pgnorm <- function(x, mu = 0, alpha = 1, beta = 1, lower.tail = TRUE,
                   log.p = FALSE) {
  if (alpha <= 0 | beta <= 0) {
    cat("Not defined for negative values of alpha and/or beta.\n")
    return(rep(NaN, length(x)))
  }
  p <- 1/2 + sign(x - mu)*(pgamma((abs(x - mu)/alpha)^beta, 1/beta))/(2*gamma(1/beta))
  if (lower.tail) {
    if (!log.p) {
      return(p)
    } else {
      return(log(p))
    }
  } else if (!lower.tail) {
    if (!log.p) {
      return(1 - p)
    } else {
      return(log(1 - p))
    }
  }
}
#' @export
qgnorm <- function(p, mu = 0, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE) {
  if (alpha <= 0 | beta <= 0) {
    cat("Not defined for negative values of alpha and/or beta.\n")
    return(rep(NaN, n))
  }
  if (lower.tail & !log.p) {
    p <- p
  } else if (lower.tail & log.p) {
    p <- exp(p)
  } else if (!lower.tail & !log.p) {
    p <- 1 - p
  } else {
    p <- log(1 - p)
  }

  lambda <- (1/alpha)^beta
  return(sign(p - 0.5)*qgamma(abs(p - 0.5)*2, shape = 1/beta, scale = 1/lambda)^(1/beta) + mu)
}
#' @export
rgnorm <- function(n, mu = 0, alpha = 1, beta = 1) {
  if (alpha <= 0 | beta <= 0) {
    cat("Not defined for negative values of alpha and/or beta.\n")
    return(rep(NaN, n))
  }
  lambda <- (1/alpha)^beta
  unifs <- runif(n)
  scales <- qgamma(unifs, shape = 1/beta, scale = 1/lambda)^(1/beta)
  return(scales*((-1)^rbinom(n, 1, 0.5)) + mu)
}

