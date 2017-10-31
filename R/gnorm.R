#' The generalized normal distribution
#'
#' @name gnorm
#' @description Density, distribution function and random generation for the generalized normal/exponential power distribution. \cr \cr
#' A generalized normal random variable \eqn{x} with parameters \eqn{\mu}, \eqn{\alpha > 0} and \eqn{\beta > 0} has density:\cr
#' \deqn{p(x) = \beta exp{-(|x - \mu|/\alpha)^\beta}/(2\alpha \Gamma(1/\beta)).} \cr
#' The mean and variance of \eqn{x} are \eqn{\mu} and \eqn{\alpha^2 \Gamma(3/\beta)/\Gamma(1/\beta)}, respectively.
#' @aliases dgnorm pgnorm qgnorm rgnorm
#' @usage
#' dgnorm(x, mu = 0, alpha = 1, beta = 1, log = FALSE)
#' pgnorm(q, mu = 0, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)
#' qgnorm(p, mu = 0, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)
#' rgnorm(n, mu = 0, alpha = 1, beta = 1)
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations
#' @param mu location parameter
#' @param alpha scale parameter
#' @param beta shape parameter
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X\leq x]}, otherwise \eqn{P[X> x]}
#' @source \code{dgnorm}, \code{pgnorm}, \code{qgnorm} and\code{rgnorm} are all parametrized as in Version 1 of the \href{https://en.wikipedia.org/wiki/Generalized_normal_distribution}{Generalized Normal Distribution Wikipedia page},
#' which uses the parametrization given by in Nadarajah (2005).
#' The same distribution was described much earlier by Subbotin (1923) and named the exponential power distribution by Box and Tiao (1973). \cr \cr
#' Box, G. E. P. and G. C. Tiao. "Bayesian inference in Statistical Analysis." Addison-Wesley Pub. Co., Reading, Mass (1973). \cr
#' Nadarajah, Saralees. "A generalized normal distribution." Journal of Applied Statistics 32.7 (2005): 685-694. \cr
#' Subbotin, M. T. "On the Law of Frequency of Error." Mat. Sb. 31.2 (1923):  206-301.
#' @importFrom stats pgamma qgamma rbinom runif
#' @export
dgnorm <- function(x, mu = 0, alpha = 1, beta = 1, log = FALSE) {
  if (alpha <= 0 | beta <= 0) {
    cat("Not defined for negative values of alpha and/or beta.\n")
    return(rep(NaN, length(x)))
  }
  if (!log) {
    return(exp(-(abs(x - mu)/alpha)^beta)*beta/(2*alpha*gamma(1/beta)))
  } else {
    return(-(abs(x - mu)/alpha)^beta + log(beta) - (log(2) + log(alpha) + log(gamma(1/beta))))
  }
}
#' @export
pgnorm <- function(q, mu = 0, alpha = 1, beta = 1, lower.tail = TRUE,
                   log.p = FALSE) {
  if (alpha <= 0 | beta <= 0) {
    cat("Not defined for negative values of alpha and/or beta.\n")
    return(rep(NaN, length(q)))
  }
  p <- 1/2 + sign(q - mu)*pgamma(abs(q - mu)^beta, shape = 1/beta, rate = (1/alpha)^beta)/2
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
    return(rep(NaN, p))
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

