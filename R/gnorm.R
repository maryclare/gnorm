#' The generalized normal distribution
#'
#' @description Density, distribution function and random generation for the generalized normal distribution
#'
#' @usage \code{dgnorm(x, mu = 0, alpha = 1, beta = 1, log = FALSE)}
#' @usage \code{x, mu = 0, alpha = 1, beta = 1, lower.tail = TRUE, log.p = FALSE)}
#'
#' @param \code{x} vector of quantiles
#' @param \code{n}n number of observations. If length(n) > 1, the length is taken to be the number required
#' @param \code{mu} location parameter
#' @param \code{alpha} scale parameter
#' @param \code{beta} shape parameter
#' @source \code{dgnorm}, \code{pgnorm} and\code{rgnorm} are all parametrized as in Version 1 of the following: https://en.wikipedia.org/wiki/Generalized_normal_distribution
#' @export
dgnorm <- function(x, mu = 0, alpha = 1, beta = 1, log = FALSE) {
  if (!log) {
    return(exp(-(abs(x - mu)/a)^beta)*beta/(2*alpha*gamma(1/beta)))
  } else {
    if (!log) {
      return(-(abs(x - mu)/a)^beta + log(beta) - (log(2) + log(alpha) + log(gamma(1/beta))))
    }
  }
}
#' @export
pgnorm <- function(x, mu = 0, alpha = 1, beta = 1, lower.tail = TRUE,
                   log.p = FALSE) {
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
rgnorm <- function(n, mu = 0, alpha = 1, beta = 1) {
  lambda <- 1/alpha
  sigma.sq.beta <- (gamma(3/beta)/gamma(1/beta))*lambda
  unifs <- runif(n)
  scales <- qgamma(unifs, shape = 1/beta, scale = 1)^(1/beta)
  a <- sqrt(gamma(1/beta)/gamma(3/beta))
  return(sqrt(sigma.sq.beta)*a*scales*((-1)^rbinom(n, 1, 0.5)) + mu)
}

