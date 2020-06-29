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
#' Subbotin, M. T. "On the Law of Frequency of Error." Matematicheskii Sbornik 31.2 (1923):  206-301.
#' @importFrom stats pgamma qgamma rbinom runif
#' @examples
#' # Plot generalized normal/exponential power density
#' # that corresponds to the standard normal density
#' xs <- seq(-1, 1, length.out = 100)
#' plot(xs, dgnorm(xs, mu = 0, alpha = sqrt(2), beta = 2), type = "l",
#'      xlab = "x", ylab = expression(p(x)))
#'
#' # Plot the generalized normal/exponential power CDF
#' # that corresponds to the standard normal CDF
#' s <- seq(-1, 1, length.out = 100)
#' plot(xs, pgnorm(xs, 0, sqrt(2), 2), type = "l", xlab = "q",
#'      ylab = expression(paste("Pr(", x<=q, ")", sep = "")))
#'
#' # Plot the generalized normal/exponential power inverse CDF
#' # that corresponds to the standard normal inverse CDF
#' xs <- seq(0, 1, length.out = 100)
#' plot(xs, qgnorm(xs, 0, sqrt(2), 2), type = "l", xlab = "p",
#'      ylab = expression(paste("q: p = Pr(", x<=q, ")", sep = "")))
#'
#' # Make a histogram of draws from the generalized normal/exponential
#' # power distribution that corresponds to a standard normal distribution
#' xs <- rgnorm(100, 0, sqrt(2), 2)
# 'hist(xs, xlab = "x", freq = FALSE, main = "Histogram of Draws")
#'
#' @export
dgnorm <- function(x, mu = 0, alpha = 1, beta = 1,
                   log = FALSE) {
  maxLength <- max(length(x),length(mu),length(alpha),length(beta))
  # Create a vector to return and a vector of logical
  gnormValues <- vector("numeric",maxLength)
  valuesToUse <- vector("logical",maxLength)
  x <- rep(x,maxLength)[1:maxLength]
  mu <- rep(mu,maxLength)[1:maxLength]
  alpha <- rep(alpha,maxLength)[1:maxLength]
  beta <- rep(beta,maxLength)[1:maxLength]
  
  # A failsafe for NaN / NAs of alpha / beta
  if(any(is.nan(alpha))){
    alpha[is.nan(alpha)] <- 0;
  }
  if(any(is.nan(beta))){
    beta[is.nan(beta)] <- 0;
  }
  if(any(is.na(alpha))){
    alpha[is.na(alpha)] <- 0;
  }
  if(any(is.na(beta))){
    beta[is.na(beta)] <- 0;
  }
  
  # Do checks and substitute the unacceptable values by -Infinity
  if (any(alpha <= 0) || any(beta <= 0)){
    valuesToUse[] <- !(alpha <= 0 | beta <= 0)
    gnormValues[!valuesToUse] <- NaN
    warning("NaNs produced")
  }
  if (!log) {
    gnormValues[valuesToUse] <- (exp(-(abs(x[valuesToUse]-mu[valuesToUse])/
                                         alpha[valuesToUse])^beta[valuesToUse])*
                                   beta[valuesToUse]/(2*alpha[valuesToUse]*gamma(1/beta[valuesToUse])))
  } else {
    gnormValues[valuesToUse] <- -((abs(x[valuesToUse] - mu[valuesToUse])/alpha[valuesToUse])^beta[valuesToUse] +
                                    log(beta[valuesToUse]) - (log(2) + log(alpha[valuesToUse]) +
                                                                log(gamma(1/beta[valuesToUse]))))
  }
  return(gnormValues)
}

#' @export
pgnorm <- function(q, mu = 0, alpha = 1, beta = 1,
                   lower.tail = TRUE, log.p = FALSE) {
  maxLength <- max(length(q),length(mu),length(alpha),length(beta))
  # Create a vector to return and a vector of logical
  p <- vector("numeric",maxLength)
  valuesToUse <- vector("logical",maxLength)
  q <- rep(q,maxLength)[1:maxLength]
  mu <- rep(mu,maxLength)[1:maxLength]
  alpha <- rep(alpha,maxLength)[1:maxLength]
  beta <- rep(beta,maxLength)[1:maxLength]
  
  # A failsafe for NaN / NAs of alpha / beta
  if(any(is.nan(alpha))){
    alpha[is.nan(alpha)] <- 0;
  }
  if(any(is.nan(beta))){
    beta[is.nan(beta)] <- 0;
  }
  if(any(is.na(alpha))){
    alpha[is.na(alpha)] <- 0;
  }
  if(any(is.na(beta))){
    beta[is.na(beta)] <- 0;
  }
  # Do checks and substitute the unacceptable values by -Infinity
  if (any(alpha <= 0) || any(beta <= 0)){
    valuesToUse[] <- !(alpha <= 0 | beta <= 0)
    p[!valuesToUse] <- NaN
    warning("NaNs produced")
  }
  
  p[valuesToUse] <- (1/2 + sign(q[valuesToUse] - mu[valuesToUse])*
                       pgamma(abs(q[valuesToUse] - mu[valuesToUse])^beta[valuesToUse],
                              shape = 1/beta[valuesToUse],
                              rate = (1/alpha[valuesToUse])^beta[valuesToUse])/2)
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
qgnorm <- function(p, mu = 0, alpha = 1, beta = 1,
                   lower.tail = TRUE, log.p = FALSE) {
  maxLength <- max(length(p),length(mu),length(alpha),length(beta))
  # Create a vector to return and a vector of logical
  gnormValues <- vector("numeric",maxLength)
  valuesToUse <- vector("logical",maxLength)
  p <- rep(p,maxLength)[1:maxLength]
  mu <- rep(mu,maxLength)[1:maxLength]
  alpha <- rep(alpha,maxLength)[1:maxLength]
  beta <- rep(beta,maxLength)[1:maxLength]
  
  # A failsafe for NaN / NAs of alpha / beta
  if(any(is.nan(alpha))){
    alpha[is.nan(alpha)] <- 0;
  }
  if(any(is.nan(beta))){
    beta[is.nan(beta)] <- 0;
  }
  if(any(is.na(alpha))){
    alpha[is.na(alpha)] <- 0;
  }
  if(any(is.na(beta))){
    beta[is.na(beta)] <- 0;
  }
  # Do checks and substitute the unacceptable values by -Infinity
  if (any(alpha <= 0) || any(beta <= 0)){
    valuesToUse[] <- !(alpha <= 0 | beta <= 0)
    gnormValues[!valuesToUse] <- NaN
    warning("NaNs produced")
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
  
  lambda <- (1/alpha[valuesToUse])^beta[valuesToUse]
  gnormValues[valuesToUse] <- (sign(p[valuesToUse]-0.5)*qgamma(abs(p[valuesToUse] - 0.5)*2,
                                                               shape = 1/beta[valuesToUse],
                                                               scale = 1/lambda)^(1/beta[valuesToUse]) +
                                 mu[valuesToUse])
  return(gnormValues)
}
#' @export
rgnorm <- function(n, mu = 0, alpha = 1, beta = 1) {
  maxLength <- max(length(mu),length(alpha),length(beta))*n;
  # Create a vector to return and a vector of logical
  gnormValues <- vector("numeric",maxLength)
  valuesToUse <- vector("logical",maxLength)
  mu <- rep(mu,maxLength)[1:maxLength]
  alpha <- rep(alpha,maxLength)[1:maxLength]
  beta <- rep(beta,maxLength)[1:maxLength]
  
  # A failsafe for NaN / NAs of alpha / beta
  if(any(is.nan(alpha))){
    alpha[is.nan(alpha)] <- 0;
  }
  if(any(is.nan(beta))){
    beta[is.nan(beta)] <- 0;
  }
  if(any(is.na(alpha))){
    alpha[is.na(alpha)] <- 0;
  }
  if(any(is.na(beta))){
    beta[is.na(beta)] <- 0;
  }
  # Do checks and substitute the unacceptable values by -Infinity
  if (any(alpha <= 0) || any(beta <= 0)){
    valuesToUse[] <- !(alpha <= 0 | beta <= 0)
    gnormValues[!valuesToUse] <- NaN
    warning("NaNs produced")
  }
  
  lambda <- (1/alpha[valuesToUse])^beta[valuesToUse]
  unifs <- runif(n)
  scales <- qgamma(unifs, shape = 1/beta[valuesToUse], scale = 1/lambda)^(1/beta[valuesToUse])
  gnormValues[valuesToUse] <- scales*((-1)^rbinom(n, 1, 0.5)) + mu[valuesToUse]
  return(gnormValues)
}

