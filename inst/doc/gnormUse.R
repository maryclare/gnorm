## ---- message = FALSE----------------------------------------------------
library(devtools)
install_github("maryclare/gnorm")
library(gnorm)

## ------------------------------------------------------------------------
set.seed(100)

## ---- fig.show='hold', fig.align = 'center', fig.width = 4, fig.height = 4----
xs <- seq(-1, 1, length.out = 100)
plot(xs, dgnorm(xs, mu = 0, alpha = sqrt(2), beta = 2), type = "l", 
     xlab = "x", ylab = expression(p(x)))

## ---- fig.show='hold', fig.align = 'center', fig.width = 4, fig.height = 4----
xs <- seq(-1, 1, length.out = 100)
plot(xs, pgnorm(xs, 0, sqrt(2), 2), type = "l", xlab = "q", ylab = expression(paste("Pr(", x<=q, ")", sep = "")))

