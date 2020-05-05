library(mvtnorm)
library(kernlab)
library(optimx)





# RBF kernel

ker <- function(x, l, sigf) {
  rbf <- rbfdot(sigma = 1/l^2)
  return(sigf^2 * kernelMatrix(rbf, x = x))
}

ker2 <- function(x, y, l, sigf) {
  rbf <- rbfdot(sigma = 1/l^2)
  return(sigf^2 * kernelMatrix(rbf, x = x, y = y))
}




mu <- function(t, a = 1, b0 = 1, b1 = 1, v = 1) {
  a / (1 + b0 * exp(- b1 * t))^(v^2)
}


#--- Multistart
parseq <- c(.01,.1,1)
parmat <- expand.grid(parseq, parseq, parseq, parseq, parseq, parseq, parseq)

Hyper.ms <- function(x, y) {
  
  marlik <- function(theta) {
    x <- as.matrix(x)
    n <- dim(x)[1]
    #theta <- theta^2
    k <- ker(x = x, l = theta[1], sigf = theta[2])
    -dmvnorm(x = y, mean = mu(x, a = theta[4], b0 = theta[5], b1 = theta[6], v = theta[7]),  sigma = k + theta[3]^2 * diag(n), log = T)
  }
  
  hyp <- multistart(parmat=parmat, fn = marlik, method = 'Nelder-Mead',
                    control=list(maxit = 10000))
  #print(hyp)
  return(hyp)
  
}


fitlogistic <- function(n, tmin = 0, tmax = 7) {
  
  #simulate data
  t <- seq(tmin, tmax, length.out = n)
  mut <- mu(t,a = 1, b0 = 2, b1 = 1, v = 1)
  covmat <- ker(t, l = 1.5, sigf = .05)
  genY <- rmvnorm(1, mut, covmat) 
  genY <- genY + rnorm(n,0,.06)
  genY <- as.numeric(genY)
  
  
  #estimate hyperparamters
  dt <- cbind(genY, t)
  dt <- as.data.frame(dt)
  nonlin <- nls(genY ~ a / (1 + b0 * exp(- b1 * t))^(v^2), data = dt, start = list(a = 1, b0 = 1, b1 = 1, v = 1))
  nonlin.par <- coef(nonlin)
  
  theta.ms <- Hyper.ms(x = t, y = genY)
  theta <- theta.ms[which(theta.ms$value == min(theta.ms$value)),]
  theta <- theta[1:7]
  theta <- as.numeric(theta)
  
  
  #fit
  tstar <- seq(-15,15, length.out = 100)
  mustar <- mu(tstar, a = 1, b0 = 2, b1 = .7, v = 1)
  mu.pred <- mu(t,a = theta[4], b0 = theta[5], b1 = theta[6], v = theta[7])
  mustar.pred <- mu(tstar,a = theta[4], b0 = theta[5], b1 = theta[6], v = theta[7])
  
  n <- length(t)
  nx <- length(tstar)
  kx <- ker2(x = tstar, y = t, l = theta[1], sigf = theta[2])
  kxx <- ker(x = tstar, l = theta[1], sigf = theta[2]) + theta[3]^2 * diag(nx)
  k <- ker(x = t, l = theta[1], sigf = theta[2]) + theta[3]^2 * diag(n)
  kinv <- chol2inv(chol(k))
  posmu <- kx %*% (kinv %*% matrix(genY - mu.pred, ncol = 1))
  posmu <- pred + posmu
  
  return(list(mu.pred = mu.pred, mustar.pred = mustar.pred, posmu = posmu, 
              theta = theta, t = t, tstar = tstar, mustar = mustar, genY = genY))
  
}


lf <- fitlogistic(50)
plot(lf$tstar, lf$mustar, type = 'l', lwd = 3, col = 1, ylim = c(0, 1.1), 
     main = 'Parametric and GP based prediction of logistic growth curve', cex.main = 1.5)
lines(lf$tstar, lf$mustar.pred, lwd = 3, col = 2)
points(lf$t, lf$genY, pch = 16, cex =1.2, col = rgb(0,0,0,.2))
lines(lf$tstar, lf$posmu, lwd = 3, col = 3)
legend('bottomright', c('True', 'Parametric fit', 'Posterior mean'), lty = 1, 
       lwd = 3, col = 1:3, bty = 'n')



