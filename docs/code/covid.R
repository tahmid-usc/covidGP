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




mu <- function(t, a = 1, b0 = 1, b1 = 1) {
  a / (1 + b0 * exp(- b1 * t))^2
}



t <- seq(-3,10, length.out = 100)
mut <- mu(t,a = 1, b0 = 2, b1 = 1)
plot(t, mut, type = 'l')

covmat <- ker(t, l = 1.5, sigf = .05)
ft <- rmvnorm(1, mut, covmat)
plot(t, ft, type = 'l', lwd = 3, ylim = c(-0.1,1.1), col = rgb(0,0,0,.3))
for (i in 1:10){
  ft <- rmvnorm(1, mut, covmat)
  lines(t, ft, lwd = 3, , col = rgb(0,0,0,.3))
}



t <- seq(0,7, length.out = 100)
mut <- mu(t,a = 1, b0 = 2, b1 = 1)
covmat <- ker(t, l = 1.5, sigf = .05)
gt <- rmvnorm(1, mean = rep(0, 100), sigma =  covmat)

plot(t, gt, type = 'l')


genY <- rmvnorm(1, mut, covmat) 
plot(t, genY, type = 'l')

genY <- genY + rnorm(100,0,.06)
genY <- as.numeric(genY)
plot(t, genY, type = 'l')





#--------------------------


Hyper <- function(x, y, init.val = c(1,1,1,1,1,1)) {
  
  marlik <- function(theta) {
    x <- as.matrix(x)
    n <- dim(x)[1]
    k <- ker(x = x, l = theta[1], sigf = theta[2])
    -dmvnorm(x = y, mean = mu(x, a = theta[4], b0 = theta[5], b1 = theta[6]),  sigma = k + theta[3]^2 * diag(n), log = T)
  }
  
  hyp <- optim(par= init.val, fn = marlik, method = 'Nelder-Mead',
               control=list(maxit = 10000))
  print(hyp)
  return(hyp$par)
  
}

#--- Multistart
parseq <- c(.01,.1,1)
parmat <- expand.grid(parseq, parseq, parseq, parseq, parseq, parseq)

Hyper.ms <- function(x, y) {
  
  marlik <- function(theta) {
    x <- as.matrix(x)
    n <- dim(x)[1]
    #theta <- theta^2
    k <- ker(x = x, l = theta[1], sigf = theta[2])
    -dmvnorm(x = y, mean = mu(x, a = theta[4], b0 = theta[5], b1 = theta[6]),  sigma = k + theta[3]^2 * diag(n), log = T)
  }
  
  hyp <- multistart(parmat=parmat, fn = marlik, method = 'Nelder-Mead',
               control=list(maxit = 10000))
  #print(hyp)
  return(hyp)
  
}

#---

dt <- cbind(genY, t)
dt <- as.data.frame(dt)
nonlin <- nls(genY ~ a / (1 + b0 * exp(- b1 * t))^2, data = dt, start = list(a = 1, b0 = 1, b1 = 1))
nonlin.par <- coef(nonlin)

theta <- Hyper(x = t, y = genY, init.val = c(1,1,1,nonlin.par))
#theta.ms <- Hyper.ms(x = t, y = genY)
#theta <- theta.ms[which(theta.ms$value == min(theta.ms$value)),]
#theta <- theta[1:6]
#theta <- as.numeric(theta)


tstar <- seq(-15,15, length.out = 100)
mutstar <- mu(tstar,a = 1, b0 = 2, b1 = .7)
mu.pred <- mu(t,a = theta[4], b0 = theta[5], b1 = theta[6])
mustar.pred <- mu(tstar,a = theta[4], b0 = theta[5], b1 = theta[6])

plot(tstar, mustar.pred, type = 'l', lwd = 2, col = 1, ylim = c(0, max(1,theta[4])), main = 'Mean function estimation')
lines(tstar, mutstar, lwd = 2, col = 2)
legend('bottomright', c( 'Fitted',  'True'), lty = 1, lwd = 3, col = 1:2)


n <- length(t)
nx <- length(tstar)
kx <- ker2(x = tstar, y = t, l = theta[1], sigf = theta[2])
kxx <- ker(x = tstar, l = theta[1], sigf = theta[2]) + theta[3]^2 * diag(nx)
k <- ker(x = t, l = theta[1], sigf = theta[2]) + theta[3]^2 * diag(n)
kinv <- chol2inv(chol(k))
posmu <- kx %*% (kinv %*% matrix(genY - mu.pred, ncol = 1))
posmu <- mustar.pred + posmu

plot(tstar, posmu, type = 'l', lwd = 2, col = 2, main = 'Posterior mean', ylim = c(0,1.1))
points(t, genY)
lines(tstar, mutstar, lwd = 2)

possigma <- kxx - kx %*% (kinv %*% t(kx))
#diag(possigma)[diag(sigma)<0] <- 0
ll <- posmu - 1.96 * sqrt(diag(possigma))
ul <- posmu + 1.96 * sqrt(diag(possigma))

plot(tstar, posmu, type = 'l', lwd = 2, col = 2, main = 'Posterior mean', ylim = c(min(ll), max(ul)))
points(t, genY)
lines(tstar, mutstar, lwd = 2)
lines(tstar, ll, lwd = 2, lty = 2)
lines(tstar, ul, lwd = 2, lty = 2)
