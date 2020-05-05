library(dplyr)
library(lubridate)



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


covid <- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv")
covid$date <- as.Date(covid$date)
covid <- covid %>% filter(state == "South Carolina") %>%  
  mutate(cumCase = cumsum(cases), t = (as.numeric(date) - 18327)/100, y = cases/ max(cases))


plot(covid$t, covid$y)

dt <- cbind(y = covid$y, t = covid$t)
dt <- as.data.frame(dt)
nonlin <- nls(y ~ a / (1 + b0 * exp(- b1 * t))^2, data = dt, start = list(a = 1, b0 = 1, b1 = .5))
nonlin.par <- coef(nonlin)

theta <- Hyper(x = covid$t, y = covid$y, init.val = c(.1,.1,.1,coef(nonlin)))


theta.ms <- Hyper.ms(x = t, y = genY)
theta <- theta.ms[which(theta.ms$value == min(theta.ms$value)),]
theta <- theta[1:6]
theta <- as.numeric(theta)

#save(theta, file = 'theta_multistart.Rdata')

tstar <- seq(0,1, length.out = 100)
mu.pred <- mu(covid$t,a = theta[4], b0 = theta[5], b1 = theta[6])
mustar.pred <- mu(tstar,a = theta[4], b0 = theta[5], b1 = theta[6])

plot(tstar, mustar.pred, type = 'l', lwd = 2, col = 1, ylim = c(0, max(1,theta[4])), main = 'Mean function estimation')
points(covid$t, covid$y, col = 2)
legend('bottomright', c( 'Fitted',  'True'), lty = 1, lwd = 3, col = 1:2)


n <- length(covid$t)
nx <- length(tstar)
kx <- ker2(x = tstar, y = covid$t, l = theta[1], sigf = theta[2])
kxx <- ker(x = tstar, l = theta[1], sigf = theta[2]) + theta[3]^2 * diag(nx)
k <- ker(x = covid$t, l = theta[1], sigf = theta[2]) + theta[3]^2 * diag(n)
kinv <- chol2inv(chol(k))
posmu <- kx %*% (kinv %*% matrix(covid$y - mu.pred, ncol = 1))
posmu <- mustar.pred + posmu

plot(tstar, posmu, type = 'l', lwd = 3, col = 2, main = 'Posterior mean')
points(covid$t, covid$y)
lines(tstar, mustar.pred, lwd = 3, cex = 1.2)
legend('bottomright', c('Parametric fit', 'Posterior Mean'), lty = 1, 
       lwd = 3, col = 1:2, bty = 'n')



#----------

datestar <- covid$date
datestar <- c(datestar, seq.Date(from = max(datestar), by = "day", length.out = 90))
tstar <- (as.numeric(datestar) - 18327)/100

mu.pred <- mu(covid$t,a = theta[4], b0 = theta[5], b1 = theta[6])
mustar.pred <- mu(tstar,a = theta[4], b0 = theta[5], b1 = theta[6])

plot(tstar, mustar.pred, type = 'l', lwd = 2, col = 1, ylim = c(0, max(1,theta[4])), main = 'Mean function estimation')
points(covid$t, covid$y, col = 2)
legend('bottomright', c( 'Fitted',  'True'), lty = 1, lwd = 3, col = 1:2)


n <- length(covid$t)
nx <- length(tstar)
kx <- ker2(x = tstar, y = covid$t, l = theta[1], sigf = theta[2])
kxx <- ker(x = tstar, l = theta[1], sigf = theta[2]) + theta[3]^2 * diag(nx)
k <- ker(x = covid$t, l = theta[1], sigf = theta[2]) + theta[3]^2 * diag(n)
kinv <- chol2inv(chol(k))
posmu <- kx %*% (kinv %*% matrix(covid$y - mu.pred, ncol = 1))
posmu <- mustar.pred + posmu


possigma <- kxx - kx %*% (kinv %*% t(kx))
#diag(possigma)[diag(sigma)<0] <- 0

posmu <- posmu * max(covid$cases)
mustar.pred <- mustar.pred * max(covid$cases)
ll <- posmu - 1.96 * sqrt(diag(possigma)) * max(covid$cases)
ul <- posmu + 1.96 * sqrt(diag(possigma)) * max(covid$cases)


preddt <- data.frame(posmu, mustar.pred, ll, ul)
save(preddt, file = 'pred_SC')



plot(datestar, posmu, lwd = 3, type = 'l', col = 2, main = 'Cumulative COVID-19 Cases Prediction for South Carolina', xaxt = 'n', ylim = c(0, max(ul)), xlab = "Day", ylab = 'Cumulative Cases')
polygon(c(datestar,rev(datestar)),c(ll,rev(ul)),col=rgb (.1,.1,.1,.2),border=NA, xaxt = 'n')
lines(datestar, posmu, lwd = 3, col = 2, xaxt = 'n')
axis.Date(x = datestar, side = 1, format = '%m-%d', at = datestar)
points(covid$date, covid$y * max(covid$cases), cex = 1.2)
lines(datestar, mustar.pred, lwd = 3)
legend('bottomright', c('Parametric fit', 'Posterior Mean'), lty = 1, 
       lwd = 3, col = 1:2, bty = 'n')


