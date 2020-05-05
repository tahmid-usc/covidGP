

library(dplyr)
library(lubridate)



covid <- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv")
covid$date <- as.Date(covid$date)
mindate <- as.numeric(min(covid$date))
covid.total <- covid %>% group_by(date) %>% summarize(total_cases = sum(cases)) %>% 
  mutate(t = (as.numeric(date) - mindate)/100, y = total_cases/ mean(total_cases))
covid <- covid.total



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
  a / (1 + b0 * exp(- b1 * t))^v
}



#--------------------------


Hyper <- function(x, y, init.val = c(1,1,1,1,1,1,1)) {
  
  marlik <- function(theta) {
    x <- as.matrix(x)
    n <- dim(x)[1]
    k <- ker(x = x, l = theta[1], sigf = theta[2])
    -dmvnorm(x = y, mean = mu(x, a = theta[4], b0 = theta[5], b1 = theta[6], v= theta[7]),  sigma = k + theta[3]^2 * diag(n), log = T)
  }
  
  hyp <- optim(par= init.val, fn = marlik, method = 'Nelder-Mead',
               control=list(maxit = 10000))
  print(hyp)
  return(hyp$par)
  
}

#--- Multistart
#parseq <- c(.01,.1,1)
#parmat <- expand.grid(parseq, parseq, parseq, parseq, parseq, parseq, parseq)

parseq <- seq(0.01, 20, by = .1)
parmat <- matrix(parseq, ncol = 7, nrow = 20)

Hyper.ms <- function(x, y) {
  
  marlik <- function(theta) {
    x <- as.matrix(x)
    n <- dim(x)[1]
    #theta <- theta^2
    k <- ker(x = x, l = theta[1], sigf = theta[2])
    -dmvnorm(x = y, mean = mu(x, a = theta[4], b0 = theta[5], b1 = theta[6], v= theta[7]),  sigma = k + theta[3]^2 * diag(n), log = T)
  }
  
  hyp <- multistart(parmat=parmat, fn = marlik, method = 'Nelder-Mead',
                    control=list(maxit = 10000))
  #print(hyp)
  return(hyp)
  
}

#---

plot(covid$t, covid$y)

dt <- cbind(y = covid$y, t = covid$t)
dt <- as.data.frame(dt)
nonlin <- nls(y ~ a / (1 + b0 * exp(- b1 * t))^v , data = dt, start = list(a = 1, b0 = 1, b1 = 1, v = 1))
nonlin.par <- coef(nonlin)

theta <- Hyper(x = covid$t, y = covid$y, init.val = rep(1,7))


#save(theta.ms, file = 'theta_total_v.Rdata')

theta.ms <- Hyper.ms(x = covid$t, y = covid$y)
theta.sort <- theta.ms[order(theta.ms$value),]


theta <- theta.sort[1,1:7]
theta <- as.numeric(theta)


tstar <- seq(0,5, length.out = 100)
mu.pred <- mu(covid$t,a = theta[4], b0 = theta[5], b1 = theta[6], v = theta[7])
mustar.pred <- mu(tstar,a = theta[4], b0 = theta[5], b1 = theta[6], v = theta[7])

plot(tstar, mustar.pred, type = 'l', lwd = 2, col = 1, ylim = c(0, max(1,theta[4])), main = 'Mean function estimation')
points(covid$t, covid$y, col = 2)
legend('bottomright', c( 'Fitted',  'True'), lty = 1, lwd = 3, col = 1:2)


#----------

datestar <- covid$date
datestar <- c(datestar, seq.Date(from = max(datestar), by = "day", length.out = 90))
tstar <- (as.numeric(datestar) - mindate)/100

mu.pred <- mu(covid$t,a = theta[4], b0 = theta[5], b1 = theta[6], v = theta[7])
mustar.pred <- mu(tstar,a = theta[4], b0 = theta[5], b1 = theta[6], v = theta[7])

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

posmu <- posmu * mean(covid$total_cases)
mustar.pred <- mustar.pred * mean(covid$total_cases)
ll <- posmu - 1.96 * sqrt(diag(possigma)) * mean(covid$total_cases)
ul <- posmu + 1.96 * sqrt(diag(possigma)) * mean(covid$total_cases)


preddt <- data.frame(posmu, mustar.pred, ll, ul)
#save(preddt, file = 'pred_total_v.csv')



plot(datestar, posmu, lwd = 3, type = 'l', col = 2, main = 'Cumulative COVID-19 Cases Prediction for South Carolina', xaxt = 'n', ylim = c(0, max(ul)), xlab = "Day", ylab = 'Cumulative Cases')
polygon(c(datestar,rev(datestar)),c(ll,rev(ul)),col=rgb (.1,.1,.1,.2),border=NA, xaxt = 'n')
lines(datestar, posmu, lwd = 3, col = 2, xaxt = 'n')
axis.Date(x = datestar, side = 1, format = '%m-%d', at = datestar)
points(covid$date, covid$y * mean(covid$total_cases), cex = 1.2)
lines(datestar, mustar.pred, lwd = 3)
legend('bottomright', c('Parametric fit', 'Posterior Mean'), lty = 1, 
       lwd = 3, col = 1:2, bty = 'n')


