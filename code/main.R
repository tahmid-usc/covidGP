library(dplyr)
library(lubridate)
library(readr)
library(optimx)
library(kernlab)
library(mvtnorm)
source('code/functions.R')

#run appropriate portion of data.R code file to load the data


# Try to estimmate the intial values using nonlinear least square
dt <- cbind(y = covid$y, t = covid$t)
dt <- as.data.frame(dt)
nonlin <- nls(y ~ a / (1 + b0 * exp(- b1 * t))^v , data = dt, start = list(a = 1, b0 = 1, b1 = 1, v = 1))

#Estimate hyperparameters (single starting value)
theta <- Hyper(x = covid$t, y = covid$y, init.val = rep(1,1,1, coef(nonlin)))


#Estimate hyperparamters (multistarts) (works better)


theta.ms <- Hyper.ms(x = covid$t, y = covid$y)
theta.sort <- theta.ms[order(theta.ms$value),]
theta <- theta.sort[1,1:7] #might change set of values if the best doesnt work well
theta <- as.numeric(theta)

#check if fitting is good
tstar <- seq(0,max(covid$t)*1.5, length.out = 100)
mu.pred <- mu(covid$t,a = theta[4], b0 = theta[5], b1 = theta[6], v = theta[7])
mustar.pred <- mu(tstar,a = theta[4], b0 = theta[5], b1 = theta[6], v = theta[7])

plot(tstar, mustar.pred, type = 'l', lwd = 2, col = 1, ylim = c(0, max(1,theta[4])), main = 'Mean function estimation')
points(covid$t, covid$y, col = 2)
legend('bottomright', c( 'Fitted',  'True'), lty = 1, lwd = 3, col = 1:2)


#save(theta.ms, file = 'theta_multistart_v.Rdata')

#load("theta_multistart_v.Rdata") # load for pretrained model

#---------- If hyperparamters are okay then do final prediction

datestar <- covid$date
datestar <- c(datestar, seq.Date(from = max(datestar), by = "day", length.out = 90))
tstar <- (as.numeric(datestar) - mindate)/100

mu.pred <- mu(covid$t,a = theta[4], b0 = theta[5], b1 = theta[6], v = theta[7])
mustar.pred <- mu(tstar,a = theta[4], b0 = theta[5], b1 = theta[6], v = theta[7])


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
#save(preddt, file = 'pred_SC_v.csv')


#plot predictions
plot(datestar, posmu, lwd = 3, type = 'l', col = 2, main = 'Cumulative COVID-19 cases Prediction for South Carolina', xaxt = 'n', ylim = c(0, max(ul)), xlab = "Day", ylab = 'Cumulative cases')
polygon(c(datestar,rev(datestar)),c(ll,rev(ul)),col=rgb (.1,.1,.1,.2),border=NA, xaxt = 'n')
lines(datestar, posmu, lwd = 3, col = 2, xaxt = 'n')
axis.Date(x = datestar, side = 1, format = '%m-%d', at = datestar)
points(covid$date, covid$y * max(covid$cases), cex = 1.2)
lines(datestar, mustar.pred, lwd = 3)
legend('bottomright', c('Parametric fit', 'Posterior Mean'), lty = 1, 
       lwd = 3, col = 1:2, bty = 'n')


