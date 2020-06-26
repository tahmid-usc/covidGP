library(dplyr)
library(lubridate)
library(readr)
library(optimx)
library(ggplot2)
library(kernlab)
library(mvtnorm)
source('code/functions.R')

covid <- covid_data("Ohio")
mindate <- as.numeric(min(covid$date))
maxNum <- max(covid$cases)
covid <- covid %>% mutate(t = (as.numeric(date) - mindate)/100, y = cases / maxNum)



# Try to estimmate the intial values using nonlinear least square
dt <- cbind(y = covid$y, t = covid$t)
dt <- as.data.frame(dt)
nonlin <- nls(y ~ a / (1 + v * exp(- b1 * (t - b0)))^(1/v) , 
              data = dt, start = list(a = 2, b0 = 2, b1 = 2, v = 2))

#Estimate hyperparameters (single starting value)
theta <- Hyper(x = covid$t, y = covid$y, init.val = rep(1,7))


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

posmu <- posmu * maxNum
mustar.pred <- mustar.pred * maxNum
ll <- posmu - 1.96 * sqrt(diag(possigma)) * maxNum
ul <- posmu + 1.96 * sqrt(diag(possigma)) * maxNum


preddt <- data.frame(pred_date = datestar, posterior_mean = posmu, parametric_fit = mustar.pred, lower = ll, upper = ul)
#save(preddt, file = 'pred_SC_v.csv')


#plot predictions
plot(datestar, posmu, lwd = 3, type = 'l', col = 2, main = 'Cumulative COVID-19 cases Prediction for Ohio', xaxt = 'n', ylim = c(0, max(ul)), xlab = "Day", ylab = 'Cumulative cases')
polygon(c(datestar,rev(datestar)),c(ll,rev(ul)),col=rgb (.1,.1,.1,.2),border=NA, xaxt = 'n')
lines(datestar, posmu, lwd = 3, col = 2, xaxt = 'n')
axis.Date(x = datestar, side = 1, format = '%m-%d', at = datestar)
points(covid$date, covid$y * max(covid$cases), cex = 1.2)
lines(datestar, mustar.pred, lwd = 3)
legend('bottomright', c('Parametric fit', 'Posterior Mean'), lty = 1, 
       lwd = 3, col = 1:2, bty = 'n')
#_----------
#ggoplot


maxVal <- round(theta[4] * maxNum,0)


ggplot(data = preddt, aes(x = pred_date, y = posterior_mean)) +
  geom_line(aes(x = pred_date, y = parametric_fit, colour = "Parametric Fit"), size = 1.2) +
  geom_smooth(aes(ymin = lower, ymax = upper), stat = "identity", fill = "skyblue") +   
  geom_line(size = 1.2, aes( colour = "Posterior Mean"), alpha = .9) +
  scale_x_date(date_breaks = "week", date_labels = "%b %d") +
  scale_y_continuous(n.breaks = 20) + 
  geom_point(data = covid, aes(x = date, y = cases), size = 1, colour = 1 , alpha = .5) + 
  geom_hline(yintercept = theta[4] * max(covid$cases), colour = 'grey50') +
  annotate("text", x = min(covid$date) , y = .99 * maxVal , label = paste(maxVal), colour = "grey30", size = 2) +
  labs(x = "Date", y = "Total Cases", colour = "") + ggtitle("Cumulative COVID-19 cases Prediction for Ohio") + 
  ggsave("plot/test.png", height = 4, width = 12)

