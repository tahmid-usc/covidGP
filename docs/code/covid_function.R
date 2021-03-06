

library(dplyr)
library(lubridate)
library(readr)
library(optimx)
library(kernlab)
library(mvtnorm)
source('code/functions.R')



covidGP <- function(State = "South Carolina") {
  # Read data and transform date and cases
  # change state name as desired
  
  covid <- covid_data(State) 
  mindate <- as.numeric(min(covid$date))
  covid <- covid %>% mutate(t = (as.numeric(date) - mindate)/100, y = cases / max(cases))
  
  print('Data Loaded')
  #Estimate hyperparamters (multistarts) (works better)
  #theta <- Hyper(x = covid$t, y = covid$y)
  theta.ms <- Hyper.ms(x = covid$t, y = covid$y)
  theta.sort <- theta.ms[order(theta.ms$value),]
  theta <- theta.sort[1,1:7] #might change set of values if the best doesnt work well
  theta <- as.numeric(theta)
  
  print('Hyperparameters estimated')
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
  
  print('Posterior distribution computed')
  

  #plot predictions
  plot(datestar, posmu, lwd = 3, type = 'l', col = 2, main = paste('Cumulative COVID-19 cases Prediction for', State), xaxt = 'n', ylim = c(0, max(ul)), xlab = "Day", ylab = 'Cumulative cases')
  polygon(c(datestar,rev(datestar)),c(ll,rev(ul)),col=rgb (.1,.1,.1,.2),border=NA, xaxt = 'n')
  lines(datestar, posmu, lwd = 3, col = 2, xaxt = 'n')
  axis.Date(x = datestar, side = 1, format = '%m-%d', at = datestar)
  points(covid$date, covid$y * max(covid$cases), cex = 1.2)
  lines(datestar, mustar.pred, lwd = 3)
  legend('bottomright', c('Parametric fit', 'Posterior Mean'), lty = 1, 
         lwd = 3, col = 1:2, bty = 'n')
  
  preddt <- data.frame(pred_date = datestar, posterior_mean = posmu, parametric_fit = mustar.pred, lower = ll, upper = ul)
  write.csv(preddt, file = paste0('prediction/', State, '.csv'), row.names = F)
  save(theta.sort, file = paste0('prediction/', 'theta_', State, '.Rdata'))

  png(file =  paste0('plot/', State, '.png'))
  #plot predictions
  plot(datestar, posmu, lwd = 3, type = 'l', col = 2, main = paste('Cumulative COVID-19 cases Prediction for', State), xaxt = 'n', ylim = c(0, max(ul)), xlab = "Day", ylab = 'Cumulative cases')
  polygon(c(datestar,rev(datestar)),c(ll,rev(ul)),col=rgb (.1,.1,.1,.2),border=NA, xaxt = 'n')
  lines(datestar, posmu, lwd = 3, col = 2, xaxt = 'n')
  axis.Date(x = datestar, side = 1, format = '%m-%d', at = datestar)
  points(covid$date, covid$y * max(covid$cases), cex = 1.2)
  lines(datestar, mustar.pred, lwd = 3)
  legend('bottomright', c('Parametric fit', 'Posterior Mean'), lty = 1, 
         lwd = 3, col = 1:2, bty = 'n')
  dev.off()
}


oh.gp <- covidGP("Ohio")
