library(dplyr)
library(lubridate)
library(readr)
library(optimx)
library(ggplot2)
library(kernlab)
library(mvtnorm)
library(forecast)
source('code/functions.R')


forecast_error <- function(state = 'US') {
  
  covid.full <- covid_data(state) 
  covid <- covid.full %>% filter(date < max(date) - 6)
  mindate <- as.numeric(min(covid$date))
  maxNum <- max(covid$cases)
  covid <- covid %>% mutate(t = (as.numeric(date) - mindate)/100, y = cases / maxNum)
  
  #Estimate hyperparamters (multistarts) (works better)
  
  theta.ms <- Hyper.ms(x = covid$t, y = covid$y)
  theta.sort <- theta.ms[order(theta.ms$value),]
  theta <- theta.sort[1,1:7] #might change set of values if the best doesnt work well
  theta <- as.numeric(theta)
  
  
  #----------------
  #Compute test error
  
  datestar <- covid.full$date
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
  posmu <- posmu * maxNum
  
  posmu <- ts(posmu, start = mindate)
  
  train.err <- accuracy(posmu, x = covid.full$cases, test = 1:n)
  test.err <- accuracy(posmu, x = covid.full$cases, test = (n+1):dim(covid.full)[1])
  
  return(cbind(train.err, test.err))
  
}

US.err <- forecast_error("US") 

covid.us <- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv")
state <-  unique(covid.us$state)
state <- sort(state[-c(51,52,53,55)])
state <- c(state, "US")


state.err <- c()
for (i in 1:length(county)) {
  county.err <- rbind(state.err, state[i], forecast_error(state[i]))
}

err.name <- c('ME',     'RMSE',      'MAE',       'MPE',     'MAPE')

#write.csv(county.err, "county_err.csv")

err.mat <- matrix(county.err, nrow = 11)
err.mat <- t(err.mat)
head(err.mat)
err.mat <- as.data.frame(err.mat)
names(err.mat) <- c('county', paste0('train_',err.name), paste0('test_',err.name))
head(err.mat)
err.df <- apply(err.mat[,-1], 2, function(x){ as.numeric(x) })
err.df <- as.data.frame(err.df, stringsAsFactors = F)
err.df <- cbind(county, err.df)
names(err.df) <- c('county', paste0('train_',err.name), paste0('test_',err.name))



ggplot(data = err.df,aes(x = county, y = train_MAPE)) +
  geom_col() +
  xlab('MAPE') +
  ylab('County') + 
  scale_y_continuous(n.breaks = 20) +
  coord_flip() 
  
#  ggsave('plot/cases/error/train_MAPE.png')

par(mar = c(12, 5, 2, 2)) # Set the margin on all sides to 2
barplot(cbind(train_MAE, test_MAE) ~ county, err.df,horiz = F, beside = T, legend.text = c('Train', 'Test'), axisnames = T, las = 2, 
        axes = T, col = c('cadetblue1','dodgerblue4'), border = F, ylab = 'Mean Absolute Error (MAE)', xlab = '', args.legend = list(bty = 'n'))

barplot(cbind(train_RMSE, test_RMSE) ~ county, err.df,horiz = F, beside = T, legend.text = c('Train', 'Test'), axisnames = T, las = 2, 
        axes = T, col = c('cadetblue1','dodgerblue4'), border = F, ylab = 'RMSE', xlab = '', args.legend = list(bty = 'n'))

barplot(cbind(train_MAPE, test_MAPE) ~ county, err.df,horiz = F, beside = T, legend.text = c('Train', 'Test'), axisnames = T, las = 2, 
        axes = T, col = c('cadetblue1','dodgerblue4'), border = F, ylab = 'MAPE', xlab = '', args.legend = list(bty = 'n'))
