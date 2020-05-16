library(brms)
library(dplyr)
library(lubridate)
library(readr)

# Read covid data

covid_data <- function (State  = 'South Carolina') {
  covid <- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv")
  if(State == "US"){
    # Total USA: transform date and cases
    covid$date <- as.Date(covid$date)
    covid.total <- covid %>% group_by(date) %>% summarize(cases = sum(cases)) 
    covid <- covid.total}
  else{
    # Read data and transform date and cases
    # change state name as desired
    print(State)
    covid$date <- as.Date(covid$date)
    covid <- covid %>% dplyr:::filter(state == State) 
    
  }
  
  return(covid)
}




covid <- covid_data("US") 
mindate <- as.numeric(min(covid$date))
covid <- covid %>% mutate(t = (as.numeric(date) - mindate)/100, y = cases / max(cases))


#prior mean function (Richards growth equation)
mu <- function(t, a = 1, b0 = 1, b1 = 1, v = 1) {
  a / (1 + v * exp(- b1 * (t - b0)))^(1/v)
}

prior1 <- prior(exponential(1), nlpar = "a") +
  prior(exponential(1), nlpar = "b1") +
  prior(exponential(1), nlpar = "b0") +
  prior(beta(1,1), nlpar = "v")


dt <- cbind(y = covid$y, t = covid$t)
dt <- as.data.frame(dt)


fit1 <- brm(bf(y ~ a / (1 + v * exp(- b1 * (t - b0)))^(1/v), a + b0 + b1 + v ~ 1, nl = TRUE),
            data = dt, prior = prior1, iter = 50000,  control = list(adapt_delta = 0.9))

summary(fit1)
plot(fit1)

plot(conditional_effects(fit1), points = TRUE)



datestar <- covid$date
datestar <- c(datestar, seq.Date(from = max(datestar), by = "day", length.out = 90))
tstar <- (as.numeric(datestar) - mindate)/100

ndt <- data.frame(t = tstar)
pred <- predict(fit1, newdata = ndt)
pred <- as.data.frame(pred) %>% cbind(ndt)



est <- pred$Estimate * max(covid$cases)
ll <- pred$Q2.5  * max(covid$cases)
ul <- pred$Q97.5  * max(covid$cases)


#plot predictions
plot(datestar, est, lwd = 3, type = 'l', col = 2, main = 'Cumulative COVID-19 cases Prediction for Ohio', xaxt = 'n', ylim = c(0, max(ul)), xlab = "Day", ylab = 'Cumulative cases')
polygon(c(datestar,rev(datestar)),c(ll,rev(ul)),col=rgb (.1,.1,.1,.2),border=NA, xaxt = 'n')
lines(datestar, est, lwd = 3, col = 2, xaxt = 'n')
axis.Date(x = datestar, side = 1, format = '%m-%d', at = datestar)
points(covid$date, covid$y * max(covid$cases), cex = 1.2)





