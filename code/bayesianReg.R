library(brms)
library(dplyr)
library(lubridate)
library(readr)
library(optimx)
library(kernlab)
library(mvtnorm)
source('code/functions.R')


covid <- covid_data("South Carolina") 
mindate <- as.numeric(min(covid$date))
covid <- covid %>% mutate(t = (as.numeric(date) - mindate)/100, y = cases / max(cases))


#prior mean function (Richards growth equation)
mu <- function(t, a = 1, b0 = 1, b1 = 1, v = 1) {
  a / (1 + b0 * exp(- b1 * t))^v
}

prior1 <- prior(exponential(10000), nlpar = "a") +
  prior(exponential(1), nlpar = "b1") +
  prior(normal(1, 10), nlpar = "b0") +
  prior(exponential(1), nlpar = "v")


dt <- cbind(y = covid$y, t = covid$t)
dt <- as.data.frame(dt)


fit1 <- brm(bf(y ~ a / (1 + b0 * exp(- b1 * t))^(1/v), a + b0 + b1 + v ~ 1, nl = TRUE),
            data = dt, prior = prior1)