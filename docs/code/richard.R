#prior mean function (Richards growth equation)
mu <- function(t, a = 1, v = 1, k = 1, tm = 1) {
  a / (1 + v * exp(- k * (t  - tm)))^(1/v)
}





covid <- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv")
covid$date <- as.Date(covid$date)
covid <- covid %>% filter(state == "Iowa") 
mindate <- as.numeric(min(covid$date))
covid <- covid %>% mutate(t = (as.numeric(date) - mindate)/100, y = cases / max(cases))





Hyper <- function(x, y, init.val = c(1,1,1,1,1,1,1)) {
  
  marlik <- function(theta) {
    x <- as.matrix(x)
    n <- dim(x)[1]
    k <- ker(x = x, l = theta[1], sigf = theta[2])
    -dmvnorm(x = y, mean = mu(x, a = theta[4], v = theta[5], k = theta[6],tm = theta[7]),  sigma = k + theta[3]^2 * diag(n), log = T)
  }
  
  hyp <- optim(par= init.val, fn = marlik, method = 'Nelder-Mead',
               control=list(maxit = 10000, trace = 0))
  #print(hyp)
  return(hyp$par)
  
}



#--- Multistart
#parseq <- c(.01,.1,1)
#parmat <- expand.grid(parseq, parseq, parseq, parseq, parseq, parseq, parseq)

parseq <- seq(0.01, 20, by = .5)
parmat <- matrix(parseq, ncol = 7, nrow = length(parseq))

Hyper.ms <- function(x, y) {
  
  marlik <- function(theta) {
    x <- as.matrix(x)
    n <- dim(x)[1]
    #theta <- theta^2
    k <- ker(x = x, l = theta[1], sigf = theta[2])
    -dmvnorm(x = y, mean = mu(x, a = theta[4], v = theta[5], k = theta[6],tm = theta[7]),  sigma = k + theta[3]^2 * diag(n), log = T)
  }
  
  hyp <- multistart(parmat=parmat, fn = marlik, method = 'Nelder-Mead',
                    control=list(maxit = 100000, trace = 0))
  #print(hyp)
  return(hyp)
  
}






theta.ms <- Hyper.ms(x = covid$t, y = covid$y)
theta.sort <- theta.ms[order(theta.ms$value),]
theta <- theta.sort[1,1:7] #might change set of values if the best doesnt work well
theta <- as.numeric(theta)

#check if fitting is good
tstar <- seq(0,max(covid$t)*3, length.out = 100)
mu.pred <- mu(covid$t,a = theta[4], v = theta[5], k = theta[6], tm = theta[7])
mustar.pred <- mu(tstar, a = theta[4], v = theta[5], k = theta[6], tm = theta[7])

plot(tstar, mustar.pred, type = 'l', lwd = 2, col = 1, ylim = c(0, max(1,theta[4])), main = 'Mean function estimation')
points(covid$t, covid$y, col = 2)
legend('bottomright', c( 'Fitted',  'True'), lty = 1, lwd = 3, col = 1:2)





datestar <- covid$date
datestar <- c(datestar, seq.Date(from = max(datestar), by = "day", length.out = 90))
tstar <- (as.numeric(datestar) - mindate)/100

mu.pred <- mu(covid$t,a = theta[4], v = theta[5], k = theta[6], tm = theta[7])
mustar.pred <- mu(tstar,a = theta[4], v = theta[5], k = theta[6], tm = theta[7])


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
plot(datestar, posmu, lwd = 3, type = 'l', col = 2, main = 'Cumulative COVID-19 cases Prediction for Ohio', xaxt = 'n', ylim = c(0, max(ul)), xlab = "Day", ylab = 'Cumulative cases')
polygon(c(datestar,rev(datestar)),c(ll,rev(ul)),col=rgb (.1,.1,.1,.2),border=NA, xaxt = 'n')
lines(datestar, posmu, lwd = 3, col = 2, xaxt = 'n')
axis.Date(x = datestar, side = 1, format = '%m-%d', at = datestar)
points(covid$date, covid$y * max(covid$cases), cex = 1.2)
lines(datestar, mustar.pred, lwd = 3)
legend('bottomright', c('Parametric fit', 'Posterior Mean'), lty = 1, 
       lwd = 3, col = 1:2, bty = 'n')


