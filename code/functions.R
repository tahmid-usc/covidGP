# RBF kernel

ker <- function(x, l, sigf) {
  rbf <- rbfdot(sigma = 1/l^2)
  return(sigf^2 * kernelMatrix(rbf, x = x))
}

ker2 <- function(x, y, l, sigf) {
  rbf <- rbfdot(sigma = 1/l^2)
  return(sigf^2 * kernelMatrix(rbf, x = x, y = y))
}


#prior mean function (Richards growth equation)
mu <- function(t, a = 1, b0 = 1, b1 = 1, v = 1) {
  a / (1 + b0 * exp(- b1 * t))^v
}


#-----------------Hyper parameter optimization


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

parseq <- seq(0.01, 10, by = .1)
parmat <- matrix(parseq, ncol = 7, nrow = length(parseq))

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