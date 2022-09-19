gev_llik <- function(par,Y){
  
  m <- length(Y)
  sigma <- exp(par[2])
  mu <- par[1]
  xi <- par[3]
  
  for (val in Y) {
    if(1 + xi * (val - mu)/ sigma <= 0){
      return(-Inf)
    }
  }
  
  A <- -m * log(sigma)
  
  if(xi > -0.05 && xi < 0.05) {
    B <- -sum(Y - mu) / sigma
    C <- -sum(exp(-((Y - mu) / sigma)))
  } else {
    B <- -(1 + 1/xi) * sum(log(pmax(1 + xi * ((Y - mu)/sigma))))
    C <- -sum(pmax(1 + xi * ((Y - mu)/sigma))^(-1/xi))
  }

  
  
  
  return(A + B + C)
  
}

neg_gev_llik <- function(par,Y){
  num <- gev_llik(par,Y)
  
  return(-num)
}


# Random test
a <- c(3,1,0)
b <- c(2,2,3,4,3,4,3,3,3,3,1)

gev_llik(a,b)


pareto_llik <- function(par,Y){
  
  k <- length(Y)
  sigma <- exp(par[1])
  xi <- par[2]
  
  for (val in Y) {
    if(1 + xi * val / sigma <= 0){
      return(-Inf)
    }
  }
  
  A <- - k * log(sigma)
  
  if(xi > -0.05 && xi < 0.05) {
    B <- -sum(Y) / sigma
    } else {
    B <- -(1 + 1/xi) * sum(log(pmax(1 + xi * (Y / sigma))))
    }
  
  
  return(A + B)
}

Y <- c(1,2,3,2,2,2,3,2,1)
paretopar <- c(2,0)
pareto_llik(paretopar,Y)

return_level <- function(par, p) {
  
  sigma <- exp(par[2])
  mu <- par[1]
  xi <- par[3]
  y_p <- -log(1-p)
  
  if(xi > -0.05 && xi < 0.05) {
    z_p <- mu - sigma * log(y_p)
  } else{
    z_p <- mu - (sigma/xi) * (1-y_p^(-xi))
  }
  
  
  
  return(z_p)
}


  
