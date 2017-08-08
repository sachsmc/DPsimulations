library(DPpackagemod)
library(splines)

gen_effects <- function(scenario, n = 20) {
  
  J <- 6
  X <- matrix(rnorm(n * J, sd = 1), nrow = n, ncol = J)
  err <- rep(0, n) #025 * rt(nrow(X), df = 15)
  
  if(scenario == "0") {
    
    Y <- 1.5 + rnorm(n, sd = 1.2)
    
  } else if(scenario == "A") {
    
    Y <- 1.5 + 1 * X[, 1] + err
    
  } else if(scenario == "B") {
    
    X.star <- cbind(bs(X[, 1], knots = c(0), degree = 1), X[, -1])
    beta <- c(-1.75, .1, 0, 0, 0, 0, 0)
    Y <-  1 - X.star %*% beta + err
    
  } else if(scenario == "C") {
    
    Y <- 3.5 - (X[, 1] - 1) * (X[, 1] - .5) + err
    
  } else if(scenario == "D") {  ## multivariate
    
    beta <- c(.25, .25, .5, 0, 0, 0)
    Y <- 1.5 + X %*% beta + err
    
  } else if(scenario == "E") { ## multivariate
  
    X.star <- cbind(bs(X[, 1], knots = c(-.25, .25)), 
                    bs(X[, 2], knots = c(-.5, .25)),
                    bs(X[, 3], knots = c(-.25, .5)),
                    X[, 4:6])
    beta <- c(rep(c(.5, .25, .5, -.2, .75), 3), 0, 0, 0)
    
    Y <- .75 + X.star %*% beta + err
      
  } else stop("No scenario")
  
  list(Y = Y, X = X)
  
}



samp.indi.data <- function(i, effects) {
  
  enns <- sample(ceiling(runif(n, 100, 500)))
  
  Y <- effects$Y
  X <- effects$X
  J <- ncol(X)
  Z <- rbinom(enns[i], 1, .5)
  mu.S <- vapply(X[i, ], function(x) 2.5 + x * Z, FUN.VALUE = rep(0, length(Z)))
  lambda.T <- -1 + Y[i] * Z
  
  S <- mu.S + matrix(rnorm(enns[i] * J), nrow = enns[i], ncol = J)
  T <- rpois(enns[i], exp(lambda.T))
  data.frame(Z, T, S)

  
}

analyze.indi.data <- function(test, i) {
  
  
  t1 <- glm(T ~ Z, data = test, family = "poisson")
  tbeta <- t1$coefficients[2]
  tse <- sqrt(diag(vcov(t1)))[2]
  
  ana.1X <- function(x) {
    
    form <- as.formula(paste0("X", x, " ~ Z"))
    s1 <- glm(form, data = test)
    c(s1$coefficients[2], sqrt(diag(vcov(s1)))[2])
    
  }
  
  J <- ncol(test) - 2
  sests <- vapply(1:J, FUN = ana.1X, FUN.VALUE = rep(0, 2))
  
  data.frame(trial = i, ests = c(tbeta, sests[1, ]), ses = c(tse, sests[2, ]), var = c("Y", paste0("X", 1:J)))
  
}




