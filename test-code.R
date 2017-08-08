library(DPpackagemod)

n <- 20
freq <- TRUE

J <- 6
X <- matrix(rnorm(n * J, sd = 1), nrow = n, ncol = J)
Y <-  2 + .25 * X[, 1] + .25 * rt(nrow(X), df = 15)

if(J > 3) {
  Y <-  2 + .125 * X[, 1] + .125 * X[, 2] + .125 * X[, 3] + .25 * rt(nrow(X), df = 15)  
  
}

enns <- sample(ceiling(runif(n, 100, 500)))

samp.indi.data <- function(i) {
  
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
  sests <- vapply(1:J, FUN = ana.1X, FUN.VALUE = rep(0, 2))
  
  data.frame(trial = i, ests = c(tbeta, sests[1, ]), ses = c(tse, sests[2, ]), var = c("Y", paste0("X", 1:J)))
  
}

indi.results <- do.call(rbind, lapply(1:n, function(i) analyze.indi.data(samp.indi.data(i), i)))

## estimate full model

y <- subset(indi.results, var == "Y")$ests
y.se <- subset(indi.results, var == "Y")$ses

x <- subset(indi.results, grepl("X", var))
x.est <- as.matrix(reshape(x, direction = "wide", drop = "ses", v.names = "ests", timevar = "var", idvar = "trial")[, -1])
x.ses <- as.matrix(reshape(x, direction = "wide", drop = "ests", v.names = "ses", timevar = "var", idvar = "trial")[, -1])

xpred <- as.matrix(x.est)

test <- sample_from_model(y, x.est, y.se, x.ses, xpred)

plot(by(test$t.hat.post, test$trial, mean) ~ y)



