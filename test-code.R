library(DPpackagemod)

n <- 20
freq <- TRUE

X <- rnorm(n, sd = 3)
Y <-  2 + .25 * X + .25 * rt(length(X), df = 15)
enns <- sample(ceiling(runif(n, 100, 500)))

samp.indi.data <- function(i) {
  
  Z <- rbinom(enns[i], 1, .5)
  mu.S <- 2.5 + X[i] * Z
  lambda.T <- -1 + Y[i] * Z
  
  S <- rnorm(enns[i], mean = mu.S)
  T <- rpois(enns[i], exp(lambda.T))
  data.frame(Z, S, T)
  
  
}

analyze.indi.data <- function(test) {
  
  s1 <- bayesglm(S ~ Z, data = test)
  t1 <- bayesglm(T ~ Z, data = test, family = "poisson")
  spost <- sim(s1, n.sims = 1000)@coef[, 2]
  tpost <- sim(t1, n.sims = 1000)@coef[, 2]
  
  c(mean(tpost),mean(spost), sd(tpost),  sd(spost))
  
}

analyze.indi.data.freq <- function(test) {
  
  s1 <- glm(S ~ Z, data = test)
  t1 <- glm(T ~ Z, data = test, family = "poisson")
  spost <- s1$coefficients[2]
  tpost <- t1$coefficients[2]
  
  c(tpost,spost, sqrt(diag(vcov(t1)))[2],  sqrt(diag(vcov(s1)))[2])
  
}

if(freq) myanalyze <- analyze.indi.data.freq else myanalyze <- analyze.indi.data

indi.results <- do.call(rbind, lapply(lapply(1:n, samp.indi.data), myanalyze))

y <- indi.results[, 1]
x1 <- cbind(indi.results[, 2], indi.results[, 2] + rnorm(nrow(indi.results)))
y.se <- indi.results[, 3]
x1.se <- indi.results[, c(4, 4)]


chidat <- data.frame(x1 = x1, y = y)

xpred <- x1


x <- seq(-3, 3, by = .1)
y <- dnorm(x)

test <- sample_from_density(500, x, y)
hist(test, freq = FALSE)
curve(dnorm(x), add = TRUE)


