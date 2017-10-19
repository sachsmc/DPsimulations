library(rjags)
library(splines)


get_prior <- function(w, j) {
  
  wbar <- apply(w, 2, mean) 
  wcov <- var(w) 
  list(a0 = 10, 
       b0 = 1, 
       nu1 = 3 + 1.5 * j, 
       nu2 = 3 + 1.5 * j, 
       s2 = 0.5 * wcov, 
       m2 = wbar, 
       psiinv2 = 2 * solve(wcov),
       tau1 = 6.01, tau2 = 3.01)
  
}



cat("
    model {
    ## Priors
    alpha ~ dnorm(0, .001)
    
    for(j in 1:k){
    beta[j] ~ dnorm(0, .001)
    }
    
    ## Likelihood
    for (i in 1:n){
    
    tauy[i] <- 1 / (sdy[i] * sdy[i])
    for(j in 1:k){
    taux[i,j] <- 1 / (sdx[i,j] * sdx[i,j])
    truex[i,j] ~ dnorm(0, 0.01)
    x[i,j] ~ dnorm(truex[i,j], taux[i,j])
    
    }
    truey[i] ~ dnorm(mu[i], 1)
    y[i] ~ dnorm(truey[i], tauy[i])
    mu[i] <- alpha + beta %*% truex[i,]
    }
    }
    ", fill=T, file="xyerror.txt")


sample_from_density <- function(n, x, y) {
  
  sample(x, size = n, prob = y, replace = TRUE)
  
}

all_sets_formula <- function(vars) {
  
  all.sets <- NULL
  for(i in 1:length(vars)) {
    
    all.sets <-  c(all.sets, combn(vars, i, FUN = function(x) c(x), simplify = FALSE))
    
  }
  all.sets
  
}



sample_from_model <- function(y, x1, y.se, x1.se, xpred, ...) {
  
  ses <- cbind(y.se, x1.se)
  w <- cbind(y, x1)
  
  
  prior <- get_prior(w, ncol(x1))
  
  mcmc <- list(nburn=2000,
               nsave=1000,
               nskip=2,
               ndisplay=1e6)
  
  mygrid <- seq(min(y) - 2 * sd(y), max(y) + 2 * sd(y), length.out = 500)
  
  sink(tempfile())
  fitLDDP <- DPcdensity(y = y, x = x1, mus=w, sds=ses, xpred = xpred,
                        grid=mygrid, 
                        compute.band=FALSE,
                        type.band="HPD",
                        prior=prior, 
                        mcmc=mcmc, 
                        state=NULL, 
                        status=TRUE, 
                        work.dir = tempdir())
  sink(NULL)
  
  y.hat.samps <- lapply(1:length(y), function(i){
    
    data.frame(trial = i, t.hat.post = sample_from_density(100, fitLDDP$grid, fitLDDP$densp.m[i, ]))
    
  })
  
  do.call(rbind, y.hat.samps)
  
}



sample_from_marginal <- function(y, x1, y.se, x1.se, ...) {
  
   # bundle data
  jags_d <- list(x = x1, y = y, 
                 sdx = x1.se, sdy = y.se, n = length(y), k = ncol(x1))
  
  # initiate model
  mod2 <- jags.model("xyerror.txt", data=jags_d,
                     n.chains=1, n.adapt=1000, quiet = TRUE)
  
  # simulate posterior
  out <- coda.samples(mod2, n.iter=10000, thin=10,
                      variable.names=c("truey"), quiet = TRUE)[[1]]
  
  colnames(out) <- paste("trial", 1:ncol(out), sep = ".")
  out <- as.data.frame(out)
  
  y.hat.samps <- reshape(out, varying = colnames(out), direction = "long")[, -3]
  colnames(y.hat.samps) <- c("trial", "y.marginal")
  
  y.hat.samps
  
}



predict_from_model <- function(y, x1, y.se, x1.se, xpred, xpred.se, band = FALSE, ...) {
  
  # bundle data
  jags_d <- list(x = rbind(x1, xpred), y = c(y, NA), 
                 sdx = rbind(x1.se, xpred.se), sdy = c(y.se, mean(y.se)), n = length(y) + 1, k = ncol(x1))
  
  # initiate model
  mod2 <- jags.model("xyerror.txt", data=jags_d,
                     n.chains=1, n.adapt=1000, quiet = TRUE)
  
  # simulate posterior
  out <- coda.samples(mod2, n.iter=10000, thin=10,
                      variable.names=c("alpha", "beta", "truey"), quiet = TRUE)
  
  yhat.dist <- out[[1]][, ncol(out[[1]])]
  
  if(band) {
    data.frame(xpred = xpred, yhat = median(yhat.dist), yhat.lower = quantile(yhat.dist, .025), yhat.upper = quantile(yhat.dist, 0.975))
  } else {
    data.frame(xpred = xpred, yhat = median(yhat.dist))
  }
  
}


leave_one_out <- function(trials, xvars) {
  
  predictions <- as.data.frame(matrix(NA, ncol = 4 + length(xvars), nrow = length(unique(trials$trial))))
  for(j in unique(trials$trial)) {
    
    y <- subset(trials, var == "Y" & trial != j)$ests
    y.se <- subset(trials, var == "Y" & trial != j)$ses
    
    x <- subset(trials, var %in% xvars & trial != j)
    x.est <- as.matrix(reshape(x, direction = "wide", drop = "ses", v.names = "ests", timevar = "var", idvar = "trial")[, -1])
    x.ses <- as.matrix(reshape(x, direction = "wide", drop = "ests", v.names = "ses", timevar = "var", idvar = "trial")[, -1])
    
    xpred.x <- subset(trials, var %in% xvars & trial == j)
    xpred <- as.matrix(reshape(xpred.x, direction = "wide", drop = "ses", v.names = "ests", timevar = "var", idvar = "trial")[, -1])
    xpred.se <- as.matrix(reshape(xpred.x, direction = "wide", drop = "ests", v.names = "ses", timevar = "var", idvar = "trial")[, -1])
    
    predfit <- predict_from_model(y, x.est, y.se, x.ses, xpred, xpred.se, band = TRUE)
    predictions[j, ] <- cbind(trial.leftout = j, predfit)
    
  }
  colnames(predictions) <- c("trial.leftout", colnames(predfit))
  predictions
  
}


pred_null <- function(trials) {
  
  prediction <- NULL
  for(j in unique(trials$trial)) {
    
    y <- subset(trials, var == "Y" & trial != j)$ests
    y.se <- subset(trials, var == "Y" & trial != j)$ses
    y.bar <- mean(rnorm(length(y) * 1000, mean = rep(y, 1000), sd = rep(y.se, 1000)))
    
    prediction <- rbind(prediction, data.frame(trial = j, yhat.0 = y.bar))
    
  }
  prediction
  
}



gen_effects <- function(scenario, n = 21) {
  
  J <- 6
  X <- matrix(rnorm(n * J, sd = 1), nrow = n, ncol = J)
  err <- rep(0, n) #025 * rt(nrow(X), df = 15)
  
  if(scenario == "0") {
    
    Y <- 1.5 + rnorm(n, sd = 1.2)
    
  } else if(scenario == "A") {
    
    Y <- 1.5 + 1 * X[, 1] + err
    
  } else if(scenario == "B") {
    
    X.star <- cbind(bs(X[, 1], knots = c(-.5), degree = 1), X[, -1])
    beta <- c(-1.75, .5, 0, 0, 0, 0, 0)
    Y <-  1 - X.star %*% beta + err
    
  } else if(scenario == "C") {
    
    Y <- 3.5 - (.45 * X[, 1] - 1) * (.35 * X[, 1] - .5) + err
    
  } else if(scenario == "D") {  ## multivariate
    
    beta <- c(.25, .25, .5, 0, 0, 0)
    Y <- 1.5 + X %*% beta + err
    
  } else if(scenario == "E") { ## multivariate
  
    X.star <- cbind(bs(X[, 1], knots = c(-.25, .25)), 
                    bs(X[, 2], knots = c(-.5, .25)),
                    bs(X[, 3], knots = c(-.25, .5)),
                    X[, 4:6])
    beta <- c(rep(c(1.5, 1.25, .5, -.2, 1.25), 3), 0, 0, 0)
    
    Y <- -1 + X.star %*% beta + err
      
  }else if(scenario == "F") {
    
    Y <- 1.5 + .5 * X[, 1] + .5 * X[, 2] + X[, 1] * X[, 2] + err
    
  } else stop("No scenario")
  
  list(Y = Y, X = X)
  
}



samp.indi.data <- function(i, effects, enns) {
  
 
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



run_one_rep <- function(scenario, n = 21, output, seed) {
  
  
  set.seed(seed)
  
  enns <- sample(ceiling(runif(n, 100, 500)))
  
  eff <- gen_effects(scenario, n)
  trials <- do.call(rbind, lapply(1:n, function(i) analyze.indi.data(samp.indi.data(i, eff, enns), i)))
  ## add pairwise interactions
  
  xes.tmp <- subset(trials, grepl("X", var))
  xes.tmp2 <- subset(merge(xes.tmp, xes.tmp, by = "trial"), var.x != var.y)
  xes.tmp2$ests <- with(xes.tmp2, ests.x * ests.y)
  xes.tmp2$ses <- with(xes.tmp2, sqrt(ests.x^2 * ses.y^2 + 
                                        ests.y^2 * ses.x^2 + 
                                        ses.x^2 * ses.y^2))
  xes.tmp2$var <- with(xes.tmp2, paste(var.x, var.y, sep = ":"))
  
  xes.tmp3 <- xes.tmp2[, c("trial", "ests", "ses", "var")]
  trials <- rbind(trials, xes.tmp3)  
  
  leftout <- subset(trials, trial == n)
  #trials <- subset(trials, trial != n)
  
  y <- subset(trials, var == "Y")$ests
  y.se <- subset(trials, var == "Y")$ses
  
  x <- subset(trials, grepl("^X", var))
  
  xvars <- as.character(sort(unique(x$var)))
  xvars <- xvars[!grepl(":", xvars)]
  
  x.est <- as.matrix(reshape(x, direction = "wide", drop = "ses", v.names = "ests", timevar = "var", idvar = "trial")[, -1])
  x.ses <- as.matrix(reshape(x, direction = "wide", drop = "ests", v.names = "ses", timevar = "var", idvar = "trial")[, -1])
  
  
  Y.full.sample <- sample_from_marginal(y, x.est, y.se, x.ses)
  Y.full.sample$seed <- seed
  Y.full.sample$scenario <- scenario
  
  ## compute D0
  
  
  pred.0 <- pred_null(trials)
  varsets <- all_sets_formula(xvars)
  
  ## all pairwise interactions
  
  wh2 <- lapply(varsets, function(x){ 
    if(length(x) == 2) {
      
      return(c(x, paste(x, collapse = ":")))
    } else NULL
  
  })
  wh2 <- wh2[sapply(wh2, function(x) !is.null(x))]
  varsets <- c(varsets, wh2)
  
  D.ests <- rep(NA, length(varsets))
  D.true <- rep(NA, length(varsets))
  P.D0 <- rep(NA, length(varsets))
  
  var.groups <- sapply(varsets, paste, collapse = "|")
  compare.sets <- vector(mode = "list", length = length(var.groups))
  
  for(i in 1:length(varsets)) {
    
    pred.ins <- leave_one_out(trials, varsets[[i]])
    compareto <- merge(pred.ins, pred.0, by.x = "trial.leftout", by.y = "trial")
    compareto <- merge(compareto, data.frame(trial.leftout = 1:n, y = eff$Y[1:n]), by = "trial.leftout")
    compareto$absdiff.true <- with(compareto, abs(y - yhat))
    
    compareto$varset <- var.groups[i]
    compareto$scenario <- scenario
    compareto$replicate.seed <- seed
    compare.sets[[i]] <- compareto[, -grep("^xpred", colnames(compareto))]
  
  }
  
  outdat <- do.call(rbind, compare.sets)
  
  write.table(Y.full.sample, file = gsub(".txt", "marginal.txt", output), sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(outdat, file = output, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  
  
}

