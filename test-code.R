library(DPpackagemod)
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



predict_from_model <- function(y, x1, y.se, x1.se, xpred, band = FALSE, ...) {
  
  ses <- cbind(y.se, x1.se)
  w <- cbind(y, x1)
  
  
  prior <- get_prior(w, ncol(x1))
  
  mcmc <- list(nburn=2000,
               nsave=1000,
               nskip=2,
               ndisplay=1e6)
  
  sink(tempfile())
  fitLDDP <- DPcdensity(y = y, x = x1, mus=w, sds=ses, xpred = xpred,
                        ngrid=1, 
                        compute.band=band,
                        type.band="HPD",
                        prior=prior, 
                        mcmc=mcmc, 
                        state=NULL, 
                        status=TRUE, 
                        work.dir = tempdir())
  sink(NULL)
  
  if(band) {
    data.frame(xpred = xpred, yhat = fitLDDP$meanfp.m, yhat.lower = fitLDDP$meanfp.l, yhat.upper = fitLDDP$meanfp.h)
  } else {
    data.frame(xpred = xpred, yhat = fitLDDP$meanfp.m)
  }
  
}


leave_one_out <- function(trials, xvars) {
  
  predictions <- as.data.frame(matrix(NA, ncol = 2 + length(xvars), nrow = length(unique(trials$trial))))
  for(j in unique(trials$trial)) {
    
    y <- subset(trials, var == "Y" & trial != j)$ests
    y.se <- subset(trials, var == "Y" & trial != j)$ses
    
    x <- subset(trials, var %in% xvars & trial != j)
    x.est <- as.matrix(reshape(x, direction = "wide", drop = "ses", v.names = "ests", timevar = "var", idvar = "trial")[, -1])
    x.ses <- as.matrix(reshape(x, direction = "wide", drop = "ests", v.names = "ses", timevar = "var", idvar = "trial")[, -1])
    
    xpred.x <- subset(trials, var %in% xvars & trial == j)
    xpred <- as.matrix(reshape(xpred.x, direction = "wide", drop = "ses", v.names = "ests", timevar = "var", idvar = "trial")[, -1])
    
    predfit <- predict_from_model(y, x.est, y.se, x.ses, xpred)
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
  
  leftout <- subset(trials, trial == n)
  trials <- subset(trials, trial != n)
  
  y <- subset(trials, var == "Y")$ests
  y.se <- subset(trials, var == "Y")$ses
  
  x <- subset(trials, grepl("^X", var))
  
  xvars <- as.character(sort(unique(x$var)))
  
  x.est <- as.matrix(reshape(x, direction = "wide", drop = "ses", v.names = "ests", timevar = "var", idvar = "trial")[, -1])
  x.ses <- as.matrix(reshape(x, direction = "wide", drop = "ests", v.names = "ses", timevar = "var", idvar = "trial")[, -1])
  
  Y.full.sample <- sample_from_model(y, x.est, y.se, x.ses, x.est)
  
  ## compute D0
  
  
  pred.0 <- pred_null(trials)
  varsets <- all_sets_formula(xvars)
  
  D.ests <- rep(NA, length(varsets))
  D.true <- rep(NA, length(varsets))
  P.D0 <- rep(NA, length(varsets))
  
  var.groups <- sapply(varsets, paste, collapse = "|")
  compare.sets <- vector(mode = "list", length = length(var.groups))
  
  for(i in 1:length(varsets)) {
    
    pred.ins <- leave_one_out(trials, varsets[[i]])
    compareto <- merge(pred.ins, Y.full.sample, by.x = "trial.leftout", by.y = "trial")
    compareto <- merge(compareto, pred.0, by.x = "trial.leftout", by.y = "trial")
    compareto$absdiff <- with(compareto, abs(t.hat.post - yhat))
    compareto$absdiff.0 <- with(compareto, abs(t.hat.post - yhat.0))
    
    P.D0[i] <- mean(with(compareto, absdiff.0 < absdiff))
    D.ests[i] <- mean(compareto$absdiff)
    D.true[i] <- mean(abs(y - pred.ins$yhat))
    compare.sets[[i]] <- compareto[, c("trial.leftout", "yhat", "t.hat.post")]
  
  }
  
  
  btr.than.0 <- P.D0 < 0.4
  if(!any(btr.than.0)) {
    
    pred0.ci <- t.test(y)$conf.int
    
    outdat <- data.frame(Scenario = scenario, 
                      D.est = mean(compareto$absdiff.0), 
                      D.true = mean(abs(y - pred.0$yhat.0)), 
                      model.sel = "Null", 
                      y.j1 = subset(leftout, var == "Y")$ests, 
                      yhat.j1 = mean(y), 
                      yhat.j1.lower = pred0.ci[1], yhat.j1.upper = pred0.ci[2])
    
  } else {
    
    
    kp <- D.ests[btr.than.0]
    kp.true <- D.true[btr.than.0]
    compare.sets <- compare.sets[btr.than.0]
    
    cand.vars <- varsets[btr.than.0]
    
    mydex <- which.min(kp)
    
    ## compare all sets to mydex
    
    P.Dkp <- rep(NA, length(kp))
    for(j in (1:length(kp))) {
      
      if(j == mydex){
        P.Dkp[j] <- 1
        next
      }
    
      tmp.compare <- merge(compare.sets[[mydex]], compare.sets[[j]], by = c("trial.leftout", "t.hat.post"))
      tmp.compare <- within(tmp.compare, {
        absdiff.x <- abs(t.hat.post - yhat.x)
        absdiff.y <- abs(t.hat.post - yhat.y)
      })
      P.Dkp[j] <- mean(with(tmp.compare, absdiff.y < absdiff.x))
      
    }
    
    ## keep the smallest subset not significantly different
    
    mydex <- which.min(sapply(cand.vars[P.Dkp > .5], length))
    
    xkeep <- as.numeric(gsub("X", "", varsets[[mydex]]))
    x.est <- as.matrix(reshape(x, direction = "wide", drop = "ses", v.names = "ests", timevar = "var", idvar = "trial")[, -1])[, xkeep, drop = FALSE]
    x.ses <- as.matrix(reshape(x, direction = "wide", drop = "ests", v.names = "ses", timevar = "var", idvar = "trial")[, -1])[, xkeep, drop = FALSE]
    
    xpred.x <- subset(leftout, var %in% varsets[[mydex]])
    xpred <- as.matrix(reshape(xpred.x, direction = "wide", drop = "ses", v.names = "ests", timevar = "var", idvar = "trial")[, -1])
    
    predfit <- predict_from_model(y, x.est, y.se, x.ses, xpred, band = TRUE)
    
    outdat <- data.frame(Scenario = scenario, 
                      D.est = kp[mydex], 
                      D.true = kp.true[mydex], 
                      model.sel = var.groups[btr.than.0][mydex], 
                      y.j1 = subset(leftout, var == "Y")$ests, yhat.j1 = predfit$yhat, 
                      yhat.j1.lower = predfit$yhat.lower, yhat.j1.upper = predfit$yhat.upper) 
    
  }
  
  
  write.table(outdat, file = output, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  
  
}

