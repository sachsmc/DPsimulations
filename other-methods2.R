library(glmnet)

run.it <- function(scenario, n = 13) {
  
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
  
  
  dat <- data.frame(y, x.est)
  
  fit.init <- lm(y ~ ., data = dat)
  #fit.null <- lm(y ~ 1, data = dat)
  
  #res.null <- step(fit.null, scope = y ~ ests.X1 + ests.X2 + ests.X3 + ests.X4 + ests.X5 + ests.X6, trace = 0)
  res.step <- step(fit.init, scope = y ~ ., trace = 0)
  
  var.seld.step <- sort(gsub("ests.", "", grep("ests", names(res.step$coefficients), value = TRUE)))
  #var.seld.null <- sort(gsub("ests.", "", grep("ests", names(res.null$coefficients), value = TRUE)))
  
  las.fit <- cv.glmnet(x.est, y)
  beta.est <- las.fit$glmnet.fit$beta[, which(las.fit$lambda == las.fit$lambda.min)]
  var.seld.lass <- paste0("X", 1:6)[which(beta.est != 0)]
  
  if(length(var.seld.step) == 0) var.seld.step <- "Null"
  if(length(var.seld.lass) == 0) var.seld.lass <- "Null"
  
  c(paste(var.seld.step, collapse = "|"), paste(var.seld.lass, collapse = "|"))
  
}


res0 <- t(replicate(200, run.it("0")))
resA <- t(replicate(200, run.it("A")))
resB <- t(replicate(200, run.it("B")))
resC <- t(replicate(200, run.it("C")))
resD <- t(replicate(200, run.it("D")))
resE <- t(replicate(200, run.it("E")))
resF <- t(replicate(200, run.it("F")))

freq.res <- cbind(c("0", LETTERS[1:6]), as.data.frame(rbind(
  c(sapply(data.frame(res0), function(x) c(mean(x == "Null"), 0))),
  c(sapply(data.frame(resA), function(x) c(mean(x == "X1"), mean(grepl("X1", x))))),
  c(sapply(data.frame(resB), function(x) c(mean(x == "X1"), mean(grepl("X1", x))))),
  c(sapply(data.frame(resC), function(x) c(mean(x == "X1"), mean(grepl("X1", x))))),
  c(sapply(data.frame(resD), function(x) c(mean(x == "X1|X2|X3"), mean(grepl("X1|X2|X3", x))))),
  c(sapply(data.frame(resE), function(x) c(mean(x == "X1|X2|X3"), mean(grepl("X1|X2|X3", x))))),
  c(sapply(data.frame(resF), function(x) c(mean(x == "X1|X2"), mean(grepl("X1|X2", x))))))))


# |scenario | cor.mean| cor.mean.plus| cor.median| cor.median.plus| cor.true| cor.true.plus|
  # |:--------|--------:|-------------:|----------:|---------------:|--------:|-------------:|
  # |0        |     0.21|          0.21|       0.09|            0.09|     0.22|          0.22|
  # |A        |     0.29|          1.00|       0.27|            1.00|     0.42|          1.00|
  # |B        |     0.29|          0.75|       0.23|            0.86|     0.36|          0.79|
  # |C        |     0.38|          1.00|       0.32|            0.99|     0.44|          1.00|
  # |D        |     0.47|          0.94|       0.41|            0.86|     0.62|          0.97|
  # |E        |     0.01|          0.03|       0.01|            0.05|     0.02|          0.05|
  # |F        |     0.19|          0.43|       0.14|            0.40|     0.21|          0.43|


load("sim-table2-dpp.RData")
load("sim-table2-linear.RData")


fin.res.table1 <- cbind(dpp.simtab2, linear.simtab2[, 2:3], freq.res[, -1])
colnames(fin.res.table1) <- c("Scenario", rep(c("Exact", "Partial"), 5))

knitr::kable(fin.res.table1, digits = 2, format = "latex")

