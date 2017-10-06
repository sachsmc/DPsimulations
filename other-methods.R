set.seed(seed)


run.it <- function(scenario) {

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


res <- replicate(200, run.it("A"))

df <- data.frame(t(res))
sapply(df, function(x) mean(x == "X1"))
