source("test-code-save-all-Ds.R")

run_one_rep("E", n = 25, seed = 5543, output = "G:/STAFF/Michael Sachs/Stats Projects/Bayesian NP Meta Analysis/Simulation Code/R0fileN25.csv")
run_one_rep("E", n = 50, seed = 5543, output = "G:/STAFF/Michael Sachs/Stats Projects/Bayesian NP Meta Analysis/Simulation Code/R0fileN50.csv")
run_one_rep("E", n = 100, seed = 5543, output = "G:/STAFF/Michael Sachs/Stats Projects/Bayesian NP Meta Analysis/Simulation Code/R0fileN100.csv")


library(data.table)
rep1 <- fread("R0fileN25.csv")

colnames(rep1) <- c("trial.leftout", "yhat", "t.hat.post", "yhat.0", "y.true",  "absdiff", "absdiff.0", "absdiff.true", "varset", "scenario", "seed")

D0.obs <- rep1[varset == "X1"]$absdiff.0

sig.y <- sd(rep1[varset == "X1"]$t.hat.post)

hist(D0.obs, freq = FALSE, breaks = 35)
curve(ifelse(x > 0, 1, 0) * (dnorm(x, sd = sig.y) + dnorm(-x, sd = sig.y)), add = TRUE, col = "red")


d0all[, predy.ideal := mean(t.hat.post), by = trial.leftout]
d0all[, absdiff.ideal := abs(predy.ideal - t.hat.post)]

hist(d0all$absdiff)
hist(d0all$absdiff.ideal)


D0.true <- abs(rnorm(length(D0.obs), sd = sig.y))
wilcox.test(D0.obs, D0.true)

h1 <- outer(D0.true, D0.obs, function(x, y) x > y)
mean(h1)



d0all <- res.dt[scenario == "0"]

library(ggplot2)

ggplot(d0all, aes(x = absdiff.0)) + geom_histogram() + facet_wrap(~ trial.leftout)

hist(d0all[trial.leftout == 1]$absdiff.0, freq = FALSE)
sig.y <- sd(d0all[trial.leftout == 1]$t.hat.post)
curve(ifelse(x > 0, 1, 0) * (dnorm(x, sd = sig.y) + dnorm(-x, sd = sig.y)), add = TRUE, col = "red")


d0all[, predy.ideal := mean(t.hat.post), by = trial.leftout]
d0all[, absdiff.ideal := abs(predy.ideal - t.hat.post)]

breaks <- seq(0, max(d0all$absdiff.0), length.out = 100)

htarg <- hist(d0all$absdiff, breaks = breaks, plot = FALSE)
hideal <- hist(d0all$absdiff.ideal, breaks = breaks, plot = FALSE)
h0 <- hist(d0all$absdiff.0, breaks = breaks, plot = FALSE)


KL.div <- function(h1, h2) {
  
  breaks <- seq(0, max(c(h1,h2)), length.out = 100)
  pp <- hist(h1, breaks = breaks, plot = FALSE)
  qq <- hist(h2, breaks = breaks, plot = FALSE)
  
  excl <- pp$density == 0 | qq$density == 0
  
  sum((pp$density * (log(pp$density) - log(qq$density)))[!excl])
  
}

KL.div(d0all$absdiff, d0all$absdiff.0)
KL.div(d0all$absdiff.ideal, d0all$absdiff)

rep1[, predy.ideal := mean(t.hat.post), by = .(seed, varset, trial.leftout)]
rep1[, absdiff.ideal := abs(predy.ideal - t.hat.post)]

kl.summary <- rep1[, .(kl.0 = KL.div(absdiff, absdiff.0), kl.ideal = KL.div(absdiff.ideal, absdiff)), by = .(seed, varset)]

