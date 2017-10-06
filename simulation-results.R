library(ggplot2)

res <- read.csv("results/sim-results-noP-2017-08-26.csv", header = FALSE)
res <- read.csv("results/sim-results_25_2017-08-24.csv", header = FALSE)
colnames(res) <- c("scenario", "D.hat", "D.true", "variables", "true.Y", "hat.Y", "lower", "upper")

ggplot(res, aes(x = D.hat, y = D.true)) + geom_point() + 
  facet_wrap(~ scenario, scales = "free") + geom_abline(intercept = 0, slope = 1, color = "grey75", size = 1.25) + 
  theme_bw() + xlab(expression(widetilde(D))) + ylab(expression(widehat(D))) + 
  theme(axis.title.y = element_text(angle = 0.0, vjust = .5))
#ggsave("simfig.pdf", width = 7.25, height = 6.5)

troof <- c("Null", "X1", "X1", "X1", "X1|X2|X3", "X1|X2|X3")
names(troof) <- c("0", LETTERS[1:5])

res$true.vars <- troof[res$scenario]
res$exact.vars <- with(res, true.vars == variables)
res$partial.vars <- NA
for(i in 1:nrow(res)) {
  
  res$partial.vars[i] <- grepl(res$true.vars[i], res$variables[i])
  
}


library(data.table)

res.dt <- data.table(res)
sim.table <- res.dt[, list(bias.D = mean(D.hat- D.true), sd.D = sd(D.hat), sel.null = 100*mean(variables == "Null"),
              sel.exact = 100*mean(exact.vars), sel.partial = 100*mean(partial.vars), 
              bias.pred.Y = mean(hat.Y - true.Y), cover.pred.Y = 100*mean(true.Y < upper + D.hat & true.Y > lower - D.hat)), by = scenario]


library(knitr)

kable(sim.table, digits = 2, format = "latex")
kable(sim.table, digits = 2)

library(ggridges)

#res.dt <- fread("results/Dsamples-simTcmodel.csv")
res.dt <- fread("results-n12-saveD.csv")
colnames(res.dt) <- c("trial.leftout", "yhat", "t.hat.post", "yhat.0", "y.true",  "absdiff", "absdiff.0", "absdiff.true", "varset", "scenario", "seed")

setkey(res.dt, seed, scenario, varset, trial.leftout)
res.dt <- res.dt[!is.na(yhat)]


toadd <- res.dt[varset == "X1"]
toadd$absdiff <- toadd$absdiff.0
toadd$varset <- "Null"

res.dt <- rbind(res.dt, toadd)


explot <- res.dt[seed == 1466 & scenario == "0"]

explot2 <- explot[, .(absdiff, varset)]

reord <- explot2[, mean(absdiff), by = varset]
setkey(reord, V1)

explot2$varset <- factor(explot2$varset, levels = reord$varset, ordered = TRUE)

ggplot(explot2, aes(x = absdiff, y = varset)) + geom_density_ridges(bandwidth = .025) + 
  geom_vline(xintercept = reord$V1[which(reord$varset == "0")])


ggplot(explot2, aes(x = absdiff)) + geom_histogram(breaks = seq(0, 3, by = .1)) + facet_wrap(~varset)

ggplot(explot2, aes(x = absdiff)) + stat_ecdf() + facet_wrap(~varset)


mymary <- function(x, i) {
  
  hx <- hist(x, breaks = c(0,1, 2, 3, 15), plot = FALSE)
  hx$density[i]
  
}

resdt.q <- res.dt[, .(mean.D.1 =  mymary(absdiff, 1), 
                      mean.D.2 =  mymary(absdiff, 2), 
                      mean.D.3 =  mymary(absdiff, 3),
                      #mean.D.4 =  mymary(absdiff, 4),
                      Z.D0 = mean(absdiff.0),
                      P.D0.lt = mean(absdiff.0 < absdiff)), by = .(seed, scenario, varset)]


selme2 <- function(varset, a1, a2, a3, P.D0) {
  
  pp <- .5
  
  if(all(P.D0 > pp)) return("Null") else {
    
    varset[P.D0 < pp][which.max(rank(a1[P.D0 < pp]) + rank(-a2[P.D0 < pp]) + rank(-a3[P.D0 < pp]))]
    
  }
    
}

selme <- function(quant, P.D0, varlist, D0) {
  
  if(all(P.D0 > 0.1)) return("Null") else {
    
    c("Null", varlist[P.D0 < 0.3])[which.min(c(D0[1], quant[P.D0 < 0.3]))]
  }
  
}

resdt.sel <- resdt.q[, .(sel.mean = selme2(varset, mean.D.1, mean.D.2, mean.D.3, P.D0.lt)), 
                     by = .(seed, scenario)]


troof <- c("Null", "X1", "X1", "X1", "X1|X2|X3", "X1|X2|X3")
names(troof) <- c("0", LETTERS[1:5])

resdt.sel[, true.vars := troof[scenario]]

resdt.sel[, ":="(corsel.mean = sel.mean == true.vars)]
resdt.sel$corsel.mean.plus <- NA
for(i in 1:nrow(resdt.sel)){
  resdt.sel$corsel.mean.plus[i] <- grepl(resdt.sel$true.vars[i], resdt.sel$sel.mean[i], fixed = TRUE)
}
knitr::kable(resdt.sel[, .(cor.mean = mean(corsel.mean), 
                           cor.mean.plus = mean(corsel.mean.plus)), by = scenario], digits = 2)
