library(ggplot2)

KL.div <- function(h1, h2) {
  
  breaks <- seq(0, max(c(h1,h2)), length.out = 50)
  pp <- hist(h1, breaks = breaks, plot = FALSE)
  qq <- hist(h2, breaks = breaks, plot = FALSE)
  
  excl <- pp$density == 0 | qq$density == 0
  
  sum((pp$density * (log(pp$density) - log(qq$density)))[!excl])
  
}


mymary <- function(x, i) {
  
  hx <- hist(x, breaks = c(0,1, 2, 3, 1e5), plot = FALSE)
  hx$density[i]
  
}

library(data.table)

marg.dt <- fread("results/marginals-n12.csv")
#marg.dt <- fread("results/linear-marginals-n12.csv")
colnames(marg.dt) <- c("trial.leftout", "t.hat.post", "seed", "scenario")

pred.dt <- fread("results/predictions-n12.csv")
#pred.dt <- fread("results/linear-pred-n12.csv")
colnames(pred.dt) <- c("trial.leftout", "yhat", "yhat.lower", "yhat.upper", "yhat.0", "y.true", "absdiff.true", "varset", "scenario", "seed")

pred.dt <- pred.dt[!grepl(":", varset)]

setkey(marg.dt, seed, scenario, trial.leftout)
setkey(pred.dt, seed, scenario, trial.leftout)

seeds <- unique(pred.dt$seed)

res.fin <- vector("list", length = length(seeds))

for(j in seeds) {
  A <- marg.dt[trial.leftout < 13 & seed == j][, .SD[sample(1:nrow(.SD), 500, replace = TRUE),], by = .(scenario, trial.leftout)]
  B <- pred.dt[trial.leftout < 13 & seed == j]
  setkey(A, scenario, trial.leftout)
  setkey(B, scenario, trial.leftout)
  
  res.dt <- A[B, allow.cartesian = TRUE]
  res.dt[, predy.ideal := mean(t.hat.post), by = .(scenario, varset, trial.leftout)]
  res.dt[, absdiff.ideal := abs(predy.ideal - t.hat.post)]
  res.dt[, absdiff := abs(yhat - t.hat.post)]
  res.dt[, absdiff.0 := abs(yhat.0 - t.hat.post)]
  res.dt[, varset.minD := .SD[, .(muD = mean(absdiff)), by = varset]$varset[which.min(.SD[, .(muD = mean(absdiff)), by = varset]$muD)], by = .(scenario)]
  res.dt[, inseed := 1:.N, by = .(scenario, varset, trial.leftout)]
  tomerge <- res.dt[varset == varset.minD]
  
  res.dt <- merge(res.dt, tomerge[, .(scenario, trial.leftout, inseed, varset.minD, absdiff.minD = absdiff)], by = c("scenario", "trial.leftout", "inseed"))
  
  resdt.q <- res.dt[, .(seed = j, 
                        mean.D = mean(absdiff),
                        med.D0 = median(absdiff.0),
                        Z.D0 = mean(absdiff.0),
                        P.D0.lt = mean(absdiff.0 < absdiff), 
                        P.Dmin.lt = mean(absdiff.minD < absdiff),
                        P.Dmin.lt0 = mean(absdiff.minD < absdiff.0),
                        median.D = median(absdiff),
                         mean.trueD = mean(absdiff.true), 
                        mean.true.D0 = mean(abs(yhat.0 - y.true))), by = .(scenario, varset)]
  
  
  
  res.fin[[which(seeds == j)]] <- resdt.q
  rm(tomerge, A, B)
  gc()
  
}

resdt.q <- do.call(rbind, res.fin)

setkey(resdt.q, seed, scenario, varset)

toadd <- resdt.q[varset == "X1"]
toadd$mean.D <- toadd$Z.D0
toadd$median.D <- toadd$med.D0
toadd$mean.trueD <- toadd$mean.true.D0
toadd$P.D0.lt <- 0
toadd$varset <- "Null"
toadd$P.Dmin.lt <- toadd$P.Dmin.lt0

resdt.q <- rbind(resdt.q, toadd)

### all pairwise comparisons including D0
p.sel <- pred.dt
p.sel[, absdiff.0 := abs(yhat.0 - y.true)]

toadd <- p.sel[varset == "X1"]
toadd$absdiff.true <- toadd$absdiff.0
toadd$varset <- "Null"

p.sel <- rbind(p.sel, toadd)

## screen sets that satisfy P(D0 < D) < 0.5

p.sel[, pass.screen := mean(absdiff.0 < absdiff.true) < 0.5, by = .(seed, scenario, varset)]
p.sel2 <- p.sel[pass.screen == TRUE]

mean.pairwise <- function(varset, absdiff.true) {
  
  ## for each k in varset, get average P(Dk > Dl) for all other l
  nl <- unique(varset)
  res <- lapply(nl, function(x) {
    
    mean(outer(absdiff.true[varset == x], absdiff.true[varset != x], function(x1, x2) x1 > x2))
    
  })
  
  nl[which.min(unlist(res))]
  
}


pair.sel <- p.sel2[, .(sel.pair.pees = mean.pairwise(varset, absdiff.true)), by = .(seed, scenario)]
troof <- c("Null", "X1", "X1", "X1", "X1|X2|X3", "X1|X2|X3")
names(troof) <- c("0", LETTERS[1:6])

pair.sel[, true.vars := troof[scenario]]

knitr::kable(pair.sel[, .(cor = mean(sel.pair.pees == true.vars)), by = .(scenario)], digits = 2)

####



selme <- function(quant, P.D0, varlist) {
  
  if(all(P.D0 > 0.5)) return("Null") else {
    
    varlist[P.D0 < 0.5][which.min(quant[P.D0 < 0.5])]
  }
  
}


selme3 <- function(quant, P.D0, P.Dide, varlist) {
  
  if(all(P.D0 > 0.4)) return("Null") else {
    
    
    sc1var <- varlist[P.D0 < 0.4]
    minD <- which.min(quant[P.D0 < 0.4])
    minVar <- sc1var[minD]
    
    newCand <- sc1var[P.Dide[P.D0 < 0.4] < 0.5]
    newCandDs <- quant[P.D0 < 0.4][P.Dide[P.D0 < 0.4] < 0.5]
    if(length(newCand) == 0) return(minVar) else {
    
    if("Null" %in% newCand) return("Null") else {
      
      newCand[order(nchar(newCand), newCandDs)][1]
      
    }
    }
    
  }
  
}

resdt.sel <- resdt.q[!is.na(mean.D), .(sel.mean = selme(mean.D, P.D0.lt, varset), 
                         sel.median = selme(median.D, P.D0.lt, varset), 
                         sel.trueD = selme(mean.trueD, P.D0.lt, varset), 
                         sel.tab2 = selme3(mean.D, P.D0.lt, P.Dmin.lt, varset)), 
                     by = .(seed, scenario)]


troof <- c("Null", "X1", "X1", "X1", "X1|X2|X3", "X1|X2|X3", "X1|X2")
names(troof) <- c("0", LETTERS[1:6])

resdt.sel[, true.vars := troof[scenario]]

resdt.sel[, ":="(corsel.mean = sel.mean == true.vars, 
                 corsel.median = sel.median == true.vars, 
                 corsel.true = sel.trueD == true.vars)]
resdt.sel$corsel.mean.plus <- NA
resdt.sel$corsel.median.plus <- NA
resdt.sel$corsel.true.plus <- NA
for(i in 1:nrow(resdt.sel)){
  resdt.sel$corsel.mean.plus[i] <- grepl(resdt.sel$true.vars[i], resdt.sel$sel.mean[i], fixed = TRUE)
  resdt.sel$corsel.median.plus[i] <- grepl(resdt.sel$true.vars[i], resdt.sel$sel.median[i], fixed = TRUE)
  resdt.sel$corsel.true.plus[i] <- grepl(resdt.sel$true.vars[i], resdt.sel$sel.trueD[i], fixed = TRUE)
}


# dpp.simtab2 <- resdt.sel[, .(cor.mean = mean(corsel.mean), 
#                            cor.mean.plus = mean(corsel.mean.plus), 
#                            cor.Dhat = mean(corsel.true), 
#                            cor.Dhat.plus = mean(corsel.true.plus)), by = scenario]
# 
# 
# save(dpp.simtab2, file = "sim-table2-dpp.RData")
# 


resdt.sel[, ":="(corsel.tab2 = sel.tab2 == true.vars,
                 nullsel.tab2 = sel.tab2 == "Null")]
resdt.sel$corsel.tab2.plus <- NA
for(i in 1:nrow(resdt.sel)){
  resdt.sel$corsel.tab2.plus[i] <- grepl(resdt.sel$true.vars[i], resdt.sel$sel.tab2[i], fixed = TRUE)
}

knitr::kable(resdt.sel[, .(cor = mean(corsel.tab2), 
              corplos = mean(corsel.tab2.plus), 
              null = mean(nullsel.tab2)), by = scenario], digits = 2, format = "latex")


## confidence band coverage

topred <- pred.dt[trial.leftout == 13]

setkey(resdt.sel, seed, scenario, sel.mean)
setkey(resdt.q, seed, scenario, varset)

selected <- merge(resdt.sel, resdt.q[, .(scenario, varset, seed, mean.D)], by.x = c("seed", "scenario", "sel.mean"), 
      by.y = c("seed", "scenario", "varset"), all.x = TRUE, all.y = FALSE)

checkcover <- merge(topred, selected, by.x = c("seed", "scenario", "varset"), by.y = c("seed", "scenario", "sel.mean"))

checkcover[, covered := yhat.lower - mean.D < y.true & yhat.upper + mean.D > y.true]

cover.linear <- checkcover[, .(coverage = mean(covered), mean.Dhat = mean(absdiff.true), 
                               mean.D = mean(mean.D), sd.D = sd(mean.D)), by = "scenario"]

save(cover.linear, file = "sim-table1-linear.RData")


tab2.res <- cbind(cover.dpp[, .(scenario, sprintf("%.2f", coverage), sprintf("%.2f", mean.Dhat), 
                                sprintf("%.2f (%.2f)", mean.D, sd.D))], 
cover.linear[, .(sprintf("%.2f", coverage), sprintf("%.2f", mean.Dhat), 
                 sprintf("%.2f (%.2f)", mean.D, sd.D))])

knitr::kable(tab2.res, format = "latex")


## scenario E for larger samples


library(ggplot2)
library(data.table)
library(gridExtra)
n12 <- fread("results/predictions-n12.csv")
n25 <- fread("R0fileN25.csv")
n50 <- fread("R0fileN50.csv")
n100 <- fread("R0fileN100.csv")

colnames(n25) <- colnames(n50) <- colnames(n100) <-  c("trial.leftout", "yhat", "t.hat.post", "yhat.0", "y.true",  "absdiff", "absdiff.0", "absdiff.true", "varset", "scenario", "seed")
colnames(n12) <- c("trial.leftout", "yhat", "yhat.lower", "yhat.upper", "yhat.0", "y.true", "absdiff.true", "varset", "scenario", "seed")

db <- stat_smooth(se = FALSE, method = "loess", span = 1)
li <- c(-1, 3)
p1 <- ggplot(n12[scenario == "E" & varset == "X1|X2|X3" & seed == 7272], aes(x = yhat, y = y.true)) + geom_point() + db +
  ylim(li) + xlim(li) + ggtitle("N = 12")
p2 <- ggplot(n25[scenario == "E" & varset == "X1|X2|X3"], aes(x = yhat, y = y.true)) + geom_point() + db + 
  ylim(li) + xlim(li)+ ggtitle("N = 25")
p3 <- ggplot(n50[scenario == "E" & varset == "X1|X2|X3"], aes(x = yhat, y = y.true)) + geom_point() + db + 
  ylim(li) + xlim(li)+ ggtitle("N = 50")
p4 <- ggplot(n100[scenario == "E" & varset == "X1|X2|X3"], aes(x = yhat, y = y.true)) + geom_point() + db + 
  ylim(li) + xlim(li)+ ggtitle("N = 100")

pdf("supple-figure.pdf", width = 6.75, height = 6.5)
grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()
