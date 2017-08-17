library(ggplot2)


res <- read.csv("results/sim-results-2017-08-17.csv", header = FALSE)
colnames(res) <- c("scenario", "D.hat", "D.true", "variables")

ggplot(res, aes(x = D.hat, y = D.true)) + geom_point() + 
  facet_wrap(~ scenario, scales = "free") + geom_abline(intercept = 0, slope = 1, color = "grey75", size = 1.25) + 
  theme_bw() + xlab(expression(widetilde(D))) + ylab(expression(widehat(D))) + 
  theme(axis.title.y = element_text(angle = 0.0, vjust = .5))
ggsave("simfig.pdf", width = 7.25, height = 6.5)

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
sim.table <- res.dt[, list(bias.D = mean(D.hat- D.true), sd.D = sd(D.hat - D.true), 
              sel.exact = 100*mean(exact.vars), sel.partial = 100*mean(partial.vars)), by = scenario]


library(knitr)

kable(sim.table, digits = 2, format = "latex")
