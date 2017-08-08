

sample_from_model <- function(y, x1, y.se, x1.se, xpred, ...) {
  
  ses <- cbind(y.se, x1.se)
  w <- cbind(y, x1)
  
  
  prior <- get_prior(w)
    
  mcmc <- list(nburn=5000,
               nsave=1000,
               nskip=2,
               ndisplay=1e6)
  
  mygrid <- seq(min(y) - 2 * sd(y), max(y) + 2 * sd(y), length.out = 500)
  
  fitLDDP <- DPcdensity(y = y, x = x1, mus=w, sds=ses, xpred = xpred,
                        grid=mygrid, 
                        compute.band=FALSE,
                        type.band="HPD",
                        prior=prior, 
                        mcmc=mcmc, 
                        state=NULL, 
                        status=TRUE)
  
  
  y.hat.samps <- lapply(1:length(y), function(i){
    
    data.frame(trial = i, t.hat.post = sample_from_density(100, fitLDDP$grid, fitLDDP$densp.m[i, ]))
    
  })
  
  do.call(rbind, y.hat.samps)
  
}



predict_from_model <- function(y, x1, y.se, x1.se, xpred, ...) {
  
  ses <- cbind(y.se, x1.se)
  w <- cbind(y, x1)
  
  
  prior <- get_prior(w)
  
  mcmc <- list(nburn=5000,
               nsave=1000,
               nskip=2,
               ndisplay=1e6)
  
  
  fitLDDP <- DPcdensity(y = y, x = x1, mus=w, sds=ses, xpred = xpred,
                        ngrid=1, 
                        compute.band=FALSE,
                        type.band="HPD",
                        prior=prior, 
                        mcmc=mcmc, 
                        state=NULL, 
                        status=TRUE)
  

  data.frame(xpred = xpred, yhat = fitLDDP$meanfp.m)
  
  
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
  
  
}


