

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

