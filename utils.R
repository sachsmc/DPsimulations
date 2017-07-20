
get_prior <- function(w) {
  
  wbar <- apply(w, 2, mean) 
  wcov <- var(w) 
  list(a0 = 10, 
       b0 = 1, 
       nu1 = 4, 
       nu2 = 4, 
       s2 = 0.5 * wcov, 
       m2 = wbar, 
       psiinv2 = 2 * solve(wcov),
       tau1 = 6.01, tau2 = 3.01)
  
}

sample_from_density <- function(n, x, y) {
  
  sample(x, size = n, prob = y, replace = TRUE)
  
}

