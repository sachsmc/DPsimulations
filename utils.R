
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
  
  unlist(sapply(1:length(vars), function(i) {
    combn(vars, i, FUN = function(x) paste(x, collapse = "|"), simplify = TRUE)
  }))
}
