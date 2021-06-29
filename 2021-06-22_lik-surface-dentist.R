# testing the likleihood surface of various model fits
source("hOUwieSimmap.R")
require(dentist)
require(parallel)

getDentResults <- function(houwie.model){
  dentFunc <- function(par){
    out <- hOUwie(houwie.model$phy, houwie.model$data, houwie.model$rate.cat, discrete_model = houwie.model$discrete_model, continuous_model = houwie.model$continuous_model, houwie.model$nSim, p = par)
    return(-out$loglik)
  }
  pars <- houwie.model$p
  names(pars) <- 1:length(pars)
  dented_results <- dent_walk(par=pars, fn=dentFunc, best_neglnL=-houwie.model$loglik, delta=2, nsteps=1000, print_freq=250)
  return(dented_results)
}

# scales results
load("empirical/scales2009/scales-results.Rsave")
dent.out <- mclapply(results, function(x) getDentResults(x), mc.cores = length(results))

