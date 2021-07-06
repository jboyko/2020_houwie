setwd("~/2020_hOUwie/")

source("hOUwieSimmap.R")
source("Utils.R")

require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)
require(geiger)

# prerequisites
nCores <- 80
nIter <- 80
nTip <- 100
nSim <- 50

# the 32 2-state models we want to fit
continuous_models_cd <- getAllContinuousModelStructures(2)
discrete_model_cd <- equateStateMatPars(getRateCatMat(2), 1:2)
continuous_models_cid <- continuous_models_cd[,c(1,1,2,2),]
discrete_model_cid <- getFullMat(list(discrete_model_cd, discrete_model_cd), discrete_model_cd)
# dist(t(apply(continuous_models_cd[,,1:5], 3, function(x) c(x[1,], x[2,] - min(x[2,]) + 1, x[3,] - min(x[3,]) + 1))))

# deciding what tree size we would like to use
phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))

# deciding which params to use. we should vary delta sigma, delta alpha and delta theta. this will be the number of replicates we have per generating model.
# a function which generates a parameter set given a continuous model
generateParameters <- function(continuous_model, minAlpha, maxAlpha, minSigma2, maxSigma2, minTheta, maxTheta, discrete_model, minRate, maxRate){
  k.alpha <- length(unique(na.omit(continuous_model[1,])))
  k.sigma <- length(unique(na.omit(continuous_model[2,])))
  k.theta <- length(unique(na.omit(continuous_model[3,])))
  k.discrete <- max(na.omit(discrete_model))
  par_alpha <- runif(k.alpha, minAlpha, maxAlpha)
  par_sigma <- runif(k.sigma, minSigma2, maxSigma2)
  par_theta <- runif(k.theta, minTheta, maxTheta)
  par_rates <- runif(k.discrete, minRate, maxRate)
  p <- c(par_rates, par_alpha, par_sigma, par_theta)
  return(p)
}

singleRun <- function(full_data, continuous_models_cd, continuous_models_cid, discrete_model_cd, discrete_model_cid){
  data <- full_data$data
  data[data[,2]==3,2] <- 1
  data[data[,2]==4,2] <- 2
  cd_out <- apply(continuous_models_cd, 3, function(x) hOUwie(phy = phy, data = data, rate.cat = 1, nSim = nSim, discrete_model = discrete_model_cd, continuous_model = x))
  cat("CD models complete.")
  cid_out <- apply(continuous_models_cid, 3, function(x) hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, discrete_model = discrete_model_cid, continuous_model = x))
  cat("CID models complete.")
  return(c(cd_out, cid_out))
}

minAlpha = 0 
maxAlpha = 2
minSigma2 = 0.1
maxSigma2 = 4
minTheta = 2
maxTheta = 10
minRate = 0.25
maxRate = 0.75
# simulate a dataset and return all simulation information
# generateParameters(continuous_models_cid[,,16], minAlpha, maxAlpha, minSigma2, maxSigma2, minTheta, maxTheta, discrete_model_cid, minRate, maxRate)

pars_cd <- apply(continuous_models_cd, 3, function(x) lapply(1:nIter, function(y) generateParameters(x, minAlpha, maxAlpha, minSigma2, maxSigma2, minTheta, maxTheta, discrete_model_cd, minRate, maxRate)))
pars_cid <- apply(continuous_models_cid, 3, function(x) lapply(1:nIter, function(y) generateParameters(x, minAlpha, maxAlpha, minSigma2, maxSigma2, minTheta, maxTheta, discrete_model_cid, minRate, maxRate)))

for(i in 1:length(pars)){
  cat("Begining parameter set", i, "...\n")
  generating_model <- continuous_models_cd[,,i]
  generating_model_pars <- pars_cd[[i]]
  all_data <- lapply(generating_model_pars, function(x) generateData(phy, discrete_model_cd, generating_model, x))
  save(all_data, file = paste0("sim_data/", "cd-data-generating_model_", i, ".Rsave"))
  out <- mclapply(all_data, function(x) singleRun(x, continuous_models_cd, continuous_models_cid, discrete_model_cd, discrete_model_cid), mc.cores = nCores)
  save(out, file = paste0("sim_fits/", "cd-fit-generating_model_", i, ".Rsave"))
}


# res <- mclapply(1:nCores, function(x) singleRun(x), mc.cores = nCores)

