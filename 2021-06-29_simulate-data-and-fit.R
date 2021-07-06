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
require(lhs)

# prerequisites
nCores <- 80
nIter <- 80
nTip <- 100

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
generateContinuousParameters <- function(continuous_model, minAlpha, maxAlpha, minSigma2, maxSigma2, minTheta, maxTheta){
  k.alpha <- length(unique(na.omit(continuous_model[1,])))
  k.sigma <- length(unique(na.omit(continuous_model[2,])))
  k.theta <- length(unique(na.omit(continuous_model[3,])))
  par_alpha <- runif(k.alpha, minAlpha, maxAlpha)
  par_sigma <- runif(k.sigma, minSigma2, maxSigma2)
  par_theta <- runif(k.theta, minTheta, maxTheta)
  p <- c(par_alpha, par_sigma, par_theta)
  return(p)
}

minAlpha = 0 
maxAlpha = 2
minSigma2 = 0.1
maxSigma2 = 4
minTheta = 2
maxTheta = 10

# simulate a dataset and return all simulation information
generateContinuousParameters(continuous_models_cd[,,1], minAlpha, maxAlpha, minSigma2, maxSigma2, minTheta, maxTheta)



pars <- apply(continuous_models_cd, 3, function(x) sapply(1:nIter, function(y) generateContinuousParameters(x, minAlpha, maxAlpha, minSigma2, maxSigma2, minTheta, maxTheta)))



nCores <- 50
nTip <- 100

CID.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
CID.cor <- getFullMat(list(CID.cor, CID.cor), CID.cor)
CID.cor <- equateStateMatPars(CID.cor, c(1,2,3))
CD.ou <- CID.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, 2)
CID.ou[3,] <- c(3,3,4,4)
CD.ou[3,] <- c(3,4,3,4)

CDBMS.ou <- CIDBMS.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(CID.cor)[1])
CIDBMS.ou[2,] <- c(1,1,2,2)
CDBMS.ou[2,] <- c(1,2,1,2)
CDBMS.ou[3,] <- CIDBMS.ou[3,] <- c(3,3,3,3)


fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.cor)[1])
pars = c(0.1, 1, 1, 5, 15)  # mk, alpha, sigma, theta1, theta2
data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
nSim = 50
data <- data.houwie$data
data[data[,2]==3,2] <- 1
data[data[,2]==4,2] <- 2
#out <- hOUwie(phy = phy, data = data, rate.cat = 1, nSim = nSim, discrete_model = fit.cor, continuous_model = fit.ou)
q0.1 <- hOUwie(phy = phy, data = data, rate.cat = 1, nSim = nSim, discrete_model = fit.cor, continuous_model = fit.ou,
       p = c(0.1, 1, 1, 5, 15))
q1.0 <- hOUwie(phy = phy, data = data, rate.cat = 1, nSim = nSim, discrete_model = fit.cor, continuous_model = fit.ou,
       p = c(1, 1, 1, 5, 15))
q10. <- hOUwie(phy = phy, data = data, rate.cat = 1, nSim = nSim, discrete_model = fit.cor, continuous_model = fit.ou,
       p = c(10, 1, 1, 5, 15))

singleRun <- function(iter){
  fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
  fit.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.cor)[1])
  pars = c(0.1, 5, 2, 5, 15)  # mk, alpha, sigma, theta1, theta2
  data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
  nSim = 50
  data <- data.houwie$data
  data[data[,2]==3,2] <- 1
  data[data[,2]==4,2] <- 2
  out <- hOUwie(phy = phy, data = data, rate.cat = 1, nSim = nSim, discrete_model = fit.cor, continuous_model = fit.ou)
  return(out)
}

res <- mclapply(1:nCores, function(x) singleRun(x), mc.cores = nCores)




