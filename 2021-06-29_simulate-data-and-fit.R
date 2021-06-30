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

nCores <- 50
nTip <- 100
phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))

CID.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
CID.cor <- getFullMat(list(CID.cor, CID.cor), CID.cor)
CID.cor <- equateStateMatPars(CID.cor, c(1,2,3))
CD.ou <- CID.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(CID.cor)[1])
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




