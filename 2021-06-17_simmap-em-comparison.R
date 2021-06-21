source("~/2020_hOUwie/hOUwieSimmap.R")
source("~/2020_hOUwie/Utils.R")

# source("/space_2/jamesboyko/2020_hOUwie/hOUwieEM.R")
# source("/space_2/jamesboyko/2020_hOUwie/Utils.R")

require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)
require(geiger)

nCores <- 80
nTip <- 100
phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))


singleRun <- function(iter){
  fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
  fit.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(fit.cor)[1])
  pars = c(1, 1, 10, 10)  # mk, alpha, sigma, theta1, theta2
  nSim = 100
  data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
  data <- data.houwie$data
  data[data[,2]==3,2] <- 1
  data[data[,2]==4,2] <- 2
  # start with simmap
  source("/space_2/jamesboyko/2020_hOUwie/hOUwieSimmap.R")
  # source("~/2020_hOUwie/hOUwieSimmap.R")
  BMS_map <- hOUwie(phy = phy, data = data, rate.cat = 1, nSim = nSim, model.cor = "ER", model.ou = "BMS", weighted = TRUE)
  OUM_map <- hOUwie(phy = phy, data = data, rate.cat = 1, nSim = nSim, model.cor = "ER", model.ou = "OUM", weighted = TRUE)
  # got to em
  source("/space_2/jamesboyko/2020_hOUwie/hOUwieEM.R")
  # source("~/2020_hOUwie/hOUwieEM.R")
  BMS_em <- hOUwie.EM(phy = phy, data = data, discrete_model = "ER", continuous_model = "BMS", rate.cat = 1, niter = 500)
  OUM_em <- hOUwie.EM(phy = phy, data = data, discrete_model = "ER", continuous_model = "OUM", rate.cat = 1, niter = 500)
  obj <- list(BMS_map = BMS_map, OUM_map = OUM_map, BMS_em = BMS_em, OUM_em)
  save(obj, file= paste0("BMS_em-map-compare_", iter, ".Rsace"))
  return(obj)
}

test <- mclapply(1:80, function(x) singleRun(x), mc.cores = 80)

singleRun <- function(iter){
  fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
  fit.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.cor)[1])
  pars = c(1, 5, 5, 2, 10)  # mk, alpha, sigma, theta1, theta2
  data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
  nSim = 50
  data <- data.houwie$data
  data[data[,2]==3,2] <- 1
  data[data[,2]==4,2] <- 2
  # start with simmap
  # source("/space_2/jamesboyko/2020_hOUwie/hOUwieSimmap.R")
  source("~/2020_hOUwie/hOUwieSimmap.R")
  BMS_map <- hOUwie(phy = phy, data = data, rate.cat = 1, nSim = nSim, model.cor = "ER", model.ou = "BMS", weighted = TRUE)
  OUM_map <- hOUwie(phy = phy, data = data, rate.cat = 1, nSim = nSim, model.cor = "ER", model.ou = "OUM", weighted = TRUE)
  # got to em
  # source("/space_2/jamesboyko/2020_hOUwie/hOUwieEM.R")
  source("~/2020_hOUwie/hOUwieEM.R")
  BMS_em <- hOUwie.EM(phy = phy, data = data, discrete_model = "ER", continuous_model = "BMS", rate.cat = 1, niter = 500)
  OUM_em <- hOUwie.EM(phy = phy, data = data, discrete_model = "ER", continuous_model = "OUM", rate.cat = 1, niter = 500)
}



