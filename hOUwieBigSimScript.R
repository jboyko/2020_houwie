source("~/2020_hOUwie/hOUwieSimmap.R")
source("~/2020_hOUwie/Utils.R")

source("/space_2/jamesboyko/2020_hOUwie/hOUwieSimmap.R")
source("/space_2/jamesboyko/2020_hOUwie/Utils.R")

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


# set up the discrete character model
Disc_model <- equateStateMatPars(getRateCatMat(2), c(1,2))
Disc_model <- getFullMat(list(Disc_model, Disc_model), Disc_model)
Disc_model <- equateStateMatPars(Disc_model, c(1,2,3))

# set up the continuous character model
HYBBMS <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(Disc_model)[1])
CDBMS <- CIDBMS <- HYBBMS
CDBMS[2,] <- c(1,2,1,2)
CDBMS[3,] <- 3
CIDBMS[2,] <- c(1,1,2,2)
CIDBMS[3,] <- 3
# OUM models set up to be character dependent, independent, and a mix of both
HYBOUM <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(Disc_model)[1])
CDOUM <- CIDOUM <- HYBOUM
CDOUM[3,] <- c(3,4,3,4)
CIDOUM[3,] <- c(3,3,4,4)

HYBOUM <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(Disc_model)[1])
CDOUMV <- CIDOUMV <- HYBOUM
CDOUMV[2,] <- c(2,3,2,3)
CDOUMV[3,] <- c(4,5,4,5)
CIDOUMV[3,] <- c(2,2,3,3)
CIDOUMV[3,] <- c(4,4,5,5)


singleRun <- function(iter, pars, nSim, sim.index_disc, sim.index_ou, out.dir){
  data.houwie <- generateData(phy, sim.index_disc, sim.index_ou, pars)
  data <- data.houwie$data
  data[data[,2]==3,2] <- 1
  data[data[,2]==4,2] <- 2
  BM1 <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = Disc_model, model.ou = "BM1", weighted = TRUE)
  OU1 <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = Disc_model, model.ou = "OU1", weighted = TRUE)
  CDOUM <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = Disc_model, index.ou = CDOUM, weighted = TRUE)
  CDBMS <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = Disc_model, index.ou = CDBMS, weighted = TRUE)
  CIDOUM <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                   index.cor = Disc_model, index.ou = CIDOUM, weighted = TRUE)
  CIDBMS <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                   index.cor = Disc_model, index.ou = CIDBMS, weighted = TRUE)
  HYBOUM <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                   index.cor = Disc_model, index.ou = HYBOUM, weighted = TRUE)
  CDOUMV <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                   index.cor = Disc_model, index.ou = CDOUMV, weighted = TRUE)
  CIDOUMV <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                   index.cor = Disc_model, index.ou = CIDOUMV, weighted = TRUE)
  out <- list("BM" = BM1, "OU1" = OU1, "CDOUM" = CDOUM, "CDBMS" = CDBMS, "CIDOUM" = CIDOUM, "CIDBMS" = CIDBMS, "HYBOUM" = HYBOUM, "CIDOUMV" = CIDOUMV, "CDOUMV" = CDOUMV)
  save(out, file = paste0(out.dir, "hOUwie-", iter, ".Rsave"))
  return(out)
}


out.dirs <- list("/space_2/jamesboyko/2020_hOUwie/ModelTesting/CDOUM/",
                 "/space_2/jamesboyko/2020_hOUwie/ModelTesting/CDBMS/",
                 "/space_2/jamesboyko/2020_hOUwie/ModelTesting/CIDOUM/",
                 "/space_2/jamesboyko/2020_hOUwie/ModelTesting/CIDBMS/",
                 "/space_2/jamesboyko/2020_hOUwie/ModelTesting/HYBOUM/",
                 "/space_2/jamesboyko/2020_hOUwie/ModelTesting/CDOUMV/",
                 "/space_2/jamesboyko/2020_hOUwie/ModelTesting/CIDOUMV/")
sim.indexes_ou <- list(CDOUM = CDOUM,
                       CDBMS = CDBMS,
                       CIDOUM = CIDOUM,
                       CIDBMS = CIDBMS,
                       HYBOUM = HYBOUM,
                       CDOUMV = CDOUMV,
                       CIDOUMV = CIDOUMV)
pars <- list(c(1, 5, 5, 2, 10),
             c(1, 1, 10, 10),
             c(1, 5, 5, 2, 10),
             c(1, 1, 10, 10),
             c(1, 5, 5, 2, 10, 10, 18),
             c(1, 5, 10, 5, 2, 10),
             c(1, 5, 10, 5, 2, 10))


for(i in 1:length(pars)){
  res <- mclapply(1:100, function(x) singleRun(x, pars[[i]], 100, Disc_model, sim.indexes_ou[[i]], out.dirs[[i]]), mc.cores = nCores)
}


res <- mclapply(1:100, function(x) singleRun(x), mc.cores = nCores)

getModelTable(list("BM" = BM1, "OU1" = OU1, "CDOUM" = OUM, "CDBMS" = BMS, "CIDOUM" = CIDOUM, "CIDBMS" = CIDBMS))
