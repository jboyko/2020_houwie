# get a phylogeny rescaled to a height of 1
source("~/2020_hOUwie/hOUwie.R")
source("~/2020_hOUwie/Utils.R")
#source("/space_2/jamesboyko/2020_hOUwie/hOUwie.R")
#source("/space_2/jamesboyko/2020_hOUwie/Utils.R")
require(OUwie)
require(corHMM)
require(parallel)

# preliminaries
n.iter <- 100
nmaps <- 250
ncores <- 30


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ER BMS  
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

tmp <- generateData(phy, index.cor, index.ou, pars)

## vary the Mk rate
nTip <- 100
phy <- sim.bdtree(b = 1, d = 0.5, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))
index.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
index.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(index.cor)[1])

pars <- c(1, 1, 10, 5)
mclapply(1:n.iter, function(x) 
  SingleModelTestRun(phy, index.cor, index.ou, "ER", "BMS", 1, pars, nmaps, x, nstarts = 10), 
  mc.cores = ncores)

pars <- c(2, 1, 10, 5)
mclapply(1:n.iter, function(x) 
  SingleModelTestRun(phy, index.cor, index.ou, "ER", "BMS", 1, pars, nmaps, x), 
  mc.cores = ncores)

pars <- c(4, 1, 10, 5)
mclapply(1:n.iter, function(x) 
  SingleModelTestRun(phy, index.cor, index.ou, "ER", "BMS", 1, pars, nmaps, x),
  mc.cores = ncores)

## vary the tree size
pars <- c(1, 1, 10, 5)
mclapply(1:n.iter, function(x) 
  SingleModelTestRun(phy, index.cor, index.ou, "ER", "BMS", 1, pars, nmaps, x), 
  mc.cores = ncores)

pars <- c(2, 1, 10, 5)
mclapply(1:n.iter, function(x) 
  SingleModelTestRun(phy, index.cor, index.ou, "ER", "BMS", 1, pars, nmaps, x), 
  mc.cores = ncores)

pars <- c(4, 1, 10, 5)
mclapply(1:n.iter, function(x) 
  SingleModelTestRun(phy, index.cor, index.ou, "ER", "BMS", 1, pars, nmaps, x),
  mc.cores = ncores)



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ER OUM - vary the sigma alpha ratio
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# nTip <- 50
phy <- sim.bdtree(b = 1, d = 0.5, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))
index.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
index.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(index.cor)[1])

# for server
pars <- c(1, 2, 1, 5, 10)
mclapply(1:n.iter, function(x) 
  SingleModelTestRun(phy, index.cor, index.ou, "ER", "OUM", 1, pars, nmaps, x, nstarts = 10), 
  mc.cores = ncores)

pars <- c(1, 8, 0.5, 5, 10)
mclapply(1:n.iter, function(x) 
  SingleModelTestRun(phy, index.cor, index.ou, "ER", "OUM", 1, pars, nmaps, x, nstarts = 10), 
         mc.cores = ncores)

pars <- c(1, 10, 0.1, 5, 10)
mclapply(1:n.iter, function(x) 
  SingleModelTestRun(phy, index.cor, index.ou, "ER", "OUM", 1, pars, nmaps, x, nstarts = 10), 
  mc.cores = ncores)

pars <- c(0.25, 2, 1, 5, 10)
mclapply(1:n.iter, function(x) 
  SingleModelTestRun(phy, index.cor, index.ou, "ER", "OUM", 1, pars, nmaps, x, nstarts = 10), 
  mc.cores = ncores)

pars <- c(1, 2, 2, 5, 10)
mclapply(1:n.iter, function(x) 
  SingleModelTestRun(phy, index.cor, index.ou, "ER", "OUM", 1, pars, nmaps, x, nstarts = 10), 
  mc.cores = ncores)

pars <- c(1, 0.1, 0.2, 5, 10)
mclapply(1:n.iter, function(x) 
  SingleModelTestRun(phy, index.cor, index.ou, "ER", "OUM", 1, pars, nmaps, x, nstarts = 10), 
  mc.cores = ncores)

pars <- c(1, 1, 2, 5, 10)
mclapply(1:n.iter, function(x) 
  SingleModelTestRun(phy, index.cor, index.ou, "ER", "OUM", 1, pars, nmaps, x, nstarts = 10), 
  mc.cores = ncores)

pars <- c(1, 10, 20, 5, 10)
mclapply(1:n.iter, function(x) 
  SingleModelTestRun(phy, index.cor, index.ou, "ER", "OUM", 1, pars, nmaps, x, nstarts = 10), 
  mc.cores = ncores)


pars <- c(1, 1, 2, 5, 10)
mclapply(1:n.iter, function(x) 
  SingleModelTestRunParsimony(phy, index.cor, index.ou, "ER", "OUM", 1, pars, nmaps, x, nstarts = 10), 
  mc.cores = ncores)


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ER OUMVA
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
phy <- sim.bdtree(b = 1, d = 0.5, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))
index.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
index.ou <- getOUParamStructure("OUMVA", "three.point", FALSE, FALSE, dim(index.cor)[1])

# for server
pars <- c(1, 4, 1, 1, 0.25, 5, 10)
mclapply(1:100, function(x) SingleModelTestRun(phy, index.cor, index.ou, "ER", "OUMVA", 1, pars, nmaps, x), mc.cores = ncores)

pars <- c(2, 4, 1, 1, 0.25, 5, 10)
mclapply(1:100, function(x) SingleModelTestRun(phy, index.cor, index.ou, "ER", "OUMVA", 1, pars, nmaps, x), mc.cores = ncores)

pars <- c(4, 4, 1, 1, 0.25, 5, 10)
mclapply(1:120, function(x) SingleModelTestRun(phy, index.cor, index.ou, "ER", "OUMVA", 1, pars, nmaps, x), mc.cores = ncores)

pars <- c(8, 4, 1, 1, 0.25, 5, 10)
mclapply(1:120, function(x) SingleModelTestRun(phy, index.cor, index.ou, "ER", "OUMVA", 1, pars, nmaps, x), mc.cores = ncores)


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## CID OUM
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
nTip <- 500
phy <- sim.bdtree(b = 1, d = 0.5, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))
sim.index.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
sim.index.cor <- getFullMat(list(sim.index.cor, sim.index.cor), sim.index.cor)
sim.index.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(sim.index.cor)[1])
sim.index.ou[3, ] <- c(3,3,4,4)

fit.index.cor1 <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.index.ou1 <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.index.cor1)[1])

fit.index.cor2 <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.index.cor2 <- getFullMat(list(fit.index.cor2, fit.index.cor2), fit.index.cor2)
fit.index.cor2 <- equateStateMatPars(fit.index.cor2, c(1,2,3))
fit.index.ou2 <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.index.cor2)[1])
fit.index.ouCD <- fit.index.ou2
fit.index.ouCD[3,] <- c(3,4,3,4)
fit.index.ouCID <- fit.index.ou2
fit.index.ouCID[3,] <- c(3,3,4,4)

fit.index.cor <- list(fit.index.cor1, fit.index.cor2, fit.index.cor2, fit.index.cor2)
fit.index.ou <- list(fit.index.ou1, fit.index.ouCD, fit.index.ouCID, fit.index.ou2)
fit.rate.cat <- c(1,2,2,2)

pars <- c(1, 1, 1,
          10, 20, 5, 10)
pars <- c(1, 1, 1,
          1, 2, 5, 10)

mclapply(1:100, function(x) ModelSetRun(phy, pars, "CID_OUM", sim.index.cor, sim.index.ou, 2, fit.index.cor, fit.index.ou, fit.rate.cat, nmaps, x), mc.cores = ncores)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## CD OUM
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
nTip <- 100
phy <- sim.bdtree(b = 1, d = 0.5, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))
sim.index.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
sim.index.cor <- getFullMat(list(sim.index.cor, sim.index.cor), sim.index.cor)
sim.index.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(sim.index.cor)[1])
sim.index.ou[3, ] <- c(3,4,3,4)

fit.index.cor1 <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.index.ou1 <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.index.cor1)[1])

fit.index.cor2 <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.index.cor2 <- getFullMat(list(fit.index.cor2, fit.index.cor2), fit.index.cor2)
fit.index.cor2 <- equateStateMatPars(fit.index.cor2, c(1,2,3))
fit.index.ou2 <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.index.cor2)[1])
fit.index.ouCD <- fit.index.ou2
fit.index.ouCD[3,] <- c(3,4,3,4)
fit.index.ouCID <- fit.index.ou2
fit.index.ouCID[3,] <- c(3,3,4,4)

fit.index.cor <- list(fit.index.cor1, fit.index.cor2, fit.index.cor2, fit.index.cor2)
fit.index.ou <- list(fit.index.ou1, fit.index.ouCD, fit.index.ouCID, fit.index.ou2)
fit.rate.cat <- c(1,2,2,2)

pars <- c(1, 1, 1,
          10, 20, 5, 10)
pars <- c(1, 1, 1,
          1, 2, 5, 10)

mclapply(1:100, function(x) ModelSetRun(phy, pars, "CD_OUM", sim.index.cor, sim.index.ou, 2, fit.index.cor, fit.index.ou, fit.rate.cat, nmaps, x), mc.cores = ncores)


