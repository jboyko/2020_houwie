# this script will generate the datasets that we are interested in comparing amongst the different hOUwie models
# what are the problems of current methods? what are the consequences of ignoring the evolution of the discrete character when modeling an OU process? what are the consequences of ignoring the continuous character when modeling an Mk? 

# imports
source("~/2020_hOUwie/hOUwieSimmap.R")
source("~/2020_hOUwie/Utils.R")

# source("/space_2/jamesboyko/2020_hOUwie/hOUwieSimmap.R")
# source("/space_2/jamesboyko/2020_hOUwie/Utils.R")

require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)
require(geiger)

# parameter definition
nTip <- 100
phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))

fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.cor)[1])
pars = c(1, 5, 5, 2, 10)  # mk, alpha, sigma, theta1, theta2
data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
data <- data.houwie$data
data[data[,2]==3,2] <- 1
data[data[,2]==4,2] <- 2

corRes <- corHMM(phy, data[,c(1,2)], 1, model = "ER")
corRes$solution

