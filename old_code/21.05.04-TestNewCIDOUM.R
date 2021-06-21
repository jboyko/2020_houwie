source("~/2020_hOUwie/hOUwie.R")
source("~/2020_hOUwie/Utils.R")
# source("/space_2/jamesboyko/2020_hOUwie/hOUwie.R")
# source("/space_2/jamesboyko/2020_hOUwie/Utils.R")

require(OUwie)
require(corHMM)
require(parallel)
require(phytools)

nTip <- 100
phy <- sim.bdtree(b = 1, d = 0.5, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))

fit.CIDOUM.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.CIDOUM.cor <- getFullMat(list(fit.CIDOUM.cor, fit.CIDOUM.cor), fit.CIDOUM.cor)
fit.CIDOUM.cor <- equateStateMatPars(fit.CIDOUM.cor, c(1,2,3))
fit.CIDOUM.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.CIDOUM.cor)[1])
fit.CIDOUM.ou[3,] <- c(3,3,4,4)

pars = c(1, 10, 5, 5, 10)  # mk, alpha, sigma, theta1, theta2

data.houwie <- generateData(phy, fit.CIDOUM.cor, fit.CIDOUM.ou, pars)
cols<-setNames(c("gold","purple","red","black"),
               c("1","2","3","4"))
plotDataSet(data.houwie); legend("bottomleft", legend = c("1A", "2A", "1B", "2B"), pch=16, col = cols)

####### ####### ####### ####### ####### ####### ####### ####### ####### 
# models to fit 
####### ####### ####### ####### ####### ####### ####### ####### ####### 
# CD BMS
fit.CDBMS.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.CDBMS.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(fit.CDBMS.cor)[1])
# CD OUM
fit.CDOUM.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.CDOUM.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.CDOUM.cor)[1])
# fit CID BMS
fit.CIDBMS.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.CIDBMS.cor <- getFullMat(list(fit.CIDBMS.cor, fit.CIDBMS.cor), fit.CIDBMS.cor)
fit.CIDBMS.cor <- equateStateMatPars(fit.CIDBMS.cor, c(1,2,3))
fit.CIDBMS.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(fit.CIDBMS.cor)[1])
fit.CIDBMS.ou[2,] <- c(1,1,2,2)
fit.CIDBMS.ou[3,] <- 3
# fit CID OUM
fit.CIDOUM.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.CIDOUM.cor <- getFullMat(list(fit.CIDOUM.cor, fit.CIDOUM.cor), fit.CIDOUM.cor)
fit.CIDOUM.cor <- equateStateMatPars(fit.CIDOUM.cor, c(1,2,3))
fit.CIDOUM.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.CIDOUM.cor)[1])
fit.CIDOUM.ou[3,] <- c(3,3,4,4)
# fit Hyb BMS
fit.HYBBMS.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.HYBBMS.cor <- getFullMat(list(fit.HYBBMS.cor, fit.HYBBMS.cor), fit.HYBBMS.cor)
fit.HYBBMS.cor <- equateStateMatPars(fit.HYBBMS.cor, c(1,2,3))
fit.HYBBMS.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(fit.HYBBMS.cor)[1])
# fit Hyb OUM
fit.HYBOUM.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.HYBOUM.cor <- getFullMat(list(fit.HYBOUM.cor, fit.HYBOUM.cor), fit.HYBOUM.cor)
fit.HYBOUM.cor <- equateStateMatPars(fit.HYBOUM.cor, c(1,2,3))
fit.HYBOUM.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.HYBOUM.cor)[1])

dat <- data.houwie$data
dat[,2] [dat[,2] == 3] <- 1
dat[,2] [dat[,2] == 4] <- 2
please <- fit.hOUwie.set(phy = phy, data = dat, rate.cats = c(2,2,2), 
                         index.cor.list = list(fit.CIDBMS.cor, fit.CIDOUM.cor, fit.HYBOUM.cor),
                         index.ou.list = list(fit.CIDBMS.ou, fit.CIDOUM.ou, fit.HYBOUM.ou),
                         nSim = 10)

#HYBFit <- hOUwie(phy = phy, data = dat, rate.cat = 2, nSim = 10, index.cor = fit.HYBOUM.cor, index.ou = fit.HYBOUM.ou)
p = c(1, 10, 5, 5, 5, 10, 10)  # mk, alpha, sigma, theta1, theta2
par = c(1, 7, 3.5, 10)

hOUwie(phy = phy, data = dat, rate.cat = 2, nSim = 10, index.cor = fit.CIDBMS.cor, index.ou = fit.CIDBMS.ou, p = par)
hOUwie(phy = phy, data = dat, rate.cat = 2, nSim = 10, index.cor = fit.CIDOUM.cor, index.ou = fit.CIDOUM.ou, p = pars)
hOUwie(phy = phy, data = dat, rate.cat = 2, nSim = 10, index.cor = fit.HYBOUM.cor, index.ou = fit.HYBOUM.ou, p = p)

HYBRecon <- hOUwieRecon(nodes = "external", type = "marginal", nMap = 10, nCores = 1, HYBFit)

#CID BMS
hOUwie(phy = data.houwie$simmap[[1]], data = dat, rate.cat = 2, nSim = 10, index.cor = fit.CIDBMS.cor, index.ou = fit.CIDBMS.ou, p = par)
# CID OUM
hOUwie(phy = data.houwie$simmap[[1]], data = dat, rate.cat = 2, nSim = 10, index.cor = fit.CIDOUM.cor, index.ou = fit.CIDOUM.ou, p = pars)
# CID HYB
hOUwie(phy = data.houwie$simmap[[1]], data = dat, rate.cat = 2, nSim = 10, index.cor = fit.HYBOUM.cor, index.ou = fit.HYBOUM.ou, p = p)

