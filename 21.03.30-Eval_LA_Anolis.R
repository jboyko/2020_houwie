source("~/2020_hOUwie/hOUwie.R")
source("~/2020_hOUwie/Utils.R")
#source("/space_2/jamesboyko/2020_hOUwie/hOUwie.R")
#source("/space_2/jamesboyko/2020_hOUwie/Utils.R")
require(OUwie)
require(corHMM)
require(parallel)

load("~/2020_hOUwie/data/LA_Anolis.Rsave")

phy <- LA_Anolis$tree
phy$node.label <- NULL
# phy$edge.length <- phy$edge.length/max(branching.times(phy))
dat <- LA_Anolis$data
dat[as.numeric(dat[,2]) > 1,2] <- 2
dat[,2] <- as.numeric(dat[,2])

ARD <- getStateMat4Dat(data = dat[,c(1,2)])$rate.mat
ER <- equateStateMatPars(ARD, c(1,2))
a2s <- dropStateMatPars(ARD, 1)

CD.Cor <- getFullMat(list(a2s, a2s), ER)
CD.Cor <- equateStateMatPars(CD.Cor, c(1,2,3)) # not really any case for one evolving faster than another
CD.Cor[2,] <- c(1,0,0,0)
CD.Cor[4,] <- c(0,0,1,0)
CD.OU <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(CD.Cor)[1])
CD.OU[3,] <- c(3,4,3,5)

TC.cor <- getFullMat(list(a2s, ER), ER)
TC.cor[1,] <- c(0,1,0,1)
TC.cor[2,] <- c(0,0,1,0)
TC.cor[3,] <- c(1,0,0,0)
TC.cor[4,] <- c(0,0,0,0)
TC.OU <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(TC.cor)[1])
TC.OU[3,] <- c(3,3,3,4)

index.cor.list <- list(
  CD = CD.Cor,
  TC = TC.cor
)

index.ou.list <- list(
  CD = CD.OU,
  TC = TC.OU
)

rate.cats <- c(2,2)

# debug(fit.hOUwie.set)
test <- fit.hOUwie.set(phy = phy, data = dat, rate.cats = rate.cats, index.cor.list = index.cor.list, index.ou.list = index.ou.list, root.p = "maddfitz", nSim = 50, lb.cor = 0.002, ub.cor = 0.2)

LA_Anolis_hOUwieSetFit <- test

save(LA_Anolis_hOUwieSetFit, file = "~/2020_hOUwie/data/LA_Anolis_hOUwieSetFit.Rsave")

