source("~/2020_hOUwie/hOUwieEM.R")
source("~/2020_hOUwie/Utils.R")
# source("hOUwieEM.R")
# source("Utils.R")

require(OUwie)
require(corHMM)
data(tworegime)
phy <- tree
phy$node.label <- NULL
# data <- trait
discrete_model <- "ARD"
continuous_model <- "OUMV"
rate.cat <- 1
shift.point <- 0.5
mserr <- "none"
dual = FALSE
collapse <- TRUE
root.station <- FALSE
get.root.theta <- FALSE
lb.disc = NULL
ub.disc = NULL
lb.cont = NULL
ub.cont = NULL
opts = NULL
# 
# CID Models
CID.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
CID.cor <- getFullMat(list(CID.cor, CID.cor), CID.cor)
CID.cor <- equateStateMatPars(CID.cor, c(1,2,3))

CDOUM.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(CID.cor)[1])
CDOUM.ou[3,] <- c(3,4,3,4)
CIDOUM.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(CID.cor)[1])
CIDOUM.ou[3,] <- c(3,3,4,4)

CDBMS.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(CID.cor)[1])
CDBMS.ou[2,] <- c(1,2,1,2)
CDBMS.ou[3,] <- c(3,3,3,3)
CIDBMS.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(CID.cor)[1])
CIDBMS.ou[2,] <- c(1,1,2,2)
CIDBMS.ou[3,] <- c(3,3,3,3)

#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BMS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(fit.cor)[1])
pars = c(0.5, 1, 10, 8)  # mk, alpha, sigma, theta1, theta2
data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
cols<-setNames(c("gold","red", "purple", "black"),
               c("1","2","3","4"))
plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)

fitBM1 <- hOUwie.EM(phy, data.houwie$data, "ER", "BM1", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 100)
fitOU1 <- hOUwie.EM(phy, data.houwie$data, "ER", "OU1", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 100)
fitBMS <- hOUwie.EM(phy, data.houwie$data, "ER", "BMS", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 100)
fitOUM <- hOUwie.EM(phy, data.houwie$data, "ER", "OUM", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 100)

#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# OUM
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.cor)[1])
pars = c(0.5, 2, 2, 3, 8)  # mk, alpha, sigma, theta1, theta2
data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
cols<-setNames(c("gold","red", "purple", "black"),
               c("1","2","3","4"))
plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)

fitBM1 <- hOUwie.EM(phy, data.houwie$data, "ER", "BM1", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 100)
fitOU1 <- hOUwie.EM(phy, data.houwie$data, "ER", "OU1", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 100)
fitBMS <- hOUwie.EM(phy, data.houwie$data, "ER", "BMS", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 100)
fitOUM <- hOUwie.EM(phy, data.houwie$data, "ER", "OUM", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 100)


require(parallel)
singleRun <- function(iter){
  fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
  fit.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.cor)[1])
  pars = c(0.5, 2, 2, 3, 8)  # mk, alpha, sigma, theta1, theta2
  data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
  data <- data.houwie$data
  data[data[,2]==3,2] <- 1
  data[data[,2]==4,2] <- 2
  fitBM1 <- hOUwie.EM(phy, data, "ER", "BM1", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 1000)
  fitOU1 <- hOUwie.EM(phy, data, "ER", "OU1", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 1000)
  fitBMS <- hOUwie.EM(phy, data, "ER", "BMS", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 1000)
  fitOUM <- hOUwie.EM(phy, data, "ER", "OUM", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 1000)
  fitCIDOUM <- hOUwie.EM(phy, data, CID.cor, CIDOUM.ou, 2, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 1000)
  fitCIDBMS <- hOUwie.EM(phy, data, CID.cor, CIDBMS.ou, 2, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 1000)
  out <- list(fitBM1, fitOU1, fitBMS=fitBMS, fitOUM=fitOUM, fitCIDBMS=fitCIDBMS, fitCIDOUM=fitCIDOUM)
  save(out, file = paste0("/space_2/jamesboyko/2020_hOUwie/ModelTests/hOUwieEM-CDOUM-", iter, ".Rsave"))
  return(out)
}

mclapply(1:200, function(x) singleRun(x), mc.cores = 40)


#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CIDOUM
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.cor <- getFullMat(list(fit.cor, fit.cor), fit.cor)
fit.cor <- equateStateMatPars(fit.cor, c(1,2,3))
fit.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.cor)[1])
fit.ou[3,] <- c(3,3,4,4)
pars = c(0.5, 5, 5, 2, 10)  # mk, alpha, sigma, theta1, theta2
data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
# cols<-setNames(c("gold","red", "purple", "black"),
#                c("1","2","3","4"))
# plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
data <- data.houwie$data
data[data[,2]==3,2] <- 1
data[data[,2]==4,2] <- 2
fitBM1 <- hOUwie.EM(phy, data, "ER", "BM1", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 500)
fitOU1 <- hOUwie.EM(phy, data, "ER", "OU1", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 500)
fitBMS <- hOUwie.EM(phy, data, "ER", "BMS", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 500)
fitOUM <- hOUwie.EM(phy, data, "ER", "OUM", rate.cat, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 500)
fitCIDOUM <- hOUwie.EM(phy, data, CID.cor, CIDOUM.ou, 2, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 500)
fitCIDBMS <- hOUwie.EM(phy, data, CID.cor, CIDBMS.ou, 2, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 500)

require(parallel)
singleRun <- function(iter){
  fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
  fit.cor <- getFullMat(list(fit.cor, fit.cor), fit.cor)
  fit.cor <- equateStateMatPars(fit.cor, c(1,2))
  fit.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.cor)[1])
  fit.ou[3,] <- c(3,3,4,4)
  pars = c(1, 0.5, 5, 5, 2, 10)  # mk, alpha, sigma, theta1, theta2
  data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
  # cols<-setNames(c("gold","red", "purple", "black"),
  #                c("1","2","3","4"))
  # plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
  data <- data.houwie$data
  data[data[,2]==3,2] <- 1
  data[data[,2]==4,2] <- 2
  fitBM1 <- hOUwie.EM(phy, data, fit.cor, "BM1", 2, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 1000)
  fitOU1 <- hOUwie.EM(phy, data, fit.cor, "OU1", 2, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 1000)
  fitBMS <- hOUwie.EM(phy, data, fit.cor, CDBMS.ou, 2, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 1000)
  fitOUM <- hOUwie.EM(phy, data, fit.cor, CDOUM.ou, 2, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 1000)
  fitCIDOUM <- hOUwie.EM(phy, data, fit.cor, CIDOUM.ou, 2, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 1000)
  fitCIDBMS <- hOUwie.EM(phy, data, fit.cor, CIDBMS.ou, 2, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 1000)
  out <- list(fitBM1=fitBM1, fitOU1=fitOU1, fitBMS=fitBMS, fitOUM=fitOUM, fitCIDBMS=fitCIDBMS, fitCIDOUM=fitCIDOUM)
  save(out, file = paste0("/space_2/jamesboyko/2020_hOUwie/ModelTests/hOUwieEM-CIDOUM-", iter, ".Rsave"))
  return(out)
}

do.call(rbind, lapply(out, function(x) x[[1]][1:3]))
test <- do.call(rbind, lapply(out, function(x) x[[1]][1:3]))
colnames(test) <- c("Total", "Discrete", "Contin")

nodelabels(pch = 16, col = cols[fitCIDOUM[[2]]])
nodelabels(pch = 16, col = cols[fitCIDBMS[[2]]])

mclapply(1:200, function(x) singleRun(x), mc.cores = 40)


BestLnLik <- c()
for(i in 1:length(dir())){
  print(i)
  load(dir()[i])
  BestLnLik <- c(which.max(unlist(lapply(out, function(x) x[[1]][1]))), BestLnLik)
}

MkLnLik <- c()
for(i in 1:length(dir())){
  print(i)
  load(dir()[i])
  MkLnLik <- c(which.max(unlist(lapply(out, function(x) x[[1]][2]))), MkLnLik)
}

OULnLik <- c()
for(i in 1:length(dir())){
  print(i)
  load(dir()[i])
  OULnLik <- c(which.max(unlist(lapply(out, function(x) x[[1]][3]))), OULnLik)
}

data.frame(
  model = rle(sort(BestLnLik))$values,
  best.no = rle(sort(BestLnLik))$lengths
)

data.frame(
  model = rle(sort(MkLnLik))$values,
  best.MK.no = rle(sort(MkLnLik))$lengths
)

data.frame(
  model = rle(sort(OULnLik))$values,
  best.OU.no = rle(sort(OULnLik))$lengths
)

rle(sort(OULnLik))

do.call(rbind, lapply(out, function(x) x[[1]][1:3]))


#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CIDBMS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
require(parallel)
singleRun <- function(iter){
  fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
  fit.cor <- getFullMat(list(fit.cor, fit.cor), fit.cor)
  fit.cor <- equateStateMatPars(fit.cor, c(1,2,3))
  fit.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(fit.cor)[1])
  fit.ou[2,] <- c(1,1,2,2)
  fit.ou[3,] <- c(3,3,3,3)
  pars = c(1, 1, 10, 10)  # mk, alpha, sigma, theta1, theta2
  data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
  data <- data.houwie$data
  data[data[,2]==3,2] <- 1
  data[data[,2]==4,2] <- 2
  fitBM1 <- hOUwie.EM(phy, data, fit.cor, "BM1", 2, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 1000)
  fitOU1 <- hOUwie.EM(phy, data, fit.cor, "OU1", 2, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 1000)
  fitBMS <- hOUwie.EM(phy, data, fit.cor, CDBMS.ou, 2, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 1000)
  fitOUM <- hOUwie.EM(phy, data, fit.cor, CDOUM.ou, 2, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 1000)
  fitCIDOUM <- hOUwie.EM(phy, data, fit.cor, CIDOUM.ou, 2, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 1000)
  fitCIDBMS <- hOUwie.EM(phy, data, fit.cor, CIDBMS.ou, 2, shift.point, mserr, dual, collapse, root.station, get.root.theta, lb.disc, ub.disc, lb.cont, ub.cont, opts, 1000)
  out <- list(fitBM1=fitBM1, fitOU1=fitOU1, fitBMS=fitBMS, fitOUM=fitOUM, fitCIDBMS=fitCIDBMS, fitCIDOUM=fitCIDOUM)
  save(out, file = paste0("/space_2/jamesboyko/2020_hOUwie/ModelTests/hOUwieEM-CIDBMS-", iter, ".Rsave"))
  return(out)
}

do.call(rbind, lapply(out, function(x) x[[1]][1:3]))
test <- do.call(rbind, lapply(out, function(x) x[[1]][1:3]))
colnames(test) <- c("Total", "Discrete", "Contin")

nodelabels(pch = 16, col = cols[fitCIDOUM[[2]]])
nodelabels(pch = 16, col = cols[fitCIDBMS[[2]]])

mclapply(1:200, function(x) singleRun(x), mc.cores = 40)


BestLnLik <- c()
for(i in 1:length(dir())){
  print(i)
  load(dir()[i])
  BestLnLik <- c(which.max(unlist(lapply(out, function(x) x[[1]][1]))), BestLnLik)
}

MkLnLik <- c()
for(i in 1:length(dir())){
  print(i)
  load(dir()[i])
  MkLnLik <- c(which.max(unlist(lapply(out, function(x) x[[1]][2]))), MkLnLik)
}

OULnLik <- c()
for(i in 1:length(dir())){
  print(i)
  load(dir()[i])
  OULnLik <- c(which.max(unlist(lapply(out, function(x) x[[1]][3]))), OULnLik)
}

data.frame(
  model = rle(sort(BestLnLik))$values,
  best.no = rle(sort(BestLnLik))$lengths
)

data.frame(
  model = rle(sort(MkLnLik))$values,
  best.MK.no = rle(sort(MkLnLik))$lengths
)

data.frame(
  model = rle(sort(OULnLik))$values,
  best.OU.no = rle(sort(OULnLik))$lengths
)

rle(sort(OULnLik))

do.call(rbind, lapply(out, function(x) x[[1]][1:3]))
