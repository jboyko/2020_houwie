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

nCores <- 80
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

#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BM1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.ou <- getOUParamStructure("BM1", "three.point", FALSE, FALSE, dim(fit.cor)[1])
pars = c(1, 1, 10)  # mk, alpha, sigma, theta1, theta2
data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
# cols<-setNames(c("gold","red", "purple", "black"),
#                c("1","2","3","4"))
# plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
data <- data.houwie$data
data[data[,2]==3,2] <- 1
data[data[,2]==4,2] <- 2
BM1 <- hOUwie(phy = phy, data = data, rate.cat = 1, nSim = 50, model.cor = "ER", model.ou = "BM1", weighted = TRUE)
OU1 <- hOUwie(phy = phy, data = data, rate.cat = 1, nSim = 2, model.cor = "ER", model.ou = "OU1", weighted = TRUE)
OUM <- hOUwie(phy = phy, data = data, rate.cat = 1, nSim = 2, model.cor = "ER", model.ou = "OUM", weighted = TRUE)
BMS <- hOUwie(phy = phy, data = data, rate.cat = 1, nSim = 2, model.cor = "ER", model.ou = "BMS", weighted = TRUE)
CIDOUM <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = 2, index.cor = CID.cor, index.ou = CID.ou, weighted = TRUE)
CIDBMS <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = 2, index.cor = CID.cor, index.ou = CIDBMS.ou, weighted = TRUE)
out <- list("BM" = BM1, "OU1" = OU1, "CDOUM" = OUM, "CDBMS" = BMS, "CIDOUM" = CIDOUM, "CIDBMS" = CIDBMS)
getModelTable(list("BM" = BM1, "OU1" = OU1, "CDOUM" = OUM, "CDBMS" = BMS, "CIDOUM" = CIDOUM, "CIDBMS" = CIDBMS))

#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# OU1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.ou <- getOUParamStructure("OU1", "three.point", FALSE, FALSE, dim(fit.cor)[1])
pars = c(1, 1, 1, 8)  # mk, alpha, sigma, theta1, theta2
data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
# cols<-setNames(c("gold","red", "purple", "black"),
#                c("1","2","3","4"))
# plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
data <- data.houwie$data
data[data[,2]==3,2] <- 1
data[data[,2]==4,2] <- 2
BM1 <- hOUwie(phy = phy, data = data, rate.cat = 1, nSim = 2, model.cor = "ER", model.ou = "BM1")
OU1 <- hOUwie(phy = phy, data = data, rate.cat = 1, nSim = 2, model.cor = "ER", model.ou = "OU1")
OUM <- hOUwie(phy = phy, data = data, rate.cat = 1, nSim = 2, model.cor = "ER", model.ou = "OUM")
BMS <- hOUwie(phy = phy, data = data, rate.cat = 1, nSim = 2, model.cor = "ER", model.ou = "BMS")
CIDOUM <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = 2, index.cor = CID.cor, index.ou = CID.ou)
CIDBMS <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = 2, index.cor = CID.cor, index.ou = CIDBMS.ou)
out <- list("BM" = BM1, "OU1" = OU1, "CDOUM" = OUM, "CDBMS" = BMS, "CIDOUM" = CIDOUM, "CIDBMS" = CIDBMS)
getModelTable(list("BM" = BM1, "OU1" = OU1, "CDOUM" = OUM, "CDBMS" = BMS, "CIDOUM" = CIDOUM, "CIDBMS" = CIDBMS))

#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CDBMS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
singleRun <- function(iter){
  fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
  fit.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(fit.cor)[1])
  pars = c(1, 1, 10, 10)  # mk, alpha, sigma, theta1, theta2
  nSim = 50
  data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
  # cols<-setNames(c("gold","red", "purple", "black"),
  #                c("1","2","3","4"))
  # plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
  data <- data.houwie$data
  data[data[,2]==3,2] <- 1
  data[data[,2]==4,2] <- 2
  BM1 <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = CID.cor, model.ou = "BM1", weighted = TRUE)
  OU1 <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = CID.cor, model.ou = "OU1", weighted = TRUE)
  OUM <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = CID.cor, index.ou = CD.ou, weighted = TRUE)
  BMS <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = CID.cor, index.ou = CDBMS.ou, weighted = TRUE)
  CIDOUM <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                   index.cor = CID.cor, index.ou = CID.ou, weighted = TRUE)
  CIDBMS <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                   index.cor = CID.cor, index.ou = CIDBMS.ou, weighted = TRUE)
  out <- list("BM" = BM1, "OU1" = OU1, "CDOUM" = OUM, "CDBMS" = BMS, "CIDOUM" = CIDOUM, "CIDBMS" = CIDBMS)
  save(out, file = paste0("hOUwie-CDBMS-", iter, ".Rsave"))
  return(out)
}
res <- mclapply(1:100, function(x) singleRun(x), mc.cores = nCores)

getModelTable(list("BM" = BM1, "OU1" = OU1, "CDOUM" = OUM, "CDBMS" = BMS, "CIDOUM" = CIDOUM, "CIDBMS" = CIDBMS))

#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CDOUM
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
singleRun <- function(iter){
  fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
  fit.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.cor)[1])
  pars = c(1, 5, 5, 2, 10)  # mk, alpha, sigma, theta1, theta2
  data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
  nSim = 50
  # cols<-setNames(c("gold","red", "purple", "black"),
  #                c("1","2","3","4"))
  # plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
  data <- data.houwie$data
  data[data[,2]==3,2] <- 1
  data[data[,2]==4,2] <- 2
  BM1 <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = CID.cor, model.ou = "BM1", weighted = TRUE)
  OU1 <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = CID.cor, model.ou = "OU1", weighted = TRUE)
  OUM <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = CID.cor, index.ou = CD.ou, weighted = TRUE)
  BMS <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = CID.cor, index.ou = CDBMS.ou, weighted = TRUE)
  CIDOUM <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                   index.cor = CID.cor, index.ou = CID.ou, weighted = TRUE)
  CIDBMS <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                   index.cor = CID.cor, index.ou = CIDBMS.ou, weighted = TRUE)
  out <- list("BM" = BM1, "OU1" = OU1, "CDOUM" = OUM, "CDBMS" = BMS, "CIDOUM" = CIDOUM, "CIDBMS" = CIDBMS)
  save(out, file = paste0("hOUwie-CDOUM-", iter, ".Rsave"))
  return(out)
}
res <- mclapply(1:100, function(x) singleRun(x), mc.cores = nCores)
# getModelTable(list("BM" = BM1, "OU1" = OU1, "CDOUM" = OUM, "CDBMS" = BMS, "CIDOUM" = CIDOUM, "CIDBMS" = CIDBMS))

#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CIDOUM
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
singleRun <- function(iter){
  fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
  fit.cor <- getFullMat(list(fit.cor, fit.cor), fit.cor)
  fit.cor <- equateStateMatPars(fit.cor, c(1,2,3))
  fit.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.cor)[1])
  fit.ou[3,] <- c(3,3,4,4)
  pars = c(1, 5, 5, 2, 10)  # mk, alpha, sigma, theta1, theta2
  data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
  nSim = 50
  # cols<-setNames(c("gold","red", "purple", "black"),
  #                c("1","2","3","4"))
  # plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
  data <- data.houwie$data
  data[data[,2]==3,2] <- 1
  data[data[,2]==4,2] <- 2
  BM1 <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = CID.cor, model.ou = "BM1", weighted = TRUE)
  OU1 <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = CID.cor, model.ou = "OU1", weighted = TRUE)
  OUM <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = CID.cor, index.ou = CD.ou, weighted = TRUE)
  BMS <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = CID.cor, index.ou = CDBMS.ou, weighted = TRUE)
  CIDOUM <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                   index.cor = CID.cor, index.ou = CID.ou, weighted = TRUE)
  CIDBMS <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                   index.cor = CID.cor, index.ou = CIDBMS.ou, weighted = TRUE)
  out <- list("BM" = BM1, "OU1" = OU1, "CDOUM" = OUM, "CDBMS" = BMS, "CIDOUM" = CIDOUM, "CIDBMS" = CIDBMS)
  # getModelTable(out)
  save(out, file = paste0("hOUwie-CIDOUM-", iter, ".Rsave"))
  return(out)
}
res <- mclapply(1:100, function(x) singleRun(x), mc.cores = nCores)


#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CIDBMS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
singleRun <- function(iter){
  fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
  fit.cor <- getFullMat(list(fit.cor, fit.cor), fit.cor)
  fit.cor <- equateStateMatPars(fit.cor, c(1,2,3))
  fit.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(fit.cor)[1])
  fit.ou[2,] <- c(1,1,2,2)
  fit.ou[3,] <- c(3,3,3,3)
  pars = c(1, 1, 10, 10)  # mk, alpha, sigma, theta1, theta2
  data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
  nSim = 50
  # cols<-setNames(c("gold","red", "purple", "black"),
  #                c("1","2","3","4"))
  # plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
  data <- data.houwie$data
  data[data[,2]==3,2] <- 1
  data[data[,2]==4,2] <- 2
  BM1 <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = CID.cor, model.ou = "BM1", weighted = TRUE)
  OU1 <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = CID.cor, model.ou = "OU1", weighted = TRUE)
  OUM <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = CID.cor, index.ou = CD.ou, weighted = TRUE)
  BMS <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = CID.cor, index.ou = CDBMS.ou, weighted = TRUE)
  CIDOUM <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                   index.cor = CID.cor, index.ou = CID.ou, weighted = TRUE)
  CIDBMS <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                   index.cor = CID.cor, index.ou = CIDBMS.ou, weighted = TRUE)
  out <- list("BM" = BM1, "OU1" = OU1, "CDOUM" = OUM, "CDBMS" = BMS, "CIDOUM" = CIDOUM, "CIDBMS" = CIDBMS)
  # getModelTable(out)
  save(out, file = paste0("hOUwie-CIDBMS-", iter, ".Rsave"))
  return(out)
}
res <- mclapply(1:100, function(x) singleRun(x), mc.cores = nCores)
getModelTable(list("BM" = BM1, "OU1" = OU1, "CDOUM" = OUM, "CDBMS" = BMS, "CIDOUM" = CIDOUM, "CIDBMS" = CIDBMS))

#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# HYBOUM
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
singleRun <- function(iter){
  fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
  fit.cor <- getFullMat(list(fit.cor, fit.cor), fit.cor)
  fit.cor <- equateStateMatPars(fit.cor, c(1,2,3))
  fit.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.cor)[1])
  pars = c(2, 5, 2.5, 2, 10, 10, 18)  # mk, alpha, sigma, theta1, theta2
  data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
  nSim = 2
  # cols<-setNames(c("gold","red", "purple", "black"),
  #                c("1","2","3","4"))
  # plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
  data <- data.houwie$data
  data[data[,2]==3,2] <- 1
  data[data[,2]==4,2] <- 2
  BM1 <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = CID.cor, model.ou = "BM1", weighted = TRUE)
  OU1 <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = CID.cor, model.ou = "OU1", weighted = TRUE)
  OUM <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = CID.cor, index.ou = CD.ou, weighted = TRUE, p = c(0.7148757,  0.1462326, 19.5353580,  6.2355249, 14.7275790))
  BMS <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                index.cor = CID.cor, index.ou = CDBMS.ou, weighted = TRUE)
  CIDOUM <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                   index.cor = CID.cor, index.ou = CID.ou, weighted = TRUE)
  CIDBMS <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                   index.cor = CID.cor, index.ou = CIDBMS.ou, weighted = TRUE)
  HYBOUM <- hOUwie(phy = phy, data = data, rate.cat = 2, nSim = nSim, 
                   index.cor = CID.cor, index.ou = fit.ou, weighted = TRUE)
  out <- list("BM" = BM1, "OU1" = OU1, "CDOUM" = OUM, "CDBMS" = BMS, "CIDOUM" = CIDOUM, "CIDBMS" = CIDBMS, "HYBOUM" = HYBOUM)
  # getModelTable(out)
  save(out, file = paste0("hOUwie-CIDBMS-", iter, ".Rsave"))
  return(out)
}

res <- mclapply(1:100, function(x) singleRun(x), mc.cores = nCores)
getModelTable(list("BM" = BM1, "OU1" = OU1, "CDOUM" = OUM, "CDBMS" = BMS, "CIDOUM" = CIDOUM, "CIDBMS" = CIDBMS))
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# SUMMARY
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
require(OUwie)
require(corHMM)
require(geiger)
# source("~/2020_hOUwie/hOUwieSimmap.R")
# source("~/2020_hOUwie/Utils.R")
source("/space_2/jamesboyko/2020_hOUwie/hOUwieSimmap.R")
source("/space_2/jamesboyko/2020_hOUwie/Utils.R")
files <- dir("HYBOUM/", full.names=TRUE)
# files <- files[-grep("old", files)]
BestAICc <- c()
OUParsB <- OUParsA <- c()
TruePars <- matrix(c(5, 5, 5, 5, 5, 5, 5, 5, 2, 2, 10, 10), 3, 4, byrow = TRUE)
for(i in sequence(length(files))){
  load(files[i])
  BestAICc <- rbind(BestAICc, unlist(lapply(out, function(x) x$AICc)))
  # OUParsA <- rbind(OUParsA, c(rowMeans((out$CIDOUM$solution.ou - TruePars))))
  # OUParsB <- rbind(OUParsB, c(rowMeans(out$CDOUM$solution.ou - TruePars)^2))
}
rowSums(apply(BestAICc, 1, function(x) x == min(x)))

colMeans(OUParsA)
colMeans(OUParsB)





