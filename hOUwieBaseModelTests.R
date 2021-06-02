source("~/2020_hOUwie/hOUwie.R")
source("~/2020_hOUwie/Utils.R")

# source("hOUwie.R")
# source("Utils.R")

require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)
require(geiger)

nTip <- 50
phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))

CID.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
CID.cor <- getFullMat(list(CID.cor, CID.cor), CID.cor)
CID.cor <- equateStateMatPars(CID.cor, c(1,2,3))
CID.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(CID.cor)[1])
CID.ou[3,] <- c(3,3,4,4)

CIDBMS.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(CID.cor)[1])
CIDBMS.ou[2,] <- c(1,1,2,2)
CIDBMS.ou[3,] <- c(3,3,3,3)

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
dataa <- data.houwie
BM1a <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BM1")
OU1a <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OU1")
OUMa <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", type = "pupko")
BMSa <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS")
CIDOUMa <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = CID.cor, index.ou = CID.ou)
CIDBMSa <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = CID.cor, index.ou = CIDBMS.ou)
# BM1 SIM
getModelTable(list("BM" = BM1a, "OU1" = OU1a, "OUM" = OUMa, "BMS" = BMSa, "CIDOUM" = CIDOUMa, "CIDBMS" = CIDBMSa))

tmp <- getModelTable(list("BM" = BM1a, "OU1" = OU1a, "OUM" = OUMa, "BMS" = BMSa))

BM1Split <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BM1", split.LnLiks = TRUE, p  = c(1.027641, 1.173652, 9.967964))
OU1Split <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OU1", split.LnLiks = TRUE, p  = c(1.027458, 0.006931472, 1.176605111, 9.968395922))
OUMSplit <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", split.LnLiks = TRUE, p  = c(1.027383, 0.01159115, 1.14646765, 9.24371352, 119.68763220))
BMSSplit <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS", split.LnLiks = TRUE, p  = c(1.027757, 0.2295847, 1.545263, 10.3921793))
cbind(tmp, rbind(BM1Split, OU1Split, OUMSplit, BMSSplit))

#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# OU1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.ou <- getOUParamStructure("OU1", "three.point", FALSE, FALSE, dim(fit.cor)[1])
pars = c(2, 1, 1, 8)  # mk, alpha, sigma, theta1, theta2
data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
# cols<-setNames(c("gold","red", "purple", "black"),
#                c("1","2","3","4"))
# plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
data <- data.houwie$data
data[data[,2]==3,2] <- 1
data[data[,2]==4,2] <- 2
datab <- data.houwie
BM1b <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BM1")
OU1b <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OU1")
OUMb <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM")
BMSb <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS")
CIDOUMb <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = CID.cor, index.ou = CID.ou)
CIDBMSb <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = CID.cor, index.ou = CIDBMS.ou)
# OU1 SIM
getModelTable(list("BM" = BM1b, "OU1" = OU1b, "OUM" = OUMb, "BMS" = BMSb, "CIDOUM" = CIDOUMb, "CIDBMS" = CIDBMSb))

#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BMS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(fit.cor)[1])
pars = c(2, 1, 10, 8)  # mk, alpha, sigma, theta1, theta2
data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
# cols<-setNames(c("gold","red", "purple", "black"),
#                c("1","2","3","4"))
# plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
data <- data.houwie$data
data[data[,2]==3,2] <- 1
data[data[,2]==4,2] <- 2
datac <- data.houwie
BM1c <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BM1")
OU1c <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OU1")
OUMc <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM")
BMSc <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS")
CIDOUMc <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = CID.cor, index.ou = CID.ou)
CIDBMSc <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = CID.cor, index.ou = CIDBMS.ou)
# BMS SIM
getModelTable(list("BM" = BM1c, "OU1" = OU1c, "OUM" = OUMc, "BMS" = BMSc, "CIDOUM" = CIDOUMc, "CIDBMS" = CIDBMSc))

#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# OUM
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.cor)[1])
pars = c(0.5, 2, 2, 3, 8)  # mk, alpha, sigma, theta1, theta2
data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
# cols<-setNames(c("gold","red", "purple", "black"),
#                c("1","2","3","4"))
# plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
data <- data.houwie$data
data[data[,2]==3,2] <- 1
data[data[,2]==4,2] <- 2
datad <- data.houwie
BM1d <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BM1", all.roots = TRUE)
OU1d <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OU1", all.roots = TRUE)
OUMd <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", type = "joint", all.roots = TRUE)
BMSd <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS", type = "joint", all.roots = TRUE)
CIDOUMd <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = CID.cor, index.ou = CID.ou, type = "joint", all.roots = TRUE)
CIDBMSd <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = CID.cor, index.ou = CIDBMS.ou, type = "joint", all.roots = TRUE)
# OUM SIM
getModelTable(list("BM" = BM1d, "OU1" = OU1d, "OUM" = OUMd, "BMS" = BMSd, "CIDOUM" = CIDOUMd, "CIDBMS" = CIDBMSd))

# out <- mclapply(OPTS, function(x) hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = CID.cor, index.ou = CID.ou, opts = x))

# AR <- hOUwieRecon(CIDOUMd)
# BR <- hOUwieRecon(OUMd)
debug(hOUwieRecon.dev)
pA = c(0.307842, 1.907731, 2.320655, 0.655545, 6.618652)
OUMPupA <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", type = "pupko", p = pA)
OUMJntA <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", type = "joint", p = pA)
pB = c(13.53353, 0.1327821, 1.441436, 0.2521498, 8.863342)
OUMPupB <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", type = "pupko", p = pB)
OUMJntB <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", type = "joint", p = pB)

print(pA)
print(c("PureJoint" = OUMPupA$loglik, "ReconMap" = OUMJntA$loglik))
print(pB)
print(c("PureJoint" = OUMPupB$loglik, "ReconMap" = OUMJntB$loglik))


CIDBMSd <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = CID.cor, index.ou = CIDBMS.ou, type = "pupko", p = c(4.050993, 0.0006931472, 8.679413, 9.702098))

# LnLik: -57.44525 pars: 0.4910104 4.782075 5.017083 3.103897 8.167658 TRUE
# LnLik: -57.39112 pars: 0.4910325 4.782075 5.015803 3.103897 8.167658 TRUE
# LnLik: -57.7103 pars: 0.4683347 5.12695 5.260873 3.067666 8.091859 FALSE
# LnLik: -77.10152 pars: 0.4464151 5.12695 5.260873 3.060585 8.091859 FALSE

A <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", p = c(0.4910104, 4.782075, 5.017083, 3.103897, 8.167658), all.roots = TRUE)
B <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", p = c(0.4910104, 4.782075, 5.017083, 3.103897, 8.167658), all.roots = FALSE)

# LnLik: -57.7103 pars: 0.4683347 5.12695 5.260873 3.067666 8.091859 
C <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", p = c(0.4683347, 5.12695, 5.260873, 3.067666, 8.091859), all.roots = TRUE, split.LnLiks = FALSE)
D <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", p = c(0.4683347, 5.12695, 5.260873, 3.067666, 8.091859), all.roots = FALSE, split.LnLiks = FALSE)

# LnLik: -77.10152 pars: 0.4464151 5.12695 5.260873 3.060585 8.091859 
E <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", p = c(0.4464151, 5.12695, 5.260873, 3.060585, 8.091859), all.roots = TRUE, split.LnLiks = FALSE)
G <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", p = c(0.4464151, 5.12695, 5.260873, 3.060585, 8.091859), all.roots = FALSE, split.LnLiks = FALSE)

# 
H <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", p = c(0.4683347, 5.12695, 5.260873, 3.060585, 8.091859), all.roots = FALSE, split.LnLiks = FALSE)
I <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", p = c(0.4464151, 5.12695, 5.260873, 3.060585, 8.091859), all.roots = FALSE, split.LnLiks = FALSE)


AR <- hOUwieRecon(D)
BR <- hOUwieRecon(G)
par(mfrow=c(1,2))
cols<-setNames(c("gold","red", "purple", "black"),
               c("1","2","3","4"))
plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
nodelabels(pch = 16, col = cols[AR$NodeStates[,1]])
tiplabels(pch = 16, col = cols[AR$TipStates[,1]])
cols<-setNames(c("gold","red", "purple", "black"),
               c("1","2","3","4"))
plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
nodelabels(pch = 16, col = cols[BR$NodeStates[,1]])
tiplabels(pch = 16, col = cols[BR$TipStates[,1]])


# A <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = CID.cor, index.ou = CID.ou, 
#        p = c(.89, 0.01, 4.09, 3, 6.53))
# AR <- hOUwieRecon(A)
# B <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = CID.cor, index.ou = CID.ou, 
#        p = c(.89, 0.01, 4.09, 4, 6.53))
# BR <- hOUwieRecon(B)
# C <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = CID.cor, index.ou = CID.ou, 
#        p = c(.89, 1, 4.09, 3, 6.53))
# D <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = CID.cor, index.ou = CID.ou, 
#        p = c(.89, 1, 4.09, 4, 6.53))
# par(mfrow=c(1,2))
# cols<-setNames(c("gold","red", "purple", "black"),
#                c("1","2","3","4"))
# plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
# nodelabels(pch = 16, col = cols[AR$NodeStates[,1]])
# tiplabels(pch = 16, col = cols[AR$TipStates[,1]])
# cols<-setNames(c("gold","red", "purple", "black"),
#                c("1","2","3","4"))
# plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
# nodelabels(pch = 16, col = cols[BR$NodeStates[,1]])
# tiplabels(pch = 16, col = cols[BR$TipStates[,1]])

#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CIDOUM
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.cor <- getFullMat(list(fit.cor, fit.cor), fit.cor)
fit.cor <- equateStateMatPars(fit.cor, c(1,2,3))
fit.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.cor)[1])
fit.ou[3,] <- c(3,3,4,4)
pars = c(1, 5, 5, 2, 10)  # mk, alpha, sigma, theta1, theta2
data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
cols<-setNames(c("gold","red", "purple", "black"),
               c("1","2","3","4"))
plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
data <- data.houwie$data
data[data[,2]==3,2] <- 1
data[data[,2]==4,2] <- 2
datae <- data.houwie
BM1e <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BM1", sample.recons = TRUE)
OU1e <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OU1", sample.recons = TRUE)
OUMe <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", sample.recons = TRUE)
BMSe <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS", sample.recons = TRUE)
CIDOUMe <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = CID.cor, index.ou = CID.ou, sample.recons = TRUE) #LnLik: -74.60732 Mk: -13.76855 OU: -60.83876 pars: 0.2638472 1.780456 5.694926 0.6037362 8.576074 
CIDBMSe <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = CID.cor, index.ou = CIDBMS.ou, sample.recons = TRUE)
# CID SIM
getModelTable(list("BM" = BM1e, "OU1" = OU1e, "OUM" = OUMe, "BMS" = BMSe, "CIDOUM" = CIDOUMe, "CIDBMS" = CIDBMSe))

#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CIDBMS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.cor <- getFullMat(list(fit.cor, fit.cor), fit.cor)
fit.cor <- equateStateMatPars(fit.cor, c(1,2,3))
fit.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(fit.cor)[1])
fit.ou[2,] <- c(1,1,2,2)
fit.ou[3,] <- c(3,3,3,3)
pars = c(1, 1, 10, 10)  # mk, alpha, sigma, theta1, theta2
data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
# cols<-setNames(c("gold","red", "purple", "black"),
#                c("1","2","3","4"))
# plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
data <- data.houwie$data
data[data[,2]==3,2] <- 1
data[data[,2]==4,2] <- 2
data <- data.houwie
BM1e <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BM1")
OU1e <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OU1")
OUMe <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM")
BMSe <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS")
CIDOUMe <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = CID.cor, index.ou = CID.ou)
CIDBMSe <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = CID.cor, index.ou = CIDBMS.ou)
# CID SIM
getModelTable(list("BM" = BM1e, "OU1" = OU1e, "OUM" = OUMe, "BMS" = BMSe, "CIDOUM" = CIDOUMe, "CIDBMS" = CIDBMSe))
