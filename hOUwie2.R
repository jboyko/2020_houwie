source("~/2020_hOUwie/hOUwie.R")
source("~/2020_hOUwie/Utils.R")


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

# fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
# fit.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.cor)[1])
# pars = c(1, 5, 5, 3, 8)  # mk, alpha, sigma, theta1, theta2

# fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
# fit.ou <- getOUParamStructure("BM1", "three.point", FALSE, FALSE, dim(fit.cor)[1])
# pars = c(1, 1, 10)  # mk, alpha, sigma, theta1, theta2

fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.cor <- getFullMat(list(fit.cor, fit.cor), fit.cor)
fit.cor <- equateStateMatPars(fit.cor, c(1,2,3))
fit.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.cor)[1])
fit.ou[3,] <- c(3,3,4,4)
pars = c(1, 3, 3, 2, 10)  # mk, alpha, sigma, theta1, theta2

# fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
# fit.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(fit.cor)[1])
# pars = c(2, 1, 10, 8)  # mk, alpha, sigma, theta1, theta2

# fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
# fit.ou <- getOUParamStructure("OU1", "three.point", FALSE, FALSE, dim(fit.cor)[1])
# pars = c(2, 1, 1, 8)  # mk, alpha, sigma, theta1, theta2

CID.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
CID.cor <- getFullMat(list(CID.cor, CID.cor), CID.cor)
CID.cor <- equateStateMatPars(CID.cor, c(1,2,3))
CIDOUM.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(CID.cor)[1])
CIDOUM.ou[3,] <- c(3,3,4,4)
CIDBMS.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(CID.cor)[1])
CIDBMS.ou[2,] <- c(1,1,2,2)
CIDBMS.ou[3,] <- c(3,3,3,3)


data.houwie <- generateData(phy, fit.cor, fit.ou, pars)

cols<-setNames(c("gold","red", "purple", "black"),
               c("1","2","3","4"))
plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
data <- data.houwie$data
data[data[,2]==3,2] <- 1
data[data[,2]==4,2] <- 2

dat2 <- data[,3]
names(dat2) <- data[,1]
dat3 <- data[,2]
names(dat3) <- data[,1]
fitContinuous(phy, dat2, model = "OU")
# fitContinuous(phy, dat2, model = "BM")
# fitDiscrete(phy, dat3, "ER")
# fitContinuous(phy, dat2, model = "BM")$opt$lnL
# fitContinuous(phy, dat2, model = "white")$opt$lnL
#undebug(hOUwie.dev)
BM1 <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BM1", split.LnLiks = TRUE)
OU1 <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OU1", split.LnLiks = TRUE)
OUM <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", split.LnLiks = TRUE)
BMS <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS", split.LnLiks = TRUE)
CIDOUM <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 20, index.cor = CID.cor, index.ou = CIDOUM.ou, split.LnLiks = TRUE)
CIDBMS <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 20, index.cor = CID.cor, index.ou = CIDBMS.ou, split.LnLiks = TRUE)

p1 <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = fit.cor, index.ou = fit.ou, p = c(1.13, 0.53, 3, 1.48, 18.03))
# p2 <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = fit.cor, index.ou = fit.ou, p = c(1, 0.01, 3, 2, 10))
p3 <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = fit.cor, index.ou = fit.ou, p = c(1, 3, 3, 10, 2))


test <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM")
p <- exp(test$solution)
test$objective
#hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", p = p)
OUMRECON <- hOUwieRecon(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", p = p)
o <- round(p, 1)
parameters <- c(paste0("Mk=",o[1], ", alp=", o[2], ", sig2=", o[3], ", thta=", o[4]),
                paste0("Mk=",o[1], ", alp=", o[2], ", sig2=", o[3], ", thta=", o[5]))
plotDataSet(data.houwie); legend("bottomleft", legend = parameters, pch=16, col = cols)
tiplabels(pch = 16, col = cols[OUMRECON$TipStates[,1]])
nodelabels(pch = 16, col = cols[OUMRECON$NodeStates[,1]])


test2 <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OU1")
p <- exp(test2$solution)
p <- c(1, 2, 2, 8.787807)
hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OU1", p = p)

BM1RECON <- hOUwieRecon(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BM1", p = p)
o <- round(p, 1)
parameters <- c(paste0("Mk=",o[1], ", sig2=", o[2], ", thta=", o[3]),
                paste0("Mk=",o[1], ", sig2=", o[2], ", thta=", o[3]))
plotDataSet(data.houwie); legend("bottomleft", legend = parameters, pch=16, col = cols)
tiplabels(pch = 16, col = cols[BM1RECON$TipStates[,1]])
nodelabels(pch = 16, col = cols[BM1RECON$NodeStates[,1]])


test4 <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = fit.cor, index.ou = fit.ou)
p <- exp(test4$solution)
test4$objective
#hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", p = p)
CIDRECON <- hOUwieRecon(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = fit.cor, index.ou = fit.ou,p = p)
o <- round(p, 1)
parameters <- c(paste0("Mk=",o[1], ", alp=", o[2], ", sig2=", o[3], ", thta=", o[4]),
                paste0("Mk=",o[1], ", alp=", o[2], ", sig2=", o[3], ", thta=", o[4]),
                paste0("Mk=",o[1], ", alp=", o[2], ", sig2=", o[3], ", thta=", o[5]),
                paste0("Mk=",o[1], ", alp=", o[2], ", sig2=", o[3], ", thta=", o[5]))
plotDataSet(data.houwie); legend("topleft", legend = parameters, pch=16, col = cols)
tiplabels(pch = 16, col = cols[CIDRECON$TipStates[,1]])
nodelabels(pch = 16, col = cols[CIDRECON$NodeStates[,1]])






test3 <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS")
p <- exp(test3$solution)
# p <- c(exp(test2$solution)[c(1,2,2,3)])
# p[1] <- 100
test3$objective
# hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS", p = p)
BMSRECON <- hOUwieRecon(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS", p = p)
o <- round(p, 1)
parameters <- c(paste0("Mk=",o[1], ", sig2=", o[2], ", thta=", o[4]),
                paste0("Mk=",o[1], ", sig2=", o[3], ", thta=", o[4]))
plotDataSet(data.houwie); legend("bottomleft", parameters, pch=16, col = cols)
tiplabels(pch = 16, col = cols[BMSRECON$TipStates[,1]])
nodelabels(pch = 16, col = cols[BMSRECON$NodeStates[,1]])


test5 <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUMV")
#p <- exp(test5$solution)
p <- c(100, 1e-10, exp(test3$solution)[c(2,3,4,4)])
# p[1] <- 100
test5$objective
# hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUMV", p = p)
OUMVRECON <- hOUwieRecon(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUMV", p = p)
o <- round(p, 1)
parameters <- c(paste0("Mk=",o[1], ", sig2=", o[2], ", thta=", o[4]),
                paste0("Mk=",o[1], ", sig2=", o[3], ", thta=", o[4]))
plotDataSet(data.houwie); legend("bottomleft", parameters, pch=16, col = cols)
tiplabels(pch = 16, col = cols[OUMVRECON$TipStates[,1]])
nodelabels(pch = 16, col = cols[OUMVRECON$NodeStates[,1]])



fit.CDOUM.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.CDOUMVA.ou <- getOUParamStructure("OUMVA", "three.point", FALSE, FALSE, dim(fit.CDOUM.cor)[1])
fit.CDBMS.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(fit.CDOUM.cor)[1])
p <- exp(test3$solution)
hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS", p = p)
p2 <- c(p[1], 1e-10, 1e-10, p[2], p[3], p[4], p[4])
hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUMVA", p = p2)







source("~/2020_hOUwie/hOUwie.R")
source("~/2020_hOUwie/Utils.R")


require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)

nTip <- 20
phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))

fit.CIDOUM.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.CIDOUM.cor <- getFullMat(list(fit.CIDOUM.cor, fit.CIDOUM.cor), fit.CIDOUM.cor)
fit.CIDOUM.cor <- equateStateMatPars(fit.CIDOUM.cor, c(1,2,3))
fit.CIDOUM.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.CIDOUM.cor)[1])
fit.CIDOUM.ou[3,] <- c(3,3,4,4)
fit.CIDBMS.ou <- getOUParamStructure("BMS", "three.point", FALSE, FALSE, dim(fit.CIDOUM.cor)[1])
fit.CIDBMS.ou[2,] <- c(1,1,2,2)
fit.CIDBMS.ou[3,] <- c(3,3,3,3)


pars = c(1, 5, 10, 3, 8)  # mk, alpha, sigma, theta1, theta2

data.houwie <- generateData(phy, fit.CIDOUM.cor, fit.CIDOUM.ou, pars)

cols<-setNames(c("gold","purple","red","black"),
               c("1","2","3","4"))
plotDataSet(data.houwie); legend("bottomleft", legend = c("1A", "2A", "1B", "2B"), pch=16, col = cols)
data <- data.houwie$data
data[data[,2]==3,2] <- 1
data[data[,2]==4,2] <- 2
data[,2] <- as.numeric(data[,2])
data$reg <- c("Y", "N")[data$reg]
data <- data.frame(sp = data$sp, X =data$reg, Y = data$reg, x = data$x)

test <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM")
# test$objective
# exp(test$solution)
p <- exp(test$solution)
hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", p = p)
OUMRECON <- hOUwieRecon(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", p = p)
tiplabels(pch = 16, col = cols[OUMRECON$TipStates[,1]])
nodelabels(pch = 16, col = cols[OUMRECON$NodeStates[,1]])

# test2 <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = fit.CIDOUM.cor, index.ou = fit.CIDOUM.ou)
# test2$objective
# exp(test2$solution)
# p <- exp(test2$solution)
# hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = fit.CIDOUM.cor, index.ou = fit.CIDOUM.ou, p = p)

test3 <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS")
test3$objective
exp(test3$solution)
p <- exp(test3$solution)
hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS", p = p)
hOUwieRecon(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS", p = p)

# test4 <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BM1")
# test4$objective
# exp(test4$solution)
# p <- exp(test4$solution)
# hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BM1", p = p)
# 
# test5 <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OU1")
# test5$objective
# exp(test5$solution)
# p <- exp(test5$solution)
# hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OU1", p = p)
# 
# test6 <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = fit.CIDOUM.cor, index.ou = fit.CIDBMS.ou)
# test6$objective
# exp(test6$solution)
# p <- exp(test6$solution)
# hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = fit.CIDOUM.cor, index.ou = fit.CIDBMS.ou, p = p)

# to do: 
# add root priors correctly!
hOUwieRecon.dev(est.pars, phy=phy, organizedData=organizedData, rate.cat=rate.cat, index.cor=index.cor, index.ou=index.ou)
