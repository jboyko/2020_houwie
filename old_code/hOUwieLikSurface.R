# devtools::install_github("bomeara/dentist")
dentFuc1 <- function(par, phy, data){
  p <- c(1.13, 0.53, 3, par)
  out <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = fit.cor, index.ou = fit.ou, p = p)
  return(-out$loglik)
}

source("~/2020_hOUwie/hOUwie.R")
source("~/2020_hOUwie/Utils.R")


require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)
require(dentist)

nTip <- 20
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

data.houwie <- generateData(phy, fit.cor, fit.ou, pars)

cols<-setNames(c("gold","red", "purple", "black"),
               c("1","2","3","4"))
plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
data <- data.houwie$data
data[data[,2]==3,2] <- 1
data[data[,2]==4,2] <- 2

# model fits
BM1 <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BM1")
OU1 <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OU1")
OUM <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM")
BMS <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS")
CID <- hOUwie(phy = phy, data = data, rate.cat = 2, nBins = 10, index.cor = fit.cor, index.ou = fit.ou, ip = pars)

p <- c(1.13, 0.53, 3, 1.48, 18.03)[c(4,5)]
names(p) <- c("theta1", "theta2")
dented_results <- dent_walk(par=p, fn=dentFuc1, best_neglnL=37.05124, delta=4, nsteps=500, print_freq=250, phy=phy, data=data)
print(dented_results)
plot(dented_results)



test3 <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS")
test3$objective
exp(test3$solution)
p <- exp(test3$solution)
hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS", p = p)
hOUwieRecon(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS", p = p)

# sig2_2 at a reasonable value
p <- exp(test3$solution)[c(1,2)]
names(p) <- c("Mk", "sigma1")
dented_results2 <- dent_walk(par=p, fn=dentFuc2, best_neglnL=test3$objective,  delta=2, nsteps=1000, print_freq=250, phy=phy, data=data)
print(dented_results2)
plot(dented_results2)

# Mk at 100
p <- exp(test3$solution)[c(2,3)]
names(p) <- c("sigma1", "sigma2")
dented_results3 <- dent_walk(par=p, fn=dentFuc3, best_neglnL=test3$objective,  delta=2, nsteps=1000, print_freq=250, phy=phy, data=data)
print(dented_results3)
plot(dented_results3)

# Mk at 1
par <- exp(test3$solution)[c(2,3)]
names(par) <- c("sigma1", "sigma2")
p <- c(1, par, 9.477262e+01)
out <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS", p = p, ub.cor = 1000)
dented_results4 <- dent_walk(par=par, fn=dentFuc4, best_neglnL=out$objective,  delta=2, nsteps=1000, print_freq=250, phy=phy, data=data)
print(dented_results4)
plot(dented_results4)

# Mk at 0.1
p <- exp(test3$solution)[c(2,3)]
names(p) <- c("sigma1", "sigma2")
dented_results5 <- dent_walk(par=p, fn=dentFuc5, best_neglnL=test3$objective,  delta=2, nsteps=1000, print_freq=250, phy=phy, data=data)
print(dented_results5)
plot(dented_results5)









source("~/2020_hOUwie/hOUwie.R")
source("~/2020_hOUwie/Utils.R")


require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)
require(dentist)

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


pars = c(0.1, 5, 10, 3, 8)  # mk, alpha, sigma, theta1, theta2

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


# model fits
Test <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM")
Test$objective
exp(Test$solution)
p <- exp(Test$solution)
hOUwieRecon(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "OUM", p = p)

Test3 <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS")
Test3$objective
exp(Test3$solution)
p <- exp(Test3$solution)
hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS", p = p)
hOUwieRecon(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS", p = p)







source("~/2020_hOUwie/hOUwie.R")
source("~/2020_hOUwie/Utils.R")


require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)
require(dentist)

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


pars = c(0.1, 1, 10, 10)  # mk, sigma1, sigma2, theta

data.houwie <- generateData(phy, fit.CIDOUM.cor, fit.CIDBMS.ou, pars)

cols<-setNames(c("gold","purple","red","black"),
               c("1","2","3","4"))
plotDataSet(data.houwie); legend("bottomleft", legend = c("1A", "2A", "1B", "2B"), pch=16, col = cols)
data <- data.houwie$data
data[data[,2]==3,2] <- 1
data[data[,2]==4,2] <- 2
data[,2] <- as.numeric(data[,2])
data$reg <- c("Y", "N")[data$reg]
data <- data.frame(sp = data$sp, X =data$reg, Y = data$reg, x = data$x)


Tast3 <- hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS")
Tast3$objective
exp(Tast3$solution)
p <- exp(Tast3$solution)
hOUwie(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS", p = p)
hOUwieRecon(phy = phy, data = data, rate.cat = 1, nBins = 10, model.cor = "ER", model.ou = "BMS", p = p)
