# proof of concept for the tip importance
source("~/2020_hOUwie/hOUwieSimmap.R")
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

Q <- fit.cor
diag(Q) <- -rowSums(Q)
pars = c(1, 5, 5, 2, 10)
alpha = rep(pars[2], 4)
sigma.sq = rep(pars[3], 4)
theta = c(pars[4], pars[4], pars[5], pars[5])
tip.paths <- lapply(1:length(phy$tip.label), function(x) OUwie:::getPathToRoot(phy, x))
nSim = 100

OU.loglikTrue <- OUwie.basic(data.houwie$simmap[[1]], data.houwie$data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, mserr="none")

simmapA <- corHMM:::makeSimmap(phy, data[c(1,2)], Q, 2, root.p = "yang", nSim = nSim)
OU.loglikA <- lapply(simmapA, function(x) OUwie.basic(x, data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, mserr="none"))
OU.loglikA <- unlist(OU.loglikA)
OU.loglikA <- max(OU.loglikA) + log(sum(exp(OU.loglikA - max(OU.loglikA)))) - log(nSim)

simmapB <- corHMM:::makeSimmap(phy, data.houwie$data[c(1,2)], Q, 1, root.p = "yang", nSim = nSim)
OU.loglikB <- lapply(simmapB, function(x) OUwie.basic(x, data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, mserr="none"))
OU.loglikB <- unlist(OU.loglikB)
OU.loglikB <- max(OU.loglikB) + log(sum(exp(OU.loglikB - max(OU.loglikB)))) - log(nSim)

data.sample <- data[c(1,2)]
TipOUSamples <- cbind("1"= dnorm(data[,3], theta[1], sqrt(sigma.sq[1]/2*alpha[1])),
                      "2" = dnorm(data[,3], theta[2], sqrt(sigma.sq[2]/2*alpha[2])),
                      "3" = dnorm(data[,3], theta[3], sqrt(sigma.sq[3]/2*alpha[3])),
                      "4" = dnorm(data[,3], theta[4], sqrt(sigma.sq[4]/2*alpha[4])))
# TipOUSamples <- cbind("1"= dnorm(data[,3], theta[1], sqrt(sigma.sq[1])),
#                       "2" = dnorm(data[,3], theta[2], sqrt(sigma.sq[2])),
#                       "3" = dnorm(data[,3], theta[3], sqrt(sigma.sq[3])),
#                       "4" = dnorm(data[,3], theta[4], sqrt(sigma.sq[4])))

for(i in 1:length(data[,2])){
  if(data[i,2] == 1){
    TipOUSamples[i,c(2,4)] <- 0
  }
  if(data[i,2] == 2){
    TipOUSamples[i,c(1,3)] <- 0
  }
}
TipOUProbs <- TipOUSamples
TipOUProbs <- TipOUProbs/rowSums(TipOUProbs)
data.sample[,2] <- apply(TipOUSamples, 1, function(x) sample(1:4, 1, prob = x))
simmapC <- corHMM:::makeSimmap(phy, data.sample, Q, 1, root.p = "yang", nSim = nSim)
OU.loglikC <- lapply(simmapC, function(x) OUwie.basic(x, data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, mserr="none"))
OU.loglikC <- unlist(OU.loglikC)
OU.loglikC <- max(OU.loglikC) + log(sum(exp(OU.loglikC - max(OU.loglikC)))) - log(nSim)

c(TrueLoglik = OU.loglikTrue, UnknownTips = OU.loglikA, KnownTips = OU.loglikB, EstimatedTips = OU.loglikC)

log(mean(exp(unlist(lapply(simmapA, function(x) getMapProbability(x$maps, Q))))))
log(mean(exp(unlist(lapply(simmapB, function(x) getMapProbability(x$maps, Q))))))
log(mean(exp(unlist(lapply(simmapC, function(x) getMapProbability(x$maps, Q))))))

# # # # #
# checking whether i can use probs for tip values
# # # # #
#debug(corHMM:::makeSimmap)

simmapD <- corHMM:::makeSimmap(phy, data[c(1,2)], Q, 2, root.p = "yang", nSim = nSim, tip.liks = TipOUProbs)
OU.loglikD <- lapply(simmapD, function(x) OUwie.basic(x, data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, mserr="none"))
OU.loglikD <- unlist(OU.loglikD)
OU.loglikD <- max(OU.loglikD) + log(sum(exp(OU.loglikD - max(OU.loglikD)))) - log(nSim)

c(UnknownTips = OU.loglikA, KnownTips = OU.loglikB, EstimatedTips = OU.loglikC, EstimatedMaps = OU.loglikD)
#
simmapE <- list()
for(i in sequence(nSim)){
  cat("\r", i)
  check <- TRUE
  while(check){
    data.sample[,2] <- apply(TipOUSamples, 1, function(x) sample(1:4, 1, prob = x))
    check <- length(unique(data.sample[,2])) != 4
  }
  simmapE[[i]] <- corHMM:::makeSimmap(phy, data.sample, Q, 1, root.p = "yang", nSim = nSim)[[1]]
}
OU.loglikE <- lapply(simmapE, function(x) OUwie.basic(x, data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, mserr="none"))
OU.loglikE <- unlist(OU.loglikE)
OU.loglikE <- max(OU.loglikE) + log(sum(exp(OU.loglikE - max(OU.loglikE)))) - log(nSim)

# what happens when there are no Mk differences
c(UnknownTips = OU.loglikA, KnownTips = OU.loglikB, EstimatedTips = OU.loglikC, EstimatedMaps = OU.loglikD, SampleTips = OU.loglikE)








# proof of concept for the tip importance
source("~/2020_hOUwie/hOUwieSimmap.R")
source("~/2020_hOUwie/Utils.R")

require(OUwie)
require(corHMM)
require(parallel)
require(phytools)
require(expm)
require(POUMM)
require(geiger)

nTip <- 100
phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTip) 
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))

fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
fit.cor <- getFullMat(list(fit.cor, fit.cor), fit.cor)
# fit.cor <- equateStateMatPars(fit.cor, c(1,2))
fit.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.cor)[1])
fit.ou[3,] <- c(3,3,4,4)
pars = c(0.1, 2, 1, 5, 5, 2, 10)  # mk, alpha, sigma, theta1, theta2
data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
cols<-setNames(c("gold","red", "purple", "black"),
               c("1","2","3","4"))
plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
data <- data.houwie$data
data[data[,2]==3,2] <- 1
data[data[,2]==4,2] <- 2

Q <- fit.cor
diag(Q) <- -rowSums(Q)
pars = c(1, 5, 5, 2, 10)
alpha = rep(pars[2], 4)
sigma.sq = rep(pars[3], 4)
theta = c(pars[4], pars[4], pars[5], pars[5])
tip.paths <- lapply(1:length(phy$tip.label), function(x) OUwie:::getPathToRoot(phy, x))
nSim = 100

simmapA <- corHMM:::makeSimmap(phy, data[c(1,2)], Q, 2, root.p = "yang", nSim = nSim)
OU.loglikA <- lapply(simmapA, function(x) OUwie.basic(x, data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, mserr="none"))
OU.loglikA <- unlist(OU.loglikA)
OU.loglikA <- max(OU.loglikA) + log(sum(exp(OU.loglikA - max(OU.loglikA)))) - log(nSim)

simmapB <- corHMM:::makeSimmap(phy, data.houwie$data[c(1,2)], Q, 1, root.p = "yang", nSim = nSim)
OU.loglikB <- lapply(simmapB, function(x) OUwie.basic(x, data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, mserr="none"))
OU.loglikB <- unlist(OU.loglikB)
OU.loglikB <- max(OU.loglikB) + log(sum(exp(OU.loglikB - max(OU.loglikB)))) - log(nSim)

data.sample <- data[c(1,2)]
TipOUSamples <- cbind("1"= dnorm(data[,3], theta[1], sqrt(sigma.sq[1]/2*alpha[1])),
                      "2" = dnorm(data[,3], theta[2], sqrt(sigma.sq[2]/2*alpha[2])),
                      "3" = dnorm(data[,3], theta[3], sqrt(sigma.sq[3]/2*alpha[3])),
                      "4" = dnorm(data[,3], theta[4], sqrt(sigma.sq[4]/2*alpha[4])))
# TipOUSamples <- cbind("1"= dnorm(data[,3], theta[1], sqrt(sigma.sq[1])),
#                       "2" = dnorm(data[,3], theta[2], sqrt(sigma.sq[2])),
#                       "3" = dnorm(data[,3], theta[3], sqrt(sigma.sq[3])),
#                       "4" = dnorm(data[,3], theta[4], sqrt(sigma.sq[4])))

for(i in 1:length(data[,2])){
  if(data[i,2] == 1){
    TipOUSamples[i,c(2,4)] <- 0 
  }
  if(data[i,2] == 2){
    TipOUSamples[i,c(1,3)] <- 0 
  }
}
TipOUProbs <- TipOUSamples
TipOUProbs <- TipOUProbs/rowSums(TipOUProbs)
data.sample[,2] <- apply(TipOUSamples, 1, function(x) sample(1:4, 1, prob = x))
simmapC <- corHMM:::makeSimmap(phy, data.sample, Q, 1, root.p = "yang", nSim = nSim)
OU.loglikC <- lapply(simmapC, function(x) OUwie.basic(x, data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, mserr="none"))
OU.loglikC <- unlist(OU.loglikC)
OU.loglikC <- max(OU.loglikC) + log(sum(exp(OU.loglikC - max(OU.loglikC)))) - log(nSim)

# c(UnknownTips = OU.loglikA, KnownTips = OU.loglikB, EstimatedTips = OU.loglikC)

# # # # #
# checking whether i can use probs for tip values
# # # # #
#debug(corHMM:::makeSimmap)

simmapD <- corHMM:::makeSimmap(phy, data[c(1,2)], Q, 2, root.p = "yang", nSim = nSim, tip.liks = TipOUProbs)
OU.loglikD <- lapply(simmapD, function(x) OUwie.basic(x, data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, mserr="none"))
OU.loglikD <- unlist(OU.loglikD)
OU.loglikD <- max(OU.loglikD) + log(sum(exp(OU.loglikD - max(OU.loglikD)))) - log(nSim)

# c(UnknownTips = OU.loglikA, KnownTips = OU.loglikB, EstimatedTips = OU.loglikC, EstimatedMaps = OU.loglikD)

simmapE <- list()
for(i in sequence(nSim)){
  cat("\r", i)
  check <- TRUE
  while(check){
    data.sample[,2] <- apply(TipOUSamples, 1, function(x) sample(1:4, 1, prob = x))
    check <- length(unique(data.sample[,2])) != 4
  }
  simmapE[[i]] <- corHMM:::makeSimmap(phy, data.sample, Q, 1, root.p = "yang", nSim = nSim)[[1]]
}
OU.loglikE <- lapply(simmapE, function(x) OUwie.basic(x, data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, mserr="none"))
OU.loglikE <- unlist(OU.loglikE)
OU.loglikE <- max(OU.loglikE) + log(sum(exp(OU.loglikE - max(OU.loglikE)))) - log(nSim)

# what happens when there are Mk differences
c(UnknownTips = OU.loglikA, KnownTips = OU.loglikB, EstimatedTips = OU.loglikC, EstimatedMaps = OU.loglikD, SampleTips = OU.loglikE)
