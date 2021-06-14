probPath <- function(path, Q){
  nTrans <- length(path)
  P <- vector("numeric", length(path))
  for(i in sequence(nTrans-1)){
    state_i <- as.numeric(names(path)[1])
    state_j <- as.numeric(names(path)[2])
    time_i <- as.numeric(path[1])
    time_j <- as.numeric(sum(path[-1]))
    rate_i <- abs(Q[state_i,state_j])
    rate_j <- abs(Q[state_j,state_j])
    P[i] <- dexp(time_i, rate_i) * (1 - pexp(rate_j, time_j))
    path <- path[-1]
  }
  state_j <- as.numeric(names(path))
  time_j <- as.numeric(path)
  rate_j <- abs(Q[state_j,state_j])
  P[nTrans] <- 1 - pexp(rate_j, time_j)
  P <- prod(P)
  if(P == 0){
    P <- 1e-100
  }
  return(P)
}

Q <- matrix(c(-1,1,1,-1), 2, 2)
path <- rep(0.01, 20)
names(path) <- rep(c(1,2), 10)
probPath(path, Q)

pathB <- sum(path)
names(pathB) <- 1
probPath(pathB, Q / 100)



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


opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)

dev.path <- function(p, fit.cor, simmap){
  Q <- fit.cor * p
  diag(Q) <- -rowSums(Q)
  return(-getMapProbability(simmap$maps, Q))
}

out <- nloptr(1, dev.path, opts = opts, fit.cor=fit.cor, simmap=data.houwie$simmap[[1]])
