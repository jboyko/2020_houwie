probPathA <- function(path, Q){
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
  P <- sum(log(P))
  return(P)
}

probPathB <- function(path, Q){
  nTrans <- length(path)
  P <- vector("numeric", length(path))
  for(i in sequence(nTrans-1)){
    state_i <- as.numeric(names(path)[1])
    state_j <- as.numeric(names(path)[2])
    time_i <- as.numeric(path[1])
    rate_i <- abs(Q[state_i,state_j])
    P[i] <- dexp(time_i, rate_i)
    path <- path[-1]
  }
  state_j <- as.numeric(names(path))
  time_j <- as.numeric(path)
  rate_j <- abs(Q[state_j,state_j])
  P[nTrans] <- 1 - pexp(rate_j, time_j)
  P <- sum(log(P))
  return(P)
}

getMapProbability <- function(maps, Q, type = "A"){
  if(type == "A"){
    BranchProbs <- lapply(maps, function(x) probPathA(x, Q))
  }
  if(type == "B"){
    BranchProbs <- lapply(maps, function(x) probPathB(x, Q))
  }
  LnLik_map <- sum(unlist(BranchProbs))
  return(LnLik_map)
}

require(corHMM)
require(nloptr)

data(primates)
phy <- primates[[1]]
phy <- multi2di(phy)
data <- primates[[2]]
dat <- data.frame(sp = data[,1], d = rowSums(data[,c(2,3)]))

##run corhmm
MK <- corHMM(phy, dat, 1, model = "ER")
##get simmap from corhmm solution
model <- MK$solution
simmap <- makeSimmap(tree=phy, data=data, model=model, rate.cat=1, nSim=1, nCores=1)

Q <- MK$solution
Q[is.na(MK$solution)] <- 0
diag(Q) <- -rowSums(Q)

typeA <- getMapProbability(simmap[[1]]$maps, Q, "A")
typeB <- getMapProbability(simmap[[1]]$maps, Q, "B")

dat <- data.frame(sp = data[,1], d = rowSums(data[,c(2,3)]))

# what if we optimize a map 
optMap <- function(p, phy, dat, nMap, type){
  model <- getRateCatMat(3)
  model[model > 0] <- p
  diag(model) <- -rowSums(model)
  simmap <- makeSimmap(tree=phy, data=dat, model=model, rate.cat=1, nSim=nMap, nCores=1)
  # MK <- corHMM(phy, dat, 1, model = "ER", p = p)
  maps <- lapply(simmap, function(x) x$maps)
  llik <- max(unlist(lapply(maps, function(x) getMapProbability(x, model, type))))
  # llik <- llik + MK$loglik
  cat("\r", llik)
  return(-llik)
}

opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)

# with 0.1 we're good
outA <- nloptr(x0 = 0.1, eval_f = optMap, lb = 1e-5, ub = 1, opts = opts, phy = phy, nMap = 10, dat = dat, type = "A")
outB <- nloptr(x0 = 0.1, eval_f = optMap, lb = 1e-5, ub = 1, opts = opts, phy = phy, nMap = 10, dat = dat, type = "B")

# with 0.1 we're good
outA <- nloptr(x0 = 0.9, eval_f = optMap, lb = 1e-5, ub = 1, opts = opts, phy = phy, nMap = 10, dat = dat, type = "A")
outB <- nloptr(x0 = 0.1, eval_f = optMap, lb = 1e-5, ub = 1, opts = opts, phy = phy, nMap = 10, dat = dat, type = "B")


# LikTypeB <- sapply(xs, function(x) optMap(x, phy, 10, "B"))
# 
# xs <- seq(from = 0.5, to =1, length.out = 100)
# 
# par(mfrow=c(1,2))
# plot(y = -LikTypeB, x = log(xs))
# plot(y = -LikTypeB, x = xs)

# 
# 
# Q <- matrix(c(-1,1,1,-1), 2, 2)
# path <- rep(0.01, 20)
# names(path) <- rep(c(1,2), 10)
# 
# probPathA(path, Q/10)
# probPathB(path, Q/10)
# 
# probPathA(path, Q)
# probPathB(path, Q)
# 
# probPathA(path, Q*10)
# probPathB(path, Q*10)
# 
# probPathA(path, Q*100)
# probPathB(path, Q*100)
# 
# 
# 
# pathB <- sum(path)
# names(pathB) <- 1
# 
# probPathA(pathB, Q)
# probPathB(pathB, Q)
# 
# probPathA(pathB, Q*10)
# probPathB(pathB, Q*10)
# 
# probPathA(pathB, Q*100)
# probPathB(pathB, Q*100)
# 
# 

# 
# # proof of concept for the tip importance
# source("~/2020_hOUwie/hOUwieSimmap.R")
# source("~/2020_hOUwie/Utils.R")
# 
# require(OUwie)
# require(corHMM)
# require(parallel)
# require(phytools)
# require(expm)
# require(POUMM)
# require(geiger)
# 
# nTip <- 50
# phy <- sim.bdtree(b = 1, d = 0, stop = "taxa", n = nTip) 
# phy <- drop.extinct(phy)
# phy$edge.length <- phy$edge.length/max(branching.times(phy))
# 
# fit.cor <- equateStateMatPars(getRateCatMat(2), c(1,2))
# fit.cor <- getFullMat(list(fit.cor, fit.cor), fit.cor)
# fit.cor <- equateStateMatPars(fit.cor, c(1,2,3))
# fit.ou <- getOUParamStructure("OUM", "three.point", FALSE, FALSE, dim(fit.cor)[1])
# fit.ou[3,] <- c(3,3,4,4)
# pars = c(1, 5, 5, 2, 10)  # mk, alpha, sigma, theta1, theta2
# data.houwie <- generateData(phy, fit.cor, fit.ou, pars)
# cols<-setNames(c("gold","red", "purple", "black"),
#                c("1","2","3","4"))
# plotDataSet(data.houwie); legend("bottomleft", legend = c("1A","2A","1B","2B"), pch=16, col = cols)
# data <- data.houwie$data
# data[data[,2]==3,2] <- 1
# data[data[,2]==4,2] <- 2
# 
# Q <- fit.cor
# diag(Q) <- -rowSums(Q)
# pars = c(1, 5, 5, 2, 10)
# alpha = rep(pars[2], 4)
# sigma.sq = rep(pars[3], 4)
# theta = c(pars[4], pars[4], pars[5], pars[5])
# tip.paths <- lapply(1:length(phy$tip.label), function(x) OUwie:::getPathToRoot(phy, x))
# nSim = 100
# 
# OU.loglikTrue <- OUwie.basic(data.houwie$simmap[[1]], data.houwie$data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths, mserr="none")
# 
# 
# opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
# 
# dev.path <- function(p, fit.cor, simmap){
#   Q <- fit.cor * p
#   diag(Q) <- -rowSums(Q)
#   return(-getMapProbability(simmap$maps, Q))
# }
# 
# out <- nloptr(1, dev.path, opts = opts, fit.cor=fit.cor, simmap=data.houwie$simmap[[1]])

likelihood_given_number_transitions <- function(ntransitions=1, rate=0.1, time=50) {
  timeinterval <- time/(ntransitions+1)
  return( ((rate*exp(-rate*timeinterval))^ntransitions) * (exp(-rate*timeinterval)))
}
all_values <- expand.grid(ntransitions=c(0,4,100,1000), rate=c(10^seq(from=-2, to=.4, length.out=100)))
all_values <- all_values[order(all_values$ntransitions),]
all_values$likelihood <- NA
for (i in sequence(nrow(all_values))) {
  all_values$likelihood[i] <- likelihood_given_number_transitions(all_values$ntransitions[i], all_values$rate[i])
}
print(all_values[which.max(all_values$likelihood),])
all_values$lnL <- log(all_values$likelihood)
all_values$ntransitions <- as.factor(all_values$ntransitions)
library(ggplot2) 
ggplot(all_values, aes(x=rate, y=lnL, group=ntransitions)) + geom_line(aes(color=ntransitions))
