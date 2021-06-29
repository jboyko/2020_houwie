# simulation tests for optimizing the optimization of hOUwie 

# get a phylogeny rescaled to a height of 1
getPhy <-function(nTip, b=1, d=0.5){
  phy <- NULL
  while(length(phy$tip.label) != nTip){
    phy <- sim.bdtree(b = b, d = d, stop = "taxa", n = nTip)
    phy <- drop.extinct(phy)
  }
  Tmax <- max(branching.times(phy))
  phy$edge.length <- phy$edge.length/Tmax
  return(phy)
}


require(corHMM)
require(OUwie)
require(parallel)
require(geiger)


nTip <- 500
phy <- getPhy(nTip)
# about 25 changes expected
p.mk <- c(0.5, 0.5)
Q = matrix(c(-p.mk[1],p.mk[2],p.mk[1],-p.mk[2]), 2, 2)
# always start in state 1
root.p = c(1, 0)
alpha = c(4, 4)
sig2= c(1, 1)
theta = c(3, 8)
theta0 = 3
rate.cat = 1
model.cor = "ER"
model.ou = "OUM"
data <- OUwie:::hOUwie.sim(phy, Q, root.p, alpha, sig2, theta0, theta)[[1]]
p <- c(p.mk[1], alpha[1], sig2[1], theta)

# find the minimum to a high precision
Opt.MinFm <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.8)
Search.MinFm <- OUwie:::hOUwie(phy, data, rate.cat, model.cor = model.cor, model.ou = model.ou, weighted = FALSE, nSim = 1000, nCores = 5, opts = Opt.MinFm)
  
StopVal <- Search.MinFm$objective + Search.MinFm$objective * .Machine$double.eps^0.5

# test out several algorithms under normal circumstance
AlgorList <- list(
  Opt.SBPLX=list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "stopval"=StopVal),
  Opt.COBYL=list("algorithm"="NLOPT_LN_COBYLA", "maxeval"="1000000", "stopval"=StopVal),
  Opt.BOBYQ=list("algorithm"="NLOPT_LN_BOBYQA", "maxeval"="1000000", "stopval"=StopVal),
  Opt.PRAXI=list("algorithm"="NLOPT_LN_PRAXIS", "maxeval"="1000000", "stopval"=StopVal),
  Opt.NELDE=list("algorithm"="NLOPT_LN_NELDERMEAD", "maxeval"="1000000", "stopval"=StopVal)
)

# AlgorList <- list(
#   Opt.SBPLX=list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5),
#   Opt.COBYL=list("algorithm"="NLOPT_LN_COBYLA", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5),
#   Opt.BOBYQ=list("algorithm"="NLOPT_LN_BOBYQA", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5),
#   Opt.PRAXI=list("algorithm"="NLOPT_LN_PRAXIS", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5),
#   Opt.NELDE=list("algorithm"="NLOPT_LN_NELDERMEAD", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5)
# )

Fit_100 <- mclapply(AlgorList, function(x) OUwie:::hOUwie(phy, data, rate.cat, model.cor = model.cor, model.ou = model.ou, weighted = FALSE, nSim = 100, nCores = 1, opts = x), mc.cores=5)

Fit_500 <- mclapply(AlgorList, function(x) OUwie:::hOUwie(phy, data, rate.cat, model.cor = model.cor, model.ou = model.ou, weighted = FALSE, nSim = 500, nCores = 1, opts = x), mc.cores=5)

Fit_1000 <- mclapply(AlgorList, function(x) OUwie:::hOUwie(phy, data, rate.cat, model.cor = model.cor, model.ou = model.ou, weighted = FALSE, nSim = 1000, nCores = 1, opts = x), mc.cores=5)


Search.SBPLX <- OUwie:::hOUwie(phy, data, rate.cat, model.cor = model.cor, model.ou = model.ou, weighted = FALSE, nSim = 100, nCores = 1, opts = Opt.SBPLX)

lapply(Fit_100, function(x) exp(x$solution) - p)
lapply(Fit_100, function(x) x$objective)
lapply(Fit_100, function(x) x$iterations)

