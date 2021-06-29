require(corHMM)
require(OUwie)
require(parallel)
require(geiger)
require(phytools)
source("~/2020_hOUwie/hOUwieSimmap.R")

nTip <- 250
q <- 1 
alpha <- c(1,1,1,2,1) 
sigma.sq <- c(1,2,3,2,1) 
theta = c(4, 5, 6, 2,2)

phy <- sim.bdtree(b = 1, d = 0.5, stop = "taxa", n = nTip) # this will be scaled to H=1
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))

Q = getRateCatMat(5)/10
diag(Q) <- -rowSums(Q)
# Q = matrix(c(-q,q,q,-q), 2, 2)

root.p = c(0,0,0,0,0) # we will sample the root with equal probability of either state
root.p[sample(c(1,2,3,4,5), 1, prob =c(0.5, 0.5, 0.5, 0.5, 0.5))] <- 1

theta0 = rnorm(1, theta[which(root.p == 1)], sqrt(sigma.sq[which(root.p == 1)]/2*alpha[which(root.p == 1)])) 
full.data <- hOUwie.sim(phy, Q, root.p, alpha, sigma.sq, theta0, theta)

plotSimmap(full.data$simmap[[1]])

OUwie.basic(full.data$simmap[[1]], full.data$data, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point")
OUwie.fixed(full.data$simmap[[1]], full.data$data, model = "OUMV", simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="invert")
OUwie.fixed(full.data$simmap[[1]], full.data$data, model = "OUMV", simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point")

# OUwie.basic(simmap[[1]], data.ou, simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm=algorithm, tip.paths=tip.paths, mserr=mserr)
# OUwie.fixed(simmap[[1]], data.ou, model = "OUM", simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="invert", tip.paths=tip.paths, mserr=mserr)
# OUwie.fixed(simmap[[1]], data.ou, model = "OUM", simmap.tree=TRUE, scaleHeight=FALSE, alpha=alpha, sigma.sq=sigma.sq, theta=theta, algorithm="three.point", tip.paths=tip.paths)
