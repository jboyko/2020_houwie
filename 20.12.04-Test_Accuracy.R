# This script will test for an improvement of bias/variance for parameters when doing corHMM then OUwie vs. hOUwie

require(corHMM)
require(OUwie)
require(parallel)
require(geiger)
require(proftools)

# generate data
phy <- sim.bdtree(b = 1, d = 0.5, stop = "taxa", n = 500)
phy <- drop.extinct(phy)
phy$edge.length <- phy$edge.length/max(branching.times(phy))
nSim <- 100

# Mk
root.p = c(1, 0)
p.mk <- c(1, 1)
Q = matrix(c(-p.mk[1],p.mk[2],p.mk[1],-p.mk[2]), 2, 2)
theta = c(3, 8)
theta0 = 3
  # OUM as true
alpha = c(2, 2)
sig2= c(1, 1)
data.mk.oum <- lapply(1:nSim, function(x) OUwie:::hOUwie.sim(phy, Q, root.p, alpha, sig2, theta0, theta, 1)[[1]])
  # OUMV as true
alpha = c(2, 2)
sig2= c(0.5, 2)
data.mk.oumv <- lapply(1:nSim, function(x) OUwie:::hOUwie.sim(phy, Q, root.p, alpha, sig2, theta0, theta)[[1]])
  # OUMVA as true
alpha = c(1, 4)
sig2= c(0.5, 2)
data.mk.oumva <- lapply(1:nSim, function(x) OUwie:::hOUwie.sim(phy, Q, root.p, alpha, sig2, theta0, theta)[[1]])

# HMM
root.p = c(1, 0, 0, 0)
p.mk <- c(4, 1, 2, 2)
Q = matrix(c(0, p.mk[1], p.mk[3], 0,
             p.mk[1], 0, 0, p.mk[3],
             p.mk[4], 0, 0, p.mk[2],
             0, p.mk[4], p.mk[2], 0), 4, 4, byrow = TRUE)
diag(Q) <- -rowSums(Q)

# State dependent as true
theta = c(3, 8, 3, 8)
theta0 = 3
# OUM as true
alpha = c(2, 2, 2, 2)
sig2= c(1, 1, 1, 1)
data.sd.oum <- lapply(1:nSim, function(x) OUwie:::hOUwie.sim(phy, Q, root.p, alpha, sig2, theta0, theta, 1)[[1]])
# OUMV as true
alpha = c(2, 2, 2, 2)
sig2= c(0.5, 1, 2, 4)
data.sd.oumv <- lapply(1:nSim, function(x) OUwie:::hOUwie.sim(phy, Q, root.p, alpha, sig2, theta0, theta)[[1]])
# OUMVA as true
alpha = c(1, 0.5, 2, 4)
sig2= c(0.5, 2, 1, 4)
data.sd.oumva <- lapply(1:nSim, function(x) OUwie:::hOUwie.sim(phy, Q, root.p, alpha, sig2, theta0, theta)[[1]])

# State independent as true
theta = c(3, 8, 3, 8)
theta0 = 3
# OUM as true
alpha = c(2, 2, 2, 2)
sig2= c(1, 1, 1, 1)
data.si.oum <- lapply(1:nSim, function(x) OUwie:::hOUwie.sim(phy, Q, root.p, alpha, sig2, theta0, theta, 1)[[1]])
# OUMV as true
alpha = c(2, 2, 2, 2)
sig2= c(0.5, 0.5, 2, 2)
data.si.oumv <- lapply(1:nSim, function(x) OUwie:::hOUwie.sim(phy, Q, root.p, alpha, sig2, theta0, theta)[[1]])
# OUMVA as true
alpha = c(1, 1, 4, 4)
sig2= c(0.5, 0.5, 2, 2)
data.si.oumva <- lapply(1:nSim, function(x) OUwie:::hOUwie.sim(phy, Q, root.p, alpha, sig2, theta0, theta)[[1]])
